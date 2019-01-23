#include "Macro.h"
#include "CUFLU.h"

#if ( defined GPU  &&  MODEL == ELBDM )



// useful macros
#define to1D1(z,y,x) ( __umul24(z, FLU_NXT*FLU_NXT) + __umul24(y, FLU_NXT) + x )
#define to1D2(z,y,x) ( __umul24(z-FLU_GHOST_SIZE, PS2*PS2) + __umul24(y-FLU_GHOST_SIZE, PS2) + x-FLU_GHOST_SIZE )

#ifdef LAPLACIAN_4TH
#  define LAP_GHOST  2
#  define LAP1(In,i)    (  real(1.0/ 12.0)*( - In[i-2] + (real)16.0*In[i-1] - (real)30.0*In[i  ] \
                                             - In[i+2] + (real)16.0*In[i+1] )  )
#  define LAP2(In,i)    (  real(1.0/144.0)*( + In[i-4] - (real)32.0*In[i-3] + (real)316.0*In[i-2] - (real)992.0*In[i-1] \
                                             + In[i+4] - (real)32.0*In[i+3] + (real)316.0*In[i+2] - (real)992.0*In[i+1] \
                                             +  (real)1414.0*In[i  ] )  )
#  ifndef CONSERVE_MASS
#  define LAP3(In,i)    (  real(1.0/1728.0)* \
       (  -In[i-6] + (real)48*In[i-5] - (real)858*In[i-4] + (real)7024*In[i-3] - (real)27279*In[i-2] + (real)58464*In[i-1] \
          -In[i+6] + (real)48*In[i+5] - (real)858*In[i+4] + (real)7024*In[i+3] - (real)27279*In[i+2] + (real)58464*In[i+1] \
          - (real)74796*In[i  ] )  )
#  endif

#else // #ifdef LAPLACIAN_4TH

#  define LAP_GHOST     1
#  define LAP1(In,i)    ( + In[i-1] - (real)2.0*In[i  ] + In[i+1] )
#  define LAP2(In,i)    ( + In[i-2] - (real)4.0*In[i-1] + (real)6.0*In[i  ] - (real)4.0*In[i+1] + In[i+2] )
#  ifndef CONSERVE_MASS
#  define LAP3(In,i)    ( + In[i-3] - (real)6.0*In[i-2] + (real)15.0*In[i-1] - (real)20.0*In[i  ] \
                          + In[i+3] - (real)6.0*In[i+2] + (real)15.0*In[i+1] )

#  endif

#endif // #ifdef LAPLACIAN_4TH ... else ...


static __device__ void CUFLU_Advance( real g_Fluid_In [][FLU_NIN ][ FLU_NXT*FLU_NXT*FLU_NXT ],
                                      real g_Fluid_Out[][FLU_NOUT][ PS2*PS2*PS2 ],
                                      real g_Flux     [][9][NFLUX_TOTAL][ PS2*PS2 ],
                                      const real dt, const real _dh, const real Eta1, const real Eta2, const bool StoreFlux,
                                      const real Taylor3_Coeff, const uint j_gap, const uint k_gap,
                                      real s_In[][FLU_BLOCK_SIZE_Y][FLU_NXT], real s_Half[][FLU_BLOCK_SIZE_Y][FLU_NXT],
                                      real s_Flux1[][PS2+1], real s_Flux2[][PS2+1], const bool FinalOut, const int XYZ, const real MinDens );




//-------------------------------------------------------------------------------------------------------
// Function    :  CUFLU_ELBDMSolver
// Description :  GPU ELBDM kinematic solver based on expanding the propagator to 3rd order
//
// Note        :  1. The three-dimensional evolution is achieved by applying x, y, and z operators successively.
//                   Since these operators commute, the order of applying them are irrelevant.
//                   --> Input pamameter "XYZ" is actually useless
//                   --> Nevertheless, the symmetry in different directions will be broken if CONSERVE_MASS is on
//                2. The implementation is very similar to the function " CUFLU_FluidSolver_RTVD"
//                4. Prefix "g" for pointers pointing to the "Global" memory space
//                   Prefix "s" for pointers pointing to the "Shared" memory space
//
// Parameter   :  g_Fluid_In     : Global memory array storing the input variables
//                g_Fluid_Out    : Global memory array to store the output variables
//                g_Flux         : Global memory array to store the output fluxes (useful only if StoreFlux == true)
//                dt             : Time interval to advance solution
//                _dh            : 1 / grid size
//                Eta            : Particle mass / Planck constant
//                StoreFlux      : true --> store the coarse-fine fluxes
//                                      --> useful only if CONSERVE_MASS is defined
//                Taylor3_Coeff  : Coefficient in front of the third term in the Taylor expansion
//                XYZ            : true  : x->y->z ( forward sweep)
//                                 false : z->y->x (backward sweep)
//                                 --> Meaningless if CONSERVE_MASS is off since the operators along different directions
//                                     commute
//                                 --> Meaningful if CONSERVE_MASS is on, in which the symmetry along different directions
//                                     are broken ...
//                MinDens        : Minimum allowed density
//-------------------------------------------------------------------------------------------------------
__global__ void CUFLU_ELBDMSolver( real g_Fluid_In [][FLU_NIN ][ FLU_NXT*FLU_NXT*FLU_NXT ],
                                   real g_Fluid_Out[][FLU_NOUT][ PS2*PS2*PS2 ],
                                   real g_Flux     [][9][NFLUX_TOTAL][ PS2*PS2 ],
                                   const real dt, const real _dh, const real Eta1, const real Eta2, const bool StoreFlux,
                                   const real Taylor3_Coeff, const bool XYZ, const real MinDens )
{

   __shared__ real s_In  [FLU_NIN][FLU_BLOCK_SIZE_Y][FLU_NXT];
#  ifdef CONSERVE_MASS
   __shared__ real s_Half[FLU_NIN][FLU_BLOCK_SIZE_Y][FLU_NXT];
   __shared__ real s_Flux1[FLU_BLOCK_SIZE_Y][PS2+1];
   __shared__ real s_Flux2[FLU_BLOCK_SIZE_Y][PS2+1];
#  else
   real (*s_Half)[FLU_BLOCK_SIZE_Y][FLU_NXT] = NULL;  // useless if CONSERVE_MASS is off
   real (*s_Flux1)[PS2+1]                     = NULL;  // useless if CONSERVE_MASS is off
   real (*s_Flux2)[PS2+1]                     = NULL;  // useless if CONSERVE_MASS is off
#  endif

   if ( XYZ )
   {
      CUFLU_Advance( g_Fluid_In, g_Fluid_Out, g_Flux, dt, _dh, Eta1, Eta2, StoreFlux, Taylor3_Coeff,
                                  0,              0, s_In, s_Half, s_Flux1, s_Flux2, false, 0, MinDens );
      CUFLU_Advance( g_Fluid_In, g_Fluid_Out, g_Flux, dt, _dh, Eta1, Eta2, StoreFlux, Taylor3_Coeff,
                     FLU_GHOST_SIZE,              0, s_In, s_Half, s_Flux1, s_Flux2, false, 3, MinDens );
      CUFLU_Advance( g_Fluid_In, g_Fluid_Out, g_Flux, dt, _dh, Eta1, Eta2, StoreFlux, Taylor3_Coeff,
                     FLU_GHOST_SIZE, FLU_GHOST_SIZE, s_In, s_Half, s_Flux1, s_Flux2,  true, 6, MinDens );
   }

   else
   {
      CUFLU_Advance( g_Fluid_In, g_Fluid_Out, g_Flux, dt, _dh, Eta1, Eta2, StoreFlux, Taylor3_Coeff,
                                  0,              0, s_In, s_Half, s_Flux1, s_Flux2, false, 6, MinDens );
      CUFLU_Advance( g_Fluid_In, g_Fluid_Out, g_Flux, dt, _dh, Eta1, Eta2, StoreFlux, Taylor3_Coeff,
                                  0, FLU_GHOST_SIZE, s_In, s_Half, s_Flux1, s_Flux2, false, 3, MinDens );
      CUFLU_Advance( g_Fluid_In, g_Fluid_Out, g_Flux, dt, _dh, Eta1, Eta2, StoreFlux, Taylor3_Coeff,
                     FLU_GHOST_SIZE, FLU_GHOST_SIZE, s_In, s_Half, s_Flux1, s_Flux2,  true, 0, MinDens );
   }

} // FUNCTION : CUFLU_ELBDMSolver



//-------------------------------------------------------------------------------------------------------
// Function    :  CUFLU_Advance
// Description :  Use GPU to advance solutions by one time-step
//
// Note        :  1. Based on expanding the kinematic propagator to 3rd order
//                2. Prefix "g" for pointers pointing to the "Global" memory space
//                   Prefix "s" for pointers pointing to the "Shared" memory space
//                3. The direction of the one dimensional sweep is determined by the input parameter "XYZ"
//
// Parameter   :  g_Fluid_In     : Global memory array storing the input variables
//                g_Fluid_Out    : Global memory array to store the output variables
//                g_Flux         : Global memory array to store the output fluxes (useful only if StoreFlux == true)
//                dt             : Time interval to advance solution
//                _dh            : 1 / grid size
//                Eta            : Particle mass / Planck constant
//                StoreFlux      : true --> store the coarse-fine fluxes
//                                   --> useful only if CONSERVE_MASS is defined
//                Taylor3_Coeff  : Coefficient in front of the third term in the Taylor expansion
//                j_gap          : Number of useless grids on each side in the j direction (j may not be equal to y)
//                k_gap          : Number of useless grids on each side in the k direction (k mya not be equal to z)
//                s_In           : Shared memory array to store the input data
//                s_Half         : Shared memory array to store the half-step solution
//                s_Flux         : Shared memory array to store the boundary fluxes
//                FinalOut       : true --> store the updated data to g_Fluid_Out
//                XYZ            : 0 : Update the solution in the x direction
//                                 3 : Update the solution in the y direction
//                                 6 : Update the solution in the z direction
//                                 --> This parameter is also used to determine the place to store the output fluxes
//                MinDens        : Minimum allowed density
//-------------------------------------------------------------------------------------------------------
__device__ void CUFLU_Advance( real g_Fluid_In [][FLU_NIN ][ FLU_NXT*FLU_NXT*FLU_NXT ],
                               real g_Fluid_Out[][FLU_NOUT][ PS2*PS2*PS2 ],
                               real g_Flux     [][9][NFLUX_TOTAL][ PS2*PS2 ],
                               const real dt, const real _dh, const real Eta1, const real Eta2, const bool StoreFlux, const real Taylor3_Coeff,
                               const uint j_gap, const uint k_gap, real s_In[][FLU_BLOCK_SIZE_Y][FLU_NXT],
                               real s_Half[][FLU_BLOCK_SIZE_Y][FLU_NXT], real s_Flux1[][PS2+1], real s_Flux2[][PS2+1], const bool FinalOut,
                               const int XYZ, const real MinDens )
{

   const real _Eta1         = (real)1.0/Eta1;
   const real _Eta2         = (real)1.0/Eta2;
   const real dT1           = (real)0.5*dt*_Eta1;
   const real dT2           = (real)0.5*dt*_Eta2;
   const real _Eta12_dh     = (real)0.5*_dh*_Eta1;
   const real _Eta22_dh     = (real)0.5*_dh*_Eta2;
   const real Coeff11       = dT1*_dh*_dh;
   const real Coeff12       = dT2*_dh*_dh;
#  ifdef CONSERVE_MASS
   const real Coeff21       = Taylor3_Coeff*SQR(Coeff11);
   const real Coeff22       = Taylor3_Coeff*SQR(Coeff12);
#  else
   const real Coeff21       = (real)0.5*SQR(Coeff11);
   const real Coeff22       = (real)0.5*SQR(Coeff12);
   const real Coeff31       = Taylor3_Coeff*CUBE(Coeff11);
   const real Coeff32       = Taylor3_Coeff*CUBE(Coeff12);
#  endif

   const uint bx           = blockIdx.x;
   const uint tx           = threadIdx.x;
   const uint ty           = threadIdx.y;
   const uint tid          = __umul24(ty,FLU_BLOCK_SIZE_X) + tx;
   const uint size_j       = FLU_NXT - (j_gap<<1);
   const uint size_k       = FLU_NXT - (k_gap<<1);
   const uint NColumnTotal = __umul24( size_j, size_k );    // total number of data columns to be updated
   const uint i            = tx + FLU_GHOST_SIZE;           // (i,j,k): array indices used in g_Fluid_In
   const uint j_end       = FLU_NXT - j_gap;
         uint j           = j_gap + ty%size_j;
         uint k           = k_gap + ty/size_j;
         uint Column0     = 0;                              // the total number of columns that have been updated
         uint NColumnOnce = MIN( NColumnTotal, FLU_BLOCK_SIZE_Y );

   double Amp1_New, Amp2_New;            // use double precision to reduce the round-off error in the mass conservation
   real   Re1_Old, Im1_Old, Re1_New, Im1_New;
   real   Re2_Old, Im2_Old, Re2_New, Im2_New;
   uint   Idx1, Idx2, Idx3, delta_k;

#  ifdef CONSERVE_MASS
   const uint NThread     = FLU_BLOCK_SIZE_X*FLU_BLOCK_SIZE_Y;
   const uint NHalf       = FLU_NXT - 4*LAP_GHOST;
   const real dT1_dh2     = dT1*_dh*_dh;
   const real dT2_dh2     = dT2*_dh*_dh;
   const uint txp         = tx + 1;

   double Amp1_Old, Amp2_Old, Amp1_Corr, Amp2_Corr;  // use double precision to reduce the round-off error in the mass conservation
   real   R1, I1, dR1, dI1;
   real   R2, I2, dR2, dI2;
   uint   Idx;
   uint   si, sj;                                           // array indices used in the shared memory array
   uint   f, fp1;                                           // array indices used in the s_Flux array
#  ifdef LAPLACIAN_4TH
   uint   fm1, fp2;
#  endif
#  endif // #ifdef CONSERVE_MASS


// determine the array indices for loading the ghost-zone data
   bool LoadGhost = false;                                  // true --> load the ghost-zone data
   uint LoadGhost_i;
   int  LoadGhost_di, LoadGhost_dIdx1;

   if ( tx < 2*FLU_GHOST_SIZE )
   {
      LoadGhost = true;

      if ( tx < FLU_GHOST_SIZE )    LoadGhost_di = -FLU_GHOST_SIZE;
      else                          LoadGhost_di = -FLU_GHOST_SIZE + PS2;

      switch ( XYZ )
      {
         case 0:  LoadGhost_dIdx1 = LoadGhost_di;                                break;
         case 3:  LoadGhost_dIdx1 = __mul24( LoadGhost_di, FLU_NXT );            break;
         case 6:  LoadGhost_dIdx1 = __mul24( LoadGhost_di, FLU_NXT*FLU_NXT );    break;
      }

      LoadGhost_i = (int)i + LoadGhost_di;
   } // if ( tx < 2*FLU_GHOST_SIZE )


// loop over all data columns
   while ( Column0 < NColumnTotal )
   {
//    1. load data into shared memory
      if ( tid < NColumnOnce*PS2 )
      {
//       1.1 determine the array indices for loading global memory data along different directions
         switch ( XYZ )
         {
            case 0:  Idx1 = to1D1( k, j, i );    break;
            case 3:  Idx1 = to1D1( k, i, j );    break;
            case 6:  Idx1 = to1D1( i, k, j );    break;
         }

//       1.2 load the interior data into shared memory
         Re1_Old = g_Fluid_In[bx][0][Idx1];
         Im1_Old = g_Fluid_In[bx][1][Idx1];
         Re2_Old = g_Fluid_In[bx][2][Idx1];
         Im2_Old = g_Fluid_In[bx][3][Idx1];

         s_In[0][ty][i] = Re1_Old;
         s_In[1][ty][i] = Im1_Old;
         s_In[2][ty][i] = Re2_Old;
         s_In[3][ty][i] = Im2_Old;

//       1.3 load the ghost-zone data into shared memory
         if ( LoadGhost )
         {
            s_In[0][ty][LoadGhost_i] = g_Fluid_In[bx][0][ (int)Idx1 + LoadGhost_dIdx1 ];
            s_In[1][ty][LoadGhost_i] = g_Fluid_In[bx][1][ (int)Idx1 + LoadGhost_dIdx1 ];
            s_In[2][ty][LoadGhost_i] = g_Fluid_In[bx][2][ (int)Idx1 + LoadGhost_dIdx1 ];
            s_In[3][ty][LoadGhost_i] = g_Fluid_In[bx][3][ (int)Idx1 + LoadGhost_dIdx1 ];
         }
      } // if ( tid < NColumnOnce*PS2 )

      __syncthreads();


#     ifdef CONSERVE_MASS


//    2. half-step solution
      Idx = tid;
      while ( Idx < NColumnOnce*NHalf )
      {
         si = Idx % NHalf + 2*LAP_GHOST;
         sj = Idx / NHalf;

         s_Half[0][sj][si] = s_In[0][sj][si] - (real)0.5*Coeff11*LAP1( s_In[1][sj], si ) - Coeff21*LAP2( s_In[0][sj], si );
         s_Half[1][sj][si] = s_In[1][sj][si] + (real)0.5*Coeff11*LAP1( s_In[0][sj], si ) - Coeff21*LAP2( s_In[1][sj], si );
         s_Half[2][sj][si] = s_In[2][sj][si] - (real)0.5*Coeff12*LAP1( s_In[3][sj], si ) - Coeff22*LAP2( s_In[2][sj], si );
         s_Half[3][sj][si] = s_In[3][sj][si] + (real)0.5*Coeff12*LAP1( s_In[2][sj], si ) - Coeff22*LAP2( s_In[3][sj], si );

         Idx += NThread;
      } // while ( Idx < NColumnOnce*NHalf )

      __syncthreads();


//    3. calculate the face-center fluxes (the coefficient _dh has been absorted into the constant dT_dh2)
      Idx = tid;
      while ( Idx < NColumnOnce*(PS2+1) )
      {
         si  = Idx % (PS2+1);
         sj  = Idx / (PS2+1);
         f   = si + FLU_GHOST_SIZE - 1;
         fp1 = f + 1;

#        ifdef LAPLACIAN_4TH
         fm1 = f - 1;
         fp2 = f + 2;

         R1  = real(1./28.)*( -s_Half[0][sj][fm1]+(real)15*s_Half[0][sj][f]+(real)15*s_Half[0][sj][fp1]-s_Half[0][sj][fp2] );
         I1  = real(1./28.)*( -s_Half[1][sj][fm1]+(real)15*s_Half[1][sj][f]+(real)15*s_Half[1][sj][fp1]-s_Half[1][sj][fp2] );
         dR1 = real(1./12.)*( +s_Half[0][sj][fm1]-(real)15*s_Half[0][sj][f]+(real)15*s_Half[0][sj][fp1]-s_Half[0][sj][fp2] );
         dI1 = real(1./12.)*( +s_Half[1][sj][fm1]-(real)15*s_Half[1][sj][f]+(real)15*s_Half[1][sj][fp1]-s_Half[1][sj][fp2] );
         R2  = real(1./28.)*( -s_Half[2][sj][fm1]+(real)15*s_Half[2][sj][f]+(real)15*s_Half[2][sj][fp1]-s_Half[2][sj][fp2] );
         I2  = real(1./28.)*( -s_Half[3][sj][fm1]+(real)15*s_Half[3][sj][f]+(real)15*s_Half[3][sj][fp1]-s_Half[3][sj][fp2] );
         dR2 = real(1./12.)*( +s_Half[2][sj][fm1]-(real)15*s_Half[2][sj][f]+(real)15*s_Half[2][sj][fp1]-s_Half[2][sj][fp2] );
         dI2 = real(1./12.)*( +s_Half[3][sj][fm1]-(real)15*s_Half[3][sj][f]+(real)15*s_Half[3][sj][fp1]-s_Half[3][sj][fp2] );

#        else

         R1  = real(0.5)*( + s_Half[0][sj][f] + s_Half[0][sj][fp1] );
         I1  = real(0.5)*( + s_Half[1][sj][f] + s_Half[1][sj][fp1] );
         dR1 =           ( - s_Half[0][sj][f] + s_Half[0][sj][fp1] );
         dI1 =           ( - s_Half[1][sj][f] + s_Half[1][sj][fp1] );
         R2  = real(0.5)*( + s_Half[2][sj][f] + s_Half[2][sj][fp1] );
         I2  = real(0.5)*( + s_Half[3][sj][f] + s_Half[3][sj][fp1] );
         dR2 =           ( - s_Half[2][sj][f] + s_Half[2][sj][fp1] );
         dI2 =           ( - s_Half[3][sj][f] + s_Half[3][sj][fp1] );
#        endif

         s_Flux1[sj][si] = (real)2.0*( R1*dI1 - I1*dR1 );
         s_Flux2[sj][si] = (real)2.0*( R2*dI2 - I2*dR2 );
         Idx += NThread;
      } // while ( Idx < NColumnOnce*(PS2+1) )

      __syncthreads();


//    4a. full-step solution (equivalent to the 3rd-order Taylor expansion)
      if ( tid < NColumnOnce*PS2 )
      {
         Re1_New   = Re1_Old - Coeff11*LAP1( s_Half[1][ty], i );
         Im1_New   = Im1_Old + Coeff11*LAP1( s_Half[0][ty], i );
         Re2_New   = Re2_Old - Coeff12*LAP1( s_Half[3][ty], i );
         Im2_New   = Im2_Old + Coeff12*LAP1( s_Half[2][ty], i );

         Amp1_Old  = SQR( Re1_Old ) + SQR( Im1_Old );
         Amp2_Old  = SQR( Re2_Old ) + SQR( Im2_Old );
         Amp1_New  = SQR( Re1_New ) + SQR( Im1_New );
         Amp2_New  = SQR( Re2_New ) + SQR( Im2_New );
         Amp1_Corr = Amp1_Old - dT1_dh2*( s_Flux1[ty][txp] - s_Flux1[ty][tx] );
         Amp2_Corr = Amp2_Old - dT2_dh2*( s_Flux2[ty][txp] - s_Flux2[ty][tx] );

//       be careful about the negative density and the vacuum (where we might have Amp_New == 0.0)
//       if ( Amp_Corr > (real)0.0  &&  Amp_New > (real)0.0 )
         if ( Amp1_Corr >       0.0  &&  Amp1_New >       0.0 )
         {
            /*
            Re_New *= SQRT( Amp_Corr / Amp_New );
            Im_New *= SQRT( Amp_Corr / Amp_New );
            */
            Re1_New *= sqrt( Amp1_Corr / Amp1_New );  // use double precision to improve the mass conservation further
            Im1_New *= sqrt( Amp1_Corr / Amp1_New );
            Amp1_New = Amp1_Corr;
         }
         if ( Amp2_Corr >       0.0  &&  Amp2_New >       0.0 )
         {
            Re2_New *= sqrt( Amp2_Corr / Amp2_New );
            Im2_New *= sqrt( Amp2_Corr / Amp2_New );
            Amp2_New = Amp2_Corr;
         }
      } // if if ( tid < NColumnOnce*PS2 )


#     else // CONSERVE_MASS


//    4b. full-step solution if CONSERVE_MASS is not defined (equivalent to the 3rd-order Taylor expansion)
      if ( tid < NColumnOnce*PS2 )
      {
         Re1_New  = Re1_Old - Coeff11*LAP1( s_In[1][ty], i ) - Coeff21*LAP2( s_In[0][ty], i ) + Coeff31*LAP3( s_In[1][ty], i );
         Im1_New  = Im1_Old + Coeff11*LAP1( s_In[0][ty], i ) - Coeff21*LAP2( s_In[1][ty], i ) - Coeff31*LAP3( s_In[0][ty], i );
         Re2_New  = Re2_Old - Coeff12*LAP1( s_In[3][ty], i ) - Coeff22*LAP2( s_In[2][ty], i ) + Coeff32*LAP3( s_In[3][ty], i );
         Im2_New  = Im2_Old + Coeff12*LAP1( s_In[2][ty], i ) - Coeff22*LAP2( s_In[3][ty], i ) - Coeff32*LAP3( s_In[2][ty], i );
         Amp1_New = SQR( Re1_New ) + SQR( Im1_New );
         Amp2_New = SQR( Re2_New ) + SQR( Im2_New );
      }


#     endif // CONSERVE_MASS ... else ...


//    5. store the updated data (and fluxes) back to the global memory
      if ( tid < NColumnOnce*PS2 )
      {
//       5.1 data
         if ( FinalOut )
         {
//          apply the the minimum density check
//          --> to be consistent with the CPU solver, we apply it just before storing the output results to g_Fluid_Out
            if ( Amp1_New < MinDens )
            {
               const real Rescale1 = SQRT( MinDens / (real)Amp1_New );

               Re1_New *= Rescale1;
               Im1_New *= Rescale1;
               Amp1_New = MinDens;
            }
            if ( Amp2_New < MinDens )
            {
               const real Rescale2 = SQRT( MinDens / (real)Amp2_New );

               Re2_New *= Rescale2;
               Im2_New *= Rescale2;
               Amp2_New = MinDens;
            }

            switch ( XYZ )
            {
               case 0:  Idx2 = to1D2( k, j, i );    break;
               case 3:  Idx2 = to1D2( k, i, j );    break;
               case 6:  Idx2 = to1D2( i, k, j );    break;
            }

            g_Fluid_Out[bx][0][Idx2] = Amp1_New;
            g_Fluid_Out[bx][1][Idx2] = Re1_New;
            g_Fluid_Out[bx][2][Idx2] = Im1_New;
            g_Fluid_Out[bx][3][Idx2] = Amp2_New;
            g_Fluid_Out[bx][4][Idx2] = Re2_New;
            g_Fluid_Out[bx][5][Idx2] = Im2_New;
         }

         else
         {
            g_Fluid_In[bx][0][Idx1] = Re1_New;
            g_Fluid_In[bx][1][Idx1] = Im1_New;
            g_Fluid_In[bx][2][Idx1] = Re2_New;
            g_Fluid_In[bx][3][Idx1] = Im2_New;
         }


//       5.2 fluxes (for the flux-correction operation)
         if ( StoreFlux  &&  tx == 0 )
         if ( k >= FLU_GHOST_SIZE  &&  k < FLU_NXT-FLU_GHOST_SIZE )
         if ( j >= FLU_GHOST_SIZE  &&  j < FLU_NXT-FLU_GHOST_SIZE )
         {
            Idx3 = __umul24( k-FLU_GHOST_SIZE, PS2 ) + (j-FLU_GHOST_SIZE);

            g_Flux[bx][XYZ+0][0][Idx3] = s_Flux1[ty][  0]*_Eta12_dh;
            g_Flux[bx][XYZ+1][0][Idx3] = s_Flux1[ty][PS1]*_Eta12_dh;
            g_Flux[bx][XYZ+2][0][Idx3] = s_Flux1[ty][PS2]*_Eta12_dh;
            g_Flux[bx][XYZ+0][1][Idx3] = s_Flux2[ty][  0]*_Eta22_dh;
            g_Flux[bx][XYZ+1][1][Idx3] = s_Flux2[ty][PS1]*_Eta22_dh;
            g_Flux[bx][XYZ+2][1][Idx3] = s_Flux2[ty][PS2]*_Eta22_dh;
         }


//       5.3 reset the target array indices
         j += NColumnOnce;

         if ( j >= j_end )
         {
            delta_k  = ( j - j_end )/size_j + 1;
            k       += delta_k;
            j       -= __umul24( size_j, delta_k );
         }
      } // if ( tid < NColumnOnce*PS2 )

      __syncthreads();

      Column0     += NColumnOnce;
      NColumnOnce  = MIN( NColumnTotal - Column0, FLU_BLOCK_SIZE_Y );

   } // while ( Column0 < NColumnTotal )

} // FUNCTION : CUFLU_Advance



#endif // #if ( defined GPU  &&  MODEL == ELBDM )
