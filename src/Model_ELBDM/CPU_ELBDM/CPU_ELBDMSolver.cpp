#include "GAMER.h"
#include "CUFLU.h"

#if ( !defined GPU  &&  MODEL == ELBDM )



// useful macros
#define to1D(z,y,x)     ( z*FLU_NXT*FLU_NXT + y*FLU_NXT + x )

#ifdef LAPLACIAN_4TH
#  define LAP1(In,t)    (  real(1.0/ 12.0)*( - In[t-2] + (real)16.0*In[t-1] - (real)30.0*In[t] + \
                                             - In[t+2] + (real)16.0*In[t+1] )  )
#  define LAP2(In,t)    (  real(1.0/144.0)*(  In[t-4] - (real)32.0*In[t-3] + (real)316.0*In[t-2] - (real)992.0*In[t-1] + \
                                              In[t+4] - (real)32.0*In[t+3] + (real)316.0*In[t+2] - (real)992.0*In[t+1] + \
                                              (real)1414.0*In[t] )  )
#else
#  define LAP1(In,t)    ( In[t-1] - (real)2.0*In[t] + In[t+1] )
#  define LAP2(In,t)    ( In[t-2] - (real)4.0*In[t-1] + (real)6.0*In[t] - (real)4.0*In[t+1] + In[t+2] )
#endif


static void CPU_AdvanceX( real u[][ FLU_NXT*FLU_NXT*FLU_NXT ], real Flux_Array[][NFLUX_TOTAL][ PS2*PS2 ],
                          const real dt, const real dh, const real Eta1, const real Eta2, const bool StoreFlux, const real Taylor3_Coeff,
                          const int j_gap, const int k_gap, const int Flux_XYZ );
static void TransposeXY( real u[][ FLU_NXT*FLU_NXT*FLU_NXT ] );
static void TransposeXZ( real u[][ FLU_NXT*FLU_NXT*FLU_NXT ] );




//-------------------------------------------------------------------------------------------------------
// Function    :  CPU_ELBDMSolver
// Description :  CPU ELBDM kinematic solver based on expanding the propagator to the 3rd order
//
// Note        :  1. The three-dimensional evolution is achieved by applying x, y, and z operators successively.
//                   Since these operators commute, the order of applying them are irrelevant.
//                   --> Input pamameter "XYZ" is actually meaningless (if CONSERVE_MASS is off)
//                   --> Nevertheless, the symmetry in different directions will be broken if CONSERVE_MASS is on
//                2. The implementation is very similar to the function "CPU_FluidSolver_RTVD"
//
// Parameter   :  Flu_Array_In   : Array storing the input variables (only REAL1/IMAG1/REAL2/IMAG2)
//                Flu_Array_Out  : Array to store the output variables (DENS/REAL1/IMAG1/REAL2/IMAG2)
//                Flux_Array     : Array to store the output flux
//                NPatchGroup    : Number of patch groups to be evaluated
//                dt             : Time interval to advance solution
//                dh             : Grid size
//                Eta1           : Particle mass 1 / Planck constant
//                Eta2           : Particle mass 2 / Planck constant
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
void CPU_ELBDMSolver( real Flu_Array_In [][FLU_NIN ][ FLU_NXT*FLU_NXT*FLU_NXT ],
                      real Flu_Array_Out[][FLU_NOUT][ PS2*PS2*PS2 ],
                      real Flux_Array[][9][NFLUX_TOTAL][ PS2*PS2 ],
                      const int NPatchGroup, const real dt, const real dh, const real Eta1, const real Eta2, const bool StoreFlux,
                      const real Taylor3_Coeff, const bool XYZ, const real MinDens )
{

   if ( XYZ )
   {
#     pragma omp parallel for schedule( runtime )
      for (int P=0; P<NPatchGroup; P++)
      {
         CPU_AdvanceX( Flu_Array_In[P], Flux_Array[P], dt, dh, Eta1, Eta2, StoreFlux, Taylor3_Coeff,
                                    0,              0, 0 );

         TransposeXY ( Flu_Array_In[P] );

         CPU_AdvanceX( Flu_Array_In[P], Flux_Array[P], dt, dh, Eta1, Eta2, StoreFlux, Taylor3_Coeff,
                       FLU_GHOST_SIZE,              0, 3 );

         TransposeXZ ( Flu_Array_In[P] );

         CPU_AdvanceX( Flu_Array_In[P], Flux_Array[P], dt, dh, Eta1, Eta2, StoreFlux, Taylor3_Coeff,
                       FLU_GHOST_SIZE, FLU_GHOST_SIZE, 6 );

         TransposeXZ ( Flu_Array_In[P] );
         TransposeXY ( Flu_Array_In[P] );
      }
   }

   else
   {
#     pragma omp parallel for schedule( runtime )
      for (int P=0; P<NPatchGroup; P++)
      {
         TransposeXY ( Flu_Array_In[P] );
         TransposeXZ ( Flu_Array_In[P] );

         CPU_AdvanceX( Flu_Array_In[P], Flux_Array[P], dt, dh, Eta1, Eta2, StoreFlux, Taylor3_Coeff,
                                    0,              0, 6 );

         TransposeXZ ( Flu_Array_In[P] );

         CPU_AdvanceX( Flu_Array_In[P], Flux_Array[P], dt, dh, Eta1, Eta2, StoreFlux, Taylor3_Coeff,
                                    0, FLU_GHOST_SIZE, 3 );

         TransposeXY ( Flu_Array_In[P] );

         CPU_AdvanceX( Flu_Array_In[P], Flux_Array[P], dt, dh, Eta1, Eta2, StoreFlux, Taylor3_Coeff,
                       FLU_GHOST_SIZE, FLU_GHOST_SIZE, 0 );
      }
   }


// copy the updated data to Flu_Array_Out
   int  Idx1, Idx2, v_m1;
   real Amp, Rescale;   // not using double precision since MinDens will break the mass conservation anyway

#  pragma omp parallel for private( Idx1, Idx2, v_m1, Amp, Rescale ) schedule( runtime )
   for (int P=0; P<NPatchGroup; P++)
   {
//    copy data
      for (int v=1; v<FLU_NOUT; v++)
      {
         v_m1 = v-1;
         Idx1 = 0;

         for (int k=FLU_GHOST_SIZE; k<FLU_GHOST_SIZE+PS2; k++)
         for (int j=FLU_GHOST_SIZE; j<FLU_GHOST_SIZE+PS2; j++)
         for (int i=FLU_GHOST_SIZE; i<FLU_GHOST_SIZE+PS2; i++)
         {
            Idx2 = to1D(k,j,i);

            Flu_Array_Out[P][v][ Idx1++ ] = Flu_Array_In[P][v_m1][Idx2];
         }
      }

//    evaluate the new density (and apply the minimum density check)
      for (int t=0; t<PS2*PS2*PS2; t++)
      {
         Amp = SQR( Flu_Array_Out[P][1][t] ) + SQR( Flu_Array_Out[P][2][t] ) + SQR( Flu_Array_Out[P][3][t] ) + SQR( Flu_Array_Out[P][4][t] );

         if ( Amp < MinDens )
         {
            Rescale                 = SQRT( MinDens / Amp );
            Flu_Array_Out[P][1][t] *= Rescale;
            Flu_Array_Out[P][2][t] *= Rescale;
            Flu_Array_Out[P][3][t] *= Rescale;
            Flu_Array_Out[P][4][t] *= Rescale;
            Amp                     = MinDens;
         }

         Flu_Array_Out[P][0][t] = Amp;
      }
   } // for (int P=0; P<NPatchGroup; P++)

} // FUNCTION : CPU_ELBDMSolver



//-------------------------------------------------------------------------------------------------------
// Function    :  CPU_AdvanceX
// Description :  Use CPU to advance a single patch group by one time-step in the x direction
//
// Note        :  Based on expanding the kinematic propagator to the 3rd order
//
// Parameter   :  u              : Array storing the input variables (only REAL/IMAG)
//                Flux_Array     : Array to store the output flux (only density)
//                dt             : Time interval to advance solution
//                dh             : Grid size
//                Eta1           : Particle mass 1 / Planck constant
//                Eta2           : Particle mass 2 / Planck constant
//                StoreFlux      : true --> store the coarse-fine fluxes
//                                      --> useful only if CONSERVE_MASS is defined
//                Taylor3_Coeff  : Coefficient in front of the third term in the Taylor expansion
//                j_gap          : Number of cells to be skipped on each side in the y direction
//                k_gap          : Number of cells to be skipped on each side in the z direction
//                Flux_XYZ       : Parameter used to determine the place to store the output fluxes
//                                 --> (0,3,6) <-> (x/y/z) fluxes
//                                 --> useful only if CONSERVE_MASS is defined
//-------------------------------------------------------------------------------------------------------
void CPU_AdvanceX( real u[][ FLU_NXT*FLU_NXT*FLU_NXT ], real Flux_Array[][NFLUX_TOTAL][ PS2*PS2 ],
                   const real dt, const real dh, const real Eta1, const real Eta2, const bool StoreFlux, const real Taylor3_Coeff,
                   const int j_gap, const int k_gap, const int Flux_XYZ )
{

   const real _dh      = (real)1.0/dh;
   const real dT1      = (real)0.5*dt/Eta1;
   const real dT2      = (real)0.5*dt/Eta2;
   const real _Eta12_dh= (real)0.5*_dh/Eta1;
   const real _Eta22_dh= (real)0.5*_dh/Eta2;
   const real Coeff11   = dT1*_dh*_dh;
   const real Coeff12   = dT2*_dh*_dh;
   const real Coeff21   = Taylor3_Coeff*Coeff11*Coeff11;
   const real Coeff22   = Taylor3_Coeff*Coeff12*Coeff12;
   const int j_start   = j_gap;
   const int k_start   = k_gap;
   const int j_end     = FLU_NXT - j_gap;
   const int k_end     = FLU_NXT - k_gap;

   real Re_Old1 [FLU_NXT];  // one column of the real      part in the input array "u"
   real Im_Old1 [FLU_NXT];  // one column of the imaginary part in the input array "u"
   real Re_Half1[FLU_NXT];  // one column of the real      part at the half time-step
   real Im_Half1[FLU_NXT];  // one column of the imaginary part at the half time-step
   real *Re_New1 = NULL;    // pointer to store the full-step real      part
   real *Im_New1 = NULL;    // pointer to store the full-step imaginary part
   real Re_Old2 [FLU_NXT];  // one column of the real      part in the input array "u"
   real Im_Old2 [FLU_NXT];  // one column of the imaginary part in the input array "u"
   real Re_Half2[FLU_NXT];  // one column of the real      part at the half time-step
   real Im_Half2[FLU_NXT];  // one column of the imaginary part at the half time-step
   real *Re_New2 = NULL;    // pointer to store the full-step real      part
   real *Im_New2 = NULL;    // pointer to store the full-step imaginary part
   int Idx;

#  ifdef CONSERVE_MASS
   const real dT1_dh2 = dT1*_dh*_dh;
   const real dT2_dh2 = dT2*_dh*_dh;
   real   R1, I1, R2, I2, dR1, dI1, dR2, dI2, Flux1[PS2+1], Flux2[PS2+1];
   double Amp_Old, Amp_New, Amp_Corr;  // use double precision to reduce the round-off error in the mass conservation
   int    Idx2, Idx3;
#  endif


// loop over all target columns
   for (int k=k_start; k<k_end; k++)
   for (int j=j_start; j<j_end; j++)
   {
//    1. backup one column of data
//    ------------------------------------------------------------------------------------------------------------
      Idx    = to1D(k,j,0);
      Re_New1 = &u[0][Idx];
      Im_New1 = &u[1][Idx];
      Re_New2 = &u[2][Idx];
      Im_New2 = &u[3][Idx];

      memcpy( Re_Old1, Re_New1, FLU_NXT*sizeof(real) );
      memcpy( Im_Old1, Im_New1, FLU_NXT*sizeof(real) );
      memcpy( Re_Old2, Re_New2, FLU_NXT*sizeof(real) );
      memcpy( Im_Old2, Im_New2, FLU_NXT*sizeof(real) );


//    2. half-step solution
//    ------------------------------------------------------------------------------------------------------------
#     ifdef LAPLACIAN_4TH
      for (int i=4; i<FLU_NXT-4; i++)
#     else
      for (int i=2; i<FLU_NXT-2; i++)
#     endif
      {
         Re_Half1[i] = Re_Old1[i] - (real)0.5*Coeff11*LAP1( Im_Old1, i ) - Coeff21*LAP2( Re_Old1, i );
         Im_Half1[i] = Im_Old1[i] + (real)0.5*Coeff11*LAP1( Re_Old1, i ) - Coeff21*LAP2( Im_Old1, i );
         Re_Half2[i] = Re_Old2[i] - (real)0.5*Coeff12*LAP1( Im_Old2, i ) - Coeff22*LAP2( Re_Old2, i );
         Im_Half2[i] = Im_Old2[i] + (real)0.5*Coeff12*LAP1( Re_Old2, i ) - Coeff22*LAP2( Im_Old2, i );
      }


//    3. full-step solution (equivalent to the 3rd-order Taylor expansion)
//    ------------------------------------------------------------------------------------------------------------
      for (int i=FLU_GHOST_SIZE; i<FLU_NXT-FLU_GHOST_SIZE; i++)
      {
         Re_New1[i] = Re_Old1[i] - Coeff11*LAP1( Im_Half1, i );
         Im_New1[i] = Im_Old1[i] + Coeff11*LAP1( Re_Half1, i );
         Re_New2[i] = Re_Old2[i] - Coeff12*LAP1( Im_Half2, i );
         Im_New2[i] = Im_Old2[i] + Coeff12*LAP1( Re_Half2, i );
      }


//    4. enforce the mass conservation by solving the continuity eq.
//    ------------------------------------------------------------------------------------------------------------
#     ifdef CONSERVE_MASS
//    4.1. calculate the face-center fluxes (the coefficient _dh has been absorted into the constant dT_dh2)
      Idx2 = 0;
      for (int i=FLU_GHOST_SIZE-1; i<FLU_NXT-FLU_GHOST_SIZE; i++)
      {
#        ifdef LAPLACIAN_4TH
         R1  = real(1.0/28.0)*( - Re_Half1[i-1] + (real)15.0*Re_Half1[i] + (real)15.0*Re_Half1[i+1] - Re_Half1[i+2] );
         I1  = real(1.0/28.0)*( - Im_Half1[i-1] + (real)15.0*Im_Half1[i] + (real)15.0*Im_Half1[i+1] - Im_Half1[i+2] );
         R2  = real(1.0/28.0)*( - Re_Half2[i-1] + (real)15.0*Re_Half2[i] + (real)15.0*Re_Half2[i+1] - Re_Half2[i+2] );
         I2  = real(1.0/28.0)*( - Im_Half2[i-1] + (real)15.0*Im_Half2[i] + (real)15.0*Im_Half2[i+1] - Im_Half2[i+2] );
         dR1 = real(1.0/12.0)*( + Re_Half1[i-1] - (real)15.0*Re_Half1[i] + (real)15.0*Re_Half1[i+1] - Re_Half1[i+2] );
         dI1 = real(1.0/12.0)*( + Im_Half1[i-1] - (real)15.0*Im_Half1[i] + (real)15.0*Im_Half1[i+1] - Im_Half1[i+2] );
         dR2 = real(1.0/12.0)*( + Re_Half2[i-1] - (real)15.0*Re_Half2[i] + (real)15.0*Re_Half2[i+1] - Re_Half2[i+2] );
         dI2 = real(1.0/12.0)*( + Im_Half2[i-1] - (real)15.0*Im_Half2[i] + (real)15.0*Im_Half2[i+1] - Im_Half2[i+2] );
#        else
         R1  = real(0.5)*( + Re_Half1[i] + Re_Half1[i+1] );
         I1  = real(0.5)*( + Im_Half1[i] + Im_Half1[i+1] );
         R2  = real(0.5)*( + Re_Half2[i] + Re_Half2[i+1] );
         I2  = real(0.5)*( + Im_Half2[i] + Im_Half2[i+1] );
         dR1 =           ( - Re_Half1[i] + Re_Half1[i+1] );
         dI1 =           ( - Im_Half1[i] + Im_Half1[i+1] );
         dR2 =           ( - Re_Half2[i] + Re_Half2[i+1] );
         dI2 =           ( - Im_Half2[i] + Im_Half2[i+1] );
#        endif

         Flux1[ Idx2 ++ ] = (real)2.0*( R1*dI1 - I1*dR1 );
         Flux2[ Idx2 ++ ] = (real)2.0*( R2*dI2 - I2*dR2 );
      }

//    4.2. correct the amplitude
      Idx2 = 0;
      for (int i=FLU_GHOST_SIZE; i<FLU_NXT-FLU_GHOST_SIZE; i++)
      {
         Amp_Old  = SQR( Re_Old1[i] ) + SQR( Im_Old1[i] ) +  SQR( Re_Old2[i] ) + SQR( Im_Old2[i] );
         Amp_New  = SQR( Re_New1[i] ) + SQR( Im_New1[i] ) +  SQR( Re_New2[i] ) + SQR( Im_New2[i] );
         Amp_Corr = Amp_Old - dT1_dh2*( Flux1[Idx2+1] - Flux1[Idx2] ) -dT2_dh2*( Flux2[Idx2+1] - Flux2[Idx2] );

//       be careful about the negative density and the vacuum (where we might have Amp_New == 0.0)
//       if ( Amp_Corr > (real)0.0  &&  Amp_New > (real)0.0 )
         if ( Amp_Corr >       0.0  &&  Amp_New >       0.0 )
         {
            /*
            Re_New[i] *= SQRT( Amp_Corr / Amp_New );
            Im_New[i] *= SQRT( Amp_Corr / Amp_New );
            */
            Re_New1[i] *= sqrt( Amp_Corr / Amp_New );  // use double precision to improve the mass conservation further
            Im_New1[i] *= sqrt( Amp_Corr / Amp_New );
            Re_New2[i] *= sqrt( Amp_Corr / Amp_New );
            Im_New2[i] *= sqrt( Amp_Corr / Amp_New );
         }


         Idx2 ++;
      }

//    4.3. save the fluxes across all patch boundaries (remeber to put the coefficient "1/(2*Eta*dh)" back)
      if ( StoreFlux )
      if (  ( j>=FLU_GHOST_SIZE && j<FLU_NXT-FLU_GHOST_SIZE )  &&  ( k>=FLU_GHOST_SIZE && k<FLU_NXT-FLU_GHOST_SIZE )  )
      {
         Idx3 = (k-FLU_GHOST_SIZE)*PS2 + (j-FLU_GHOST_SIZE);

         Flux_Array[Flux_XYZ+0][0][Idx3] = Flux1[  0]*_Eta12_dh + Flux2[  0]*_Eta22_dh;
         Flux_Array[Flux_XYZ+1][0][Idx3] = Flux1[PS1]*_Eta12_dh + Flux2[PS1]*_Eta22_dh;
         Flux_Array[Flux_XYZ+2][0][Idx3] = Flux1[PS2]*_Eta12_dh + Flux2[Ps2]*_Eta22_dh;
      }
#     endif // #ifdef CONSERVE_MASS

   } // for j,k

} // FUNCTION : CPU_AdvanceX



//-------------------------------------------------------------------------------------------------------
// Function    :  TrasposeXY
// Description :  Transpose the x and y directions
//
// Parameter   :  u : Input wave function (density, real, imaginary)
//-------------------------------------------------------------------------------------------------------
void TransposeXY( real u[][ FLU_NXT*FLU_NXT*FLU_NXT ] )
{

   real (*u_xy)[FLU_NXT*FLU_NXT] = new real [FLU_NIN][FLU_NXT*FLU_NXT];
   int Idx1, Idx2;

   for (int k=0; k<FLU_NXT; k++)
   {
      for (int j=0; j<FLU_NXT; j++)
      for (int i=0; i<FLU_NXT; i++)
      {
         Idx1 = to1D(k,j,i);
         Idx2 = j + i*FLU_NXT;

         u_xy[0][Idx2] = u[0][Idx1];
         u_xy[1][Idx2] = u[1][Idx1];
         u_xy[2][Idx2] = u[2][Idx1];
         u_xy[3][Idx2] = u[3][Idx1];
      }

      for (int v=0; v<FLU_NIN; v++)    memcpy( &u[v][to1D(k,0,0)], u_xy[v], FLU_NXT*FLU_NXT*sizeof(real) );
   }

   delete [] u_xy;

} // FUNCTION : TrasposeXY



//-------------------------------------------------------------------------------------------------------
// Function    :  TrasposeXZ
// Description :  Transpose the x and z directions
//
// Parameter   :  u : Input wave function (density, real, imaginary)
//-------------------------------------------------------------------------------------------------------
void TransposeXZ( real u[][ FLU_NXT*FLU_NXT*FLU_NXT ] )
{

   real u_temp[FLU_NIN];
   int Idx1, Idx2;

   for (int j=0; j<FLU_NXT; j++)
   for (int k=0; k<FLU_NXT; k++)
   {
      for (int i=0; i<k; i++)
      {
         Idx1 = to1D(k,j,i);
         Idx2 = to1D(i,j,k);

         u_temp[0] = u[0][Idx1];
         u_temp[1] = u[1][Idx1];
         u_temp[2] = u[2][Idx1];
         u_temp[3] = u[3][Idx1];

         u[0][Idx1] = u[0][Idx2];
         u[1][Idx1] = u[1][Idx2];
         u[2][Idx1] = u[2][Idx2];
         u[3][Idx1] = u[3][Idx2];

         u[0][Idx2] = u_temp[0];
         u[1][Idx2] = u_temp[1];
         u[2][Idx2] = u_temp[2];
         u[3][Idx2] = u_temp[3];
      }

      Idx1 = to1D(k,j,k);
   } // j,k

} // FUNCTION : TrasposeXZ



#endif // #if ( !defined GPU  &&  MODEL == ELBDM )
