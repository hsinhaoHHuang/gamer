#include "Macro.h"
#include "CUPOT.h"

#if ( defined GPU  &&  MODEL == ELBDM  &&  defined GRAVITY )



#include "../../SelfGravity/GPU_Poisson/CUPOT_ExternalPot.cu"

// variables reside in constant memory
__constant__ double ExtPot_AuxArray_d[EXT_POT_NAUX_MAX];




//-------------------------------------------------------------------------------------------------------
// Function    :  CUPOT_ELBDMGravitySolver_SetConstMem
// Description :  Set the constant memory used by CUPOT_ELBDMGravitySolver
//
// Note        :  Adopt the suggested approach for CUDA version >= 5.0
//
// Parameter   :  None
//
// Return      :  0/-1 : successful/failed
//---------------------------------------------------------------------------------------------------
int CUPOT_ELBDMGravitySolver_SetConstMem( double ExtPot_AuxArray_h[] )
{

   if (  cudaSuccess != cudaMemcpyToSymbol( ExtPot_AuxArray_d, ExtPot_AuxArray_h, EXT_POT_NAUX_MAX*sizeof(double),
                                            0, cudaMemcpyHostToDevice)  )
      return -1;

   else
      return 0;

} // FUNCTION : CUPOT_ELBDMGravitySolver_SetConstMem



//-------------------------------------------------------------------------------------------------------
// Function    :  CUPOT_ELBDMGravitySolver
// Description :  GPU ELBDM gravity solver
//                --> Use GPU to advance wave function by exp( -i*Eta*(Phi+Lambda*Rho)*dt )
//
// Note        :  1. ELBDM gravity solver requires NO potential and fluid ghost zone
//                   --> Optimized performance can be achieved if GRA_GHOST_SIZE == 0, GRA_NXT == PATCH_SIZE
//                   --> But the code supports GRA_GHOST_SIZE > 0 as well (mainly for the STORE_POT_GHOST option) 
//                2. ELBDM gravity solver does NOT need the density information (if QUARTIC_SELF_INTERACTION is off)
//                   --> DENS component will NOT be sent in and out in this solver
//                   --> GRA_NIN == 2 (only store the real and imaginary parts)
//                   --> If QUARTIC_SELF_INTERACTION is on, the density is *calculated* here to be REAL^2+IMAG^2
//                3. Prefix "g" for pointers pointing to the "Global" memory space
//                   Prefix "s" for pointers pointing to the "Shared" memory space
//                4. No shared memory is used in this kernel since no computational stencil is required
//                   and hence no data needed to be shared
//
// Parameter   :  g_Flu_Array    : Global memory array to store the input and output data
//                g_Pot_Array    : Global memory array storing the input potential for evaluating the
//                                 gravitational acceleration
//                g_Corner_Array : Global memory array storing the physical corner coordinates of each patch
//                EtaDt          : Particle mass / Planck constant * dt
//                dh             : Cell size
//                Lambda         : Quartic self-interaction coefficient in ELBDM
//                ExtPot         : Add the external potential
//                Time           : Physical time (may be used by CUPOT_ExternalPot)
//---------------------------------------------------------------------------------------------------
__global__ void CUPOT_ELBDMGravitySolver(       real g_Flu_Array[][GRA_NIN][ PS1*PS1*PS1 ],
                                          const real g_Pot_Array[][ GRA_NXT*GRA_NXT*GRA_NXT ],
                                          const double g_Corner_Array[][3],
                                          const real Eta1Dt, const real Eta2Dt, const real dh, const real Lambda, const bool ExtPot,
                                          const double Time )
{

   const uint bx      = blockIdx.x;
   const uint tx      = threadIdx.x; 
   const uint ty      = threadIdx.y; 
   const uint tz      = threadIdx.z; 
         uint Idx_Flu =  __umul24( tz, PS1*PS1 ) + __umul24( ty, PS1 ) + tx;
         uint Idx_Pot =  __umul24( tz+GRA_GHOST_SIZE, GRA_NXT*GRA_NXT )
                       + __umul24( ty+GRA_GHOST_SIZE, GRA_NXT )
                       +           tx+GRA_GHOST_SIZE;

   real   Re1, Im1, Re2, Im2, Phase1, Phase2, Cos_Phase1, Sin_Phase1, Cos_Phase2, Sin_Phase2, Pot;
   double x, y, z;


   if ( ExtPot )
   {
      x = g_Corner_Array[bx][0] + (double)(tx*dh);
      y = g_Corner_Array[bx][1] + (double)(ty*dh);
   }

   for (uint k=tz; k<PS1; k+=GRA_BLOCK_SIZE_Z)
   {
      Re1        = g_Flu_Array[bx][0][Idx_Flu];
      Im1        = g_Flu_Array[bx][1][Idx_Flu];
      Re2        = g_Flu_Array[bx][2][Idx_Flu];
      Im2        = g_Flu_Array[bx][3][Idx_Flu];
      Pot       = g_Pot_Array[bx]   [Idx_Pot];

#     ifdef QUARTIC_SELF_INTERACTION
      Pot      += Lambda*( SQR(Re1) + SQR(Im1) + SQR(Re2) + SQR(Im2) );
#     endif

      if ( ExtPot ) {
      z         = g_Corner_Array[bx][2] + (double)(k*dh);
      Pot      += CUPOT_ExternalPot( x, y, z, Time, ExtPot_AuxArray_d ); }

      Phase1     = Eta1Dt * Pot;
      Phase2     = Eta2Dt * Pot;
      Cos_Phase1 = COS( Phase1 );
      Sin_Phase1 = SIN( Phase1 );
      Cos_Phase2 = COS( Phase2 );
      Sin_Phase2 = SIN( Phase2 );

      g_Flu_Array[bx][0][Idx_Flu] = Cos_Phase1*Re1 + Sin_Phase1*Im1;
      g_Flu_Array[bx][1][Idx_Flu] = Cos_Phase1*Im1 - Sin_Phase1*Re1;
      g_Flu_Array[bx][2][Idx_Flu] = Cos_Phase2*Re2 + Sin_Phase2*Im2;
      g_Flu_Array[bx][3][Idx_Flu] = Cos_Phase2*Im2 - Sin_Phase2*Re2;

      Idx_Flu += GRA_BLOCK_SIZE_Z*PS1*PS1;
      Idx_Pot += GRA_BLOCK_SIZE_Z*GRA_NXT*GRA_NXT;

   } // for (uint k=tz; k<PS1; k+=GRA_BLOCK_SIZE_Z)

} // FUNCTION : CUPOT_ELBDMGravitySolver



#endif // #if ( defined GPU  &&  MODEL == ELBDM  &&  defined GRAVITY )
