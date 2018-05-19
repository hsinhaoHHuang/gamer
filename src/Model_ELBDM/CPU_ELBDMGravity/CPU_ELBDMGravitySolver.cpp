#include "GAMER.h"

#if ( !defined GPU  &&  MODEL == ELBDM  &&  defined GRAVITY )




//-----------------------------------------------------------------------------------------
// Function    :  CPU_ELBDMGravitySolver
// Description :  CPU ELBDM gravity solver
//                --> Use CPU to advance wave function by exp( -i*Eta*(Phi+Lambda*Rho)*dt )
//
// Note        :  1. ELBDM gravity solver requires NO potential and fluid ghost zone
//                   --> Optimized performance can be achieved if GRA_GHOST_SIZE == 0, GRA_NXT == PATCH_SIZE
//                   --> But the code supports GRA_GHOST_SIZE > 0 as well (mainly for the STORE_POT_GHOST option)
//                2. ELBDM gravity solver does NOT need the density information (if QUARTIC_SELF_INTERACTION is off)
//                   --> DENS component will NOT be sent in and out in this solver
//                   --> GRA_NIN == 2 (only store the real and imaginary parts)
//                   --> If QUARTIC_SELF_INTERACTION is on, the density is *calculated* here to be REAL^2+IMAG^2
//
// Parameter   :  Flu_Array      : Array to store the input and output data
//                Pot_Array      : Array storing the input potential for evaluating the gravitational acceleration
//                Corner_Array   : Array storing the physical corner coordinates of each patch
//                NPatchGroup    : Number of patch groups to be evaluated
//                EtaDt          : Particle mass / Planck constant * dt
//                dh             : cell size
//                Lambda         : Quartic self-interaction coefficient in ELBDM
//                ExtPot         : Add the external potential
//                Time           : Physical time (may be used by CPU_ExternalPot)
//                ExtPot_AuxArray: Auxiliary array for adding external potential
//-----------------------------------------------------------------------------------------
void CPU_ELBDMGravitySolver(       real Flu_Array[][GRA_NIN][PATCH_SIZE][PATCH_SIZE][PATCH_SIZE],
                             const real Pot_Array[][GRA_NXT][GRA_NXT][GRA_NXT],
                             const double Corner_Array[][3],
                             const int NPatchGroup, const real Eta1Dt, const real Eta2Dt, const real dh, const real Lambda,
                             const bool ExtPot, const double Time, const double ExtPot_AuxArray[] )
{

// check
#  ifdef GAMER_DEBUG
   if ( ExtPot  &&  Corner_Array == NULL )   Aux_Error( ERROR_INFO, "Corner_Array == NULL for ExtPot !!\n" );
#  endif


   const int NPatch = NPatchGroup*8;
   real   Re1, Im1, Re2, Im2, Phase1, Phase2, Cos_Phase1, Sin_Phase1, Cos_Phase2, Sin_Phase2, Pot;
   double x, y, z;
   int    ii, jj, kk;


// loop over all patches
#  pragma omp parallel for private( Re1, Im1, Re2, Im2, Phase1, Phase2, Cos_Phase1, Sin_Phase1, Cos_Phase2, Sin_Phase2, Pot, x, y, z, ii, jj, kk ) schedule( runtime )
   for (int P=0; P<NPatch; P++)
   {
      for (int k=0; k<PATCH_SIZE; k++)    { kk = k + GRA_GHOST_SIZE;
      for (int j=0; j<PATCH_SIZE; j++)    { jj = j + GRA_GHOST_SIZE;
      for (int i=0; i<PATCH_SIZE; i++)    { ii = i + GRA_GHOST_SIZE;

         Re1        = Flu_Array[P][0][ k][ j][ i];
         Im1        = Flu_Array[P][1][ k][ j][ i];
         Re2        = Flu_Array[P][2][ k][ j][ i];
         Im2        = Flu_Array[P][3][ k][ j][ i];
         Pot       = Pot_Array[P]   [kk][jj][ii];

#        ifdef QUARTIC_SELF_INTERACTION
         Pot      += Lambda*( SQR(Re1) + SQR(Im1) +SQR(Re2) + SQR(Im2) );
#        endif

         if ( ExtPot ) {
         x         = Corner_Array[P][0] + (double)(i*dh);
         y         = Corner_Array[P][1] + (double)(j*dh);
         z         = Corner_Array[P][2] + (double)(k*dh);
         Pot      += CPU_ExternalPot( x, y, z, Time, ExtPot_AuxArray ); }

         Phase1    = Eta1Dt * Pot;
         Phase2    = Eta2Dt * Pot;
         Cos_Phase1 = COS( Phase1 );
         Sin_Phase1 = SIN( Phase1 );
         Cos_Phase2 = COS( Phase2 );
         Sin_Phase2 = SIN( Phase2 );

         Flu_Array[P][0][k][j][i] = Cos_Phase1*Re1 + Sin_Phase1*Im1;
         Flu_Array[P][1][k][j][i] = Cos_Phase1*Im1 - Sin_Phase1*Re1;
         Flu_Array[P][2][k][j][i] = Cos_Phase2*Re2 + Sin_Phase2*Im2;
         Flu_Array[P][3][k][j][i] = Cos_Phase2*Im2 - Sin_Phase2*Re2;

      }}} // i,j,k
   } // for (int P=0; P<NPatch; P++)

} // FUNCTION : CPU_ELBDMGravitySolver



#endif // #if ( !defined GPU  &&  MODEL == ELBDM  &&  defined GRAVITY )
