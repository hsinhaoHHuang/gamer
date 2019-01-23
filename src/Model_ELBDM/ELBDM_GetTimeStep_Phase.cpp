#include "GAMER.h"

#if ( MODEL == ELBDM )

static real GetMaxPhaseDerivative( const int lv );




//-------------------------------------------------------------------------------------------------------
// Function    :  ELBDM_GetTimeStep_Phase
// Description :  Estimate the evolution time-step by restricting the maximum phase rotation
//
// Note        :  1. This function should be applied to both physical and comoving coordinates and always
//                   return the evolution time-step (dt) actually used in various solvers
//                   --> Physical coordinates : dt = physical time interval
//                       Comoving coordinates : dt = delta(scale_factor) / ( Hubble_parameter*scale_factor^3 )
//                   --> We convert dt back to the physical time interval, which equals "delta(scale_factor)"
//                       in the comoving coordinates, in Mis_GetTimeStep()
//                2. dt = DT__PHASE*2*pi/Max(first temporal derivative of phase)
//
// Parameter   :  lv : Target refinement level
//
// Return      :  dt
//-------------------------------------------------------------------------------------------------------
double ELBDM_GetTimeStep_Phase( const int lv )
{

   real   MaxdS_dt;
   double dt;

// get the maximum first temporal derivative of phase
   MaxdS_dt = GetMaxPhaseDerivative( lv );

// get the time-step
   dt = DT__PHASE*2.0*M_PI/MaxdS_dt;

   return dt;

} // FUNCTION : ELBDM_GetTimeStep_Phase



//-------------------------------------------------------------------------------------------------------
// Function    :  GetMaxPhaseDerivative
// Description :  Evaluate the maximum first temporal derivative of phase for the time-step estimation
//
// Note        :  1. Invoked by ELBDM_GetTimeStep_Phase()
//                2. Currently this function does not take into account the self-interaction potential
//
// Parameter   :  lv : Target refinement level
//
// Return      :  MaxdS_dt
//-------------------------------------------------------------------------------------------------------
real GetMaxPhaseDerivative( const int lv )
{

   const real Eps               = 1.0e-2;          // soften in the denominator for calculating dS_dt
   const bool IntPhase_No       = false;
   const bool DE_Consistency_No = false;
   const real MinDens_No        = -1.0;
   const real MinPres_No        = -1.0;
   const int  NGhost            = 1;               // number of ghost zones to calculate dS_dt
   const int  Size_Flu          = PS2 + 2*NGhost;  // size of the array Flu_Array
   const int  NPG               = 1;               // number of patch groups (must be ONE here)
   const int  AMP1              = DENS1;           // array index to store amplitude1
   const int  AMP2              = DENS2;           // array index to store amplitude2

   real (*Flu_Array)[NCOMP_FLUID][Size_Flu][Size_Flu][Size_Flu] = NULL;

   real dS1_dt, dS2_dt, MaxdS_dt, _dh, _dh2, _dh_sqr, Re1, Im1, Re2, Im2, GradS1[3], GradS2[3];
   real LapAmp1_Amp1, LapAmp2_Amp2, Vel_Sqr1, Vel_Sqr2;                 // laplacian(amplitude)/amplitude, -grad(phase)^2
   real _Dens1, _Dens2, _Amp1, _Amp2;                         // 0.5/Dens/dh, 1.0/amplitude
   int  im, ip, jm, jp, km, kp, I, J, K;

#  ifdef GRAVITY
   const real PotCoeff1 = -2.0*SQR(ELBDM_ETA1);
   const real PotCoeff2 = -2.0*SQR(ELBDM_ETA2);
   const int  Size_Pot = PS2;                // size of the array Pot_Array
   real Pot;                                 // -2.0*ELBDM_ETA^2*potential
   real (*Pot_Array)[Size_Pot][Size_Pot][Size_Pot] = NULL;
#  endif

   MaxdS_dt = 0.0;
   _dh      = (real)1.0/amr->dh[lv];
   _dh2     = (real)0.5*_dh;
   _dh_sqr  = _dh*_dh;


#  ifdef GRAVITY
#  pragma omp parallel private( Flu_Array, dS1_dt, dS2_dt, Re1, Im1, Re2, Im2, GradS1, GradS2, LapAmp1_Amp1, LapAmp2_Amp2, Vel_Sqr1, Vel_Sqr2, _Dens1, _Dens2, _Amp1, _Amp2, Pot_Array, Pot, \
                                im, ip, jm, jp, km, kp, I, J, K )
#  else
#  pragma omp parallel private( Flu_Array, dS1_dt, dS2_dt, Re1, Im1, Re2, Im2, GradS1, GradS2, LapAmp1_Amp1, LapAmp2_Amp2, Vel_Sqr1, Vel_Sqr2, _Dens1, _Dens2, _Amp1, _Amp2, \
                                im, ip, jm, jp, km, kp, I, J, K )
#  endif
   {
      Flu_Array = new real [NPG][NCOMP_FLUID][Size_Flu][Size_Flu][Size_Flu];
#     ifdef GRAVITY
      Pot_Array = new real [NPG]             [Size_Pot][Size_Pot][Size_Pot];
#     endif

//    loop over all patches
#     pragma omp for reduction( max:MaxdS_dt ) schedule( runtime )
      for (int PID0=0; PID0<amr->NPatchComma[lv][1]; PID0+=NPG*8)
      {
//       prepare real part with NGhost ghost zone on each side (any interpolation scheme can be used)
         Prepare_PatchData( lv, Time[lv], &Flu_Array[0][REAL1][0][0][0], NGhost, NPG, &PID0, _REAL1,
                            INT_MINMOD1D, UNIT_PATCHGROUP, NSIDE_06, IntPhase_No, OPT__BC_FLU, BC_POT_NONE,
                            MinDens_No, MinPres_No, DE_Consistency_No );

//       prepare imag part with NGhost ghost zone on each side (any interpolation scheme can be used)
         Prepare_PatchData( lv, Time[lv], &Flu_Array[0][IMAG1][0][0][0], NGhost, NPG, &PID0, _IMAG1,
                            INT_MINMOD1D, UNIT_PATCHGROUP, NSIDE_06, IntPhase_No, OPT__BC_FLU, BC_POT_NONE,
                            MinDens_No, MinPres_No, DE_Consistency_No );

//       prepare real part with NGhost ghost zone on each side (any interpolation scheme can be used)
         Prepare_PatchData( lv, Time[lv], &Flu_Array[0][REAL2][0][0][0], NGhost, NPG, &PID0, _REAL2,
                            INT_MINMOD1D, UNIT_PATCHGROUP, NSIDE_06, IntPhase_No, OPT__BC_FLU, BC_POT_NONE,
                            MinDens_No, MinPres_No, DE_Consistency_No );

//       prepare imag part with NGhost ghost zone on each side (any interpolation scheme can be used)
         Prepare_PatchData( lv, Time[lv], &Flu_Array[0][IMAG2][0][0][0], NGhost, NPG, &PID0, _IMAG2,
                            INT_MINMOD1D, UNIT_PATCHGROUP, NSIDE_06, IntPhase_No, OPT__BC_FLU, BC_POT_NONE,
                            MinDens_No, MinPres_No, DE_Consistency_No );
//       prepare potential with no ghost zone
#        ifdef GRAVITY
         Prepare_PatchData( lv, Time[lv], &Pot_Array[0][0][0][0],            0, NPG, &PID0, _POTE,
                            (IntScheme_t)NULL_INT, UNIT_PATCHGROUP, NSIDE_00, IntPhase_No, OPT__BC_FLU, BC_POT_NONE,
                            MinDens_No, MinPres_No, DE_Consistency_No );
#        endif

//       evaluate amplitude with soften
         for (int k=0; k<Size_Flu; k++)
         for (int j=0; j<Size_Flu; j++)
         for (int i=0; i<Size_Flu; i++)
         {
            Re1 = Flu_Array[0][REAL1][k][j][i];
            Im1 = Flu_Array[0][IMAG1][k][j][i];
            Re2 = Flu_Array[0][REAL2][k][j][i];
            Im2 = Flu_Array[0][IMAG2][k][j][i];

            Flu_Array[0][AMP1][k][j][i] = SQRT( Re1*Re1 + Im1*Im1 + Eps );
            Flu_Array[0][AMP2][k][j][i] = SQRT( Re2*Re2 + Im2*Im2 + Eps );
         }

//       evaluate dS_dt
         for (int k=NGhost; k<Size_Flu-NGhost; k++)    {  km = k - 1;    kp = k + 1;   K = k - NGhost;
         for (int j=NGhost; j<Size_Flu-NGhost; j++)    {  jm = j - 1;    jp = j + 1;   J = j - NGhost;
         for (int i=NGhost; i<Size_Flu-NGhost; i++)    {  im = i - 1;    ip = i + 1;   I = i - NGhost;

            _Amp1       = (real)1.0 / Flu_Array[0][AMP1][k][j][i];
            _Amp2       = (real)1.0 / Flu_Array[0][AMP2][k][j][i];
            _Dens1      = _dh2*_Amp1*_Amp1;
            _Dens2      = _dh2*_Amp2*_Amp2;

            GradS1[0]   = _Dens1*(  Flu_Array[0][REAL1][k][j][i]*( Flu_Array[0][IMAG1][k ][j ][ip] -
                                                                Flu_Array[0][IMAG1][k ][j ][im] )
                                 -Flu_Array[0][IMAG1][k][j][i]*( Flu_Array[0][REAL1][k ][j ][ip] -
                                                                Flu_Array[0][REAL1][k ][j ][im] )  );
            GradS1[1]   = _Dens1*(  Flu_Array[0][REAL1][k][j][i]*( Flu_Array[0][IMAG1][k ][jp][i ] -
                                                                Flu_Array[0][IMAG1][k ][jm][i ] )
                                 -Flu_Array[0][IMAG1][k][j][i]*( Flu_Array[0][REAL1][k ][jp][i ] -
                                                                Flu_Array[0][REAL1][k ][jm][i ] )  );
            GradS1[2]   = _Dens1*(  Flu_Array[0][REAL1][k][j][i]*( Flu_Array[0][IMAG1][kp][j ][i ] -
                                                                Flu_Array[0][IMAG1][km][j ][i ] )
                                 -Flu_Array[0][IMAG1][k][j][i]*( Flu_Array[0][REAL1][kp][j ][i ] -
                                                              Flu_Array[0][REAL1][km][j ][i ] )  );
            GradS2[0]   = _Dens2*(  Flu_Array[0][REAL2][k][j][i]*( Flu_Array[0][IMAG2][k ][j ][ip] -
                                                                Flu_Array[0][IMAG2][k ][j ][im] )
                                 -Flu_Array[0][IMAG2][k][j][i]*( Flu_Array[0][REAL2][k ][j ][ip] -
                                                                Flu_Array[0][REAL2][k ][j ][im] )  );
            GradS2[1]   = _Dens2*(  Flu_Array[0][REAL2][k][j][i]*( Flu_Array[0][IMAG2][k ][jp][i ] -
                                                                Flu_Array[0][IMAG2][k ][jm][i ] )
                                 -Flu_Array[0][IMAG2][k][j][i]*( Flu_Array[0][REAL2][k ][jp][i ] -
                                                                Flu_Array[0][REAL2][k ][jm][i ] )  );
            GradS2[2]   = _Dens2*(  Flu_Array[0][REAL2][k][j][i]*( Flu_Array[0][IMAG2][kp][j ][i ] -
                                                                Flu_Array[0][IMAG2][km][j ][i ] )
                                 -Flu_Array[0][IMAG2][k][j][i]*( Flu_Array[0][REAL2][kp][j ][i ] -
                                                              Flu_Array[0][REAL2][km][j ][i ] )  );

            LapAmp1_Amp1 = _dh_sqr*( Flu_Array[0][AMP1][kp][j ][i ] + Flu_Array[0][AMP1][km][j ][i ] +
                                   Flu_Array[0][AMP1][k ][jp][i ] + Flu_Array[0][AMP1][k ][jm][i ] +
                                   Flu_Array[0][AMP1][k ][j ][ip] + Flu_Array[0][AMP1][k ][j ][im] -
                                   (real)6.0*Flu_Array[0][AMP1][k][j][i] ) * _Amp1;
            LapAmp2_Amp2 = _dh_sqr*( Flu_Array[0][AMP2][kp][j ][i ] + Flu_Array[0][AMP2][km][j ][i ] +
                                   Flu_Array[0][AMP2][k ][jp][i ] + Flu_Array[0][AMP2][k ][jm][i ] +
                                   Flu_Array[0][AMP2][k ][j ][ip] + Flu_Array[0][AMP2][k ][j ][im] -
                                   (real)6.0*Flu_Array[0][AMP2][k][j][i] ) * _Amp2;
            Vel_Sqr1    = -( GradS1[0]*GradS1[0] + GradS1[1]*GradS1[1] + GradS1[2]*GradS1[2] );
            Vel_Sqr2    = -( GradS2[0]*GradS2[0] + GradS2[1]*GradS2[1] + GradS2[2]*GradS2[2] );
            dS1_dt      = LapAmp1_Amp1 + Vel_Sqr1;
            dS2_dt      = LapAmp2_Amp2 + Vel_Sqr2;

#           ifdef GRAVITY
//            Pot        = PotCoeff*Pot_Array[0][K][J][I];
            dS1_dt     += PotCoeff1*Pot_Array[0][K][J][I];
            dS2_dt     += PotCoeff2*Pot_Array[0][K][J][I];
#           endif

            dS1_dt      = FABS( dS1_dt );
            dS2_dt      = FABS( dS2_dt );
            MaxdS_dt   = MAX( MaxdS_dt, MAX( dS1_dt, dS2_dt ) );

         }}} // i,j,k
      } // for (int PID0=0; PID0<amr->NPatchComma[lv][1]; PID0+=NPG*8)

      delete [] Flu_Array;
#     ifdef GRAVITY
      delete [] Pot_Array;
#     endif
   } // OpenMP parallel region


// get the maximum potential in all ranks
   real MaxdS_dt_AllRank;
#  ifdef FLOAT8
   MPI_Allreduce( &MaxdS_dt, &MaxdS_dt_AllRank, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD );
#  else
   MPI_Allreduce( &MaxdS_dt, &MaxdS_dt_AllRank, 1, MPI_FLOAT,  MPI_MAX, MPI_COMM_WORLD );
#  endif


// check
   if ( MaxdS_dt_AllRank == 0.0  &&  MPI_Rank == 0 )
      Aux_Message( stderr, "WARNING : MaxdS_dt == 0.0 at lv %d !!\n", lv );


   return MaxdS_dt_AllRank*0.5/MIN( ELBDM_ETA1, ELBDM_ETA2 );

} // FUNCTION : GetMaxPhaseDerivative



#endif // #if ( MODEL == ELBDM )
