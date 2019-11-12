#include "GAMER.h"

#ifdef PARTICLE

extern int  Plummer_RSeed;
extern real Plummer_MaxR;
extern real Plummer_Rho0;
extern real Plummer_R0;
extern int  Plummer_NBinR;
extern bool Plummer_Collision;
extern real Plummer_Collision_D;
extern real Plummer_Center[3];
extern real Plummer_BulkVel[3];
#if ( MODEL == HYDRO )
extern real Plummer_GasMFrac;
#endif

double MassProf_Plummer( const double r );
static void RanVec_FixRadius( const double r, double RanVec[] );




//-------------------------------------------------------------------------------------------------------
// Function    :  Par_Init_ByFunction
// Description :  Initialize the particle position and velocity
//
// Note        :  1. Invoked by Init_GAMER() using the function pointer "Par_Init_ByFunction_Ptr"
//                   --> This function pointer may be reset by various test problem initializers, in which case
//                       this funtion will become useless
//                2. Periodicity should be taken care of in this function
//                   --> No particles should lie outside the simulation box when the periodic BC is adopted
//                   --> However, if the non-periodic BC is adopted, particles are allowed to lie outside the box
//                       (more specifically, outside the "active" region defined by amr->Par->RemoveCell)
//                       in this function. They will later be removed automatically when calling Par_Aux_InitCheck()
//                       in Init_GAMER().
//                3. Particles set by this function are only temporarily stored in this MPI rank
//                   --> They will later be redistributed when calling Par_FindHomePatch_UniformGrid()
//                       and LB_Init_LoadBalance()
//                   --> Therefore, there is no constraint on which particles should be set by this function
//
// Parameter   :  NPar_ThisRank : Number of particles to be set by this MPI rank
//                NPar_AllRank  : Total Number of particles in all MPI ranks
//                ParMass       : Particle mass     array with the size of NPar_ThisRank
//                ParPosX/Y/Z   : Particle position array with the size of NPar_ThisRank
//                ParVelX/Y/Z   : Particle velocity array with the size of NPar_ThisRank
//                ParTime       : Particle time     array with the size of NPar_ThisRank
//                AllAttribute  : Pointer array for all particle attributes
//                                --> Dimension = [PAR_NATT_TOTAL][NPar_ThisRank]
//                                --> Use the attribute indices defined in Field.h (e.g., Idx_ParCreTime)
//                                    to access the data
//
// Return      :  ParMass, ParPosX/Y/Z, ParVelX/Y/Z, ParTime, AllAttribute
//-------------------------------------------------------------------------------------------------------
void Par_Init_ByFunction_Plummer( const long NPar_ThisRank, const long NPar_AllRank,
                                  real *ParMass, real *ParPosX, real *ParPosY, real *ParPosZ,
                                  real *ParVelX, real *ParVelY, real *ParVelZ, real *ParTime,
                                  real *AllAttribute[PAR_NATT_TOTAL] )
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );


   real *Mass_AllRank   = NULL;
   real *Pos_AllRank[3] = { NULL, NULL, NULL };
   real *Vel_AllRank[3] = { NULL, NULL, NULL };

// only the master rank will construct the initial condition
   if ( MPI_Rank == 0 )
   {
      const double TotM_Inf    = 4.0/3.0*M_PI*CUBE(Plummer_R0)*Plummer_Rho0;
      const double Vmax_Fac    = sqrt( 2.0*NEWTON_G*TotM_Inf );
      const double Coll_Offset = 0.5*Plummer_Collision_D/sqrt(3.0);

      double *Table_MassProf_r = NULL;
      double *Table_MassProf_M = NULL;
      double  TotM, ParM, dr, RanM, RanR, rL, rR, ML, MR, EstM, ErrM, ErrM_Max=-1.0, RanVec[3];
      double  Vmax, RanV, RanProb, Prob;
      int     IdxL, IdxR;

      Mass_AllRank = new real [amr->Par->NPar_Active_AllRank];
      for (int d=0; d<3; d++)
      {
         Pos_AllRank[d] = new real [amr->Par->NPar_Active_AllRank];
         Vel_AllRank[d] = new real [amr->Par->NPar_Active_AllRank];
      }


//    set the random seed
      srand( Plummer_RSeed );


//    determine the total enclosed mass within the maximum radius
      TotM = MassProf_Plummer( Plummer_MaxR );
      ParM = TotM / amr->Par->NPar_Active_AllRank;

      if ( Plummer_Collision )   ParM *= 2.0;

//    rescale particle mass to account for the gas contribution
#     if ( MODEL == HYDRO )
      ParM *= 1.0 - Plummer_GasMFrac;
#     endif


//    construct the mass profile table
      Table_MassProf_r = new double [Plummer_NBinR];
      Table_MassProf_M = new double [Plummer_NBinR];

      dr = Plummer_MaxR / (Plummer_NBinR-1);

      for (int b=0; b<Plummer_NBinR; b++)
      {
         Table_MassProf_r[b] = dr*b;
         Table_MassProf_M[b] = MassProf_Plummer( Table_MassProf_r[b] );
      }


//    set particle attributes
      for (int p=0; p<amr->Par->NPar_Active_AllRank; p++)
      {
//       mass
         Mass_AllRank[p] = ParM;


//       position (sample from the cumulative mass profile and perform linear interpolation)
         RanM = ( (double)rand()/RAND_MAX )*TotM;
         IdxL = Mis_BinarySearch_Real( Table_MassProf_M, 0, Plummer_NBinR-1, RanM );
         IdxR = IdxL + 1;
         rL   = Table_MassProf_r[IdxL];
         rR   = Table_MassProf_r[IdxR];
         ML   = Table_MassProf_M[IdxL];
         MR   = Table_MassProf_M[IdxR];

//       linear interpolation
         RanR = rL + (rR-rL)/(MR-ML)*(RanM-ML);

//       record the maximum error
         EstM     = MassProf_Plummer( RanR );
         ErrM     = fabs( (EstM-RanM)/RanM );
         ErrM_Max = fmax( ErrM, ErrM_Max );

//       randomly set the position vector with a given radius
         RanVec_FixRadius( RanR, RanVec );
         for (int d=0; d<3; d++)    Pos_AllRank[d][p] = RanVec[d] + Plummer_Center[d];

//       set position offset for the Plummer collision test
         if ( Plummer_Collision )
         for (int d=0; d<3; d++)    Pos_AllRank[d][p] += Coll_Offset*( (p<amr->Par->NPar_Active_AllRank/2)?-1.0:+1.0 );

//       check periodicity
         for (int d=0; d<3; d++)
         {
            if ( OPT__BC_FLU[d*2] == BC_FLU_PERIODIC )
               Pos_AllRank[d][p] = FMOD( Pos_AllRank[d][p]+(real)amr->BoxSize[d], (real)amr->BoxSize[d] );
         }


//       velocity
//       determine the maximum velocity (the escaping velocity)
         Vmax = Vmax_Fac*pow( SQR(Plummer_R0) + SQR(RanR), -0.25 );

//       randomly determine the velocity amplitude (ref: Aarseth, S. et al. 1974, A&A, 37, 183: Eq. [A4,A5])
         do
         {
            RanV    = ( (double)rand()/RAND_MAX );          // 0.0 ... 1.0
            RanProb = ( (double)rand()/RAND_MAX )*0.1;      // 0.0 ... 0.1
            Prob    = SQR(RanV)*pow( 1.0-SQR(RanV), 3.5 );  // < 0.1
         }
         while ( RanProb > Prob );

//       randomly set the velocity vector with the given amplitude (RanV*Vmax)
         RanVec_FixRadius( RanV*Vmax, RanVec );
         for (int d=0; d<3; d++)    Vel_AllRank[d][p] = RanVec[d] + Plummer_BulkVel[d];

      } // for (int p=0; p<amr->Par->NPar_Active_AllRank; p++)

      Aux_Message( stdout, "   Total enclosed mass within MaxR  = %13.7e\n",  TotM );
      Aux_Message( stdout, "   Total enclosed mass to inifinity = %13.7e\n",  TotM_Inf );
      Aux_Message( stdout, "   Enclosed ratio                   = %6.2f%%\n", 100.0*TotM/TotM_Inf );
      Aux_Message( stdout, "   Particle mass                    = %13.7e\n",  ParM );
      Aux_Message( stdout, "   Maximum mass error               = %13.7e\n",  ErrM_Max );


//    free memory
      delete [] Table_MassProf_r;
      delete [] Table_MassProf_M;
   } // if ( MPI_Rank == 0 )


// synchronize all particles to the physical time at the base level
   for (long p=0; p<amr->Par->NPar_AcPlusInac; p++)   amr->Par->Time[p] = Time[0];


// get the number of particles in each rank and set the corresponding offsets
   if ( amr->Par->NPar_Active_AllRank > (long)__INT_MAX__ )
      Aux_Error( ERROR_INFO, "NPar_Active_AllRank (%ld) exceeds the maximum integer (%ld) --> MPI will likely fail !!\n",
                 amr->Par->NPar_Active_AllRank, (long)__INT_MAX__ );

   int NSend[MPI_NRank], SendDisp[MPI_NRank], NPar_MyRank=(int)amr->Par->NPar_AcPlusInac;

   MPI_Gather( &NPar_MyRank, 1, MPI_INT, NSend, 1, MPI_INT, 0, MPI_COMM_WORLD );

   if ( MPI_Rank == 0 )
   {
      SendDisp[0] = 0;
      for (int r=1; r<MPI_NRank; r++)  SendDisp[r] = SendDisp[r-1] + NSend[r-1];
   }


// send particle attributes from the master rank to all ranks
   real *Mass   =   amr->Par->Mass;
   real *Pos[3] = { amr->Par->PosX, amr->Par->PosY, amr->Par->PosZ };
   real *Vel[3] = { amr->Par->VelX, amr->Par->VelY, amr->Par->VelZ };

#  ifdef FLOAT8
   MPI_Scatterv( Mass_AllRank, NSend, SendDisp, MPI_DOUBLE, Mass, amr->Par->NPar_AcPlusInac, MPI_DOUBLE, 0, MPI_COMM_WORLD );

   for (int d=0; d<3; d++)
   {
      MPI_Scatterv( Pos_AllRank[d], NSend, SendDisp, MPI_DOUBLE, Pos[d], amr->Par->NPar_AcPlusInac, MPI_DOUBLE, 0, MPI_COMM_WORLD );
      MPI_Scatterv( Vel_AllRank[d], NSend, SendDisp, MPI_DOUBLE, Vel[d], amr->Par->NPar_AcPlusInac, MPI_DOUBLE, 0, MPI_COMM_WORLD );
   }

#  else
   MPI_Scatterv( Mass_AllRank, NSend, SendDisp, MPI_FLOAT,  Mass, amr->Par->NPar_AcPlusInac, MPI_FLOAT,  0, MPI_COMM_WORLD );

   for (int d=0; d<3; d++)
   {
      MPI_Scatterv( Pos_AllRank[d], NSend, SendDisp, MPI_FLOAT,  Pos[d], amr->Par->NPar_AcPlusInac, MPI_FLOAT,  0, MPI_COMM_WORLD );
      MPI_Scatterv( Vel_AllRank[d], NSend, SendDisp, MPI_FLOAT,  Vel[d], amr->Par->NPar_AcPlusInac, MPI_FLOAT,  0, MPI_COMM_WORLD );
   }
#  endif


   if ( MPI_Rank == 0 )
   {
      delete [] Mass_AllRank;

      for (int d=0; d<3; d++)
      {
         delete [] Pos_AllRank[d];
         delete [] Vel_AllRank[d];
      }
   }


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : Par_Init_ByFunction



//-------------------------------------------------------------------------------------------------------
// Function    :  MassProf_Plummer
// Description :  Mass profile of the Plummer model
//
// Note        :  Calculate the enclosed mass for a given radius
//
// Parameter   :  r  : Input radius
//
// Return      :  Enclosed mass
//-------------------------------------------------------------------------------------------------------
double MassProf_Plummer( const double r )
{

   const double x = r / Plummer_R0;

   return 4.0/3.0*M_PI*Plummer_Rho0*CUBE(r)*pow( 1.0+x*x, -1.5 );

} // FUNCTION : MassProf_Plummer



//-------------------------------------------------------------------------------------------------------
// Function    :  RanVec_FixRadius
// Description :  Compute a random 3D vector with a fixed radius
//
// Note        :  Uniformly random sample in theta and phi does NOT give a uniformly random sample in 3D space
//                --> Uniformly random sample in a 3D sphere and then normalize all vectors to the given radius
//
// Parameter   :  r        : Input radius
//                RanVec   : Array to store the random 3D vector
//
// Return      :  RanVec
//-------------------------------------------------------------------------------------------------------
void RanVec_FixRadius( const double r, double RanVec[] )
{

   double Norm, RanR2;

   do
   {
      RanR2 = 0.0;

      for (int d=0; d<3; d++)
      {
         RanVec[d]  = ( (double)rand()/RAND_MAX )*2.0 - 1.0;
         RanR2     += SQR( RanVec[d] );
      }
   } while ( RanR2 > 1.0 );

   Norm = r / sqrt( RanR2 );

   for (int d=0; d<3; d++)    RanVec[d] *= Norm;

} // FUNCTION : RanVec_FixRadius



#endif // #ifdef PARTICLE
