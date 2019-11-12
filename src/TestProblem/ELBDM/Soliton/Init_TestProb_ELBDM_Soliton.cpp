#include "GAMER.h"
#include "TestProb.h"



// problem-specific global variables
// =======================================================================================
static int      Soliton_N;                               // total number of solitons
static int      Soliton_RSeed;                           // random seed for setting Soliton_Center (<0 --> set manually)
static double   Soliton_CoreRadiusAll;                   // core radius for all solitons (<=0.0 --> set manually)
static double   Soliton_CoreRadiusAll1;                  // core radius for all soliton1s (<=0.0 --> set manually)
static double   Soliton_CoreRadiusAll2;                  // core radius for all soliton2s (<=0.0 --> set manually)
static double   Soliton_EmptyRegion;                     // soliton-free region from the boundary (useful only when Soliton_RSeed>=0)
static double   Soliton_Distance;                        // normalized distance between solitons ( useful only when Soliton_N == 2 )
static char     Soliton_DensProf_Filename[MAX_STRING];   // filename of the soliton density profile

static int      Soliton_DensProf_NBin;                   // number of radial bins of the soliton density profile
static double  *Soliton_DensProf   = NULL;               // soliton density profile [radius/density]
static double  *Soliton_CoreRadius = NULL;               // core radius of each soliton
static double (*Soliton_Center)[3] = NULL;               // center coordinates of each soliton
static double (*Soliton_Velocity)[3] = NULL;             // velocity of each soliton
static double  *Soliton_ScaleL     = NULL;               // L/D: length/density scale factors of each soliton
                                                         //      (defined as the ratio between the core radii/peak
                                                         //      density of the target and reference soliton profiles
static double  *Soliton_ScaleD      = NULL;
static int     *Soliton_ParMass     = NULL;              // particle mass of each soliton
static double   Velocity1X;                              // velocity of the soliton with ELBDM_MASS1 ( code_length / code_time )
static double   Velocity1Y;
static double   Velocity1Z;
static double   Velocity2X;                              // velocity of the soliton with ELBDM_MASS2 ( code_length / code_time )
static double   Velocity2Y;
static double   Velocity2Z;

static double ExtPot_M;          // point source mass
static double ExtPot_Cen[3];     // point source position
// =======================================================================================




//-------------------------------------------------------------------------------------------------------
// Function    :  Validate
// Description :  Validate the compilation flags and runtime parameters for this test problem
//
// Note        :  None
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Validate()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Validating test problem %d ...\n", TESTPROB_ID );


// errors
#  if ( MODEL != ELBDM )
   Aux_Error( ERROR_INFO, "MODEL != ELBDM !!\n" );
#  endif

#  ifndef GRAVITY
   Aux_Error( ERROR_INFO, "GRAVITY must be enabled !!\n" );
#  endif

#  ifdef COMOVING
   Aux_Error( ERROR_INFO, "COMOVING must be disabled !!\n" );
#  endif

#  ifdef PARTICLE
   Aux_Error( ERROR_INFO, "PARTICLE must be disabled !!\n" );
#  endif

#  ifdef GRAVITY
//   if ( OPT__BC_POT != BC_POT_ISOLATED )
//      Aux_Error( ERROR_INFO, "must adopt isolated BC for gravity --> reset OPT__BC_POT !!\n" );
#  endif


// warnings
   if ( MPI_Rank == 0 )
   {
      if ( !OPT__INIT_RESTRICT )
         Aux_Message( stderr, "WARNING : it's recommended to enable OPT__INIT_RESTRICT !!\n" );
   }


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Validating test problem %d ... done\n", TESTPROB_ID );

} // FUNCTION : Validate



#if ( MODEL == ELBDM  &&  defined GRAVITY )
//-------------------------------------------------------------------------------------------------------
// Function    :  SetParameter
// Description :  Load and set the problem-specific runtime parameters
//
// Note        :  1. Filename is set to "Input__TestProb" by default
//                2. Major tasks in this function:
//                   (1) load the problem-specific runtime parameters
//                   (2) set the problem-specific derived parameters
//                   (3) reset other general-purpose parameters if necessary
//                   (4) make a note of the problem-specific parameters
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void SetParameter()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Setting runtime parameters ...\n" );


// (1) load the problem-specific runtime parameters
   const char FileName[] = "Input__TestProb";
   ReadPara_t *ReadPara  = new ReadPara_t;

// (1-1) add parameters in the following format:
// --> note that VARIABLE, DEFAULT, MIN, and MAX must have the same data type
// --> some handy constants (e.g., Useless_bool, Eps_double, NoMin_int, ...) are defined in "include/ReadPara.h"
// ********************************************************************************************************************************
// ReadPara->Add( "KEY_IN_THE_FILE",           &VARIABLE,                   DEFAULT,       MIN,              MAX               );
// ********************************************************************************************************************************
   ReadPara->Add( "Soliton_N",                 &Soliton_N,                 -1,             1,                NoMax_int         );
   ReadPara->Add( "Soliton_RSeed",             &Soliton_RSeed,              0,             NoMin_int,        NoMax_int         );
   ReadPara->Add( "Soliton_CoreRadiusAll",     &Soliton_CoreRadiusAll,      NoDef_double,  NoMin_double,     NoMax_double      );
   ReadPara->Add( "Soliton_CoreRadiusAll1",    &Soliton_CoreRadiusAll1,     NoDef_double,  NoMin_double,     NoMax_double      );
   ReadPara->Add( "Soliton_CoreRadiusAll2",    &Soliton_CoreRadiusAll2,     NoDef_double,  NoMin_double,     NoMax_double      );
   ReadPara->Add( "Soliton_EmptyRegion",       &Soliton_EmptyRegion,        0.0,           NoMin_double,     NoMax_double      );
   ReadPara->Add( "Soliton_Distance",          &Soliton_Distance,           0.5,          0.0,              1.0               );
   ReadPara->Add( "Soliton_DensProf_Filename",  Soliton_DensProf_Filename,  Useless_str,   Useless_str,      Useless_str       );
   ReadPara->Add( "Velocity1X",                &Velocity1X,                 0.0,           NoMin_double,     NoMax_double      );
   ReadPara->Add( "Velocity1Y",                &Velocity1Y,                 0.0,           NoMin_double,     NoMax_double      );
   ReadPara->Add( "Velocity1Z",                &Velocity1Z,                 0.0,           NoMin_double,     NoMax_double      );
   ReadPara->Add( "Velocity2X",                &Velocity2X,                 0.0,           NoMin_double,     NoMax_double      );
   ReadPara->Add( "Velocity2Y",                &Velocity2Y,                 0.0,           NoMin_double,     NoMax_double      );
   ReadPara->Add( "Velocity2Z",                &Velocity2Z,                 0.0,           NoMin_double,     NoMax_double      );
   ReadPara->Add( "ExtPot_M",          &ExtPot_M,               0.0,          0.0,       NoMax_double      );
   ReadPara->Add( "ExtPot_Cen_X",      &ExtPot_Cen[0],         -1.0,          NoMin_double,     NoMax_double      );
   ReadPara->Add( "ExtPot_Cen_Y",      &ExtPot_Cen[1],         -1.0,          NoMin_double,     NoMax_double      );
   ReadPara->Add( "ExtPot_Cen_Z",      &ExtPot_Cen[2],         -1.0,          NoMin_double,     NoMax_double      );

   ReadPara->Read( FileName );

   delete ReadPara;

// set the default center
   for (int d=0; d<3; d++)
      if ( ExtPot_Cen[d] < 0.0 )    ExtPot_Cen[d] = 0.5*amr->BoxSize[d];

// (1-2) set the default values

// (1-3) check the runtime parameters
   if ( Soliton_RSeed >= 0  &&  Soliton_EmptyRegion < 0.0 )
      Aux_Error( ERROR_INFO, "Soliton_EmptyRegion (%14.7e) < 0.0 !!\n", Soliton_EmptyRegion );
   if ( Soliton_CoreRadiusAll == NoDef_double )
      Aux_Error( ERROR_INFO, "Runtime parameter \"Soliton_CoreRadiusAll\" is not set !!\n" );

// (2) set the problem-specific derived parameters
// (2-1) allocate memory
   Soliton_CoreRadius  = new double [Soliton_N];
   Soliton_Center      = new double [Soliton_N][3];
   Soliton_Velocity    = new double [Soliton_N][3];
   Soliton_ScaleL      = new double [Soliton_N];
   Soliton_ScaleD      = new double [Soliton_N];
   Soliton_ParMass     = new int    [Soliton_N];

// (2-2) soliton scale factors
   if ( Soliton_CoreRadiusAll > 0.0 )
   {
      for (int t=0; t<Soliton_N; t++) 
      { 
         Soliton_CoreRadius[t] = Soliton_CoreRadiusAll;
         Soliton_ParMass[t] = t%2;
      }
   }

   else
   {
//    for Soliton_FixedScale <= 0.0, comment out the following line and hard code the scale factor of each soliton
//      if ( Soliton_N != 2 ) Aux_Error( ERROR_INFO, "for Soliton_CoreRadiusAll <= 0.0, please comment out this error check and hard code "
//                             "the core radius of each soliton !!\n" );
//    for (int t=0; t<Soliton_N; t++)  Soliton_CoreRadius[t] = XXX;
      const double  CoreRadius[2] = {Soliton_CoreRadiusAll1, Soliton_CoreRadiusAll2 };//{ 10.0, 166.66666667 };

      if ( Soliton_N == 1 )
      {
         Soliton_CoreRadius[0] = CoreRadius[0];
         Soliton_ParMass[0] = 0;
      }
      else if ( Soliton_N == 2 )
      {
          Soliton_CoreRadius[0] = CoreRadius[0];
          Soliton_CoreRadius[1] = CoreRadius[1];
          Soliton_ParMass[0] = 0;
          Soliton_ParMass[1] = 1;
      }
      else
      {
         if( Soliton_RSeed >= 0 )
         {
            srand( Soliton_RSeed );
            int RandomParMass;
            
            for(int t=0; t<Soliton_N; t++)
            {
               RandomParMass = rand()%2;
               if ( RandomParMass == 0 )
               {
                   Soliton_CoreRadius[t] = CoreRadius[0];
                   Soliton_ParMass[t] = 0;
               }
               if ( RandomParMass == 1 )
               {
                   Soliton_CoreRadius[t] = CoreRadius[1];
                   Soliton_ParMass[t] = 1;
               }
            }
         }
         else
         {
            //Aux_Error( ERROR_INFO," for Soliton_CoreRadius <= 0.0 and Soliton_N != 2 and Soliton_RSeed < 0, please hard code the core radius of each soliton \n ");
            Soliton_CoreRadius[0] = CoreRadius[0];
            Soliton_ParMass[0] = 0;
            Soliton_CoreRadius[1] = CoreRadius[0];
            Soliton_ParMass[1] = 0;
            Soliton_CoreRadius[2] = CoreRadius[0];
            Soliton_ParMass[2] = 0;
            Soliton_CoreRadius[3] = CoreRadius[0];
            Soliton_ParMass[3] = 0;

            Soliton_CoreRadius[4] = CoreRadius[1];
            Soliton_ParMass[4] = 1;
            Soliton_CoreRadius[5] = CoreRadius[1];
            Soliton_ParMass[5] = 1;
            Soliton_CoreRadius[6] = CoreRadius[1];
            Soliton_ParMass[6] = 1;
            Soliton_CoreRadius[7] = CoreRadius[1];
            Soliton_ParMass[7] = 1;
            Soliton_CoreRadius[8] = CoreRadius[1];
            Soliton_ParMass[8] = 1;
            Soliton_CoreRadius[9] = CoreRadius[1];
            Soliton_ParMass[9] = 1;

            for(int t=10;t<Soliton_N; t++){
               Soliton_CoreRadius[t] = CoreRadius[0];
               Soliton_ParMass[t] = 0;
            }


         }
      }
   }

// (2-3) soliton centers
   if ( Soliton_RSeed >= 0 )
   {
      const double Coord_Min[3] = { Soliton_EmptyRegion, Soliton_EmptyRegion, Soliton_EmptyRegion };
      const double Coord_Max[3] = { amr->BoxSize[0] - Soliton_EmptyRegion,
                                    amr->BoxSize[1] - Soliton_EmptyRegion,
                                    amr->BoxSize[2] - Soliton_EmptyRegion };
      srand( Soliton_RSeed );

      for (int t=0; t<Soliton_N; t++)
      for (int d=0; d<3; d++)
         Soliton_Center[t][d] = ( (double)rand()/RAND_MAX )*(Coord_Max[d]-Coord_Min[d]) + Coord_Min[d];
   }

   else
   {
      if ( Soliton_N == 1 )
      {
         for (int d=0; d<3; d++)    Soliton_Center[0][d] = 0.5*amr->BoxSize[d];
      }
      else if ( Soliton_N == 2)
      {
         Soliton_Center[0][0] = (0.5)*amr->BoxSize[0];//(0.5-0.5*Soliton_Distance)*amr->BoxSize[0];
         Soliton_Center[0][1] = (0.5)*amr->BoxSize[1];
         Soliton_Center[0][2] = (0.5)*amr->BoxSize[2];
         Soliton_Center[1][0] = (0.5)*amr->BoxSize[0];//(0.5+0.5*Soliton_Distance)*amr->BoxSize[0];
         Soliton_Center[1][1] = (0.5)*amr->BoxSize[1];
         Soliton_Center[1][2] = (0.5)*amr->BoxSize[2];
      }
      else
      {
//       for Soliton_RSeed<0, comment out the following line and hard code the center of each soliton
         //Aux_Error( ERROR_INFO, "for Soliton_RSeed < 0 and Soliton_N > 2, please comment out this error check and hard code "
         //                       "the center of each soliton !!\n" );
         srand( 1234 );
         Soliton_Center[0][0] = ((0.25+0.5*( (double)rand()/RAND_MAX ) ))*amr->BoxSize[0];//(0.25+0.5*( (double)rand()/RAND_MAX ) )*amr->BoxSize[0];
         Soliton_Center[0][1] = ((0.25+0.5*( (double)rand()/RAND_MAX ) ))*amr->BoxSize[1];
         Soliton_Center[0][2] = ((0.25+0.5*( (double)rand()/RAND_MAX ) ))*amr->BoxSize[2];

         Soliton_Center[1][0] = ((0.25+0.5*( (double)rand()/RAND_MAX ) ))*amr->BoxSize[0];
         Soliton_Center[1][1] = ((0.25+0.5*( (double)rand()/RAND_MAX ) ))*amr->BoxSize[1];
         Soliton_Center[1][2] = ((0.25+0.5*( (double)rand()/RAND_MAX ) ))*amr->BoxSize[2];

         Soliton_Center[2][0] = ((0.25+0.5*( (double)rand()/RAND_MAX ) ))*amr->BoxSize[0];
         Soliton_Center[2][1] = ((0.25+0.5*( (double)rand()/RAND_MAX ) ))*amr->BoxSize[1];
         Soliton_Center[2][2] = ((0.25+0.5*( (double)rand()/RAND_MAX ) ))*amr->BoxSize[2];

         Soliton_Center[3][0] = ((0.25+0.5*( (double)rand()/RAND_MAX ) ))*amr->BoxSize[0];
         Soliton_Center[3][1] = ((0.25+0.5*( (double)rand()/RAND_MAX ) ))*amr->BoxSize[1];
         Soliton_Center[3][2] = ((0.25+0.5*( (double)rand()/RAND_MAX ) ))*amr->BoxSize[2];
//
         Soliton_Center[4][0] = ((0.25+0.5*0.75))*amr->BoxSize[0];
         Soliton_Center[4][1] = ((0.25+0.5*0.75))*amr->BoxSize[1];
         Soliton_Center[4][2] = ((0.25+0.5*0.75))*amr->BoxSize[2];

         Soliton_Center[5][0] = ((0.25+0.5*0.25))*amr->BoxSize[0];
         Soliton_Center[5][1] = ((0.25+0.5*0.25))*amr->BoxSize[1];
         Soliton_Center[5][2] = ((0.25+0.5*0.25))*amr->BoxSize[2];
         Soliton_Center[6][0] = ((0.25+0.5*0.67))*amr->BoxSize[0];
         Soliton_Center[6][1] = ((0.25+0.5*0.33))*amr->BoxSize[1];
         Soliton_Center[6][2] = ((0.25+0.5*0.33))*amr->BoxSize[2];
         Soliton_Center[7][0] = ((0.25+0.5*0.33))*amr->BoxSize[0];
         Soliton_Center[7][1] = ((0.25+0.5*0.67))*amr->BoxSize[1];
         Soliton_Center[7][2] = ((0.25+0.5*0.67))*amr->BoxSize[2];
         Soliton_Center[8][0] = ((0.25+0.5*(0.33+0.34*( (double)rand()/RAND_MAX ) )))*amr->BoxSize[0];
         Soliton_Center[8][1] = ((0.25+0.5*(0.33+0.34*( (double)rand()/RAND_MAX ) )))*amr->BoxSize[1];
         Soliton_Center[8][2] = ((0.25+0.5*(0.33+0.34*( (double)rand()/RAND_MAX ) )))*amr->BoxSize[2];
         Soliton_Center[9][0] = ((0.25+0.5*(0.33+0.34*( (double)rand()/RAND_MAX ) )))*amr->BoxSize[0];
         Soliton_Center[9][1] = ((0.25+0.5*(0.33+0.34*( (double)rand()/RAND_MAX ) )))*amr->BoxSize[1];
         Soliton_Center[9][2] = ((0.25+0.5*(0.33+0.34*( (double)rand()/RAND_MAX ) )))*amr->BoxSize[2];
         
         
         /*Soliton_Center[10][0] = Soliton_Center[0][0]+30.0;
         Soliton_Center[10][1] = Soliton_Center[0][1]+30.0;
         Soliton_Center[10][2] = Soliton_Center[0][2]+30.0;
         Soliton_Center[11][0] = Soliton_Center[0][0]+21.2;
         Soliton_Center[11][1] = Soliton_Center[0][1]-42.4;
         Soliton_Center[11][2] = Soliton_Center[0][2]+21.2;
         Soliton_Center[12][0] = Soliton_Center[0][0]-21.2;
         Soliton_Center[12][1] = Soliton_Center[0][1]+42.4;
         Soliton_Center[12][2] = Soliton_Center[0][2]-21.2;
         Soliton_Center[0][0] = Soliton_Center[0][0]-30.0;
         Soliton_Center[0][1] = Soliton_Center[0][1]-30.0;
         Soliton_Center[0][2] = Soliton_Center[0][2]-30.0;
         Soliton_Center[3][0] = Soliton_Center[7][0];
         Soliton_Center[3][1] = Soliton_Center[7][1];
         Soliton_Center[3][2] = Soliton_Center[7][2];
         Soliton_Center[10][0] = Soliton_Center[8][0];
         Soliton_Center[10][1] = Soliton_Center[8][1];
         Soliton_Center[10][2] = Soliton_Center[8][2];
         Soliton_Center[11][0] = Soliton_Center[9][0];
         Soliton_Center[11][1] = Soliton_Center[9][1];
         Soliton_Center[11][2] = Soliton_Center[9][2];*/
         
         /*for (int t=0; t<Soliton_N; t++)
         for (int d=0; d<3; d++)          Soliton_Center[t][d] = XXX;
         */
      }
   } // if ( Soliton_RSeed >= 0 ) ... else ...

   
// (2-4) soliton velocity
   srand( 1234 );
   for (int t=0; t<Soliton_N; t++)
   {
      if ( Soliton_ParMass[t] == 0)
      {
         Soliton_Velocity[t][0] = Velocity1X;//(2.0*( (double)rand()/RAND_MAX - 0.5 ) )*Velocity1X;
         Soliton_Velocity[t][1] = Velocity1Y;//(2.0*( (double)rand()/RAND_MAX - 0.5 ) )*Velocity1Y;
         Soliton_Velocity[t][2] = Velocity1Z;//(2.0*( (double)rand()/RAND_MAX - 0.5 ) )*Velocity1Z;
      }
      else
      {
         Soliton_Velocity[t][0] = (2.0*( (double)rand()/RAND_MAX - 0.5 ) )*Velocity2X;
         Soliton_Velocity[t][1] = (2.0*( (double)rand()/RAND_MAX - 0.5 ) )*Velocity2Y;
         Soliton_Velocity[t][2] = (2.0*( (double)rand()/RAND_MAX - 0.5 ) )*Velocity2Z;
      }

   }



// (3) load the reference soliton density profile and evaluate the scale factors
   if ( OPT__INIT != INIT_BY_RESTART )
   {
//    load the reference profile
      const bool RowMajor_No  = false;    // load data into the column-major order
      const bool AllocMem_Yes = true;     // allocate memory for Soliton_DensProf
      const int  NCol         = 2;        // total number of columns to load
      const int  Col[NCol]    = {0, 1};   // target columns: (radius, density)

      Soliton_DensProf_NBin = Aux_LoadTable( Soliton_DensProf, Soliton_DensProf_Filename, NCol, Col, RowMajor_No, AllocMem_Yes );


//    get the core radius of the reference profile
      const double *RadiusRef = Soliton_DensProf + 0*Soliton_DensProf_NBin;
      const double *DensRef   = Soliton_DensProf + 1*Soliton_DensProf_NBin;
      const double  DensCore  = 0.5*DensRef[0];   // define core radius as the half-density radius

      double CoreRadiusRef = NULL_REAL;

      for (int b=1; b<Soliton_DensProf_NBin-1; b++)
      {
         if ( DensRef[b] >= DensCore  &&  DensRef[b+1] <= DensCore )
         {
            CoreRadiusRef = 0.5*( RadiusRef[b] + RadiusRef[b+1] );
            break;
         }
      }

      if ( CoreRadiusRef == NULL_REAL )
         Aux_Error( ERROR_INFO, "cannot determine the reference core radius !!\n" );


//    evaluate the scale factors of each soliton
      for (int t=0; t<Soliton_N; t++)
      {
         Soliton_ScaleL[t] = Soliton_CoreRadius[t] / CoreRadiusRef;
         if ( Soliton_ParMass[t] == 0)
         {
            Soliton_ScaleD[t] = 1.0 / ( 4.0*M_PI*NEWTON_G*SQR(ELBDM_ETA1)*POW4(Soliton_ScaleL[t]) );
         }
         else
         {
            Soliton_ScaleD[t] = 1.0 / ( 4.0*M_PI*NEWTON_G*SQR(ELBDM_ETA2)*POW4(Soliton_ScaleL[t]) );
         }
      }
   } // if ( OPT__INIT != INIT_BY_RESTART )


// (4) reset other general-purpose parameters
//     --> a helper macro PRINT_WARNING is defined in TestProb.h
   const long   End_Step_Default = __INT_MAX__;
   const double End_T_Default    = 1.0e3;

   if ( END_STEP < 0 ) {
      END_STEP = End_Step_Default;
      PRINT_WARNING( "END_STEP", END_STEP, FORMAT_LONG );
   }

   if ( END_T < 0.0 ) {
      END_T = End_T_Default;
      PRINT_WARNING( "END_T", END_T, FORMAT_REAL );
   }


// (5) make a note
   if ( MPI_Rank == 0 )
   {
      Aux_Message( stdout, "======================================================================================\n" );
      Aux_Message( stdout, "  test problem ID                           = %d\n",     TESTPROB_ID                );
      Aux_Message( stdout, "  total number of solitons                  = %d\n",     Soliton_N                  );
      Aux_Message( stdout, "  random seed for setting the center coord. = %d\n",     Soliton_RSeed              );
      Aux_Message( stdout, "  size of the soliton-free zone             = %13.7e\n", Soliton_EmptyRegion        );
      Aux_Message( stdout, "  normalized distance                       = %13.7e\n", Soliton_Distance           );
      Aux_Message( stdout, "  density profile filename                  = %s\n",     Soliton_DensProf_Filename  );
      Aux_Message( stdout, "  number of bins of the density profile     = %d\n",     Soliton_DensProf_NBin      );
      Aux_Message( stdout, "\n" );
      Aux_Message( stdout, "  Soliton info:\n" );
      Aux_Message( stdout, "  %7s  %13s  %13s  %13s  %13s  %13s  %13s  %13s  %13s  %13s %13s\n",
                   "ID", "CoreRadius", "ScaleL", "ScaleD", "Center_X", "Center_Y", "Center_Z", "Velocity_X", "Velocity_Y", "Velocity_Z", "ParticleMass" );
      for (int t=0; t<Soliton_N; t++)
      Aux_Message( stdout, "  %7d  %13.6e  %13.6e  %13.6e  %13.6e  %13.6e  %13.6e  %13.6e  %13.6e  %13.6e   ELBDM_MASS%d\n",
                   t, Soliton_CoreRadius[t], Soliton_ScaleL[t], Soliton_ScaleD[t],
                   Soliton_Center[t][0], Soliton_Center[t][1], Soliton_Center[t][2], Soliton_Velocity[t][0], Soliton_Velocity[t][1], Soliton_Velocity[t][2], Soliton_ParMass[t]+1 );
      Aux_Message( stdout, "  point source mass       = %13.7e\n",                   ExtPot_M );
      Aux_Message( stdout, "  point source position   = (%13.7e, %13.7e, %13.7e)\n", ExtPot_Cen[0],
                                                                                     ExtPot_Cen[1],
                                                                                     ExtPot_Cen[2] );
      Aux_Message( stdout, "======================================================================================\n" );
   }


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Setting runtime parameters ... done\n" );

} // FUNCTION : SetParameter



//-------------------------------------------------------------------------------------------------------
// Function    :  SetGridIC
// Description :  Set the problem-specific initial condition on grids
//
// Note        :  1. This function may also be used to estimate the numerical errors when OPT__OUTPUT_USER is enabled
//                   --> In this case, it should provide the analytical solution at the given "Time"
//                2. This function will be invoked by multiple OpenMP threads when OPENMP is enabled
//                   --> Please ensure that everything here is thread-safe
//
// Parameter   :  fluid    : Fluid field to be initialized
//                x/y/z    : Physical coordinates
//                Time     : Physical time
//                lv       : Target refinement level
//                AuxArray : Auxiliary array
//
// Return      :  fluid
//-------------------------------------------------------------------------------------------------------
void SetGridIC( real fluid[], const double x, const double y, const double z, const double Time,
                const int lv, double AuxArray[] )
{

   const double *Table_Radius  = Soliton_DensProf + 0*Soliton_DensProf_NBin;  // radius
   const double *Table_Density = Soliton_DensProf + 1*Soliton_DensProf_NBin;  // density

   double r_tar, r_ref, dens_ref, dens1, dens2, real1, imag1, real2, imag2;


// initialize density as zero since there may be multiple solitons
   fluid[DENS1] = 0.0;
   fluid[DENS2] = 0.0;
   dens1 = 0.0;
   dens2 = 0.0;
   real1 = 0.0;
   imag1 = 0.0;
   real2 = 0.0;
   imag2 = 0.0;

// loop over all solitons to get the total density
   for (int t=0; t<Soliton_N; t++)
   {
      r_tar = sqrt( SQR(x-Soliton_Center[t][0]) + SQR(y-Soliton_Center[t][1]) + SQR(z-Soliton_Center[t][2]) );

//    rescale radius (target radius --> reference radius)
      r_ref = r_tar / Soliton_ScaleL[t];

//    linear interpolation
      dens_ref = Mis_InterpolateFromTable( Soliton_DensProf_NBin, Table_Radius, Table_Density, r_ref );

      if ( dens_ref == NULL_REAL )
      {
         if      ( r_ref <  Table_Radius[0] )                        
            dens_ref = Table_Density[0];
         else if ( r_ref >= Table_Radius[Soliton_DensProf_NBin-1] )  
            dens_ref = Table_Density[Soliton_DensProf_NBin-1];
         else
            Aux_Error( ERROR_INFO, "interpolation failed at radius %13.7e (min/max radius = %13.7e/%13.7e) !!\n",
                       r_ref, Table_Radius[0], Table_Radius[Soliton_DensProf_NBin-1] );
      }

//    rescale density (reference density --> target density) and add to the fluid array
      if ( Soliton_ParMass[t] == 0 )
      {
         dens1 += dens_ref*Soliton_ScaleD[t];
         real1 += sqrt(dens_ref*Soliton_ScaleD[t])*cos( ELBDM_ETA1*( Soliton_Velocity[t][0]*x + Soliton_Velocity[t][1]*y + Soliton_Velocity[t][2]*z ) );
         imag1 += sqrt(dens_ref*Soliton_ScaleD[t])*sin( ELBDM_ETA1*( Soliton_Velocity[t][0]*x + Soliton_Velocity[t][1]*y + Soliton_Velocity[t][2]*z ) );
      }
      else
      {
         dens2 += dens_ref*Soliton_ScaleD[t]; 
         real2 += sqrt(dens_ref*Soliton_ScaleD[t])*cos( ELBDM_ETA2*( Soliton_Velocity[t][0]*x + Soliton_Velocity[t][1]*y + Soliton_Velocity[t][2]*z ) );
         imag2 += sqrt(dens_ref*Soliton_ScaleD[t])*sin( ELBDM_ETA2*( Soliton_Velocity[t][0]*x + Soliton_Velocity[t][1]*y + Soliton_Velocity[t][2]*z ) );
      }
   } // for (int t=0; t<Soliton_N; t++)


// set the real and imaginary parts
   if ( dens1 == 0.0 )
   {
     fluid[REAL1] = (real)0.0;
     fluid[IMAG1] = (real)0.0;
     fluid[DENS1] = (real)0.0;
   }
   else
   {
     fluid[REAL1] = real1*sqrt( dens1/(real1*real1+imag1*imag1) );//sqrt( dens1 )*cos( ELBDM_ETA1*( Velocity1X*x + Velocity1Y*y + Velocity1Z*z ) );
     fluid[IMAG1] = imag1*sqrt( dens1/(real1*real1+imag1*imag1) );//sqrt( dens1 )*sin( ELBDM_ETA1*( Velocity1X*x + Velocity1Y*y + Velocity1Z*z ) );
     fluid[DENS1] = SQR( fluid[REAL1] ) + SQR( fluid[IMAG1] );
   }

   if ( dens2 == 0.0 )
   {
     fluid[REAL2] = (real)0.0;
     fluid[IMAG2] = (real)0.0;
     fluid[DENS2] = (real)0.0;
   }
   else
   {
     fluid[REAL2] = real2*sqrt( dens2/(real2*real2+imag2*imag2) );//sqrt( dens2 )*cos( ELBDM_ETA2*( Velocity2X*x + Velocity2Y*y + Velocity2Z*z ) );
     fluid[IMAG2] = imag2*sqrt( dens2/(real2*real2+imag2*imag2) );//sqrt( dens2 )*sin( ELBDM_ETA2*( Velocity2X*x + Velocity2Y*y + Velocity2Z*z ) );
     fluid[DENS2] = SQR( fluid[REAL2] ) + SQR( fluid[IMAG2] );
   }


} // FUNCTION : SetGridIC



//-------------------------------------------------------------------------------------------------------
// Function    :  End_Soliton
// Description :  Free memory before terminating the program
//
// Note        :  1. Linked to the function pointer "End_User_Ptr" to replace "End_User()"
//
// Parameter   :  None
//-------------------------------------------------------------------------------------------------------
void End_Soliton()
{

   delete [] Soliton_DensProf;
   delete [] Soliton_Center;
   delete [] Soliton_Velocity;
   delete [] Soliton_ScaleL;
   delete [] Soliton_ScaleD;
   delete [] Soliton_ParMass;

} // FUNCTION : End_Soliton



//-------------------------------------------------------------------------------------------------------
// Function    :  BC
// Description :  Set the extenral boundary condition
//
// Note        :  1. Linked to the function pointer "BC_User_Ptr"
//
// Parameter   :  fluid    : Fluid field to be set
//                x/y/z    : Physical coordinates
//                Time     : Physical time
//                lv       : Refinement level
//                AuxArray : Auxiliary array
//
// Return      :  fluid
//-------------------------------------------------------------------------------------------------------
void BC( real fluid[], const double x, const double y, const double z, const double Time,
         const int lv, double AuxArray[] )
{

   fluid[REAL1] = (real)0.0;
   fluid[IMAG1] = (real)0.0;
   fluid[REAL2] = (real)0.0;
   fluid[IMAG2] = (real)0.0;
   fluid[DENS1] = (real)0.0;
   fluid[DENS2] = (real)0.0;

} // FUNCTION : BC

//-------------------------------------------------------------------------------------------------------
// Function    :  Init_ExtPot
// Description :  Set the array "ExtPot_AuxArray" used by the external potential routines
//                "CUPOT_ExternalPot.cu / CPU_ExternalPot.cpp"
//
// Note        :  1. Linked to the function pointer "Init_ExternalPot_Ptr"
//                2. Enabled by the runtime option "OPT__EXTERNAL_POT"
//
// Parameter   :  None
//-------------------------------------------------------------------------------------------------------
void Init_ExtPot()
{

// ExtPot_AuxArray has the size of EXT_POT_NAUX_MAX (default = 10)
   ExtPot_AuxArray[0] = ExtPot_Cen[0];
   ExtPot_AuxArray[1] = ExtPot_Cen[1];
   ExtPot_AuxArray[2] = ExtPot_Cen[2];
   ExtPot_AuxArray[3] = ExtPot_M*NEWTON_G;

} // FUNCTION : Init_ExtPot
#endif // #if ( MODEL == ELBDM  &&  defined GRAVITY )



//-------------------------------------------------------------------------------------------------------
// Function    :  Init_TestProb_ELBDM_Soliton
// Description :  Test problem initializer
//
// Note        :  None
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Init_TestProb_ELBDM_Soliton()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );


// validate the compilation flags and runtime parameters
   Validate();


#  if ( MODEL == ELBDM  &&  defined GRAVITY )
// set the problem-specific runtime parameters
   SetParameter();


   Init_Function_User_Ptr   = SetGridIC;
   Flag_User_Ptr            = NULL;
   Mis_GetTimeStep_User_Ptr = NULL;
   BC_User_Ptr              = BC;
   Flu_ResetByUser_Func_Ptr = NULL;
   Output_User_Ptr          = NULL;
   Aux_Record_User_Ptr      = NULL;
   End_User_Ptr             = End_Soliton;
   Init_ExternalAcc_Ptr     = NULL;
   Init_ExternalPot_Ptr     = Init_ExtPot;
#  endif // #if ( MODEL == ELBDM  &&  defined GRAVITY )


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : Init_TestProb_ELBDM_Soliton
