#include "GAMER.h"
#include "TestProb.h"



// problem-specific global variables
// =======================================================================================

static double Height1;
static double Width1;
static double Height2;
static double Width2;
static int    Center_RSeed;
static double EmptyRegion;
static double (*Center)[3] = NULL;
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

//#  ifndef GRAVITY
//   Aux_Error( ERROR_INFO, "GRAVITY must be enabled !!\n" );
//#  endif

#  ifdef COMOVING
   Aux_Error( ERROR_INFO, "COMOVING must be disabled !!\n" );
#  endif

#  ifdef GRAVITY
   if ( OPT__BC_POT != BC_POT_ISOLATED )
      Aux_Error( ERROR_INFO, "must adopt isolated BC for gravity --> reset OPT__BC_POT !!\n" );
#  endif


// warnings
   if ( MPI_Rank == 0 )
   {
      if(!OPT__INIT_RESTRICT)
      Aux_Message( stderr, "WARNING : it's recommended to enable OPT__INIT_RESTRICT !!\n" );

   } // if ( MPI_Rank == 0 )



   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Validating test problem %d ... done\n", TESTPROB_ID );

} // FUNCTION : Validate



// replace HYDRO by the target model (e.g., MHD/ELBDM) and also check other compilation flags if necessary (e.g., GRAVITY/PARTICLE)
#if ( MODEL == ELBDM ) // && defined GRAVITY )
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
// ReadPara->Add( "KEY_IN_THE_FILE",   &VARIABLE,              DEFAULT,       MIN,              MAX               );
// ********************************************************************************************************************************
   ReadPara->Add( "Height1",           &Height1,               1.0,           0.0,              NoMax_double      );
   ReadPara->Add( "Width1",            &Width1,                1.0,           0.0,              NoMax_double      );
   ReadPara->Add( "Height2",           &Height2,               1.0,           0.0,              NoMax_double      );
   ReadPara->Add( "Width2",            &Width2,                1.0,           0.0,              NoMax_double      );
   ReadPara->Add( "Center_RSeed",      &Center_RSeed,          0,             NoMin_int,        NoMax_int         );
   ReadPara->Add( "EmptyRegion",       &EmptyRegion,           0.0,           NoMin_double,     NoMax_double      );

   ReadPara->Read( FileName );

   delete ReadPara;

// (1-2) set the default values

// (1-3) check the runtime parameters
   if ( Center_RSeed >= 0 && EmptyRegion <0.0 )
      Aux_Error(ERROR_INFO, "EmptyRegion(%14.7e)<0.0 !!\n",EmptyRegion);

// (2) set the problem-specific derived parameters
   Center = new double[2][3];


   if(Center_RSeed >= 0)
   {
      const double Coord_Min[3]={ EmptyRegion, EmptyRegion, EmptyRegion };
      const double Coord_Max[3]={amr->BoxSize[0]-EmptyRegion, amr->BoxSize[1]-EmptyRegion, amr->BoxSize[2]-EmptyRegion };
      srand(Center_RSeed);

      for(int t=0;t<2;t++)
      for(int d=0;d<3;d++)
         Center[t][d] = ( (double)rand()/RAND_MAX)*(Coord_Max[d]-Coord_Min[d]) + Coord_Min[d];

   }
   else
   {
      
         Center[0][0]=0.25*amr->BoxSize[0];   
         Center[0][1]=0.25*amr->BoxSize[1];   
         Center[0][2]=0.25*amr->BoxSize[2];

         Center[1][0]=0.75*amr->BoxSize[0];   
         Center[1][1]=0.75*amr->BoxSize[1];   
         Center[1][2]=0.75*amr->BoxSize[2];   
   
 //      Aux_Error(ERROR_INFO,"please hard code the center position !!\n");
               
       
   }

// (3) reset other general-purpose parameters
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


// (4) make a note
   if ( MPI_Rank == 0 )
   {
      Aux_Message( stdout, "=============================================================================\n" );
      Aux_Message( stdout, "  test problem ID           = %d\n",     TESTPROB_ID );
      Aux_Message( stdout, "  Height of the first wavepacket  = %f\n",     Height1 );
      Aux_Message( stdout, "  Width of the first wavepacket  = %f\n",     Width1 );
      Aux_Message( stdout, "  Height of the second wavepacket  = %f\n",     Height2 );
      Aux_Message( stdout, "  Width of the second wavepacket  = %f\n",    Width2 );
      Aux_Message( stdout, "  random seed for setting the center coord. = %d\n",   Center_RSeed   );
      Aux_Message( stdout, "  size of the center-free zone     = %13.7e\n", EmptyRegion   );
      Aux_Message( stdout, "=============================================================================\n" );
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
//                3. Even when DUAL_ENERGY is adopted for HYDRO, one does NOT need to set the dual-energy variable here
//                   --> It will be calculated automatically
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
   
   fluid[REAL1] = Height1* exp(-( SQR(x-Center[0][0]) + SQR(y-Center[0][1]) + SQR(z-Center[0][2]) )/(2.0*Width1*Width1) ) ;
   fluid[IMAG1] = 0.0;
   fluid[REAL2] = Height2* exp(-( SQR(x-Center[1][0]) + SQR(y-Center[1][1]) + SQR(z-Center[1][2]) )/(2.0*Width2*Width2) ) ;
   fluid[IMAG2] = 0.0;
   fluid[DENS]  = SQR( fluid[REAL1] ) + SQR( fluid[IMAG1] ) + SQR( fluid[REAL2] ) + SQR( fluid[IMAG2] );


} // FUNCTION : SetGridIC


void End()
{
   delete []Center;
}

void BCo( real fluid[], const double x, const double y, const double z, const double Time, const int lv, double AuxArray[] )
{
   fluid[REAL1] = (real)0.0;
   fluid[IMAG1] = (real)0.0;
   fluid[REAL2] = (real)0.0;
   fluid[IMAG2] = (real)0.0;
}


#endif // #if ( MODEL == ELBDM && defined GRAVITY )
//-------------------------------------------------------------------------------------------------------
// Function    :  Init_TestProb_ELBDM_TwoMass
// Description :  Test problem initializer
//
// Note        :  None
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Init_TestProb_ELBDM_TwoMass()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );


// validate the compilation flags and runtime parameters
   Validate();


// replace HYDRO by the target model (e.g., MHD/ELBDM) and also check other compilation flags if necessary (e.g., GRAVITY/PARTICLE)
#  if ( MODEL == ELBDM ) //&& defined GRAVITY )
// set the problem-specific runtime parameters
   SetParameter();


// procedure to enable a problem-specific function:
// 1. define a user-specified function (example functions are given below)
// 2. declare its function prototype on the top of this file
// 3. set the corresponding function pointer below to the new problem-specific function
// 4. enable the corresponding runtime option in "Input__Parameter"
//    --> for instance, enable OPT__OUTPUT_USER for Output_User_Ptr
   Init_Function_User_Ptr   = SetGridIC;
   Flag_User_Ptr            = NULL;       // option: OPT__FLAG_USER;        example: Refine/Flag_User.cpp
   Mis_GetTimeStep_User_Ptr = NULL;       // option: OPT__DT_USER;          example: Miscellaneous/Mis_GetTimeStep_User.cpp
   BC_User_Ptr              = BCo;       // option: OPT__BC_FLU_*=4;       example: TestProblem/ELBDM/ExtPot/Init_TestProb_ELBDM_ExtPot.cpp --> BC()
   Flu_ResetByUser_Func_Ptr = NULL;       // option: OPT__RESET_FLUID;      example: Fluid/Flu_ResetByUser.cpp
   Output_User_Ptr          = NULL;       // option: OPT__OUTPUT_USER;      example: TestProblem/Hydro/AcousticWave/Init_TestProb_Hydro_AcousticWave.cpp --> OutputError()
   Aux_Record_User_Ptr      = NULL;       // option: OPT__RECORD_USER;      example: Auxiliary/Aux_Record_User.cpp
   End_User_Ptr             = End;       // option: none;                  example: TestProblem/Hydro/ClusterMerger_vs_Flash/Init_TestProb_ClusterMerger_vs_Flash.cpp --> End_ClusterMerger()
#  ifdef GRAVITY
   Init_ExternalAcc_Ptr     = NULL;       // option: OPT__GRAVITY_TYPE=2/3; example: SelfGravity/Init_ExternalAcc.cpp
   Init_ExternalPot_Ptr     = NULL;       // option: OPT__EXTERNAL_POT;     example: TestProblem/ELBDM/ExtPot/Init_TestProb_ELBDM_ExtPot.cpp --> Init_ExtPot()
#  endif
#  ifdef PARTICLE
   Par_Init_ByFunction_Ptr  = NULL;       // option: PAR_INIT=1;            example: Particle/Par_Init_ByFunction.cpp
#  endif
#  endif // #if ( MODEL == HYDRO )


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : Init_TestProb_Template
