#include "GAMER.h"
#include "TestProb.h"



// problem-specific global variables
// =======================================================================================
static double TotalMassRatio;     // (TotalMass1/TotalMass2)
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

#  ifndef COMOVING
   Aux_Error( ERROR_INFO, "COMOVING must be enabled !!\n" );
#  endif

#  ifdef PARTICLE
   Aux_Error( ERROR_INFO, "PARTICLE must be disabled !!\n" );
#  endif


   for (int f=0; f<6; f++)
   if ( OPT__BC_FLU[f] != BC_FLU_PERIODIC )
      Aux_Error( ERROR_INFO, "must adopt periodic BC for fluid --> reset OPT__BC_FLU* !!\n" );

#  ifdef GRAVITY
   if ( OPT__BC_POT != BC_POT_PERIODIC )
      Aux_Error( ERROR_INFO, "must adopt periodic BC for gravity --> reset OPT__BC_POT !!\n" );
#  endif


// warnings
   if ( MPI_Rank == 0 )
   {
      if ( !OPT__OUTPUT_BASEPS )
         Aux_Message( stderr, "WARNING : it's recommended to enable OPT__OUTPUT_BASEPS !!\n" );
   } // if ( MPI_Rank == 0 )



   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Validating test problem %d ... done\n", TESTPROB_ID );

} // FUNCTION : Validate



// replace HYDRO by the target model (e.g., MHD/ELBDM) and also check other compilation flags if necessary (e.g., GRAVITY/PARTICLE)
#if ( MODEL == ELBDM )
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

   ReadPara->Add( "TotalMassRatio",    &TotalMassRatio,           0.1,        0.0,              NoMax_double      );
   ReadPara->Read( FileName );

   delete ReadPara;

// (1-2) set the default values

// (1-3) check the runtime parameters


// (2) set the problem-specific derived parameters


// (3) reset other general-purpose parameters
//     --> a helper macro PRINT_WARNING is defined in TestProb.h
   const long   End_Step_Default = __INT_MAX__;
   const double End_T_Default    = __FLT_MAX__;

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
      Aux_Message( stdout, "  total mass ratio          = %13.7e\n", TotalMassRatio  );
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
void SetGridIC( real fluid_out[], const real fluid_in[], const int nvar_in,
        const double x, const double y, const double z, const double Time,
        const int lv, double AuxArray[] )
{

   if ( nvar_in != 1 )
      Aux_Error( ERROR_INFO, "nvar_in (%d) != 1 !!\n", nvar_in );

   double Rho0 = 1.0;
   

   fluid_out[DENS1] = (TotalMassRatio/(1.0+TotalMassRatio)) * ( Rho0 + (real)5.0/3.0*( fluid_in[0] - Rho0 ) );
   if ( fluid_out[DENS1] < 0.0 ) fluid_out[DENS1] = 0.0;
   fluid_out[REAL1] = SQRT( fluid_out[DENS1] );
   fluid_out[IMAG1] = (real)0.0;

   fluid_out[DENS2] = (1.0/(1.0+TotalMassRatio)) * ( Rho0 + (real)5.0/3.0*( fluid_in[0] - Rho0 ) );
   if ( fluid_out[DENS2] < 0.0 ) fluid_out[DENS2] = 0.0;
   fluid_out[REAL2] = SQRT( fluid_out[DENS2] );
   fluid_out[IMAG2] = (real)0.0;

} // FUNCTION : SetGridIC
#endif // #if ( MODEL == ELBDM )



//-------------------------------------------------------------------------------------------------------
// Function    :  Init_TestProb_Cosmo
// Description :  Test problem initializer
//
// Note        :  None
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Init_TestProb_ELBDM_Cosmo()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );


// validate the compilation flags and runtime parameters
   Validate();


// replace HYDRO by the target model (e.g., MHD/ELBDM) and also check other compilation flags if necessary (e.g., GRAVITY/PARTICLE)
#  if ( MODEL == ELBDM && defined GRAVITY && defined COMOVING )
// set the problem-specific runtime parameters
   SetParameter();


// procedure to enable a problem-specific function:
// 1. define a user-specified function (example functions are given below)
// 2. declare its function prototype on the top of this file
// 3. set the corresponding function pointer below to the new problem-specific function
// 4. enable the corresponding runtime option in "Input__Parameter"
//    --> for instance, enable OPT__OUTPUT_USER for Output_User_Ptr
   Init_Function_User_Ptr         = NULL;
   Init_ByFile_User_Ptr           = SetGridIC;
   Init_Field_User_Ptr            = NULL; // set NCOMP_PASSIVE_USER;          example: TestProblem/Hydro/Plummer/Init_TestProb_Hydro_Plummer.cpp --> AddNewField()
   Flag_User_Ptr                  = NULL; // option: OPT__FLAG_USER;          example: Refine/Flag_User.cpp
   Mis_GetTimeStep_User_Ptr       = NULL; // option: OPT__DT_USER;            example: Miscellaneous/Mis_GetTimeStep_User.cpp
   BC_User_Ptr                    = NULL; // option: OPT__BC_FLU_*=4;         example: TestProblem/ELBDM/ExtPot/Init_TestProb_ELBDM_ExtPot.cpp --> BC()
   Flu_ResetByUser_Func_Ptr       = NULL; // option: OPT__RESET_FLUID;        example: Fluid/Flu_ResetByUser.cpp
   Output_User_Ptr                = NULL; // option: OPT__OUTPUT_USER;        example: TestProblem/Hydro/AcousticWave/Init_TestProb_Hydro_AcousticWave.cpp --> OutputError()
   Aux_Record_User_Ptr            = NULL; // option: OPT__RECORD_USER;        example: Auxiliary/Aux_Record_User.cpp
   Init_User_Ptr                  = NULL; // option: none;                    example: none
   End_User_Ptr                   = NULL; // option: none;                    example: TestProblem/Hydro/ClusterMerger_vs_Flash/Init_TestProb_ClusterMerger_vs_Flash.cpp --> End_ClusterMerger()
#  ifdef GRAVITY
   Init_ExternalAcc_Ptr           = NULL; // option: OPT__GRAVITY_TYPE=2/3;   example: SelfGravity/Init_ExternalAcc.cpp
   Init_ExternalPot_Ptr           = NULL; // option: OPT__EXTERNAL_POT;       example: TestProblem/ELBDM/ExtPot/Init_TestProb_ELBDM_ExtPot.cpp --> Init_ExtPot()
   Poi_AddExtraMassForGravity_Ptr = NULL; // option: OPT__GRAVITY_EXTRA_MASS; example: none
#  endif
#  endif // #if ( MODEL == ELBDM && defined GRAVITY && defined COMOVING )


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : Init_TestProb_ELBDM_Cosmo
