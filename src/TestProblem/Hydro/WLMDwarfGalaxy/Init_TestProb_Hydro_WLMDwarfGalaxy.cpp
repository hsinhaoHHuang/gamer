#include "GAMER.h"



// problem-specific global variables
// =======================================================================================
static bool WLMDwarfGalaxy_UseMetal = false; // add and advect a metal density field
                                             // --> to enable this option, one must
                                             //     set NCOMP_PASSIVE_USER>=1 and PAR_NATT_FLT_USER>=1 in the Makefile
                                             // --> necessary if one wants to enable metal_cooling in Grackle
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


#  if ( MODEL != HYDRO )
   Aux_Error( ERROR_INFO, "MODEL != HYDRO !!\n" );
#  endif

#  ifndef PARTICLE
   Aux_Error( ERROR_INFO, "PARTICLE must be enabled !!\n" );
#  endif

#  ifndef GRAVITY
   Aux_Error( ERROR_INFO, "GRAVITY must be enabled !!\n" );
#  endif

#  ifdef COMOVING
   Aux_Error( ERROR_INFO, "COMOVING must be disabled !!\n" );
#  endif

#  ifndef STAR_FORMATION
   Aux_Error( ERROR_INFO, "STAR_FORMATION must be enabled !!\n" );
#  endif

#  ifndef FEEDBACK
   Aux_Error( ERROR_INFO, "FEEDBACK must be enabled !!\n" );
#  endif

   if ( !OPT__UNIT )
      Aux_Error( ERROR_INFO, "OPT__UNIT must be enabled !!\n" );

   for (int f=0; f<6; f++)
   if ( OPT__BC_FLU[f] == BC_FLU_PERIODIC )
      Aux_Error( ERROR_INFO, "do not use periodic BC (OPT__BC_FLU* = 1) for this test !!\n" );

#  ifdef GRAVITY
   if ( OPT__BC_POT == BC_POT_PERIODIC )
      Aux_Error( ERROR_INFO, "do not use periodic BC (OPT__BC_POT = 1) for this test !!\n" );
#  endif

#  ifdef PARTICLE
   if ( OPT__INIT == INIT_BY_FILE  &&  amr->Par->Init != PAR_INIT_BY_FILE )
      Aux_Error( ERROR_INFO, "please set PAR_INIT = 3 (by FILE) !!\n" );
#  endif

#  ifdef FEEDBACK
   if ( !FB_RESOLVED_SNEII )
      Aux_Error( ERROR_INFO, "FB_RESOLVED_SNEII must be enabled !!\n" );
#  endif


   if ( MPI_Rank == 0 )
   {
#     if ( FLU_SCHEME != MHM )
         Aux_Message( stderr, "WARNING : it's recommended to adopt the MHM scheme for this test !!\n" );
#     endif

#     ifndef DUAL_ENERGY
         Aux_Message( stderr, "WARNING : it's recommended to enable DUAL_ENERGY for this test !!\n" );
#     endif

#     if ( !defined MHD  &&  RSOLVER != HLLC )
         Aux_Message( stderr, "WARNING : it's recommended to adopt the HLLC Riemann solver for this test !!\n" );
#     endif

#     if ( MODEL == HYDRO )
      if ( MINMOD_COEFF > 1.5 )
         Aux_Message( stderr, "WARNING : it's recommended to set MINMOD_COEFF <= 1.5 for this test !!\n" );
#     endif
   } // if ( MPI_Rank == 0 )


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Validating test problem %d ... done\n", TESTPROB_ID );

} // FUNCTION : Validate



#if ( MODEL == HYDRO  &&  defined MASSIVE_PARTICLES )
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
//                3. Must call EoS_Init() before calling any other EoS routine
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

// add parameters in the following format:
// --> note that VARIABLE, DEFAULT, MIN, and MAX must have the same data type
// --> some handy constants (e.g., NoMin_int, Eps_float, ...) are defined in "include/ReadPara.h"
// ********************************************************************************************************************************
// ReadPara->Add( "KEY_IN_THE_FILE",                  &VARIABLE,                 DEFAULT,       MIN,              MAX               );
// ********************************************************************************************************************************
   ReadPara->Add( "WLMDwarfGalaxy_UseMetal",          &WLMDwarfGalaxy_UseMetal,  true,          Useless_bool,     Useless_bool      );

   ReadPara->Read( FileName );

   delete ReadPara;

// check the runtime parameters
#  ifdef SUPPORT_GRACKLE
   if ( GRACKLE_METAL  &&  !WLMDwarfGalaxy_UseMetal )
      Aux_Error( ERROR_INFO, "please enable \"WLMDwarfGalaxy_UseMetal\" for \"GRACKLE_METAL\" !!\n" );
#  endif


// (2) set the problem-specific derived parameters


// (3) reset other general-purpose parameters
//     --> a helper macro PRINT_RESET_PARA is defined in Macro.h
   const long   End_Step_Default = __INT_MAX__;
   const double End_T_Default    =  1000.0*Const_Myr/UNIT_T;

   if ( END_STEP < 0 ) {
      END_STEP = End_Step_Default;
      PRINT_RESET_PARA( END_STEP, FORMAT_LONG, "" );
   }

   if ( END_T < 0.0 ) {
      END_T = End_T_Default;
      PRINT_RESET_PARA( END_T, FORMAT_REAL, "" );
   }


// (4) make a note
   if ( MPI_Rank == 0 )
   {
      Aux_Message( stdout, "=============================================================================\n" );
      Aux_Message( stdout, "  test problem ID   = %d\n",             TESTPROB_ID                             );
      Aux_Message( stdout, "  UseMetal          = %d\n",             WLMDwarfGalaxy_UseMetal                 );
      Aux_Message( stdout, "=============================================================================\n" );
   }


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Setting runtime parameters ... done\n" );

} // FUNCTION : SetParameter



//-------------------------------------------------------------------------------------------------------
// Function    :  SetGridIC
// Description :  Initialize the fluid for the WLM dwarf galaxy test
//
// Note        :  1. NOT needed beacuse the fluid is initialized by file
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

   Aux_Error( ERROR_INFO, "OPT__INIT = 1 is not supported for this test problem !!\n" );

} // FUNCTION : SetGridIC



//-------------------------------------------------------------------------------------------------------
// Function    :  AddNewField_WLMDwarfGalaxy
// Description :  Add the problem-specific grid fields
//
// Note        :  1. Ref: https://github.com/gamer-project/gamer/wiki/Adding-New-Simulations#v-add-problem-specific-grid-fields-and-particle-attributes
//                2. Invoke AddField() for each of the problem-specific field:
//                   --> Field label sent to AddField() will be used as the output name of the field
//                   --> Field index returned by AddField() can be used to access the field data
//                3. Pre-declared field indices are put in Field.h
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void AddNewField_WLMDwarfGalaxy()
{

// add the metallicity field only if it has not been done
// --> since Grackle may already add this field automatically when GRACKLE_METAL is enabled
// --> also note that "Idx_Metal" has been predefined in Field.h
   if ( WLMDwarfGalaxy_UseMetal  &&  Idx_Metal == Idx_Undefined )
      Idx_Metal = AddField( "Metal", FIXUP_FLUX_YES, FIXUP_REST_YES, NORMALIZE_NO, INTERP_FRAC_YES );

} // FUNCTION : AddNewField_WLMDwarfGalaxy



//-------------------------------------------------------------------------------------------------------
// Function    :  AddNewParticleAttribute_WLMDwarfGalaxy
// Description :  Add the problem-specific particle attributes
//
// Note        :  1. Ref: https://github.com/gamer-project/gamer/wiki/Adding-New-Simulations#v-add-problem-specific-grid-fields-and-particle-attributes
//                2. Invoke AddParticleField() for each of the problem-specific particle attribute:
//                   --> Attribute label sent to AddParticleField() will be used as the output name of the attribute
//                   --> Attribute index returned by AddParticleField() can be used to access the particle attribute data
//                3. Pre-declared attribute indices are put in Field.h
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void AddNewParticleAttribute_WLMDwarfGalaxy()
{

// "Idx_ParMetalFrac" has been predefined in Field.h
   if ( WLMDwarfGalaxy_UseMetal  &&  Idx_ParMetalFrac == Idx_Undefined )
      Idx_ParMetalFrac = AddParticleAttributeFlt( "ParMetalFrac" );

} // FUNCTION : AddNewParticleAttribute_WLMDwarfGalaxy
#endif // #if ( MODEL == HYDRO  &&  defined MASSIVE_PARTICLES )



//-------------------------------------------------------------------------------------------------------
// Function    :  Init_TestProb_Hydro_WLMDwarfGalaxy_IsolatedGalaxy
// Description :  Test problem initializer
//
// Note        :  None
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Init_TestProb_Hydro_WLMDwarfGalaxy()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );


// validate the compilation flags and runtime parameters
   Validate();


#  if ( MODEL == HYDRO  &&  defined MASSIVE_PARTICLES )
// set the problem-specific runtime parameters
   SetParameter();


// set the function pointers of various problem-specific routines
   Init_Function_User_Ptr      = SetGridIC;
   Init_Field_User_Ptr         = AddNewField_WLMDwarfGalaxy;
   Par_Init_Attribute_User_Ptr = AddNewParticleAttribute_WLMDwarfGalaxy;
#  endif // if ( MODEL == HYDRO  &&  defined MASSIVE_PARTICLES )


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : Init_TestProb_Hydro_WLMDwarfGalaxy
