#include "GAMER.h"



// problem-specific global variables
// =======================================================================================
static int      HaloMerger_Halo_Num;                                       // total number of halos [2]
static int      HaloMerger_Halo_InitMode;                                  // halo initialization mode (1=single-level UM_IC of real and imaginary parts) [1]
static int      HaloMerger_Soliton_Num;                                    // total number of solitons [0]
static int      HaloMerger_Soliton_InitMode;                               // soliton initialization mode (1=table of the density profile, 2=analytical function of the density profile) [1]
       int      HaloMerger_ParCloud_Num;                                   // total number of particle clouds [0]
       int      HaloMerger_ParCloud_InitMode;                              // particle cloud initialization mode (1=table of the density profile) [1]

// External potential-related parameters to read from the input
       double   HaloMerger_ExtPot_UniDenSph_M;                             // mass of the uniform-density sphere for the external potential (must >= 0.0) [0.0]
       double   HaloMerger_ExtPot_UniDenSph_R;                             // radius of the uniform-density sphere for the external potential (must > 0.0) [1.0]
       double   HaloMerger_ExtPot_UniDenSph_CenCoordX;                     // x/y/z-coordinate of the center of the uniform-density sphere for the external potential (<0.0=auto -> box center) [-1.0]
       double   HaloMerger_ExtPot_UniDenSph_CenCoordY;
       double   HaloMerger_ExtPot_UniDenSph_CenCoordZ;
       double   HaloMerger_ExtPot_UniDenSph_VelocityX;                     // x/y/z-component of the velocity of the uniform-density sphere for the external potential [0.0]
       double   HaloMerger_ExtPot_UniDenSph_VelocityY;
       double   HaloMerger_ExtPot_UniDenSph_VelocityZ;

// Halo-related internal variables
static double (*HaloMerger_Halo_CenCoord)[3]                       = NULL; // center coordinates of each halo
static double (*HaloMerger_Halo_Velocity)[3]                       = NULL; // bulk velocity of each halo
static char   (*HaloMerger_Halo_UM_IC_Filename)[MAX_STRING]        = NULL; // UM_IC filename of each halo
static double (*HaloMerger_Halo_UM_IC_BoxLen)[3]                   = NULL; // physical length of UM_IC box of each halo
static int    (*HaloMerger_Halo_UM_IC_NCells)[3]                   = NULL; // number of cells of UM_IC of each halo
static int     *HaloMerger_Halo_UM_IC_Float8                       = NULL; // data precision of UM_IC of each halo
static double (*HaloMerger_Halo_UM_IC_dh)[3]                       = NULL; // grid size of each halo
static double (*HaloMerger_Halo_UM_IC_Range_EdgeL)[3]              = NULL; // left edge of the range of each halo
static double (*HaloMerger_Halo_UM_IC_Range_EdgeR)[3]              = NULL; // right edge of the range of each halo
static char   **HaloMerger_Halo_UM_IC_Data                         = NULL; // array to store the data read from UM_IC

// Soliton-related internal variables
static double  *HaloMerger_Soliton_CoreRadius                      = NULL; // core radius of each soliton
static double  *HaloMerger_Soliton_CoreRho                         = NULL; // peak density of each soliton
static double (*HaloMerger_Soliton_CenCoord)[3]                    = NULL; // center coordinates of each soliton
static double (*HaloMerger_Soliton_Velocity)[3]                    = NULL; // bulk velocity of each soliton
static double  *HaloMerger_Soliton_OuterSlope                      = NULL; // outer slope of the analytical density profile of each soliton
static char   (*HaloMerger_Soliton_DensProf_Filename)[MAX_STRING]  = NULL; // filename of the density profile table of each soliton
static int     *HaloMerger_Soliton_DensProf_NBin                   = NULL; // number of bins of the density profile table
static int     *HaloMerger_Soliton_DensProf_Rescale                = NULL; // whether to scale the density profile table of each soliton
static double  *HaloMerger_Soliton_DensProf_ScaleL                 = NULL; // L/D: length/density scale factors of each soliton
static double  *HaloMerger_Soliton_DensProf_ScaleD                 = NULL; //      (defined as the ratio between the core radii/peak
                                                                           //      density of the target and reference soliton profiles)
static double **HaloMerger_Soliton_DensProf                        = NULL; // array to store the density profile read from table

// ParCloud-related internal variables
       double (*HaloMerger_ParCloud_CenCoord)[3]                   = NULL; // center coordinates of each particle cloud
       double (*HaloMerger_ParCloud_Velocity)[3]                   = NULL; // bulk velocity of each particle cloud
       char   (*HaloMerger_ParCloud_DensProf_Filename)[MAX_STRING] = NULL; // filename of the density profile table of each particle cloud
       double  *HaloMerger_ParCloud_DensProf_MaxR                  = NULL; // maximum radius for particles of each particle cloud
       int     *HaloMerger_ParCloud_RSeed                          = NULL; // random seed for particles of each particle cloud
       long    *HaloMerger_ParCloud_NPar                           = NULL; // number of particles of each particle cloud
// =======================================================================================

#if ( MODEL == ELBDM  &&  defined GRAVITY )
// external potential routine
// =======================================================================================
void Init_ExtPot_ELBDM_HaloMerger();

// external particle initialization routine
// =======================================================================================
#ifdef MASSIVE_PARTICLES
void Par_Init_ByFunction_HaloMerger( const long NPar_ThisRank, const long NPar_AllRank,
                                     real *ParMass, real *ParPosX, real *ParPosY, real *ParPosZ,
                                     real *ParVelX, real *ParVelY, real *ParVelZ, real *ParTime,
                                     real *ParType, real *AllAttribute[PAR_NATT_TOTAL] );
#endif

// problem-specific functions
// =======================================================================================
static void HaloMerger_Add_Velocity( double *RealPart, double *ImagPart,
                                     const double Velocity_X, const double Velocity_Y, const double Velocity_Z,
                                     const double Position_X, const double Position_Y, const double Position_Z );

static double HaloMerger_Trilinear_Interpolation( const double Target_X, const double Target_Y, const double Target_Z,
                                                  const double Ref_Value[2][2][2],
                                                  const double Ref_X[2], const double Ref_Y[2], const double Ref_Z[2] );

static double HaloMerger_Get_Value_From_Halo_UM_IC_Data( const double x, const double y, const double z, const int v, const int index_halo );
static double LinearDensityShellMass      ( const double r0, const double r1, const double rho0, const double rho1 );
static double LogLinearDensityShellMass   ( const double r0, const double r1, const double rho0, const double rho1 );
// =======================================================================================
#endif // #if ( MODEL == ELBDM  &&  defined GRAVITY )




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
   if ( OPT__INIT == INIT_BY_FUNCTION  &&  amr->Par->Init != PAR_INIT_BY_FUNCTION )
      Aux_Error( ERROR_INFO, "please set PAR_INIT = 1 (by FUNCTION) !!\n" );

#  ifndef SUPPORT_GSL
   Aux_Error( ERROR_INFO, "SUPPORT_GSL must be enabled !!\n" );
#  endif
#  endif

   if ( !OPT__UNIT )
      Aux_Error( ERROR_INFO, "OPT__UNIT must be enabled !!\n" );

   for (int f=0; f<6; f++)
   if ( OPT__BC_FLU[f] != BC_FLU_PERIODIC )
      Aux_Error( ERROR_INFO, "must adopt periodic BC for fluid --> reset OPT__BC_FLU* !!\n" );

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
// --> some handy constants (e.g., NoMin_int, Eps_float, ...) are defined in "include/ReadPara.h"
// ********************************************************************************************************************************
// ReadPara->Add( "KEY_IN_THE_FILE",                       &VARIABLE,                                   DEFAULT,          MIN,           MAX              );
// ********************************************************************************************************************************
   ReadPara->Add( "HaloMerger_Halo_Num",                   &HaloMerger_Halo_Num,                        2,                0,             NoMax_int        );
   ReadPara->Add( "HaloMerger_Halo_InitMode",              &HaloMerger_Halo_InitMode,                   1,                1,             1                );
   ReadPara->Add( "HaloMerger_Soliton_Num",                &HaloMerger_Soliton_Num,                     0,                0,             NoMax_int        );
   ReadPara->Add( "HaloMerger_Soliton_InitMode",           &HaloMerger_Soliton_InitMode,                1,                1,             2                );
   ReadPara->Add( "HaloMerger_ParCloud_Num",               &HaloMerger_ParCloud_Num,                    0,                0,             NoMax_int        );
   ReadPara->Add( "HaloMerger_ParCloud_InitMode",          &HaloMerger_ParCloud_InitMode,               1,                1,             1                );
   ReadPara->Add( "HaloMerger_ExtPot_UniDenSph_M",         &HaloMerger_ExtPot_UniDenSph_M,              0.0,              0.0,           NoMax_double     );
   ReadPara->Add( "HaloMerger_ExtPot_UniDenSph_R",         &HaloMerger_ExtPot_UniDenSph_R,              1.0,              Eps_double,    NoMax_double     );
   ReadPara->Add( "HaloMerger_ExtPot_UniDenSph_CenCoordX", &HaloMerger_ExtPot_UniDenSph_CenCoordX,     -1.0,              NoMin_double,  amr->BoxEdgeR[0] );
   ReadPara->Add( "HaloMerger_ExtPot_UniDenSph_CenCoordY", &HaloMerger_ExtPot_UniDenSph_CenCoordY,     -1.0,              NoMin_double,  amr->BoxEdgeR[1] );
   ReadPara->Add( "HaloMerger_ExtPot_UniDenSph_CenCoordZ", &HaloMerger_ExtPot_UniDenSph_CenCoordZ,     -1.0,              NoMin_double,  amr->BoxEdgeR[2] );
   ReadPara->Add( "HaloMerger_ExtPot_UniDenSph_VelocityX", &HaloMerger_ExtPot_UniDenSph_VelocityX,      0.0,              NoMin_double,  NoMax_double     );
   ReadPara->Add( "HaloMerger_ExtPot_UniDenSph_VelocityY", &HaloMerger_ExtPot_UniDenSph_VelocityY,      0.0,              NoMin_double,  NoMax_double     );
   ReadPara->Add( "HaloMerger_ExtPot_UniDenSph_VelocityZ", &HaloMerger_ExtPot_UniDenSph_VelocityZ,      0.0,              NoMin_double,  NoMax_double     );

   ReadPara->Read( FileName );

   delete ReadPara;

// (1-2) load the runtime parameters for the halos
   if ( OPT__INIT != INIT_BY_RESTART  &&  HaloMerger_Halo_Num > 0 )
   {
      // (1-2-1) allocate the memory
      HaloMerger_Halo_CenCoord              = new double [HaloMerger_Halo_Num][3];
      HaloMerger_Halo_Velocity              = new double [HaloMerger_Halo_Num][3];

      if ( HaloMerger_Halo_InitMode == 1 )
      {
      HaloMerger_Halo_UM_IC_Filename        = new char   [HaloMerger_Halo_Num][MAX_STRING];
      HaloMerger_Halo_UM_IC_BoxLen          = new double [HaloMerger_Halo_Num][3];
      HaloMerger_Halo_UM_IC_NCells          = new int    [HaloMerger_Halo_Num][3];
      HaloMerger_Halo_UM_IC_Float8          = new int    [HaloMerger_Halo_Num];
      HaloMerger_Halo_UM_IC_dh              = new double [HaloMerger_Halo_Num][3];
      HaloMerger_Halo_UM_IC_Range_EdgeL     = new double [HaloMerger_Halo_Num][3];
      HaloMerger_Halo_UM_IC_Range_EdgeR     = new double [HaloMerger_Halo_Num][3];
      HaloMerger_Halo_UM_IC_Data            = new char*  [HaloMerger_Halo_Num];
      } // if ( HaloMerger_Halo_InitMode == 1 )
      else
         Aux_Error( ERROR_INFO, "unsupported initialization mode (%s = %d) !!\n",
                    "HaloMerger_Halo_InitMode", HaloMerger_Halo_InitMode );

      // (1-2-2) read the parameters for the halos
      const char FileName_Halo[] = "Input__TestProb_Halo";
      ReadPara_t *ReadPara_Halo  = new ReadPara_t;

      // Halo-related parameters to read from the input
      char HaloMerger_Halo_i_CenCoordX[MAX_STRING];      // x/y/z-coordinate of the center of the i-th halo (<0.0=auto -> box center) [-1.0]
      char HaloMerger_Halo_i_CenCoordY[MAX_STRING];      // (Note that CenCoordX/Y/Z denotes the UM_IC box center, not the exact halo center, when HaloMerger_Halo_InitMode == 1)
      char HaloMerger_Halo_i_CenCoordZ[MAX_STRING];
      char HaloMerger_Halo_i_VelocityX[MAX_STRING];      // x/y/z-component of the bulk velocity of the i-th halo [0.0]
      char HaloMerger_Halo_i_VelocityY[MAX_STRING];
      char HaloMerger_Halo_i_VelocityZ[MAX_STRING];

      char HaloMerger_Halo_i_UM_IC_Filename[MAX_STRING]; // filename of UM_IC (binary file in vzyx format; row-major and v=field) (single AMR level) for the i-th halo (HaloMerger_Halo_InitMode == 1 only)
      char HaloMerger_Halo_i_UM_IC_BoxLenX [MAX_STRING]; // physical length in the x/y/z-direction of UM_IC box for the i-th halo (must > 0.0) (HaloMerger_Halo_InitMode == 1 only)
      char HaloMerger_Halo_i_UM_IC_BoxLenY [MAX_STRING];
      char HaloMerger_Halo_i_UM_IC_BoxLenZ [MAX_STRING];
      char HaloMerger_Halo_i_UM_IC_NCellsX [MAX_STRING]; // number of cells in the x/y/z-direction of UM_IC box for the i-th halo (must > 0) (HaloMerger_Halo_InitMode == 1 only)
      char HaloMerger_Halo_i_UM_IC_NCellsY [MAX_STRING];
      char HaloMerger_Halo_i_UM_IC_NCellsZ [MAX_STRING];
      char HaloMerger_Halo_i_UM_IC_Float8  [MAX_STRING]; // data precision of UM_IC for the i-th halo (0=float, 1=double) (HaloMerger_Halo_InitMode == 1 only) [0]

      for (int index_halo=0; index_halo<HaloMerger_Halo_Num; index_halo++)
      {
         const int index_halo_input = index_halo+1; // index of halo in the input file

         sprintf( HaloMerger_Halo_i_CenCoordX,      "HaloMerger_Halo_%d_CenCoordX",      index_halo_input );
         sprintf( HaloMerger_Halo_i_CenCoordY,      "HaloMerger_Halo_%d_CenCoordY",      index_halo_input );
         sprintf( HaloMerger_Halo_i_CenCoordZ,      "HaloMerger_Halo_%d_CenCoordZ",      index_halo_input );
         sprintf( HaloMerger_Halo_i_VelocityX,      "HaloMerger_Halo_%d_VelocityX",      index_halo_input );
         sprintf( HaloMerger_Halo_i_VelocityY,      "HaloMerger_Halo_%d_VelocityY",      index_halo_input );
         sprintf( HaloMerger_Halo_i_VelocityZ,      "HaloMerger_Halo_%d_VelocityZ",      index_halo_input );

         if ( HaloMerger_Halo_InitMode == 1 )
         {
         sprintf( HaloMerger_Halo_i_UM_IC_Filename, "HaloMerger_Halo_%d_UM_IC_Filename", index_halo_input );
         sprintf( HaloMerger_Halo_i_UM_IC_BoxLenX,  "HaloMerger_Halo_%d_UM_IC_BoxLenX",  index_halo_input );
         sprintf( HaloMerger_Halo_i_UM_IC_BoxLenY,  "HaloMerger_Halo_%d_UM_IC_BoxLenY",  index_halo_input );
         sprintf( HaloMerger_Halo_i_UM_IC_BoxLenZ,  "HaloMerger_Halo_%d_UM_IC_BoxLenZ",  index_halo_input );
         sprintf( HaloMerger_Halo_i_UM_IC_NCellsX,  "HaloMerger_Halo_%d_UM_IC_NCellsX",  index_halo_input );
         sprintf( HaloMerger_Halo_i_UM_IC_NCellsY,  "HaloMerger_Halo_%d_UM_IC_NCellsY",  index_halo_input );
         sprintf( HaloMerger_Halo_i_UM_IC_NCellsZ,  "HaloMerger_Halo_%d_UM_IC_NCellsZ",  index_halo_input );
         sprintf( HaloMerger_Halo_i_UM_IC_Float8,   "HaloMerger_Halo_%d_UM_IC_Float8",   index_halo_input );
         } // if ( HaloMerger_Halo_InitMode == 1 )
         else
            Aux_Error( ERROR_INFO, "unsupported initialization mode (%s = %d) !!\n",
                       "HaloMerger_Halo_InitMode", HaloMerger_Halo_InitMode );

      // (1-2-3) add parameters in the following format:
      // --> note that VARIABLE, DEFAULT, MIN, and MAX must have the same data type
      // --> some handy constants (e.g., NoMin_int, Eps_float, ...) are defined in "include/ReadPara.h"
      // ********************************************************************************************************************************
      // ReadPara_Halo->Add( "KEY_IN_THE_FILE",                       &VARIABLE,                                      DEFAULT,          MIN,           MAX              );
      // ********************************************************************************************************************************
         ReadPara_Halo->Add( HaloMerger_Halo_i_CenCoordX,             &HaloMerger_Halo_CenCoord[index_halo][0],      -1.0,              NoMin_double,  amr->BoxEdgeR[0] );
         ReadPara_Halo->Add( HaloMerger_Halo_i_CenCoordY,             &HaloMerger_Halo_CenCoord[index_halo][1],      -1.0,              NoMin_double,  amr->BoxEdgeR[1] );
         ReadPara_Halo->Add( HaloMerger_Halo_i_CenCoordZ,             &HaloMerger_Halo_CenCoord[index_halo][2],      -1.0,              NoMin_double,  amr->BoxEdgeR[2] );
         ReadPara_Halo->Add( HaloMerger_Halo_i_VelocityX,             &HaloMerger_Halo_Velocity[index_halo][0],       0.0,              NoMin_double,  NoMax_double     );
         ReadPara_Halo->Add( HaloMerger_Halo_i_VelocityY,             &HaloMerger_Halo_Velocity[index_halo][1],       0.0,              NoMin_double,  NoMax_double     );
         ReadPara_Halo->Add( HaloMerger_Halo_i_VelocityZ,             &HaloMerger_Halo_Velocity[index_halo][2],       0.0,              NoMin_double,  NoMax_double     );

         if ( HaloMerger_Halo_InitMode == 1 )
         {
         ReadPara_Halo->Add( HaloMerger_Halo_i_UM_IC_Filename,         HaloMerger_Halo_UM_IC_Filename[index_halo],    NoDef_str,        Useless_str,   Useless_str      );
         ReadPara_Halo->Add( HaloMerger_Halo_i_UM_IC_BoxLenX,         &HaloMerger_Halo_UM_IC_BoxLen[index_halo][0],  -1.0,              Eps_double,    NoMax_double     );
         ReadPara_Halo->Add( HaloMerger_Halo_i_UM_IC_BoxLenY,         &HaloMerger_Halo_UM_IC_BoxLen[index_halo][1],  -1.0,              Eps_double,    NoMax_double     );
         ReadPara_Halo->Add( HaloMerger_Halo_i_UM_IC_BoxLenZ,         &HaloMerger_Halo_UM_IC_BoxLen[index_halo][2],  -1.0,              Eps_double,    NoMax_double     );
         ReadPara_Halo->Add( HaloMerger_Halo_i_UM_IC_NCellsX,         &HaloMerger_Halo_UM_IC_NCells[index_halo][0],  -1,                1,             NoMax_int        );
         ReadPara_Halo->Add( HaloMerger_Halo_i_UM_IC_NCellsY,         &HaloMerger_Halo_UM_IC_NCells[index_halo][1],  -1,                1,             NoMax_int        );
         ReadPara_Halo->Add( HaloMerger_Halo_i_UM_IC_NCellsZ,         &HaloMerger_Halo_UM_IC_NCells[index_halo][2],  -1,                1,             NoMax_int        );
         ReadPara_Halo->Add( HaloMerger_Halo_i_UM_IC_Float8,          &HaloMerger_Halo_UM_IC_Float8[index_halo],      0,                0,             1                );
         } // if ( HaloMerger_Halo_InitMode == 1 )
         else
            Aux_Error( ERROR_INFO, "unsupported initialization mode (%s = %d) !!\n",
                       "HaloMerger_Halo_InitMode", HaloMerger_Halo_InitMode );

      } // for (int index_halo=0; index_halo<HaloMerger_Halo_Num; index_halo++)

      ReadPara_Halo->Read( FileName_Halo );

      delete ReadPara_Halo;

   } // if ( OPT__INIT != INIT_BY_RESTART  &&  HaloMerger_Halo_Num > 0 )

// (1-3) load the runtime parameters for the solitons
   if ( OPT__INIT != INIT_BY_RESTART  &&  HaloMerger_Soliton_Num > 0 )
   {
      // (1-3-1) allocate the memory
      HaloMerger_Soliton_CenCoord          = new double [HaloMerger_Soliton_Num][3];
      HaloMerger_Soliton_Velocity          = new double [HaloMerger_Soliton_Num][3];
      HaloMerger_Soliton_CoreRadius        = new double [HaloMerger_Soliton_Num];
      HaloMerger_Soliton_CoreRho           = new double [HaloMerger_Soliton_Num];

      if ( HaloMerger_Soliton_InitMode == 1 )
      {
      HaloMerger_Soliton_DensProf_Filename = new char    [HaloMerger_Soliton_Num][MAX_STRING];
      HaloMerger_Soliton_DensProf          = new double* [HaloMerger_Soliton_Num];
      HaloMerger_Soliton_DensProf_NBin     = new int     [HaloMerger_Soliton_Num];
      HaloMerger_Soliton_DensProf_Rescale  = new int     [HaloMerger_Soliton_Num];
      HaloMerger_Soliton_DensProf_ScaleL   = new double  [HaloMerger_Soliton_Num];
      HaloMerger_Soliton_DensProf_ScaleD   = new double  [HaloMerger_Soliton_Num];
      } // if ( HaloMerger_Soliton_InitMode == 1 )
      else if ( HaloMerger_Soliton_InitMode == 2 )
      {
      HaloMerger_Soliton_OuterSlope        = new double  [HaloMerger_Soliton_Num];
      } // else if ( HaloMerger_Soliton_InitMode == 2 )
      else
         Aux_Error( ERROR_INFO, "unsupported initialization mode (%s = %d) !!\n",
                    "HaloMerger_Soliton_InitMode", HaloMerger_Soliton_InitMode );

      // (1-3-2) read the parameters for the solitons
      const char FileName_Soliton[] = "Input__TestProb_Soliton";
      ReadPara_t *ReadPara_Soliton  = new ReadPara_t;

      // Soliton-related parameters to read from the input
      char HaloMerger_Soliton_i_CoreRadius[MAX_STRING];        // core radius of the i-th soliton (<=0.0=set by HaloMerger_Soliton_i_CoreRho) (will be overwritten if HaloMerger_Soliton_i_DensProf_Rescale == 0) [-1.0]
      char HaloMerger_Soliton_i_CoreRho   [MAX_STRING];        // peak density of the i-th soliton (will be overwritten if HaloMerger_Soliton_i_CoreRadius > 0.0) (will be overwritten if HaloMerger_Soliton_i_DensProf_Rescale == 0) [-1.0]
      char HaloMerger_Soliton_i_CenCoordX [MAX_STRING];        // x/y/z-coordinate of the center of the i-th soliton (<0.0=auto -> box center) [-1.0]
      char HaloMerger_Soliton_i_CenCoordY [MAX_STRING];
      char HaloMerger_Soliton_i_CenCoordZ [MAX_STRING];
      char HaloMerger_Soliton_i_VelocityX [MAX_STRING];        // x/y/z-component of the bulk velocity of the i-th soliton [0.0]
      char HaloMerger_Soliton_i_VelocityY [MAX_STRING];
      char HaloMerger_Soliton_i_VelocityZ [MAX_STRING];

      char HaloMerger_Soliton_i_DensProf_Filename[MAX_STRING]; // filename of the density profile table for the i-th soliton (HaloMerger_Soliton_InitMode == 1 only)
      char HaloMerger_Soliton_i_DensProf_Rescale[MAX_STRING];  // whether to scale the density profile table for the i-th soliton (HaloMerger_Soliton_InitMode == 1 only) [1]

      char HaloMerger_Soliton_i_OuterSlope[MAX_STRING];        // outer slope of the analytical density profile of the i-th soliton (HaloMerger_Soliton_InitMode == 2 only) [-8.0]

      for (int index_soliton=0; index_soliton<HaloMerger_Soliton_Num; index_soliton++)
      {
         const int index_soliton_input = index_soliton+1; // index of soliton in the input file

         sprintf( HaloMerger_Soliton_i_CoreRadius,        "HaloMerger_Soliton_%d_CoreRadius",        index_soliton_input );
         sprintf( HaloMerger_Soliton_i_CoreRho,           "HaloMerger_Soliton_%d_CoreRho",           index_soliton_input );
         sprintf( HaloMerger_Soliton_i_CenCoordX,         "HaloMerger_Soliton_%d_CenCoordX",         index_soliton_input );
         sprintf( HaloMerger_Soliton_i_CenCoordY,         "HaloMerger_Soliton_%d_CenCoordY",         index_soliton_input );
         sprintf( HaloMerger_Soliton_i_CenCoordZ,         "HaloMerger_Soliton_%d_CenCoordZ",         index_soliton_input );
         sprintf( HaloMerger_Soliton_i_VelocityX,         "HaloMerger_Soliton_%d_VelocityX",         index_soliton_input );
         sprintf( HaloMerger_Soliton_i_VelocityY,         "HaloMerger_Soliton_%d_VelocityY",         index_soliton_input );
         sprintf( HaloMerger_Soliton_i_VelocityZ,         "HaloMerger_Soliton_%d_VelocityZ",         index_soliton_input );

         if ( HaloMerger_Soliton_InitMode == 1 )
         {
         sprintf( HaloMerger_Soliton_i_DensProf_Filename, "HaloMerger_Soliton_%d_DensProf_Filename", index_soliton_input );
         sprintf( HaloMerger_Soliton_i_DensProf_Rescale,  "HaloMerger_Soliton_%d_DensProf_Rescale",  index_soliton_input );
         } // if ( HaloMerger_Soliton_InitMode == 1 )
         else if ( HaloMerger_Soliton_InitMode == 2 )
         {
         sprintf( HaloMerger_Soliton_i_OuterSlope,        "HaloMerger_Soliton_%d_OuterSlope",        index_soliton_input );
         } // else if ( HaloMerger_Soliton_InitMode == 2 )
         else
            Aux_Error( ERROR_INFO, "unsupported initialization mode (%s = %d) !!\n",
                       "HaloMerger_Soliton_InitMode", HaloMerger_Soliton_InitMode );

      // (1-3-3) add parameters in the following format:
      // --> note that VARIABLE, DEFAULT, MIN, and MAX must have the same data type
      // --> some handy constants (e.g., NoMin_int, Eps_float, ...) are defined in "include/ReadPara.h"
      // ********************************************************************************************************************************
      // ReadPara_Soliton->Add( "KEY_IN_THE_FILE",                          &VARIABLE,                                               DEFAULT,          MIN,           MAX              );
      // ********************************************************************************************************************************
         ReadPara_Soliton->Add( HaloMerger_Soliton_i_CoreRadius,            &HaloMerger_Soliton_CoreRadius[index_soliton],          -1.0,              NoMin_double,  NoMax_double     );
         ReadPara_Soliton->Add( HaloMerger_Soliton_i_CoreRho,               &HaloMerger_Soliton_CoreRho[index_soliton],             -1.0,              NoMin_double,  NoMax_double     );
         ReadPara_Soliton->Add( HaloMerger_Soliton_i_CenCoordX,             &HaloMerger_Soliton_CenCoord[index_soliton][0],         -1.0,              NoMin_double,  amr->BoxEdgeR[0] );
         ReadPara_Soliton->Add( HaloMerger_Soliton_i_CenCoordY,             &HaloMerger_Soliton_CenCoord[index_soliton][1],         -1.0,              NoMin_double,  amr->BoxEdgeR[1] );
         ReadPara_Soliton->Add( HaloMerger_Soliton_i_CenCoordZ,             &HaloMerger_Soliton_CenCoord[index_soliton][2],         -1.0,              NoMin_double,  amr->BoxEdgeR[2] );
         ReadPara_Soliton->Add( HaloMerger_Soliton_i_VelocityX,             &HaloMerger_Soliton_Velocity[index_soliton][0],          0.0,              NoMin_double,  NoMax_double     );
         ReadPara_Soliton->Add( HaloMerger_Soliton_i_VelocityY,             &HaloMerger_Soliton_Velocity[index_soliton][1],          0.0,              NoMin_double,  NoMax_double     );
         ReadPara_Soliton->Add( HaloMerger_Soliton_i_VelocityZ,             &HaloMerger_Soliton_Velocity[index_soliton][2],          0.0,              NoMin_double,  NoMax_double     );

         if ( HaloMerger_Soliton_InitMode == 1 )
         {
         ReadPara_Soliton->Add( HaloMerger_Soliton_i_DensProf_Filename,      HaloMerger_Soliton_DensProf_Filename[index_soliton],    NoDef_str,        Useless_str,   Useless_str      );
         ReadPara_Soliton->Add( HaloMerger_Soliton_i_DensProf_Rescale,      &HaloMerger_Soliton_DensProf_Rescale[index_soliton],     1,                0,             1                );
         } // if ( HaloMerger_Soliton_InitMode == 1 )
         else if ( HaloMerger_Soliton_InitMode == 2 )
         {
         ReadPara_Soliton->Add( HaloMerger_Soliton_i_OuterSlope,            &HaloMerger_Soliton_OuterSlope[index_soliton],          -8.0,              NoMin_double,  NoMax_double     );
         } // else if ( HaloMerger_Soliton_InitMode == 2 )
         else
            Aux_Error( ERROR_INFO, "unsupported initialization mode (%s = %d) !!\n",
                       "HaloMerger_Soliton_InitMode", HaloMerger_Soliton_InitMode );

      } // for (int index_soliton=0; index_soliton<HaloMerger_Soliton_Num; index_soliton++)

      ReadPara_Soliton->Read( FileName_Soliton );

      delete ReadPara_Soliton;

   } // if ( OPT__INIT != INIT_BY_RESTART  &&  HaloMerger_Soliton_Num > 0 )

// (1-4) load the runtime parameters for the particle clouds
   if ( OPT__INIT != INIT_BY_RESTART  &&  HaloMerger_ParCloud_Num > 0 )
   {
      // (1-4-1) allocate the memory
      HaloMerger_ParCloud_CenCoord          = new double [HaloMerger_ParCloud_Num][3];
      HaloMerger_ParCloud_Velocity          = new double [HaloMerger_ParCloud_Num][3];

      if ( HaloMerger_ParCloud_InitMode == 1 )
      {
      HaloMerger_ParCloud_DensProf_Filename = new char   [HaloMerger_ParCloud_Num][MAX_STRING];
      HaloMerger_ParCloud_DensProf_MaxR     = new double [HaloMerger_ParCloud_Num];
      HaloMerger_ParCloud_RSeed             = new int    [HaloMerger_ParCloud_Num];
      HaloMerger_ParCloud_NPar              = new long   [HaloMerger_ParCloud_Num];
      } // if ( HaloMerger_ParCloud_InitMode == 1 )
      else
         Aux_Error( ERROR_INFO, "unsupported initialization mode (%s = %d) !!\n",
                    "HaloMerger_ParCloud_InitMode", HaloMerger_ParCloud_InitMode );

      // (1-4-2) read the parameters for the halos
      const char FileName_ParCloud[] = "Input__TestProb_ParCloud";
      ReadPara_t *ReadPara_ParCloud  = new ReadPara_t;

      // ParCloud-related parameters to read from the input
      char HaloMerger_ParCloud_i_CenCoordX[MAX_STRING];          // x/y/z-coordinate of the center of the i-th particle cloud (<0.0=auto -> box center) [-1.0]
      char HaloMerger_ParCloud_i_CenCoordY[MAX_STRING];
      char HaloMerger_ParCloud_i_CenCoordZ[MAX_STRING];
      char HaloMerger_ParCloud_i_VelocityX[MAX_STRING];          // x/y/z-component of the bulk velocity of the i-th particle cloud [0.0]
      char HaloMerger_ParCloud_i_VelocityY[MAX_STRING];
      char HaloMerger_ParCloud_i_VelocityZ[MAX_STRING];
      char HaloMerger_ParCloud_i_DensProf_Filename[MAX_STRING];  // filename of the density profile table for the i-th particle cloud (HaloMerger_ParCloud_InitMode == 1 only)
      char HaloMerger_ParCloud_i_DensProf_MaxR[MAX_STRING];      // maximum radius for particles for the i-th particle cloud (must > 0.0) (HaloMerger_ParCloud_InitMode == 1 only) [0.5*amr->BoxSize[0]]
      char HaloMerger_ParCloud_i_RSeed[MAX_STRING];              // random seed for setting particle position and velocity for the i-th particle cloud (must >= 0) (HaloMerger_ParCloud_InitMode == 1 only) [123]
      char HaloMerger_ParCloud_i_NPar[MAX_STRING];               // number of particles for the i-th particle cloud (must >= 0) (HaloMerger_ParCloud_InitMode == 1 only) [0]

      for (int index_parcloud=0; index_parcloud<HaloMerger_ParCloud_Num; index_parcloud++)
      {
         const int index_parcloud_input = index_parcloud+1; // index of particle cloud in the input file

         sprintf( HaloMerger_ParCloud_i_CenCoordX,         "HaloMerger_ParCloud_%d_CenCoordX",         index_parcloud_input );
         sprintf( HaloMerger_ParCloud_i_CenCoordY,         "HaloMerger_ParCloud_%d_CenCoordY",         index_parcloud_input );
         sprintf( HaloMerger_ParCloud_i_CenCoordZ,         "HaloMerger_ParCloud_%d_CenCoordZ",         index_parcloud_input );
         sprintf( HaloMerger_ParCloud_i_VelocityX,         "HaloMerger_ParCloud_%d_VelocityX",         index_parcloud_input );
         sprintf( HaloMerger_ParCloud_i_VelocityY,         "HaloMerger_ParCloud_%d_VelocityY",         index_parcloud_input );
         sprintf( HaloMerger_ParCloud_i_VelocityZ,         "HaloMerger_ParCloud_%d_VelocityZ",         index_parcloud_input );

         if ( HaloMerger_ParCloud_InitMode == 1 )
         {
         sprintf( HaloMerger_ParCloud_i_DensProf_Filename, "HaloMerger_ParCloud_%d_DensProf_Filename", index_parcloud_input );
         sprintf( HaloMerger_ParCloud_i_DensProf_MaxR,     "HaloMerger_ParCloud_%d_DensProf_MaxR",     index_parcloud_input );
         sprintf( HaloMerger_ParCloud_i_RSeed,             "HaloMerger_ParCloud_%d_RSeed",             index_parcloud_input );
         sprintf( HaloMerger_ParCloud_i_NPar,              "HaloMerger_ParCloud_%d_NPar",              index_parcloud_input );
         } // if ( HaloMerger_ParCloud_InitMode == 1 )
         else
            Aux_Error( ERROR_INFO, "unsupported initialization mode (%s = %d) !!\n",
                       "HaloMerger_ParCloud_InitMode", HaloMerger_ParCloud_InitMode );

      // (1-4-3) add parameters in the following format:
      // --> note that VARIABLE, DEFAULT, MIN, and MAX must have the same data type
      // --> some handy constants (e.g., NoMin_int, Eps_float, ...) are defined in "include/ReadPara.h"
      // ********************************************************************************************************************************
      // ReadPara_ParCloud->Add( "KEY_IN_THE_FILE",                       &VARIABLE,                                               DEFAULT,              MIN,           MAX              );
      // ********************************************************************************************************************************
         ReadPara_ParCloud->Add( HaloMerger_ParCloud_i_CenCoordX,         &HaloMerger_ParCloud_CenCoord[index_parcloud][0],       -1.0,                  NoMin_double,  amr->BoxEdgeR[0] );
         ReadPara_ParCloud->Add( HaloMerger_ParCloud_i_CenCoordY,         &HaloMerger_ParCloud_CenCoord[index_parcloud][1],       -1.0,                  NoMin_double,  amr->BoxEdgeR[1] );
         ReadPara_ParCloud->Add( HaloMerger_ParCloud_i_CenCoordZ,         &HaloMerger_ParCloud_CenCoord[index_parcloud][2],       -1.0,                  NoMin_double,  amr->BoxEdgeR[2] );
         ReadPara_ParCloud->Add( HaloMerger_ParCloud_i_VelocityX,         &HaloMerger_ParCloud_Velocity[index_parcloud][0],        0.0,                  NoMin_double,  NoMax_double     );
         ReadPara_ParCloud->Add( HaloMerger_ParCloud_i_VelocityY,         &HaloMerger_ParCloud_Velocity[index_parcloud][1],        0.0,                  NoMin_double,  NoMax_double     );
         ReadPara_ParCloud->Add( HaloMerger_ParCloud_i_VelocityZ,         &HaloMerger_ParCloud_Velocity[index_parcloud][2],        0.0,                  NoMin_double,  NoMax_double     );

         if ( HaloMerger_ParCloud_InitMode == 1 )
         {
         ReadPara_ParCloud->Add( HaloMerger_ParCloud_i_DensProf_Filename,  HaloMerger_ParCloud_DensProf_Filename[index_parcloud],  NoDef_str,            Useless_str,   Useless_str      );
         ReadPara_ParCloud->Add( HaloMerger_ParCloud_i_DensProf_MaxR,     &HaloMerger_ParCloud_DensProf_MaxR[index_parcloud],      0.5*amr->BoxSize[0],  Eps_double,    NoMax_double     );
         ReadPara_ParCloud->Add( HaloMerger_ParCloud_i_RSeed,             &HaloMerger_ParCloud_RSeed[index_parcloud],              123,                  0,             NoMax_int        );
         ReadPara_ParCloud->Add( HaloMerger_ParCloud_i_NPar,              &HaloMerger_ParCloud_NPar[index_parcloud],               (long)0,              (long)0,       NoMax_long       );
         } // if ( HaloMerger_ParCloud_InitMode == 1 )
         else
            Aux_Error( ERROR_INFO, "unsupported initialization mode (%s = %d) !!\n",
                       "HaloMerger_ParCloud_InitMode", HaloMerger_ParCloud_InitMode );

      } // for (int index_parcloud=0; index_parcloud<HaloMerger_ParCloud_Num; index_parcloud++)

      ReadPara_ParCloud->Read( FileName_ParCloud );

      delete ReadPara_ParCloud;

   } // if ( OPT__INIT != INIT_BY_RESTART  &&  HaloMerger_ParCloud_Num > 0 )


// (2) check the runtime parameters
// (2-1) check the parameters for the external potential
   if ( HaloMerger_ExtPot_UniDenSph_M != 0.0  &&  OPT__EXT_POT != EXT_POT_FUNC )
      Aux_Error( ERROR_INFO, "OPT__EXT_POT must be EXT_POT_FUNC (%d) to add the external potential of a uniform-density sphere !!\n", EXT_POT_FUNC );

   double HaloMerger_ExtPot_UniDenSph_Rho;  // density of the uniform-density sphere for the external potential

   if ( OPT__EXT_POT == EXT_POT_FUNC )
   {
      // set the center
      // put at the box ceneter by default
      if ( HaloMerger_ExtPot_UniDenSph_CenCoordX < 0.0 )    HaloMerger_ExtPot_UniDenSph_CenCoordX = amr->BoxCenter[0];
      if ( HaloMerger_ExtPot_UniDenSph_CenCoordY < 0.0 )    HaloMerger_ExtPot_UniDenSph_CenCoordY = amr->BoxCenter[1];
      if ( HaloMerger_ExtPot_UniDenSph_CenCoordZ < 0.0 )    HaloMerger_ExtPot_UniDenSph_CenCoordZ = amr->BoxCenter[2];

      // set the rho as mass/volume
      HaloMerger_ExtPot_UniDenSph_Rho = HaloMerger_ExtPot_UniDenSph_M/( 4.0*M_PI*CUBE(HaloMerger_ExtPot_UniDenSph_R)/3.0 );

   } // if ( OPT__EXT_POT == EXT_POT_FUNC )

// (2-2) check the parameters for the halos
   if ( OPT__INIT != INIT_BY_RESTART  &&  HaloMerger_Halo_Num > 0 )
   {
      for (int index_halo=0; index_halo<HaloMerger_Halo_Num; index_halo++)
      {
         const int index_halo_input = index_halo+1; // index of halo in the input file

         // check the center for the halos
         for (int d=0; d<3; d++)
         if ( HaloMerger_Halo_CenCoord[index_halo][d] < 0.0 )
            HaloMerger_Halo_CenCoord[index_halo][d] = amr->BoxCenter[d]; // put at the box ceneter by default

         // check and set the UM_IC-related parameters
         if ( HaloMerger_Halo_InitMode == 1 )
         {
            // check the UM_IC file exists
            if ( !Aux_CheckFileExist(HaloMerger_Halo_UM_IC_Filename[index_halo]) )
               Aux_Error( ERROR_INFO, "Halo_%d UM_IC file \"%s\" does not exist !!\n",
                          index_halo_input, HaloMerger_Halo_UM_IC_Filename[index_halo] );

            // check the range of the halos
            for (int d=0; d<3; d++)
            {
               // derived parameters from the input for the halos
               HaloMerger_Halo_UM_IC_dh         [index_halo][d] = HaloMerger_Halo_UM_IC_BoxLen[index_halo][d]/HaloMerger_Halo_UM_IC_NCells[index_halo][d];
               HaloMerger_Halo_UM_IC_Range_EdgeL[index_halo][d] = HaloMerger_Halo_CenCoord[index_halo][d] - 0.5*HaloMerger_Halo_UM_IC_BoxLen[index_halo][d];
               HaloMerger_Halo_UM_IC_Range_EdgeR[index_halo][d] = HaloMerger_Halo_CenCoord[index_halo][d] + 0.5*HaloMerger_Halo_UM_IC_BoxLen[index_halo][d];

               // check wheter the UM_IC dh are the same along different directions
               for (int d2=0; d2<d; d2++)
                  if ( ! Mis_CompareRealValue( HaloMerger_Halo_UM_IC_dh[index_halo][d], HaloMerger_Halo_UM_IC_dh[index_halo][d2], NULL, false )  )
                     if ( MPI_Rank == 0 )
                        Aux_Message( stderr, "WARNING : The cell size (dh = %13.6e) along direction-%d is different from the cell size (dh = %13.6e) along direction-%d for Halo_%d UM_IC !!\n",
                                             HaloMerger_Halo_UM_IC_dh[index_halo][d], d, HaloMerger_Halo_UM_IC_dh[index_halo][d2], d2, index_halo_input );

               // check whether the input halos cross the boundary
               if ( HaloMerger_Halo_UM_IC_Range_EdgeR[index_halo][d] > amr->BoxEdgeR[d]  ||
                    HaloMerger_Halo_UM_IC_Range_EdgeL[index_halo][d] < amr->BoxEdgeL[d] )
                  Aux_Error( ERROR_INFO, "The edge in direction-%d [%13.6e, %13.6e] of Halo_%d UM_IC range is outside of the simulation box !!\n",
                             d, HaloMerger_Halo_UM_IC_Range_EdgeL[index_halo][d],
                             HaloMerger_Halo_UM_IC_Range_EdgeR[index_halo][d], index_halo_input );

            } // for (int d=0; d<3; d++)

            // check whether the input halos overlap with each other
            for (int index2_halo=0; index2_halo<index_halo; index2_halo++)
            {
               const int index2_halo_input = index2_halo+1; // index of halo in the input file

               bool isOverlap = true;

               // three directions
               for (int d=0; d<3; d++)
               {
                  if ( HaloMerger_Halo_UM_IC_Range_EdgeR[index_halo][d]  <= HaloMerger_Halo_UM_IC_Range_EdgeL[index2_halo][d] ||
                       HaloMerger_Halo_UM_IC_Range_EdgeR[index2_halo][d] <= HaloMerger_Halo_UM_IC_Range_EdgeL[index_halo][d] )
                  {
                     isOverlap = false;
                     break;
                  }
               } // for (int d=0; d<3; d++)

               if ( MPI_Rank == 0 )
               {
                  if ( isOverlap )
                     Aux_Message( stderr, "WARNING : Halo_%d UM_IC range overlaps with the Halo_%d UM_IC range !!\n",
                                  index_halo_input, index2_halo_input );
               }

            } // for (int index2_halo=0; index2_halo<index_halo; index2_halo++)

         } // if ( HaloMerger_Halo_InitMode == 1 )
         else
            Aux_Error( ERROR_INFO, "unsupported initialization mode (%s = %d) !!\n",
                       "HaloMerger_Halo_InitMode", HaloMerger_Halo_InitMode );

      } // for (int index_halo=0; index_halo<HaloMerger_Halo_Num; index_halo++)

      // load the UM_IC data for the halos when HaloMerger_Halo_InitMode == 1
      if ( HaloMerger_Halo_InitMode == 1 )
      {
         const long UM_IC_NVar = 2; // (Real part & Imag part)

         for (int index_halo=0; index_halo<HaloMerger_Halo_Num; index_halo++)
         {
             // UM_IC information
             const long UM_IC_NCells3D = (long)HaloMerger_Halo_UM_IC_NCells[index_halo][0]*HaloMerger_Halo_UM_IC_NCells[index_halo][1]*HaloMerger_Halo_UM_IC_NCells[index_halo][2];
             size_t load_data_size     = ( HaloMerger_Halo_UM_IC_Float8[index_halo] ) ? sizeof(double) : sizeof(float);
             const long UM_IC_TotSize  = UM_IC_NVar*UM_IC_NCells3D*load_data_size;

             // allocate the memory for the array to read the data
             HaloMerger_Halo_UM_IC_Data[index_halo] = new char [ UM_IC_TotSize ];

             // open the file
             FILE *File = fopen( HaloMerger_Halo_UM_IC_Filename[index_halo], "rb" );

             // check the file size
             fseek( File, 0, SEEK_END );
             const long ExpectSize = UM_IC_TotSize;
             const long FileSize   = ftell( File );
             if ( FileSize != ExpectSize )
                Aux_Error( ERROR_INFO, "size of the UM_IC <%s> (%ld) != expect (%ld) !!\n", HaloMerger_Halo_UM_IC_Filename[index_halo], FileSize, ExpectSize );

             // load data from the file
             fseek( File, 0, SEEK_SET );
             fread( HaloMerger_Halo_UM_IC_Data[index_halo], 1, UM_IC_TotSize, File );

             // close the file
             fclose( File );

         } // for (int index_halo=0; index_halo<HaloMerger_Halo_Num; index_halo++)

      } // if ( HaloMerger_Halo_InitMode == 1 )
      else
         Aux_Error( ERROR_INFO, "unsupported initialization mode (%s = %d) !!\n",
                    "HaloMerger_Halo_InitMode", HaloMerger_Halo_InitMode );

   } // if ( OPT__INIT != INIT_BY_RESTART  &&  HaloMerger_Halo_Num > 0 )

// (2-3) check the parameters for the solitons
   if ( OPT__INIT != INIT_BY_RESTART  &&  HaloMerger_Soliton_Num > 0 )
   {
      for (int index_soliton=0; index_soliton<HaloMerger_Soliton_Num; index_soliton++)
      {
         const int index_soliton_input = index_soliton+1; // index of soliton in the input file

         // set the parameters for the density profile
         if ( HaloMerger_Soliton_InitMode == 1 ) // the density profile from table
         {
            // check the file exist
            if ( !Aux_CheckFileExist(HaloMerger_Soliton_DensProf_Filename[index_soliton]) )
               Aux_Error( ERROR_INFO, "Soliton_%d DensProf file \"%s\" does not exist !!\n",
                       index_soliton_input, HaloMerger_Soliton_DensProf_Filename[index_soliton] );

            // load the reference profile
            const bool RowMajor_No  = false;    // load data into the column-major order
            const bool AllocMem_Yes = true;     // allocate memory for HaloMerger_Soliton_DensProf
            const int  NCol         = 2;        // total number of columns to load
            const int  Col[NCol]    = {0, 1};   // target columns: (radius, density)

            HaloMerger_Soliton_DensProf_NBin[index_soliton] = Aux_LoadTable( HaloMerger_Soliton_DensProf[index_soliton],
                                                                             HaloMerger_Soliton_DensProf_Filename[index_soliton],
                                                                             NCol, Col, RowMajor_No, AllocMem_Yes );

            // get the core radius of the reference profile
            const double *RadiusRef = HaloMerger_Soliton_DensProf[index_soliton] + 0*HaloMerger_Soliton_DensProf_NBin[index_soliton];
            const double *DensRef   = HaloMerger_Soliton_DensProf[index_soliton] + 1*HaloMerger_Soliton_DensProf_NBin[index_soliton];
            const double  DensCore  = 0.5*DensRef[0];   // define core radius as the half-density radius

            double CoreRadiusRef = NULL_REAL;

            for (int b=1; b<HaloMerger_Soliton_DensProf_NBin[index_soliton]-1; b++)
            {
               if ( DensRef[b] >= DensCore  &&  DensRef[b+1] <= DensCore )
               {
                  CoreRadiusRef = 0.5*( RadiusRef[b] + RadiusRef[b+1] );
                  break;
               }
            } // for (int b=1; b<HaloMerger_Soliton_DensProf_NBin[index_soliton]-1; b++)

            if ( CoreRadiusRef == NULL_REAL )
               Aux_Error( ERROR_INFO, "cannot determine the reference core radius !!\n" );

            if ( HaloMerger_Soliton_DensProf_Rescale[index_soliton] )
            {
               // evaluate the scale factors of each soliton
               if ( HaloMerger_Soliton_CoreRadius[index_soliton] <= 0.0 )
               {
                  if ( HaloMerger_Soliton_CoreRho[index_soliton] > 0.0 )
                  {
                     // overwrite the core radius by the value calculated from the peak density
                     HaloMerger_Soliton_DensProf_ScaleD[index_soliton] = HaloMerger_Soliton_CoreRho[index_soliton] / DensRef[0];
                     HaloMerger_Soliton_DensProf_ScaleL[index_soliton] = sqrt( sqrt( 1.0 / (4.0*M_PI*NEWTON_G*SQR(ELBDM_ETA)*HaloMerger_Soliton_DensProf_ScaleD[index_soliton]) ) );
                     HaloMerger_Soliton_CoreRadius     [index_soliton] = CoreRadiusRef*HaloMerger_Soliton_DensProf_ScaleL[index_soliton];
                  }
                  else // if ( HaloMerger_Soliton_CoreRho[index_soliton] > 0.0 )
                  {
                     Aux_Error( ERROR_INFO, "HaloMerger_Soliton_%d_CoreRadius (%13.6e) is not set properly !!\n", index_soliton_input, HaloMerger_Soliton_CoreRadius[index_soliton] );
                  } // if ( HaloMerger_Soliton_CoreRho[index_soliton] > 0.0 ) ... else
               }
               else // if ( HaloMerger_Soliton_CoreRadius[index_soliton] <= 0.0 )
               {
                  // overwrite the peak density by the value calculated from the core radius
                  HaloMerger_Soliton_DensProf_ScaleL[index_soliton] = HaloMerger_Soliton_CoreRadius[index_soliton] / CoreRadiusRef;
                  HaloMerger_Soliton_DensProf_ScaleD[index_soliton] = 1.0 / ( 4.0*M_PI*NEWTON_G*SQR(ELBDM_ETA)*POW4(HaloMerger_Soliton_DensProf_ScaleL[index_soliton]) );
                  HaloMerger_Soliton_CoreRho        [index_soliton] = HaloMerger_Soliton_DensProf_ScaleD[index_soliton]*DensRef[0];
               } // if ( HaloMerger_Soliton_CoreRadius[index_soliton] <= 0.0 ) ... else

            } // if ( HaloMerger_Soliton_DensProf_Rescale[index_soliton] )
            else
            {
               // overwrite the peak density and core radius from the table
               HaloMerger_Soliton_DensProf_ScaleL[index_soliton] = 1.0;
               HaloMerger_Soliton_DensProf_ScaleD[index_soliton] = 1.0;
               HaloMerger_Soliton_CoreRho        [index_soliton] = DensRef[0];
               HaloMerger_Soliton_CoreRadius     [index_soliton] = CoreRadiusRef;
            } // if ( HaloMerger_Soliton_DensProf_Rescale[index_soliton] ) ... else

         } // if ( HaloMerger_Soliton_InitMode == 1 )
         else if ( HaloMerger_Soliton_InitMode == 2 ) // the density profile from the analytical function
         {
            // set the core radius and the peak density
            if ( HaloMerger_Soliton_CoreRadius[index_soliton] <= 0.0 )
            {
               if ( HaloMerger_Soliton_CoreRho[index_soliton] > 0.0 )
               {
                  // overwrite the core radius by the value calculated from the peak density
                  const double m22                             = ELBDM_MASS*UNIT_M/(Const_eV/SQR(Const_c))/1.0e-22;
                  const double rc_kpc                          = sqrt( sqrt( 1.945e7/(HaloMerger_Soliton_CoreRho[index_soliton]/Const_Msun*CUBE(Const_kpc)*(UNIT_M/CUBE(UNIT_L))) )/m22 );
                  HaloMerger_Soliton_CoreRadius[index_soliton] = rc_kpc*Const_kpc/UNIT_L;
               }
               else // if ( HaloMerger_Soliton_CoreRho[index_soliton] > 0.0 )
               {
                  Aux_Error( ERROR_INFO, "HaloMerger_Soliton_%d_CoreRadius (%13.6e) is not set properly !!\n", index_soliton_input, HaloMerger_Soliton_CoreRadius[index_soliton] );
               } // if ( HaloMerger_Soliton_CoreRho[index_soliton] > 0.0 ) ... else
            }
            else // if ( HaloMerger_Soliton_CoreRadius[index_soliton] <= 0.0 )
            {
               // overwrite the peak density by the value calculated from the core radius
               const double m22                          = ELBDM_MASS*UNIT_M/(Const_eV/SQR(Const_c))/1.0e-22;
               const double rc_kpc                       = HaloMerger_Soliton_CoreRadius[index_soliton]*UNIT_L/Const_kpc;
               HaloMerger_Soliton_CoreRho[index_soliton] = 1.945e7/SQR( m22*rc_kpc*rc_kpc )*Const_Msun/CUBE(Const_kpc)/(UNIT_M/CUBE(UNIT_L));
            } // if ( HaloMerger_Soliton_CoreRadius[index_soliton] <= 0.0 ) ... else

         } // else if ( HaloMerger_Soliton_InitMode == 2 )
         else
            Aux_Error( ERROR_INFO, "unsupported initialization mode (%s = %d) !!\n",
                       "HaloMerger_Soliton_InitMode", HaloMerger_Soliton_InitMode );

         // check the runtime parameters and set the problem-specific derived parameters for the solitons
         // set the center
         for (int d=0; d<3; d++)
         if ( HaloMerger_Soliton_CenCoord[index_soliton][d] < 0.0 )
            HaloMerger_Soliton_CenCoord[index_soliton][d] = amr->BoxCenter[d]; // put at the box center by default

         // check whether the soliton touches the boundary of box
         for (int d=0; d<3; d++)
         {
            // check whether the input solitons cross the boundary
            if ( HaloMerger_Soliton_CenCoord[index_soliton][d] + 3.0*HaloMerger_Soliton_CoreRadius[index_soliton] > amr->BoxEdgeR[d]  ||
                 HaloMerger_Soliton_CenCoord[index_soliton][d] - 3.0*HaloMerger_Soliton_CoreRadius[index_soliton] < amr->BoxEdgeL[d] )
               Aux_Error( ERROR_INFO, "The Soliton_%d 3r_c-range [%13.6e, %13.6e] is outside of the simulation box in the direction-%d  !!\n",
                          index_soliton_input,
                          HaloMerger_Soliton_CenCoord[index_soliton][d] - 3.0*HaloMerger_Soliton_CoreRadius[index_soliton],
                          HaloMerger_Soliton_CenCoord[index_soliton][d] + 3.0*HaloMerger_Soliton_CoreRadius[index_soliton],
                          d );
         } // for (int d=0; d<3; d++)

         // check whether the input solitons overlap with each other
         for (int index2_soliton=0; index2_soliton<index_soliton; index2_soliton++)
         {
           const int index2_soliton_input = index2_soliton+1; // index of soliton in the input file

            bool isOverlap = true;

            if ( sqrt( SQR( HaloMerger_Soliton_CenCoord[index_soliton][0] - HaloMerger_Soliton_CenCoord[index2_soliton][0] )
                     + SQR( HaloMerger_Soliton_CenCoord[index_soliton][1] - HaloMerger_Soliton_CenCoord[index2_soliton][1] )
                     + SQR( HaloMerger_Soliton_CenCoord[index_soliton][2] - HaloMerger_Soliton_CenCoord[index2_soliton][2] ) )
                 > 3.0*( HaloMerger_Soliton_CoreRadius[index_soliton] + HaloMerger_Soliton_CoreRadius[index2_soliton] ) )
               isOverlap = false;

            if ( MPI_Rank == 0 )
            {
               if ( isOverlap )
                  Aux_Message( stderr, "WARNING : Soliton_%d 3r_c-range (center: [%13.6e, %13.6e, %13.6e], r_c: %13.6e) overlaps with the Soliton_%d 3r_c-range (center: [%13.6e, %13.6e, %13.6e], r_c: %13.6e) !!\n",
                               index_soliton_input,
                               HaloMerger_Soliton_CenCoord[index_soliton][0], HaloMerger_Soliton_CenCoord[index_soliton][1],
                               HaloMerger_Soliton_CenCoord[index_soliton][2], HaloMerger_Soliton_CoreRadius[index_soliton],
                               index2_soliton_input,
                               HaloMerger_Soliton_CenCoord[index2_soliton][0], HaloMerger_Soliton_CenCoord[index2_soliton][1],
                               HaloMerger_Soliton_CenCoord[index2_soliton][2], HaloMerger_Soliton_CoreRadius[index2_soliton] );
            }

         } // for (int index2_soliton=0; index2_soliton<index_soliton; index2_soliton++)

         // check whether the input soliton overlaps with the input halos
         for (int index_halo=0; index_halo<HaloMerger_Halo_Num; index_halo++)
         {
            const int index_halo_input = index_halo+1; // index of halo in the input file

            if ( HaloMerger_Halo_InitMode == 1 )
            {
               bool isOverlap = true;

               for (int d=0; d<3; d++)
               {
                  if ( HaloMerger_Soliton_CenCoord[index_soliton][d] - 3.0*HaloMerger_Soliton_CoreRadius[index_soliton] >= HaloMerger_Halo_UM_IC_Range_EdgeR[index_halo][d] ||
                       HaloMerger_Soliton_CenCoord[index_soliton][d] + 3.0*HaloMerger_Soliton_CoreRadius[index_soliton] <= HaloMerger_Halo_UM_IC_Range_EdgeL[index_halo][d] )
                  {
                     isOverlap = false;
                     break;
                  }
               } // for (int d=0; d<3; d++)

               if ( MPI_Rank == 0 )
               {
                  if ( isOverlap )
                     Aux_Message( stderr, "WARNING : Soliton_%d 3r_c-range (center: [%13.6e, %13.6e, %13.6e], r_c: %13.6e) overlaps with the Halo_%d UM_IC range !!\n",
                                  index_soliton_input,
                                  HaloMerger_Soliton_CenCoord[index_soliton][0], HaloMerger_Soliton_CenCoord[index_soliton][1],
                                  HaloMerger_Soliton_CenCoord[index_soliton][2], HaloMerger_Soliton_CoreRadius[index_soliton],
                                  index_halo_input );
               }

            } // if ( HaloMerger_Halo_InitMode == 1 )
            else
               Aux_Error( ERROR_INFO, "unsupported initialization mode (%s = %d) !!\n",
                          "HaloMerger_Halo_InitMode", HaloMerger_Halo_InitMode );

         } // for (int index_halo=0; index_halo<HaloMerger_Halo_Num; index_halo++)

      } // for (int index_soliton=0; index_soliton<HaloMerger_Soliton_Num; index_soliton++)

   } // if ( OPT__INIT != INIT_BY_RESTART  &&  HaloMerger_Soliton_Num > 0 )

// (2-4) check the parameters for the particle clouds
   // set the total number of particles in all particle clouds
   long HaloMerger_ParCloud_NPar_Total = (long)0;

   if ( OPT__INIT != INIT_BY_RESTART  &&  HaloMerger_ParCloud_Num > 0 )
   {
      // check the particle is enabled
#     ifndef MASSIVE_PARTICLES
      Aux_Error( ERROR_INFO, "MASSIVE_PARTICLES must be enabled for HaloMerger_ParCloud_Num > 0 !!\n" );
#     endif

      // check there fluid is freezed if it is supposed to be a particle-only case
      if ( (HaloMerger_Halo_Num + HaloMerger_Soliton_Num) == 0  &&  OPT__FREEZE_FLUID != 1 )
         Aux_Error( ERROR_INFO, "OPT__FREEZE_FLUID should be 1 for particle-only simulations !!\n" );

      for (int index_parcloud=0; index_parcloud<HaloMerger_ParCloud_Num; index_parcloud++)
      {
         const int index_parcloud_input = index_parcloud+1; // index of particle cloud in the input file

         // check the center for the particle clouds
         for (int d=0; d<3; d++)
         if ( HaloMerger_ParCloud_CenCoord[index_parcloud][d] < 0.0 )
            HaloMerger_ParCloud_CenCoord[index_parcloud][d] = amr->BoxCenter[d]; // put at the box ceneter by default

         // check and set the parameters for the particle density profile
         if ( HaloMerger_ParCloud_InitMode == 1 ) // the density profile from table
         {
            // count the total number of particles for the particle clouds
            HaloMerger_ParCloud_NPar_Total += HaloMerger_ParCloud_NPar[index_parcloud];

            // check the density profile file exists
            if ( !Aux_CheckFileExist(HaloMerger_ParCloud_DensProf_Filename[index_parcloud]) )
               Aux_Error( ERROR_INFO, "ParCloud_%d density profile file \"%s\" does not exist !!\n",
                          index_parcloud_input, HaloMerger_ParCloud_DensProf_Filename[index_parcloud] );

            // check whether the particle cloud touches the boundary of box
            for (int d=0; d<3; d++)
            {
               // check whether the input particle clouds cross the boundary
               if ( HaloMerger_ParCloud_CenCoord[index_parcloud][d] + HaloMerger_ParCloud_DensProf_MaxR[index_parcloud] > amr->BoxEdgeR[d]  ||
                    HaloMerger_ParCloud_CenCoord[index_parcloud][d] - HaloMerger_ParCloud_DensProf_MaxR[index_parcloud] < amr->BoxEdgeL[d] )
                  Aux_Error( ERROR_INFO, "The ParCloud_%d range is outside of the simulation box in the direction-%d !!\n",
                             index_parcloud_input, d );
            } // for (int d=0; d<3; d++)

            // check whether the input particle clouds overlap with each other
            for (int index2_parcloud=0; index2_parcloud<index_parcloud; index2_parcloud++)
            {
               const int index2_parcloud_input = index2_parcloud+1; // index of particle cloud in the input file

               bool isOverlap = true;

               if ( sqrt( SQR( HaloMerger_ParCloud_CenCoord[index_parcloud][0] - HaloMerger_ParCloud_CenCoord[index2_parcloud][0] )
                        + SQR( HaloMerger_ParCloud_CenCoord[index_parcloud][1] - HaloMerger_ParCloud_CenCoord[index2_parcloud][1] )
                        + SQR( HaloMerger_ParCloud_CenCoord[index_parcloud][2] - HaloMerger_ParCloud_CenCoord[index2_parcloud][2] ) )
                    > ( HaloMerger_ParCloud_DensProf_MaxR[index_parcloud] + HaloMerger_ParCloud_DensProf_MaxR[index2_parcloud] ) )
                  isOverlap = false;

               if ( MPI_Rank == 0 )
               {
                  if ( isOverlap )
                     Aux_Message( stderr, "WARNING : ParCloud_%d range (center: [%13.6e, %13.6e, %13.6e], MaxR: %13.6e) overlaps with the ParCloud_%d range (center: [%13.6e, %13.6e, %13.6e], MaxR: %13.6e) !!\n",
                                  index_parcloud_input,
                                  HaloMerger_ParCloud_CenCoord[index_parcloud][0], HaloMerger_ParCloud_CenCoord[index_parcloud][1],
                                  HaloMerger_ParCloud_CenCoord[index_parcloud][2], HaloMerger_ParCloud_DensProf_MaxR[index_parcloud],
                                  index2_parcloud_input,
                                  HaloMerger_ParCloud_CenCoord[index2_parcloud][0], HaloMerger_ParCloud_CenCoord[index2_parcloud][1],
                                  HaloMerger_ParCloud_CenCoord[index2_parcloud][2], HaloMerger_ParCloud_DensProf_MaxR[index2_parcloud] );
               }

            } // for (int index2_parcloud=0; index2_parcloud<index_parcloud; index2_parcloud++)

            // check whether the input particle cloud overlaps with the input halos
            for (int index_halo=0; index_halo<HaloMerger_Halo_Num; index_halo++)
            {
               const int index_halo_input = index_halo+1; // index of halo in the input file

               if ( HaloMerger_Halo_InitMode == 1 )
               {

                  bool isOverlap = true;

                  for (int d=0; d<3; d++)
                  {
                     if ( HaloMerger_ParCloud_CenCoord[index_parcloud][d] - HaloMerger_ParCloud_DensProf_MaxR[index_parcloud] >= HaloMerger_Halo_UM_IC_Range_EdgeR[index_halo][d] ||
                          HaloMerger_ParCloud_CenCoord[index_parcloud][d] + HaloMerger_ParCloud_DensProf_MaxR[index_parcloud] <= HaloMerger_Halo_UM_IC_Range_EdgeL[index_halo][d] )
                     {
                        isOverlap = false;
                        break;
                     }
                  } // for (int d=0; d<3; d++)

                  if ( MPI_Rank == 0 )
                  {
                     if ( isOverlap )
                        Aux_Message( stderr, "WARNING : ParCloud_%d (center: [%13.6e, %13.6e, %13.6e], MaxR: %13.6e) overlaps with the Halo_%d UM_IC range !!\n",
                                     index_parcloud_input,
                                     HaloMerger_ParCloud_CenCoord[index_parcloud][0], HaloMerger_ParCloud_CenCoord[index_parcloud][1],
                                     HaloMerger_ParCloud_CenCoord[index_parcloud][2], HaloMerger_ParCloud_DensProf_MaxR[index_parcloud],
                                     index_halo_input );
                  }

               } // if ( HaloMerger_Halo_InitMode == 1 )
               else
                  Aux_Error( ERROR_INFO, "unsupported initialization mode (%s = %d) !!\n",
                             "HaloMerger_Halo_InitMode", HaloMerger_Halo_InitMode );

            } // for (int index_halo=0; index_halo<HaloMerger_Halo_Num; index_halo++)

            // check whether the input particle cloud overlaps with the input solitons
            for (int index_soliton=0; index_soliton<HaloMerger_Soliton_Num; index_soliton++)
            {
               const int index_soliton_input = index_soliton+1; // index of soliton in the input file

               bool isOverlap = true;

               if ( sqrt( SQR( HaloMerger_ParCloud_CenCoord[index_parcloud][0] - HaloMerger_Soliton_CenCoord[index_soliton][0] )
                        + SQR( HaloMerger_ParCloud_CenCoord[index_parcloud][1] - HaloMerger_Soliton_CenCoord[index_soliton][1] )
                        + SQR( HaloMerger_ParCloud_CenCoord[index_parcloud][2] - HaloMerger_Soliton_CenCoord[index_soliton][2] ) )
                    > ( HaloMerger_ParCloud_DensProf_MaxR[index_parcloud] + 3.0*HaloMerger_Soliton_CoreRadius[index_soliton] ) )
                  isOverlap = false;

               if ( MPI_Rank == 0 )
               {
                  if ( isOverlap )
                     Aux_Message( stderr, "WARNING : ParCloud_%d (center: [%13.6e, %13.6e, %13.6e], MaxR: %13.6e) overlaps with the Soliton_%d 3r_c-range (center: [%13.6e, %13.6e, %13.6e], r_c: %13.6e) !!\n",
                                  index_parcloud_input,
                                  HaloMerger_ParCloud_CenCoord[index_parcloud][0], HaloMerger_ParCloud_CenCoord[index_parcloud][1],
                                  HaloMerger_ParCloud_CenCoord[index_parcloud][2], HaloMerger_ParCloud_DensProf_MaxR[index_parcloud],
                                  index_soliton_input,
                                  HaloMerger_Soliton_CenCoord[index_soliton][0], HaloMerger_Soliton_CenCoord[index_soliton][1],
                                  HaloMerger_Soliton_CenCoord[index_soliton][2], HaloMerger_Soliton_CoreRadius[index_soliton] );
               }

            } // for (int index_soliton=0; index_soliton<HaloMerger_Soliton_Num; index_soliton++)

         } // if ( HaloMerger_ParCloud_InitMode == 1 )
         else
            Aux_Error( ERROR_INFO, "unsupported initialization mode (%s = %d) !!\n",
                       "HaloMerger_ParCloud_InitMode", HaloMerger_ParCloud_InitMode );

      } // for (int index_parcloud=0; index_parcloud<HaloMerger_ParCloud_Num; index_parcloud++)

      // check there are particles
      if ( HaloMerger_ParCloud_NPar_Total <= 0 )
         Aux_Error( ERROR_INFO, "Total number of particles (sum of HaloMerger_ParCloud_i_NPar) must be >0 for HaloMerger_ParCloud_Num > 0 !!\n" );

   } // if ( OPT__INIT != INIT_BY_RESTART  &&  HaloMerger_ParCloud_Num > 0 )

   // overwrite the total number of particles
   if ( OPT__INIT != INIT_BY_RESTART )
   {
#     ifdef MASSIVE_PARTICLES
      amr->Par->NPar_Active_AllRank = HaloMerger_ParCloud_NPar_Total;

      PRINT_RESET_PARA( amr->Par->NPar_Active_AllRank, FORMAT_LONG,
                        "(amr->Par->NPar_Active_AllRank is originally set by PAR_NPAR in Input__Parameter)" );
#     endif
   }


// (3) reset other general-purpose parameters
//     --> a helper macro PRINT_RESET_PARA is defined in Macro.h
   const long   End_Step_Default = __INT_MAX__;
   const double End_T_Default    = 0.25;

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
      Aux_Message( stdout, "  test problem ID                    = %d\n",    TESTPROB_ID                     );
      Aux_Message( stdout, "  total number of halos              = %d\n",    HaloMerger_Halo_Num             );
      Aux_Message( stdout, "  halo initialization mode           = %d\n",    HaloMerger_Halo_InitMode        );
      Aux_Message( stdout, "  total number of solitons           = %d\n",    HaloMerger_Soliton_Num          );
      Aux_Message( stdout, "  soliton initialization mode        = %d\n",    HaloMerger_Soliton_InitMode     );
      Aux_Message( stdout, "  total number of particle clouds    = %d\n",    HaloMerger_ParCloud_Num         );
      Aux_Message( stdout, "  particle cloud initialization mode = %d\n",    HaloMerger_ParCloud_InitMode    );
      if ( OPT__EXT_POT == EXT_POT_FUNC )
      {
         Aux_Message( stdout, "  external potential of uniform-density sphere information:\n" );
         Aux_Message( stdout, "           %14s  %14s  %14s  %14s  %14s  %14s  %14s  %14s  %14s\n",
                      "Mass", "Radius", "Rho", "CenCoord_X", "CenCoord_Y", "CenCoord_Z", "Velocity_X", "Velocity_Y", "Velocity_Z"  );
         Aux_Message( stdout, "           %14.6e  %14.6e  %14.6e  %14.6e  %14.6e  %14.6e  %14.6e  %14.6e  %14.6e\n",
                      HaloMerger_ExtPot_UniDenSph_M, HaloMerger_ExtPot_UniDenSph_R, HaloMerger_ExtPot_UniDenSph_Rho,
                      HaloMerger_ExtPot_UniDenSph_CenCoordX, HaloMerger_ExtPot_UniDenSph_CenCoordY, HaloMerger_ExtPot_UniDenSph_CenCoordZ,
                      HaloMerger_ExtPot_UniDenSph_VelocityX, HaloMerger_ExtPot_UniDenSph_VelocityY, HaloMerger_ExtPot_UniDenSph_VelocityZ );
      }

      if ( OPT__INIT != INIT_BY_RESTART  &&  HaloMerger_Halo_Num > 0 )
      {
         Aux_Message( stdout, "\n  halo information:\n" );
         Aux_Message( stdout, "  %7s  %14s  %14s  %14s  %14s  %14s  %14s\n",
                      "ID", "CenCoord_X", "CenCoord_Y", "CenCoord_Z", "Velocity_X", "Velocity_Y", "Velocity_Z" );

         for (int index_halo=0; index_halo<HaloMerger_Halo_Num; index_halo++)
            Aux_Message( stdout, "  %7d  %14.6e  %14.6e  %14.6e  %14.6e  %14.6e  %14.6e\n",
                         index_halo+1,
                         HaloMerger_Halo_CenCoord[index_halo][0], HaloMerger_Halo_CenCoord[index_halo][1], HaloMerger_Halo_CenCoord[index_halo][2],
                         HaloMerger_Halo_Velocity[index_halo][0], HaloMerger_Halo_Velocity[index_halo][1], HaloMerger_Halo_Velocity[index_halo][2] );

         if ( HaloMerger_Halo_InitMode == 1 )
         {
            Aux_Message( stdout, "\n  halo UM_IC information:\n" );
            Aux_Message( stdout, "  %7s  %14s  %14s  %14s  %14s  %14s  %14s  %14s  %14s  %14s  %14s  %14s  %20s  %20s  %20s  %20s  %20s  %20s\n",
                         "ID", "UM_IC_Filename", "UM_IC_Float8",
                               "UM_IC_BoxLen_X", "UM_IC_BoxLen_Y", "UM_IC_BoxLen_Z",
                               "UM_IC_NCells_X", "UM_IC_NCells_Y", "UM_IC_NCells_Z",
                               "UM_IC_dh_X", "UM_IC_dh_Y", "UM_IC_dh_Z",
                               "UM_IC_Range_EdgeL_X", "UM_IC_Range_EdgeL_Y", "UM_IC_Range_EdgeL_Z",
                               "UM_IC_Range_EdgeR_X", "UM_IC_Range_EdgeR_Y", "UM_IC_Range_EdgeR_Z" );

            for (int index_halo=0; index_halo<HaloMerger_Halo_Num; index_halo++)
               Aux_Message( stdout, "  %7d  %14s  %14d  %14.6e  %14.6e  %14.6e  %14d  %14d  %14d  %14.6e  %14.6e  %14.6e  %20.6e  %20.6e  %20.6e  %20.6e  %20.6e  %20.6e\n",
                            index_halo+1, HaloMerger_Halo_UM_IC_Filename[index_halo], HaloMerger_Halo_UM_IC_Float8[index_halo],
                            HaloMerger_Halo_UM_IC_BoxLen[index_halo][0], HaloMerger_Halo_UM_IC_BoxLen[index_halo][1], HaloMerger_Halo_UM_IC_BoxLen[index_halo][2],
                            HaloMerger_Halo_UM_IC_NCells[index_halo][0], HaloMerger_Halo_UM_IC_NCells[index_halo][1], HaloMerger_Halo_UM_IC_NCells[index_halo][2],
                            HaloMerger_Halo_UM_IC_dh[index_halo][0], HaloMerger_Halo_UM_IC_dh[index_halo][1], HaloMerger_Halo_UM_IC_dh[index_halo][2],
                            HaloMerger_Halo_UM_IC_Range_EdgeL[index_halo][0], HaloMerger_Halo_UM_IC_Range_EdgeL[index_halo][1], HaloMerger_Halo_UM_IC_Range_EdgeL[index_halo][2],
                            HaloMerger_Halo_UM_IC_Range_EdgeR[index_halo][0], HaloMerger_Halo_UM_IC_Range_EdgeR[index_halo][1], HaloMerger_Halo_UM_IC_Range_EdgeR[index_halo][2]);

         } // if ( HaloMerger_Halo_InitMode == 1 )
         else
            Aux_Error( ERROR_INFO, "unsupported initialization mode (%s = %d) !!\n",
                       "HaloMerger_Halo_InitMode", HaloMerger_Halo_InitMode );

      } // if ( OPT__INIT != INIT_BY_RESTART  &&  HaloMerger_Halo_Num > 0 )

      if ( OPT__INIT != INIT_BY_RESTART  &&  HaloMerger_Soliton_Num > 0 )
      {
         Aux_Message( stdout, "\n  soliton information:\n" );
         Aux_Message( stdout, "  %7s  %14s  %14s  %14s  %14s  %14s  %14s  %14s  %14s\n",
                      "ID", "CoreRadius", "CoreRho", "CenCoord_X", "CenCoord_Y", "CenCoord_Z", "Velocity_X", "Velocity_Y", "Velocity_Z" );

         for (int index_soliton=0; index_soliton<HaloMerger_Soliton_Num; index_soliton++)
            Aux_Message( stdout, "  %7d  %14.6e  %14.6e  %14.6e  %14.6e  %14.6e  %14.6e  %14.6e  %14.6e\n",
                         index_soliton+1,
                         HaloMerger_Soliton_CoreRadius[index_soliton], HaloMerger_Soliton_CoreRho[index_soliton],
                         HaloMerger_Soliton_CenCoord[index_soliton][0], HaloMerger_Soliton_CenCoord[index_soliton][1], HaloMerger_Soliton_CenCoord[index_soliton][2],
                         HaloMerger_Soliton_Velocity[index_soliton][0], HaloMerger_Soliton_Velocity[index_soliton][1], HaloMerger_Soliton_Velocity[index_soliton][2] );

         if ( HaloMerger_Soliton_InitMode == 1 )
         {
            Aux_Message( stdout, "\n  soliton density profile information:\n" );
            Aux_Message( stdout, "  %7s %31s  %14s  %14s  %14s  %14s\n",
                         "ID", "DensProf_Filename", "NBin",
                               "Rescale", "ScaleL", "ScaleD"  );

            for (int index_soliton=0; index_soliton<HaloMerger_Soliton_Num; index_soliton++)
               Aux_Message( stdout, "  %7d %31s  %14d  %14d  %14.6e  %14.6e\n",
                            index_soliton+1, HaloMerger_Soliton_DensProf_Filename[index_soliton], HaloMerger_Soliton_DensProf_NBin[index_soliton],
                            HaloMerger_Soliton_DensProf_Rescale[index_soliton], HaloMerger_Soliton_DensProf_ScaleL[index_soliton], HaloMerger_Soliton_DensProf_ScaleD[index_soliton] );

         } // if ( HaloMerger_Soliton_InitMode == 1 )
         else if ( HaloMerger_Soliton_InitMode == 2 )
         {
            Aux_Message( stdout, "\n  soliton density profile information:\n" );
            Aux_Message( stdout, "  %7s  %14s\n",
                         "ID", "OuterSlope" );

            for (int index_soliton=0; index_soliton<HaloMerger_Soliton_Num; index_soliton++)
               Aux_Message( stdout, "  %7d  %14.6e\n",
                            index_soliton+1, HaloMerger_Soliton_OuterSlope[index_soliton] );

         } // else if ( HaloMerger_Soliton_InitMode == 2 )
         else
            Aux_Error( ERROR_INFO, "unsupported initialization mode (%s = %d) !!\n",
                       "HaloMerger_Soliton_InitMode", HaloMerger_Soliton_InitMode );

      } // if ( OPT__INIT != INIT_BY_RESTART  &&  HaloMerger_Soliton_Num > 0 )

      if ( OPT__INIT != INIT_BY_RESTART  &&  HaloMerger_ParCloud_Num > 0 )
      {
         Aux_Message( stdout, "\n  particle cloud information:\n" );
         Aux_Message( stdout, "  %7s  %14s  %14s  %14s  %14s  %14s  %14s\n",
                      "ID", "CenCoord_X", "CenCoord_Y", "CenCoord_Z", "Velocity_X", "Velocity_Y", "Velocity_Z" );

         for (int index_parcloud=0; index_parcloud<HaloMerger_ParCloud_Num; index_parcloud++)
            Aux_Message( stdout, "  %7d  %14.6e  %14.6e  %14.6e  %14.6e  %14.6e  %14.6e\n",
                         index_parcloud+1,
                         HaloMerger_ParCloud_CenCoord[index_parcloud][0], HaloMerger_ParCloud_CenCoord[index_parcloud][1], HaloMerger_ParCloud_CenCoord[index_parcloud][2],
                         HaloMerger_ParCloud_Velocity[index_parcloud][0], HaloMerger_ParCloud_Velocity[index_parcloud][1], HaloMerger_ParCloud_Velocity[index_parcloud][2] );

         if ( HaloMerger_ParCloud_InitMode == 1 )
         {
            Aux_Message( stdout, "\n  particle cloud density profile information:\n" );
            Aux_Message( stdout, "  %7s  %30s  %14s  %14s  %14s\n",
                         "ID", "DensProf_Filename", "DensProf_MaxR", "RSeed", "NPar" );

            for (int index_parcloud=0; index_parcloud<HaloMerger_ParCloud_Num; index_parcloud++)
               Aux_Message( stdout, "  %7d  %30s  %14.6e  %14d  %14ld\n",
                            index_parcloud+1, HaloMerger_ParCloud_DensProf_Filename[index_parcloud], HaloMerger_ParCloud_DensProf_MaxR[index_parcloud],
                            HaloMerger_ParCloud_RSeed[index_parcloud], HaloMerger_ParCloud_NPar[index_parcloud] );

            Aux_Message( stdout, "  %18s-> Total number of particles in all particle clouds =  %14ld\n", "", HaloMerger_ParCloud_NPar_Total );

         } // if ( HaloMerger_ParCloud_InitMode == 1 )
         else
            Aux_Error( ERROR_INFO, "unsupported initialization mode (%s = %d) !!\n",
                       "HaloMerger_ParCloud_InitMode", HaloMerger_ParCloud_InitMode );

      } // if ( OPT__INIT != INIT_BY_RESTART  &&  HaloMerger_ParCloud_Num > 0 )

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
   // (1) set the background
   double Real = 0.0;
   double Imag = 0.0;


   // (2) set the halos
   for (int index_halo=0; index_halo<HaloMerger_Halo_Num; index_halo++)
   {
      // wavefunction of each halo
      double Real_halo = 0.0;
      double Imag_halo = 0.0;

      // get the wave function of the halo
      switch ( HaloMerger_Halo_InitMode )
      {
         // read the value interpolated from the UM_IC data
         case 1:
         {
            // when the (x,y,z) is inside the range of this halo
            if ( x >= HaloMerger_Halo_UM_IC_Range_EdgeL[index_halo][0]  &&
                 x <= HaloMerger_Halo_UM_IC_Range_EdgeR[index_halo][0]  &&
                 y >= HaloMerger_Halo_UM_IC_Range_EdgeL[index_halo][1]  &&
                 y <= HaloMerger_Halo_UM_IC_Range_EdgeR[index_halo][1]  &&
                 z >= HaloMerger_Halo_UM_IC_Range_EdgeL[index_halo][2]  &&
                 z <= HaloMerger_Halo_UM_IC_Range_EdgeR[index_halo][2] )
            {
               // assume real part first and then imaginary part
               Real_halo = HaloMerger_Get_Value_From_Halo_UM_IC_Data( x, y, z, 0, index_halo );
               Imag_halo = HaloMerger_Get_Value_From_Halo_UM_IC_Data( x, y, z, 1, index_halo );
            }

            break;
         } // case 1

         default:
            Aux_Error( ERROR_INFO, "unsupported initialization mode (%s = %d) !!\n",
                       "HaloMerger_Halo_InitMode", HaloMerger_Halo_InitMode );

      } // switch ( HaloMerger_Halo_InitMode )

      // add the velocity
      HaloMerger_Add_Velocity( &Real_halo, &Imag_halo,
                               HaloMerger_Halo_Velocity[index_halo][0],
                               HaloMerger_Halo_Velocity[index_halo][1],
                               HaloMerger_Halo_Velocity[index_halo][2],
                               x, y, z );

      // add the wavefunction to the box
      Real += Real_halo;
      Imag += Imag_halo;

   } // for (int index_halo=0; index_halo<HaloMerger_Halo_Num; index_halo++)


   // (3) set the solitons
   for (int index_soliton=0; index_soliton<HaloMerger_Soliton_Num; index_soliton++)
   {
      // density and wave function of each soliton
      double Dens_soliton = 0.0;
      double Real_soliton = 0.0;
      double Imag_soliton = 0.0;

      // get the density and wave function of the soliton
      switch ( HaloMerger_Soliton_InitMode )
      {
         // read from the density profile table
         case 1:
         {
            // density profile
            const double *Table_Radius  = HaloMerger_Soliton_DensProf[index_soliton] + 0*HaloMerger_Soliton_DensProf_NBin[index_soliton];  // radius
            const double *Table_Density = HaloMerger_Soliton_DensProf[index_soliton] + 1*HaloMerger_Soliton_DensProf_NBin[index_soliton];  // density

            // target radius
            const double r_tar = sqrt( SQR(x - HaloMerger_Soliton_CenCoord[index_soliton][0]) +
                                       SQR(y - HaloMerger_Soliton_CenCoord[index_soliton][1]) +
                                       SQR(z - HaloMerger_Soliton_CenCoord[index_soliton][2]) );

            // rescale radius (target radius --> reference radius)
            const double r_ref = r_tar / HaloMerger_Soliton_DensProf_ScaleL[index_soliton];

            // linear interpolation
            double dens_ref = Mis_InterpolateFromTable( HaloMerger_Soliton_DensProf_NBin[index_soliton], Table_Radius, Table_Density, r_ref );

            if ( dens_ref == NULL_REAL )
            {
               if      ( r_ref <  Table_Radius[0] )
                  dens_ref = Table_Density[0];
               else if ( r_ref >= Table_Radius[HaloMerger_Soliton_DensProf_NBin[index_soliton]-1] )
                  dens_ref = 0.0;
               else
                  Aux_Error( ERROR_INFO, "interpolation failed at radius %13.7e (min/max radius = %13.7e/%13.7e) !!\n",
                             r_ref, Table_Radius[0], Table_Radius[HaloMerger_Soliton_DensProf_NBin[index_soliton]-1] );
            }

            // get the density of soliton
            Dens_soliton = dens_ref*HaloMerger_Soliton_DensProf_ScaleD[index_soliton];
            Real_soliton = sqrt( Dens_soliton );
            Imag_soliton = 0.0;

            break;
         } // case 1

         // set from the analytical soliton density profile
         case 2:
         {
            // parameters for the analytical soliton density profile
            const double m22      = ELBDM_MASS*UNIT_M/(Const_eV/SQR(Const_c))/1.0e-22;
            const double rc_kpc   = HaloMerger_Soliton_CoreRadius[index_soliton]*UNIT_L/Const_kpc;
            const double peak_rho = 1.945e7/SQR( m22*rc_kpc*rc_kpc )*Const_Msun/CUBE(Const_kpc)/(UNIT_M/CUBE(UNIT_L));

            // target radius
            const double r_tar = sqrt( SQR(x - HaloMerger_Soliton_CenCoord[index_soliton][0]) +
                                       SQR(y - HaloMerger_Soliton_CenCoord[index_soliton][1]) +
                                       SQR(z - HaloMerger_Soliton_CenCoord[index_soliton][2]) );

            // get the density of soliton
            Dens_soliton  = peak_rho*pow( 1.0+9.06e-2*SQR(r_tar/HaloMerger_Soliton_CoreRadius[index_soliton]),
                                          HaloMerger_Soliton_OuterSlope[index_soliton] );
            Real_soliton  = sqrt( Dens_soliton );
            Imag_soliton  = 0.0;

            break;
         } // case 2

         default:
            Aux_Error( ERROR_INFO, "unsupported initialization mode (%s = %d) !!\n",
                       "HaloMerger_Soliton_InitMode", HaloMerger_Soliton_InitMode );

      } // switch ( HaloMerger_Soliton_InitMode )

      // add the velocity
      HaloMerger_Add_Velocity( &Real_soliton, &Imag_soliton,
                               HaloMerger_Soliton_Velocity[index_soliton][0],
                               HaloMerger_Soliton_Velocity[index_soliton][1],
                               HaloMerger_Soliton_Velocity[index_soliton][2],
                               x, y, z );

      // add the wavefunction to the box
      Real += Real_soliton;
      Imag += Imag_soliton;

   } // for (int index_soliton=0; index_soliton<HaloMerger_Soliton_Num; index_soliton++)



#  if ( ELBDM_SCHEME == ELBDM_HYBRID )
   if ( amr->use_wave_flag[lv] ) {
#  endif
   fluid[REAL] = Real;
   fluid[IMAG] = Imag;
   fluid[DENS] = SQR(fluid[REAL]) + SQR(fluid[IMAG]);
#  if ( ELBDM_SCHEME == ELBDM_HYBRID )
   } else { // if ( amr->use_wave_flag[lv] )
   fluid[DENS] = SQR(Real) + SQR(Imag);
   fluid[PHAS] = SATAN2( Imag, Real );
   fluid[STUB] = 0.0;
   } // if ( amr->use_wave_flag[lv] ) ... else
#  endif

} // FUNCTION : SetGridIC



//-------------------------------------------------------------------------------------------------------
// Function    :  End_HaloMerger
// Description :  Free memory before terminating the program
//
// Note        :  1. Linked to the function pointer "End_User_Ptr" to replace "End_User()"
//
// Parameter   :  None
//-------------------------------------------------------------------------------------------------------
void End_HaloMerger()
{

   // Halos
   if ( OPT__INIT != INIT_BY_RESTART  &&  HaloMerger_Halo_Num > 0 )
   {
      delete [] HaloMerger_Halo_CenCoord;
      delete [] HaloMerger_Halo_Velocity;

      if ( HaloMerger_Halo_InitMode == 1 )
      {
         for (int index_halo=0; index_halo<HaloMerger_Halo_Num; index_halo++)
         {
            delete [] HaloMerger_Halo_UM_IC_Data[index_halo];
         }

         delete [] HaloMerger_Halo_UM_IC_Filename;
         delete [] HaloMerger_Halo_UM_IC_BoxLen;
         delete [] HaloMerger_Halo_UM_IC_NCells;
         delete [] HaloMerger_Halo_UM_IC_Float8;
         delete [] HaloMerger_Halo_UM_IC_dh;
         delete [] HaloMerger_Halo_UM_IC_Range_EdgeL;
         delete [] HaloMerger_Halo_UM_IC_Range_EdgeR;
         delete [] HaloMerger_Halo_UM_IC_Data;

      } // if ( HaloMerger_Halo_InitMode == 1 )
      else
         Aux_Error( ERROR_INFO, "unsupported initialization mode (%s = %d) !!\n",
                    "HaloMerger_Halo_InitMode", HaloMerger_Halo_InitMode );

   } // if ( OPT__INIT != INIT_BY_RESTART  &&  HaloMerger_Halo_Num > 0 )

   // Solitons
   if ( OPT__INIT != INIT_BY_RESTART  &&  HaloMerger_Soliton_Num > 0 )
   {
      delete [] HaloMerger_Soliton_CenCoord;
      delete [] HaloMerger_Soliton_Velocity;
      delete [] HaloMerger_Soliton_CoreRadius;
      delete [] HaloMerger_Soliton_CoreRho;

      if ( HaloMerger_Soliton_InitMode == 1 )
      {
         for (int index_soliton=0; index_soliton<HaloMerger_Soliton_Num; index_soliton++)
         {
            delete [] HaloMerger_Soliton_DensProf[index_soliton];
         }

         delete [] HaloMerger_Soliton_DensProf_Filename;
         delete [] HaloMerger_Soliton_DensProf;
         delete [] HaloMerger_Soliton_DensProf_NBin;
         delete [] HaloMerger_Soliton_DensProf_Rescale;
         delete [] HaloMerger_Soliton_DensProf_ScaleL;
         delete [] HaloMerger_Soliton_DensProf_ScaleD;

      } // if ( HaloMerger_Soliton_InitMode == 1 )
      else if ( HaloMerger_Soliton_InitMode == 2 )
      {
         delete [] HaloMerger_Soliton_OuterSlope;

      } // else if ( HaloMerger_Soliton_InitMode == 2 )
      else
         Aux_Error( ERROR_INFO, "unsupported initialization mode (%s = %d) !!\n",
                    "HaloMerger_Soliton_InitMode", HaloMerger_Soliton_InitMode );

   } // if ( OPT__INIT != INIT_BY_RESTART  &&  HaloMerger_Soliton_Num > 0 )

   // ParClouds
   if ( OPT__INIT != INIT_BY_RESTART  &&  HaloMerger_ParCloud_Num > 0 )
   {
      delete [] HaloMerger_ParCloud_CenCoord;
      delete [] HaloMerger_ParCloud_Velocity;

      if ( HaloMerger_ParCloud_InitMode == 1 )
      {
         delete [] HaloMerger_ParCloud_DensProf_Filename;
         delete [] HaloMerger_ParCloud_DensProf_MaxR;
         delete [] HaloMerger_ParCloud_RSeed;
         delete [] HaloMerger_ParCloud_NPar;
      } // if ( HaloMerger_ParCloud_InitMode == 1 )
      else
         Aux_Error( ERROR_INFO, "unsupported initialization mode (%s = %d) !!\n",
                    "HaloMerger_ParCloud_InitMode", HaloMerger_ParCloud_InitMode );

   } // if ( OPT__INIT != INIT_BY_RESTART  &&  HaloMerger_ParCloud_Num > 0 )

} // FUNCTION : End_HaloMerger



//-------------------------------------------------------------------------------------------------------
// Function    :  HaloMerger_Add_Velocity
// Description :  Multiply the wave function by a plane wave wave function with a given velocity
//
// Note        :  None
//
// Parameter   :  RealPart    : real part of the wavefunction
//                ImagPart    : imaginary part of the wavefunction
//                Velocity_X  : velocity in the x-direction
//                Velocity_Y  : velocity in the y-direction
//                Velocity_Z  : velocity in the z-direction
//                Position_X  : coordinate x of the point
//                Position_Y  : coordinate y of the point
//                Position_Z  : coordinate z of the point
// Return      :  RealPart
//                ImagPart
//-------------------------------------------------------------------------------------------------------
void HaloMerger_Add_Velocity( double *RealPart, double *ImagPart,
                              const double Velocity_X, const double Velocity_Y, const double Velocity_Z,
                              const double Position_X, const double Position_Y, const double Position_Z )
{
   // before adding the velocity
   const double Real_Old = *RealPart;
   const double Imag_Old = *ImagPart;

   // Phase = kx = (m*v/hbar)*x = eta*v*x
   const double Phase = ELBDM_ETA*( Velocity_X*Position_X + Velocity_Y*Position_Y + Velocity_Z*Position_Z );

   // psi_new = psi_old * exp(i*Phase) = (R_old + i*I_old)*(cos(Phase) + i*sin(Phase))
   const double Real_New = Real_Old * cos(Phase) - Imag_Old * sin(Phase);
   const double Imag_New = Real_Old * sin(Phase) + Imag_Old * cos(Phase);

   // return the updated wave function
   *RealPart = Real_New;
   *ImagPart = Imag_New;

} // FUNCTION : HaloMerger_Add_Velocity



//-------------------------------------------------------------------------------------------------------
// Function    :  HaloMerger_Trilinear_Interpolation
// Description :  Apply linear interpolation to get the desired value from the data at the eight corners of a 3D cube
//
// Note        :  1. Ref_Value is in the order zyx (i.e. Ref_Value[z][y][x])
//
// Parameter   :  Target_X   : target x-coordinate
//                Target_Y   : target y-coordinate
//                Target_Z   : target z-coordinate
//                Ref_Value  : reference values at the eight corners
//                Ref_X      : reference x-coordinates at the eight corners
//                Ref_Y      : reference y-coordinates at the eight corners
//                Ref_Z      : reference z-coordinates at the eight corners
// Return      :  Value_ZYX  : values interpolated from the eight corners
//-------------------------------------------------------------------------------------------------------
double HaloMerger_Trilinear_Interpolation( const double Target_X, const double Target_Y, const double Target_Z,
                                           const double Ref_Value[2][2][2],
                                           const double Ref_X[2], const double Ref_Y[2], const double Ref_Z[2] )
{
    // linear interpolation in z-direction
    const double Value_Z00 = Ref_Value[0][0][0] + ( Ref_Value[1][0][0] - Ref_Value[0][0][0])/(Ref_Z[1] - Ref_Z[0])*(Target_Z - Ref_Z[0]);
    const double Value_Z01 = Ref_Value[0][0][1] + ( Ref_Value[1][0][1] - Ref_Value[0][0][1])/(Ref_Z[1] - Ref_Z[0])*(Target_Z - Ref_Z[0]);
    const double Value_Z10 = Ref_Value[0][1][0] + ( Ref_Value[1][1][0] - Ref_Value[0][1][0])/(Ref_Z[1] - Ref_Z[0])*(Target_Z - Ref_Z[0]);
    const double Value_Z11 = Ref_Value[0][1][1] + ( Ref_Value[1][1][1] - Ref_Value[0][1][1])/(Ref_Z[1] - Ref_Z[0])*(Target_Z - Ref_Z[0]);

    // linear interpolation in y-direction
    const double Value_ZY0 = Value_Z00          + ( Value_Z10          - Value_Z00         )/(Ref_Y[1] - Ref_Y[0])*(Target_Y - Ref_Y[0]);
    const double Value_ZY1 = Value_Z01          + ( Value_Z11          - Value_Z01         )/(Ref_Y[1] - Ref_Y[0])*(Target_Y - Ref_Y[0]);

    // linear interpolation in x-direction
    const double Value_ZYX = Value_ZY0          + ( Value_ZY1          - Value_ZY0         )/(Ref_X[1] - Ref_X[0])*(Target_X - Ref_X[0]);

    return Value_ZYX;

} // FUNCTION : HaloMerger_Trilinear_Interpolation



//-------------------------------------------------------------------------------------------------------
// Function    :  HaloMerger_Get_Value_From_Halo_UM_IC_Data
// Description :  Get the target value of field v at (x,y,z) by reading the UM_IC of halo and performing linear interpolation
//
// Note        :  None
//
// Parameter   :  x                  : target x-coordinate
//                y                  : target y-coordinate
//                z                  : target z-coordinate
//                v                  : index of the target field
//                index_halo         : index of target halo
// Return      :  Interpolated_Value : value interpolated from the UM_IC of halo
//-------------------------------------------------------------------------------------------------------
double HaloMerger_Get_Value_From_Halo_UM_IC_Data( const double x, const double y, const double z, const int v, const int index_halo )
{
    // 1. UM_IC information
    const int UM_IC_Nx          = HaloMerger_Halo_UM_IC_NCells[index_halo][0];
    const int UM_IC_Ny          = HaloMerger_Halo_UM_IC_NCells[index_halo][1];
    const int UM_IC_Nz          = HaloMerger_Halo_UM_IC_NCells[index_halo][2];
    const int UM_IC_NVar        = 2; // (Real part & Imag part)
    const int UM_IC_Float8      = HaloMerger_Halo_UM_IC_Float8[index_halo];
    size_t load_data_size       = ( UM_IC_Float8 ) ? sizeof(double) : sizeof(float);

    double UM_IC_Range_EdgeL[3];
    double UM_IC_dh[3];
    for (int d=0; d<3; d++) UM_IC_Range_EdgeL[d] = HaloMerger_Halo_UM_IC_Range_EdgeL[index_halo][d];
    for (int d=0; d<3; d++) UM_IC_dh[d]          = HaloMerger_Halo_UM_IC_dh[index_halo][d];

    // Use the eight corners from the UM_IC to do the linear interpolation

    // 2. index of the bottom-left corner in the UM_IC
    const int IntCorner000_UMICIndex[3] = {(int)floor( (x - UM_IC_Range_EdgeL[0])/UM_IC_dh[0] - 0.5 ),
                                           (int)floor( (y - UM_IC_Range_EdgeL[1])/UM_IC_dh[1] - 0.5 ),
                                           (int)floor( (z - UM_IC_Range_EdgeL[2])/UM_IC_dh[2] - 0.5 )};

    // 3. physical coordinates of the eight corners
    double IntCorner_Coord[3][2];
    for (int d=0; d<3; d++)
    for (int IntCorner_ID=0; IntCorner_ID<2; IntCorner_ID++)
       IntCorner_Coord[d][IntCorner_ID] = UM_IC_Range_EdgeL[d] + (IntCorner000_UMICIndex[d] + IntCorner_ID + 0.5)*UM_IC_dh[d];

    // 4. values at the eight corners
    double IntCorner_Value[2][2][2];
    for (int IntCorner_ID_k=0; IntCorner_ID_k<2; IntCorner_ID_k++){ const int IntCorner_UMICIndex_k = IntCorner000_UMICIndex[2] + IntCorner_ID_k;
    for (int IntCorner_ID_j=0; IntCorner_ID_j<2; IntCorner_ID_j++){ const int IntCorner_UMICIndex_j = IntCorner000_UMICIndex[1] + IntCorner_ID_j;
    for (int IntCorner_ID_i=0; IntCorner_ID_i<2; IntCorner_ID_i++){ const int IntCorner_UMICIndex_i = IntCorner000_UMICIndex[0] + IntCorner_ID_i;

       if ( IntCorner_UMICIndex_i < 0  ||  IntCorner_UMICIndex_i >= UM_IC_Nx  ||
            IntCorner_UMICIndex_j < 0  ||  IntCorner_UMICIndex_j >= UM_IC_Ny  ||
            IntCorner_UMICIndex_k < 0  ||  IntCorner_UMICIndex_k >= UM_IC_Nz )
       {
          // set the value as zero when the corner is outside the UM_IC file
          IntCorner_Value[IntCorner_ID_k][IntCorner_ID_j][IntCorner_ID_i] = 0.0;
       }
       else
       {
          const long IdxIn_UM_IC = (long)v*UM_IC_Nz*UM_IC_Ny*UM_IC_Nx
                                 + (long)IntCorner_UMICIndex_k*UM_IC_Ny*UM_IC_Nx
                                 + (long)IntCorner_UMICIndex_j*UM_IC_Nx
                                 + (long)IntCorner_UMICIndex_i;

          IntCorner_Value[IntCorner_ID_k][IntCorner_ID_j][IntCorner_ID_i] = ( UM_IC_Float8 ) ?
             (double)(*((double*)&HaloMerger_Halo_UM_IC_Data[index_halo][IdxIn_UM_IC*load_data_size])):
             (double)(*( (float*)&HaloMerger_Halo_UM_IC_Data[index_halo][IdxIn_UM_IC*load_data_size]));
       }

    }}} // for (int IntCorner_ID_kji=0; IntCorner_ID_kji<2; IntCorner_ID_kji++)

    // 5. 3D linear interpolation
    const double Interpolated_Value = HaloMerger_Trilinear_Interpolation( x, y, z, IntCorner_Value, IntCorner_Coord[0], IntCorner_Coord[1], IntCorner_Coord[2] );

    return Interpolated_Value;

} // FUNCTION : HaloMerger_Get_Value_From_Halo_UM_IC_Data



//-------------------------------------------------------------------------------------------------------
// Function    :  Flag_User_UMICAMR_HaloMerger
// Description :  Template of user-defined flag criteria
//
// Note        :  1. Invoked by Flag_Check() using the function pointer "Flag_User_Ptr",
//                   which must be set by a test problem initializer
//                2. Enabled by the runtime option "OPT__FLAG_USER"
//
// Parameter   :  i,j,k     : Indices of the target element in the patch ptr[ amr->FluSg[lv] ][lv][PID]
//                lv        : Refinement level of the target patch
//                PID       : ID of the target patch
//                Threshold : User-provided threshold for the flag operation, which is loaded from the
//                            file "Input__Flag_User"
//
// Return      :  "true"  if the flag criteria are satisfied
//                "false" if the flag criteria are not satisfied
//-------------------------------------------------------------------------------------------------------
bool Flag_User_UMICAMR_HaloMerger( const int i, const int j, const int k, const int lv, const int PID, const double *Threshold )
{

// Define the AMR refinement similiar to the reconstructed halo UM_IC
   const double dh     = amr->dh[lv];                                                  // grid size
   const double Pos[3] = { amr->patch[0][lv][PID]->EdgeL[0] + (i+0.5)*dh,              // x,y,z position
                           amr->patch[0][lv][PID]->EdgeL[1] + (j+0.5)*dh,
                           amr->patch[0][lv][PID]->EdgeL[2] + (k+0.5)*dh  };

// flag cells within the target region [Threshold ... BoxSize-Threshold]
   const double EdgeL = Threshold[0];
   const double EdgeR = amr->BoxSize[0]-Threshold[0];    // here we have assumed a cubic box

   bool Flag;

   if (  Pos[0] >= (EdgeL+0.5*dh)  &&  Pos[0] < (EdgeR-0.5*dh)  &&
         Pos[1] >= (EdgeL+0.5*dh)  &&  Pos[1] < (EdgeR-0.5*dh)  &&
         Pos[2] >= (EdgeL+0.5*dh)  &&  Pos[2] < (EdgeR-0.5*dh)     )
      Flag = true;

   else
      Flag = false;

// ##########################################################################################################

   return Flag;

} // FUNCTION : Flag_User_UMICAMR_HaloMerger



//-------------------------------------------------------------------------------------------------------
// Function    :  OutputDensityProfile_HaloMerger
// Description :  Output the density profile
//
// Note        :  1.
//                2.
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void OutputDensityProfile_HaloMerger()
{

} // FUNCTION : OutputDensityProfile_HaloMerger



//-------------------------------------------------------------------------------------------------------
// Function    :  Aux_Record_HaloMerger
// Description :
//
// Note        :  1. Invoked by main() using the function pointer "Aux_Record_User_Ptr",
//                   which must be set by a test problem initializer
//                2. Enabled by the runtime option "OPT__RECORD_USER"
//
// Parameter   :  None
//-------------------------------------------------------------------------------------------------------
void Aux_Record_HaloMerger()
{
   Extrema_t Max_Dens;
   Max_Dens.Field     = _DENS;
   Max_Dens.Radius    = __FLT_MAX__; // entire domain
   Max_Dens.Center[0] = amr->BoxCenter[0];
   Max_Dens.Center[1] = amr->BoxCenter[1];
   Max_Dens.Center[2] = amr->BoxCenter[2];

   Aux_FindExtrema( &Max_Dens, EXTREMA_MAX, 0, TOP_LEVEL, PATCH_LEAF );

   const double      Center[3]      = { Max_Dens.Coord[0], Max_Dens.Coord[1], Max_Dens.Coord[2] };
   const double      MaxRadius      = 0.5*amr->BoxSize[0];
   const double      MinBinSize     = amr->dh[MAX_LEVEL];
   const bool        LogBin         = true;
   const double      LogBinRatio    = 1.025;
   const bool        RemoveEmptyBin = true;
   const long        TVar[]         = { _DENS };
   const int         NProf          = 1;
   const int         MinLv          = 0;
   const int         MaxLv          = MAX_LEVEL;
   const PatchType_t PatchType      = PATCH_LEAF_PLUS_MAXNONLEAF;
   const double      PrepTime       = -1.0;

   Profile_t Prof_Dens;
   Profile_t *Prof[] = { &Prof_Dens };

   Aux_ComputeProfile( Prof, Center, MaxRadius, MinBinSize, LogBin, LogBinRatio, RemoveEmptyBin,
                       TVar, NProf, MinLv, MaxLv, PatchType, PrepTime );

   if ( MPI_Rank == 0 )
   {
      for (int p=0; p<NProf; p++)
      {
         char Filename[MAX_STRING];
         sprintf( Filename, "DensityProfile_Step_%09ld", Step );
         FILE *File = fopen( Filename, "w" );

         fprintf( File, "#%19s  %21s  %21s  %10s\n", "Radius", "Data", "Weight", "Cells" );
         for (int b=0; b<Prof[p]->NBin; b++)
            fprintf( File, "%20.14e  %21.14e  %21.14e  %10ld\n",
                     Prof[p]->Radius[b], Prof[p]->Data[b], Prof[p]->Weight[b], Prof[p]->NCell[b] );

         fclose( File );
      }
   }

// Enclosed Mass
   const int NBin = Prof[0]->NBin;

   double* Profile_EnclosedMass = new double [NBin];
   double* Profile_NegDens      = new double [NBin];

   Profile_EnclosedMass[0]      = 4.0*M_PI/3.0*CUBE(Prof[0]->Radius[0])*Prof[0]->Data[0];
   for (int b=1; b<NBin; b++)   Profile_EnclosedMass[b] = Profile_EnclosedMass[b-1] + LogLinearDensityShellMass( Prof[0]->Radius[b-1], Prof[0]->Radius[b], Prof[0]->Data[b-1], Prof[0]->Data[b] );
   
   for (int b=0; b<NBin; b++)   Profile_NegDens[b] = -Prof[0]->Data[b];

// Soliton
   const double CoreDensity = Max_Dens.Value;
   const double CoreRadius  = Mis_InterpolateFromTable( NBin, Profile_NegDens,  Prof[0]->Radius,      -0.5*CoreDensity );
   const double CoreMass    = Mis_InterpolateFromTable( NBin, Prof[0]->Radius,  Profile_EnclosedMass,  CoreRadius      );

   delete [] Profile_EnclosedMass;
   delete [] Profile_NegDens;

   if ( MPI_Rank == 0 )
   {
      FILE *File = fopen( "Record__CoreProperties", "a" );

      static bool FirstTime = true;

      if ( FirstTime )
      {
         fprintf( File, "#%11s  %21s  %21s  %21s  %21s  %21s\n", "Step", "Time", "CoreDensity", "CoreRadius", "CoreMass", "CoreRadiusCoreMass" );

         FirstTime = false;
      } // if ( FirstTime )

      fprintf( File, "  %10ld  %21.14e  %21.14e  %21.14e  %21.14e  %21.14e\n",
                       Step, Time[0], CoreDensity, CoreRadius, CoreMass, CoreRadius*CoreMass );

      fclose( File );
   }

} // FUNCTION : Aux_Record_HaloMerger



//-------------------------------------------------------------------------------------------------------
// Function    :  LinearDensityShellMass
// Description :  Get the shell mass between to radii according to a linear density profile
//
// Note        :  1. Assume the density profile \rho(r) is linear between two end points, r0 and r1:
//                   \rho(r) = \rho_0 + (\frac{(\rho_1-\rho_0)}{r_1-r_0})(r-r_0)
//                   Then, the integrated shell mass between the two end points is
//                   M_{shell} = \int_{r_0}^{r_1} \rho(r) 4\pi r^2 dr
//                             = 4\pi r_0^2 \Delta r   [ \frac{1}{2} \rho_0 + \frac{1}{2} \rho_1 ] +
//                               4\pi r_0   \Delta r^2 [ \frac{1}{3} \rho_0 + \frac{2}{3} \rho_1 ] +
//                               4\pi       \Delta r^3 [ \frac{1}{12}\rho_0 + \frac{3}{12}\rho_1 ]
//
// Parameter   :  r0   : Radius of the start point
//                r1   : Radius of the end point
//                rho0 : Density of the start point
//                rho1 : Density of the end point
//
// Return      :  Shell mass integrated from the linear density profile
//-------------------------------------------------------------------------------------------------------
double LinearDensityShellMass( const double r0, const double r1, const double rho0, const double rho1 )
{
   const double dr  = r1 - r0;

   return M_PI*dr*( r0*r0*( 6.0*rho0 + 6.0*rho1 ) + r0*dr*( 4.0*rho0 + 8.0*rho1 ) + dr*dr*( rho0 + 3.0*rho1 ) )/3.0;

} // FUNCTION : LinearDensityShellMass



//-------------------------------------------------------------------------------------------------------
// Function    :  LogLinearDensityShellMass
// Description :  Get the shell mass between to radii according to a linear density profile in log-log space
//
// Note        :  1. Assume the density profile \rho(r) in log-log space is linear between two end points, r0 and r1:
//                   log(\rho(r)) = log(\rho_0) + (\frac{(log(\rho_1)-log(\rho_0))}{log(r_1)-log(r_0)})(log(r)-log(r_0))
//                   or equivalently,
//                   \rho(r) = \rho_0 (\frac{r}{r_0})^p,
//                   where p = log( \rho_1/\rho_0 ) / log( r_1/r_0 )
//                   Then, the integrated shell mass between the two end points is
//                   M_{shell} = \int_{r_0}^{r_1} \rho(r) 4\pi r^2 dr
//                             = 4\pi \rho_0 \frac{ r_1^{p+3} - r_0^{p+3} }{ (p+3)r_0^p }, if p != -3
//                           ( = 4\pi \rho_0 r_0^{3} log(r_1/r_0), if p = -3 )
//
// Parameter   :  r0   : Radius of the start point
//                r1   : Radius of the end point
//                rho0 : Density of the start point
//                rho1 : Density of the end point
//
// Return      :  Shell mass integrated from the linear density profile in log-log space
//-------------------------------------------------------------------------------------------------------
double LogLinearDensityShellMass( const double r0, const double r1, const double rho0, const double rho1 )
{
   const double p = log( rho1/rho0 )/log( r1/r0 );

   if ( p == -3.0 )   return 4.0*M_PI*rho0*CUBE(r0)*log( r1/r0 );

   return 4.0*M_PI*rho0*( pow( r1, p+3.0 ) - pow( r0, p+3.0 ) )/( (p+3.0)*pow( r0, p ) );

} // FUNCTION : LinearDensityShellMass
#endif // #if ( MODEL == ELBDM  &&  defined GRAVITY )



//-------------------------------------------------------------------------------------------------------
// Function    :  Init_TestProb_ELBDM_HaloMerger
// Description :  Test problem initializer
//
// Note        :  None
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Init_TestProb_ELBDM_HaloMerger()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );


// validate the compilation flags and runtime parameters
   Validate();


#  if ( MODEL == ELBDM  &&  defined GRAVITY )
// set the problem-specific runtime parameters
   SetParameter();


// set the function pointers of various problem-specific routines
   Init_Function_User_Ptr  = SetGridIC;
   End_User_Ptr            = End_HaloMerger;
   Init_ExtPot_Ptr         = Init_ExtPot_ELBDM_HaloMerger;
   Flag_User_Ptr           = Flag_User_UMICAMR_HaloMerger;
   Output_User_Ptr         = OutputDensityProfile_HaloMerger;
   Aux_Record_User_Ptr     = Aux_Record_HaloMerger;
#  ifdef MASSIVE_PARTICLES
   Par_Init_ByFunction_Ptr = Par_Init_ByFunction_HaloMerger;
#  endif // ifdef MASSIVE_PARTICLES

#  endif // if ( MODEL == ELBDM  &&  defined GRAVITY )

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : Init_TestProb_ELBDM_HaloMerger
