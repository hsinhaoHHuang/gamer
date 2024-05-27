#include "GAMER.h"



// problem-specific global variables
// =======================================================================================
static int      Halo_InputMode;                          // Halo input mode: 1/0 -> UM_IC/none
static int      Soliton_InputMode;                       // soliton input mode: 1/2 -> table/approximate analytical form
static double   Soliton_OuterSlope;                      // soliton outer slope (only used by Soliton_InputMode=2)
static int      Soliton_N;                               // total number of solitons
static double   Soliton_CoreRadiusAll;                   // core radius for all solitons
                                                         //    (<=0.0 --> hard coding each soliton)
static double   Soliton_CenterX;                         // center coordinates for soliton
static double   Soliton_CenterY;                         // center coordinates for soliton
static double   Soliton_CenterZ;                         // center coordinates for soliton
static double   Soliton_VelocityX;                       // velocity for soliton
static double   Soliton_VelocityY;                       // velocity for soliton
static double   Soliton_VelocityZ;                       // velocity for soliton
static double   Soliton2_CenterX;                        // center coordinates for soliton 2
static double   Soliton2_CenterY;                        // center coordinates for soliton 2
static double   Soliton2_CenterZ;                        // center coordinates for soliton 2
static double   Soliton2_VelocityX;                      // velocity for soliton 2
static double   Soliton2_VelocityY;                      // velocity for soliton 2
static double   Soliton2_VelocityZ;                      // velocity for soliton 2
static double   Soliton2_CoreRadius;                     // core radius for soliton 2
static char     Soliton_DensProf_Filename[MAX_STRING];   // filename of the reference soliton density profile

static int      Soliton_DensProf_NBin;                   // number of radial bins of the soliton density profile
static double  *Soliton_DensProf    = NULL;              // soliton density profile [radius/density]
static double  *Soliton_CoreRadius  = NULL;              // core radius of each soliton
static double (*Soliton_Center)[3]  = NULL;              // center coordinates of each soliton
static double (*Soliton_Velocity)[3]= NULL;              // velocity of each soliton
static double  *Soliton_ScaleL      = NULL;              // L/D: length/density scale factors of each soliton
                                                         //      (defined as the ratio between the core radii/peak
                                                         //      density of the target and reference soliton profiles)
static double  *Soliton_ScaleD      = NULL;
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
// ReadPara->Add( "KEY_IN_THE_FILE",      &VARIABLE,              DEFAULT,       MIN,              MAX               );
// ********************************************************************************************************************************
   ReadPara->Add( "Halo_InputMode",            &Halo_InputMode,             1,             0,                1                 );
   ReadPara->Add( "Soliton_InputMode",         &Soliton_InputMode,          1,             1,                2                 );
   ReadPara->Add( "Soliton_OuterSlope",        &Soliton_OuterSlope,        -8.0,           NoMin_double,     NoMax_double      );
   ReadPara->Add( "Soliton_N",                 &Soliton_N,                  0,             0,                NoMax_int         );
   ReadPara->Add( "Soliton_CoreRadiusAll",     &Soliton_CoreRadiusAll,      NoDef_double,  NoMin_double,     NoMax_double      );
   ReadPara->Add( "Soliton_CenterX",           &Soliton_CenterX,           -1.0,           NoMin_double,     NoMax_double      );
   ReadPara->Add( "Soliton_CenterY",           &Soliton_CenterY,           -1.0,           NoMin_double,     NoMax_double      );
   ReadPara->Add( "Soliton_CenterZ",           &Soliton_CenterZ,           -1.0,           NoMin_double,     NoMax_double      );
   ReadPara->Add( "Soliton_VelocityX",         &Soliton_VelocityX,          0.0,           NoMin_double,     NoMax_double      );
   ReadPara->Add( "Soliton_VelocityY",         &Soliton_VelocityY,          0.0,           NoMin_double,     NoMax_double      );
   ReadPara->Add( "Soliton_VelocityZ",         &Soliton_VelocityZ,          0.0,           NoMin_double,     NoMax_double      );
   ReadPara->Add( "Soliton2_CenterX",          &Soliton2_CenterX,          -1.0,           NoMin_double,     NoMax_double      );
   ReadPara->Add( "Soliton2_CenterY",          &Soliton2_CenterY,          -1.0,           NoMin_double,     NoMax_double      );
   ReadPara->Add( "Soliton2_CenterZ",          &Soliton2_CenterZ,          -1.0,           NoMin_double,     NoMax_double      );
   ReadPara->Add( "Soliton2_VelocityX",        &Soliton2_VelocityX,         0.0,           NoMin_double,     NoMax_double      );
   ReadPara->Add( "Soliton2_VelocityY",        &Soliton2_VelocityY,         0.0,           NoMin_double,     NoMax_double      );
   ReadPara->Add( "Soliton2_VelocityZ",        &Soliton2_VelocityZ,         0.0,           NoMin_double,     NoMax_double      );
   ReadPara->Add( "Soliton2_CoreRadius",       &Soliton2_CoreRadius,       -1.0,           NoMin_double,     NoMax_double      );
   ReadPara->Add( "Soliton_DensProf_Filename",  Soliton_DensProf_Filename,  NoDef_str,     Useless_str,      Useless_str       );

   ReadPara->Read( FileName );

   delete ReadPara;

// (1-2) set the default values

// (1-3) check the runtime parameters
   if ( OPT__INIT == INIT_BY_FUNCTION )
      Aux_Error( ERROR_INFO, "OPT__INIT=1 is not supported for this test problem !!\n" );

   if ( Soliton_N > 0 ){
   if ( Soliton_CoreRadiusAll == NoDef_double )
      Aux_Error( ERROR_INFO, "Runtime parameter \"Soliton_CoreRadiusAll\" is not set !!\n" );

// (2) reset other general-purpose parameters
//     --> a helper macro PRINT_RESET_PARA is defined in Macro.h
// (2-1) allocate memory
   Soliton_CoreRadius = new double [Soliton_N];
   Soliton_Center     = new double [Soliton_N][3];
   Soliton_Velocity   = new double [Soliton_N][3];
   Soliton_ScaleL     = new double [Soliton_N];
   Soliton_ScaleD     = new double [Soliton_N];
   if ( Soliton_CoreRadiusAll > 0.0 )
   {
      for (int t=0; t<Soliton_N; t++)  Soliton_CoreRadius[t] = Soliton_CoreRadiusAll;
   }

   else
   {
//    for Soliton_CoreRadiusAll <= 0.0, comment out the following line and hard code the core radius of each soliton
      Aux_Error( ERROR_INFO, "for Soliton_CoreRadiusAll <= 0.0, please comment out this error check and hard code "
                             "the core radius of each soliton !!\n" );
//    for (int t=0; t<Soliton_N; t++)  Soliton_CoreRadius[t] = XXX;
   }

   if ( Soliton_N == 2  &&  Soliton2_CoreRadius > 0.0 )
   {
      Soliton_CoreRadius[1] = Soliton2_CoreRadius;
   }

   if ( Soliton_CenterX  < 0.0 ) Soliton_CenterX = 0.5;
   if ( Soliton_CenterY  < 0.0 ) Soliton_CenterY = 0.5;
   if ( Soliton_CenterZ  < 0.0 ) Soliton_CenterZ = 0.5;
   if ( Soliton2_CenterX < 0.0 ) Soliton2_CenterX = 0.5;
   if ( Soliton2_CenterY < 0.0 ) Soliton2_CenterY = 0.5;
   if ( Soliton2_CenterZ < 0.0 ) Soliton2_CenterZ = 0.5;

   if ( Soliton_N == 1 )
   {
      Soliton_Center[0][0] = Soliton_CenterX*amr->BoxSize[0];
      Soliton_Center[0][1] = Soliton_CenterY*amr->BoxSize[1];
      Soliton_Center[0][2] = Soliton_CenterZ*amr->BoxSize[2];

      Soliton_Velocity[0][0] = Soliton_VelocityX;
      Soliton_Velocity[0][1] = Soliton_VelocityY;
      Soliton_Velocity[0][2] = Soliton_VelocityZ;
   }
   else if ( Soliton_N == 2)
   {
      Soliton_Center[0][0] = Soliton_CenterX*amr->BoxSize[0];
      Soliton_Center[0][1] = Soliton_CenterY*amr->BoxSize[1];
      Soliton_Center[0][2] = Soliton_CenterZ*amr->BoxSize[2];

      Soliton_Velocity[0][0] = Soliton_VelocityX;
      Soliton_Velocity[0][1] = Soliton_VelocityY;
      Soliton_Velocity[0][2] = Soliton_VelocityZ;

      Soliton_Center[1][0] = Soliton2_CenterX*amr->BoxSize[0];
      Soliton_Center[1][1] = Soliton2_CenterY*amr->BoxSize[1];
      Soliton_Center[1][2] = Soliton2_CenterZ*amr->BoxSize[2];

      Soliton_Velocity[1][0] = Soliton2_VelocityX;
      Soliton_Velocity[1][1] = Soliton2_VelocityY;
      Soliton_Velocity[1][2] = Soliton2_VelocityZ;

   }

   else
   {
//    for Soliton_RSeed<0, comment out the following line and hard code the center of each soliton
      Aux_Error( ERROR_INFO, "for Soliton_RSeed < 0 and Soliton_N > 1, please comment out this error check and hard code "
                             "the center of each soliton !!\n" );

      /*
      for (int t=0; t<Soliton_N; t++)
      for (int d=0; d<3; d++)          Soliton_Center[t][d] = XXX;
      */
   }

   if ( OPT__INIT != INIT_BY_RESTART  &&  Soliton_InputMode == 1)
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
         Soliton_ScaleD[t] = 1.0 / ( 4.0*M_PI*NEWTON_G*SQR(ELBDM_ETA)*POW4(Soliton_ScaleL[t]) );
      }
   } // if ( OPT__INIT != INIT_BY_RESTART )
   }

   const long   End_Step_Default = __INT_MAX__;
   const double End_T_Default    = 0.5;   // ~7 Gyr

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
      Aux_Message( stdout, "  test problem ID = %d\n", TESTPROB_ID         );
      Aux_Message( stdout, "  halo input mode                           = %d\n",     Halo_InputMode             );
      Aux_Message( stdout, "  total number of solitons                  = %d\n",     Soliton_N                  );
      if ( Soliton_N > 0){
      Aux_Message( stdout, "  soliton input mode                        = %d\n",     Soliton_InputMode          );
      if      ( Soliton_InputMode == 2 )
      Aux_Message( stdout, "  soliton outer slope                       = %13.6e\n", Soliton_OuterSlope         );
      else if ( Soliton_InputMode == 1 ) {
      Aux_Message( stdout, "  density profile filename                  = %s\n",     Soliton_DensProf_Filename  );
      Aux_Message( stdout, "  number of bins of the density profile     = %d\n",     Soliton_DensProf_NBin      );
      }
      Aux_Message( stdout, "\n" );
      Aux_Message( stdout, "  Soliton info:\n" );
      Aux_Message( stdout, "  %7s  %13s  %13s  %13s  %13s  %13s  %13s  %13s  %13s  %13s\n",
                   "ID", "CoreRadius", "ScaleL", "ScaleD", "Center_X", "Center_Y", "Center_Z", "Velocity_X", "Velocity_Y", "Velocity_Z" );
      for (int t=0; t<Soliton_N; t++)
      Aux_Message( stdout, "  %7d  %13.6e  %13.6e  %13.6e  %13.6e  %13.6e  %13.6e  %13.6e  %13.6e  %13.6e\n",
                   t, Soliton_CoreRadius[t], Soliton_ScaleL[t], Soliton_ScaleD[t],
                   Soliton_Center[t][0], Soliton_Center[t][1], Soliton_Center[t][2],
                   Soliton_Velocity[t][0], Soliton_Velocity[t][1], Soliton_Velocity[t][2] );
      }
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

   Aux_Error( ERROR_INFO, "OPT__INIT=1 is not supported for this test problem !!\n" );

} // FUNCTION : SetGridIC



//-------------------------------------------------------------------------------------------------------
// Function    :  Init_ByFile_ELBDM_IsolatedHalo
// Description :  Function to actually set the fluid field from the input uniform-mesh array
//
// Note        :  1. Invoked by Init_ByFile_AssignData() using the function pointer Init_ByFile_User_Ptr()
//                   --> The function pointer may be reset by various test problem initializers, in which case
//                       this funtion will become useless
//
// Parameter   :  fluid_out : Fluid field to be set
//                fluid_in  : Fluid field loaded from the uniform-mesh array (UM_IC)
//                nvar_in   : Number of variables in fluid_in
//                x/y/z     : Target physical coordinates
//                Time      : Target physical time
//                lv        : Target AMR level
//                AuxArray  : Auxiliary array
//
// Return      :  fluid_out
//-------------------------------------------------------------------------------------------------------
void Init_ByFile_ELBDM_IsolatedHalo( real fluid_out[], const real fluid_in[], const int nvar_in,
                            const double x, const double y, const double z, const double Time,
                            const int lv, double AuxArray[] )
{
   double Re, Im;

   if ( Halo_InputMode == 1 )
   {
      if ( nvar_in != 2 )  Aux_Error( ERROR_INFO, "nvar_in (%d) != 2 !!\n", nvar_in );

      Re = fluid_in[0];
      Im = fluid_in[1];
   }
   else if ( Halo_InputMode == 0 )
   {
      Re = 0.0;
      Im = 0.0;
   }
   else
      Aux_Error( ERROR_INFO, "Unsupported Halo_InputMode (%d) !!\n", Halo_InputMode );

   if ( Soliton_N > 0  ){
      double Phase;


//    initialize density as zero since there may be multiple solitons
      double Dens_soliton, Re_soliton, Im_soliton;
      Dens_soliton = 0.0;
      Re_soliton = 0.0;
      Im_soliton = 0.0;
      double r_tar;

      if ( Soliton_InputMode == 1 )
      {
         const double *Table_Radius  = Soliton_DensProf + 0*Soliton_DensProf_NBin;  // radius
         const double *Table_Density = Soliton_DensProf + 1*Soliton_DensProf_NBin;  // density

         double r_ref, dens_ref;

//       loop over all solitons to get the total density
         for (int t=0; t<Soliton_N; t++)
         {
            r_tar = sqrt( SQR(x-Soliton_Center[t][0]) + SQR(y-Soliton_Center[t][1]) + SQR(z-Soliton_Center[t][2]) );

//          rescale radius (target radius --> reference radius)
            r_ref = r_tar / Soliton_ScaleL[t];

//          linear interpolation
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

//          rescale density (reference density --> target density) and add to the fluid array
            Phase = ELBDM_ETA*( Soliton_Velocity[t][0]*x + Soliton_Velocity[t][1]*y + Soliton_Velocity[t][2]*z );
            Dens_soliton  = dens_ref*Soliton_ScaleD[t];
            Re_soliton   += sqrt( Dens_soliton ) * cos(Phase);
            Im_soliton   += sqrt( Dens_soliton ) * sin(Phase);
         } // for (int t=0; t<Soliton_N; t++)

         Re += Re_soliton;
         Im += Im_soliton;
      }
      else if ( Soliton_InputMode == 2 )
      {
         for (int t=0; t<Soliton_N; t++)
         {
            r_tar = sqrt( SQR(x-Soliton_Center[t][0]) + SQR(y-Soliton_Center[t][1]) + SQR(z-Soliton_Center[t][2]) );

            const double m22      = ELBDM_MASS*UNIT_M/(Const_eV/SQR(Const_c))/1.0e-22;
            const double rc_kpc   = Soliton_CoreRadius[t]*UNIT_L/Const_kpc;
            const double peak_rho = 1.945e7/SQR( m22*rc_kpc*rc_kpc )*Const_Msun/CUBE(Const_kpc)/(UNIT_M/CUBE(UNIT_L));

            Phase = ELBDM_ETA*( Soliton_Velocity[t][0]*x + Soliton_Velocity[t][1]*y + Soliton_Velocity[t][2]*z );
            Dens_soliton = peak_rho*pow( 1.0+9.06e-2*SQR(r_tar/Soliton_CoreRadius[t]), Soliton_OuterSlope );
            Re_soliton   += sqrt( Dens_soliton ) * cos(Phase);
            Im_soliton   += sqrt( Dens_soliton ) * sin(Phase);
         } // for (int t=0; t<Soliton_N; t++)
         Re += Re_soliton;
         Im += Im_soliton;
      }

      else
         Aux_Error( ERROR_INFO, "Unsupported Soliton_InputMode (%d) !!\n", Soliton_InputMode );
   }

#  if ( ELBDM_SCHEME == ELBDM_HYBRID )
   if ( amr->use_wave_flag[lv] ) {
#  endif
   fluid_out[DENS] = SQR( Re ) + SQR( Im );
   fluid_out[REAL] = Re;
   fluid_out[IMAG] = Im;
#  if ( ELBDM_SCHEME == ELBDM_HYBRID )
   } else { // if ( amr->use_wave_flag[lv] )
   fluid_out[DENS] = Re;
   fluid_out[PHAS] = Im;
   fluid_out[STUB] = 0.0;
   } // if ( amr->use_wave_flag[lv] ) ... else
#  endif

} // Init_ByFile_ELBDM_IsolatedHalo



//-------------------------------------------------------------------------------------------------------
// Function    :  End_IsolatedHalo
// Description :  Free memory before terminating the program
//
// Note        :  1. Linked to the function pointer "End_User_Ptr" to replace "End_User()"
//
// Parameter   :  None
//-------------------------------------------------------------------------------------------------------
void End_IsolatedHalo()
{

   if ( Soliton_N > 0 ){
      delete [] Soliton_DensProf;
      delete [] Soliton_CoreRadius;
      delete [] Soliton_Center;
      delete [] Soliton_ScaleL;
      delete [] Soliton_ScaleD;
   }

} // FUNCTION : End_IsolatedHalo



//-------------------------------------------------------------------------------------------------------
// Function    :  Flag_User_UMICAMR
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
bool Flag_User_UMICAMR( const int i, const int j, const int k, const int lv, const int PID, const double *Threshold )
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

} // FUNCTION : Flag_User_UMICAMR
#endif // #if ( MODEL == ELBDM  &&  defined GRAVITY )



//-------------------------------------------------------------------------------------------------------
// Function    :  Init_TestProb_ELBDM_IsolatedHalo
// Description :  Test problem initializer
//
// Note        :  None
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Init_TestProb_ELBDM_IsolatedHalo()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );


// validate the compilation flags and runtime parameters
   Validate();


#  if ( MODEL == ELBDM  &&  defined GRAVITY )
// set the problem-specific runtime parameters
   SetParameter();


   Init_Function_User_Ptr = SetGridIC;
   Init_ByFile_User_Ptr   = Init_ByFile_ELBDM_IsolatedHalo;
   Flag_User_Ptr          = Flag_User_UMICAMR;
   End_User_Ptr           = End_IsolatedHalo;
#  endif // if ( MODEL == ELBDM  &&  defined GRAVITY )


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : Init_TestProb_ELBDM_IsolatedHalo
