#include "GAMER.h"
#include "TestProb.h"



// problem-specific global variables
// =======================================================================================
static double Jeans_Coeff0_1;       // proportional coefficient of the growing mode of the wave function 1
static double Jeans_Coeff0_2;       // proportional coefficient of the growing mode of the wave function 2
static double Jeans_Coeff1;         // proportional coefficient 1 for two-component small y approximation solution
static double Jeans_Coeff2;         // proportional coefficient 2 for two-component small y approximation solution
static double Jeans_Coeff3;         // proportional coefficient 3 for two-component small y approximation solution
static double Jeans_Coeff4;         // proportional coefficient 4 for two-component small y approximation solution
static double Jeans_Coeff5;         // proportional coefficient 5 for two-component large y approximation solution
static double Jeans_Coeff6;         // proportional coefficient 6 for two-component large y approximation solution
static double Jeans_Coeff7;         // proportional coefficient 7 for two-component large y approximation solution
static double Jeans_Coeff8;         // proportional coefficient 8 for two-component large y approximation solution
static double Jeans_Phase0;         // initial phase shift
static double Jeans_RhoBG;          // Background density
static double Jeans_RhoBG_frac_1;   // Background density fraction 1
static double Jeans_RhoBG_frac_2;   // Background density fraction 2
static int    Jeans_AnalyticalForm; // Form of analytical solution (1:Single Exact, 2:Two Large Y, 3:Two Small Y)

static double Jeans_Wavelength;     // wavelength
static double Jeans_WaveK;          // wavenumber
static double Jeans_WaveKj_1;       // critical wavenumber 1
static double Jeans_WaveKj_2;       // critical wavenumber 2
static bool   Jeans_Stable_1;       // true/false --> Jeans stable/unstable for 1
static bool   Jeans_Stable_2;       // true/false --> Jeans stable/unstable for 2
// =======================================================================================



// inline functions to calculate the amplitudes of the real and imaginary parts
// =======================================================================================
inline double Jeans_Single_deltaK( const double Coeff0, const double y )
{
   const double y2 = y*y;

   return Coeff0 * ( 3.0*cos(y) + 3.0*y*sin(y) - y2*cos(y) ) / y2; // single exact
}

inline double Jeans_Single_MddydeltaK( const double Coeff0, const double y )
{
   const double y2 = y*y;
   const double y3 = y*y2;
   const double y4 = y2*y2;

   return -Coeff0 * ( -6.0*y*cos(y) - 6.0*y2*sin(y) + 3.0*y3*cos(y) + y4*sin(y) ) / y4; // single exact
}

inline double Jeans_Single_deltaK_LargeY( const double Coeff0, const double y )
{

   return -Coeff0 * cos(y); //approx when y>>1
}


inline double Jeans_Single_MddydeltaK_LargeY( const double Coeff0, const double y )
{
   return -Coeff0 * sin(y); //approx when y>>1
}


inline double Jeans_Single_deltaK_SmallY( const double Coeff0, const double y )
{
   const double y2 = y*y;

   return Coeff0 * 3.0/y2;  //approx when y<<1
}

inline double Jeans_Single_MddydeltaK_SmallY( const double Coeff0, const double y )
{
   const double y2 = y*y;
   const double y3 = y*y2;

   return -Coeff0 * ( -6.0 ) / y3; //approx when y<<1
}

// =======================================================================================
inline double Jeans_Two_deltaK1_SmallY( const double Coeff1, const double Coeff2, const double RhoBG_frac_1, const double RhoBG_frac_2 , const double y )
{
   const double y2 = y*y;

   return Coeff1 * 1.0/y2 - Coeff2*(RhoBG_frac_2);  //approx when y<<1
}

inline double Jeans_Two_deltaK2_SmallY( const double Coeff1, const double Coeff2, const double RhoBG_frac_1, const double RhoBG_frac_2 , const double y )
{
   const double y2 = y*y;

   return Coeff1 * 1.0/y2 + Coeff2*(RhoBG_frac_1);  //approx when y<<1
}

inline double Jeans_Two_MddydeltaK1_SmallY( const double Coeff1, const double m1_m, const double y )
{
   const double y2 = y*y;
   const double y3 = y*y2;

   return -Coeff1 * m1_m *( -2.0 ) / y3; //approx when y<<1
}

inline double Jeans_Two_MddydeltaK2_SmallY( const double Coeff1, const double m2_m, const double y )
{
   const double y2 = y*y;
   const double y3 = y*y2;

   return -Coeff1 * m2_m *( -2.0 ) / y3; //approx when y<<1
}

// =======================================================================================
inline double Jeans_Two_deltaK1_LargeY( const double Coeff5, const double m_m1, const double y )
{
   return Coeff5 * cos(m_m1*y);  //approx when y>>1
}

inline double Jeans_Two_deltaK2_LargeY( const double Coeff7, const double m_m2, const double y )
{
   return Coeff7 * cos(m_m2*y);  //approx when y>>1
}

inline double Jeans_Two_MddydeltaK1_LargeY( const double Coeff5, const double m_m1, const double y )
{
   return Coeff5 * sin(m_m1*y);  //approx when y>>1
}

inline double Jeans_Two_MddydeltaK2_LargeY( const double Coeff7, const double m_m2, const double y )
{
   return Coeff7 * sin(m_m2*y);  //approx when y>>1
}

// =======================================================================================

static void OutputError();




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

#  ifndef FLOAT8
   Aux_Error( ERROR_INFO, "FLOAT8 must be enabled !!\n" );
#  endif

   for (int f=0; f<6; f++)
   if ( OPT__BC_FLU[f] != BC_FLU_PERIODIC )
      Aux_Error( ERROR_INFO, "must adopt periodic BC for fluid --> reset OPT__BC_FLU* !!\n" );

#  ifdef GRAVITY
   if ( OPT__BC_POT != BC_POT_PERIODIC )
      Aux_Error( ERROR_INFO, "must adopt periodic BC for gravity --> reset OPT__BC_POT !!\n" );
#  endif

#  ifdef COMOVING
   if ( OMEGA_M0 != 1.0 )
      Aux_Error( ERROR_INFO, "must adopt \"OMEGA_M0 = 1.0\" !!\n" );
#  endif

   if ( amr->BoxSize[0] != amr->BoxSize[1]  ||  amr->BoxSize[0] != amr->BoxSize[2] )
      Aux_Error( ERROR_INFO, "simulation domain must be cubic --> reset BOX_SIZE !!\n" );


// warnings
   if ( MPI_Rank == 0 )
   {
      if ( !OPT__OUTPUT_USER )
         Aux_Message( stderr, "WARNING : it's recommended to enable OPT__OUTPUT_USER !!\n" );
   }


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Validating test problem %d ... done\n", TESTPROB_ID );

} // FUNCTION : Validate



#if ( MODEL == ELBDM  &&  defined GRAVITY  &&  defined COMOVING )
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
   ReadPara->Add( "Jeans_Coeff0_1",    &Jeans_Coeff0_1,       -1.0,           0.0,              NoMax_double      );
   ReadPara->Add( "Jeans_Coeff0_2",    &Jeans_Coeff0_2,       -1.0,           0.0,              NoMax_double      );
   ReadPara->Add( "Jeans_Coeff1",      &Jeans_Coeff1,          0.0,           NoMin_double,     NoMax_double      );
   ReadPara->Add( "Jeans_Coeff2",      &Jeans_Coeff2,          0.0,           NoMin_double,     NoMax_double      );
   ReadPara->Add( "Jeans_Coeff3",      &Jeans_Coeff3,          0.0,           NoMin_double,     NoMax_double      );
   ReadPara->Add( "Jeans_Coeff4",      &Jeans_Coeff4,          0.0,           NoMin_double,     NoMax_double      );
   ReadPara->Add( "Jeans_Coeff5",      &Jeans_Coeff5,          0.0,           NoMin_double,     NoMax_double      );
   ReadPara->Add( "Jeans_Coeff6",      &Jeans_Coeff6,          0.0,           NoMin_double,     NoMax_double      );
   ReadPara->Add( "Jeans_Coeff7",      &Jeans_Coeff7,          0.0,           NoMin_double,     NoMax_double      );
   ReadPara->Add( "Jeans_Coeff8",      &Jeans_Coeff8,          0.0,           NoMin_double,     NoMax_double      );
   ReadPara->Add( "Jeans_Phase0",      &Jeans_Phase0,          0.0,           NoMin_double,     NoMax_double      );
   ReadPara->Add( "Jeans_RhoBG_frac_1",&Jeans_RhoBG_frac_1,    1.0,           0.0,                       1.0      );
   ReadPara->Add( "Jeans_RhoBG_frac_2",&Jeans_RhoBG_frac_2,    0.0,           0.0,                       1.0      );
   ReadPara->Add( "Jeans_AnalyticalForm",&Jeans_AnalyticalForm,  1,             1,                         5      );

   ReadPara->Read( FileName );

   delete ReadPara;

// (1-2) set the default values

// (1-3) check the runtime parameters


// (2) set the problem-specific derived parameters
   Jeans_RhoBG      = 1.0;
   Jeans_Wavelength = amr->BoxSize[0]/sqrt(3.0);   // assuming cubic simulation domain
   Jeans_WaveK      = 2.0*M_PI/Jeans_Wavelength;
   Jeans_WaveKj_1   = POW( 6.0*Time[0]*SQR(ELBDM_ETA1), 0.25 );
   Jeans_WaveKj_2   = POW( 6.0*Time[0]*SQR(ELBDM_ETA2), 0.25 );
   Jeans_Stable_1   = ( Jeans_WaveK > Jeans_WaveKj_1 ) ? true : false;
   Jeans_Stable_2   = ( Jeans_WaveK > Jeans_WaveKj_2 ) ? true : false;


// (3) reset other general-purpose parameters
//     --> a helper macro PRINT_WARNING is defined in TestProb.h
   const long   End_Step_Default = __INT_MAX__;
// End_T : (stable/unstable) --> (1 period in the high-k limit / grow by a factor of 50 in the low-k limit)
   const double End_T_Default    = ( Jeans_Stable_1) ?
                                    A_INIT + pow( 0.5*Jeans_WaveK*Jeans_WaveK/M_PI/ELBDM_ETA1, 2.0 ) :
                                    A_INIT*50;
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
      const double y  = SQR(Jeans_WaveK)/ELBDM_ETA1*pow( Time[0], -0.5 );
      const double y1 = SQR(Jeans_WaveK)/ELBDM_ETA1*pow( Time[0], -0.5 );
      const double y2 = SQR(Jeans_WaveK)/ELBDM_ETA2*pow( Time[0], -0.5 );
      const double y1_end = SQR(Jeans_WaveK)/ELBDM_ETA1*pow( END_T, -0.5 );
      const double y2_end = SQR(Jeans_WaveK)/ELBDM_ETA2*pow( END_T, -0.5 );

      Aux_Message( stdout, "=============================================================================\n" );
      Aux_Message( stdout, "  test problem ID       = %d\n",     TESTPROB_ID                      );
      Aux_Message( stdout, "  coefficient 0 for 1   = %13.7e\n", Jeans_Coeff0_1                   );
      Aux_Message( stdout, "  coefficient 0 for 2   = %13.7e\n", Jeans_Coeff0_2                   );
      Aux_Message( stdout, "  coefficient 1         = %13.7e\n", Jeans_Coeff1                     );
      Aux_Message( stdout, "  coefficient 2         = %13.7e\n", Jeans_Coeff2                     );
      Aux_Message( stdout, "  coefficient 3         = %13.7e\n", Jeans_Coeff3                     );
      Aux_Message( stdout, "  coefficient 4         = %13.7e\n", Jeans_Coeff4                     );
      Aux_Message( stdout, "  coefficient 5         = %13.7e\n", Jeans_Coeff5                     );
      Aux_Message( stdout, "  coefficient 6         = %13.7e\n", Jeans_Coeff6                     );
      Aux_Message( stdout, "  coefficient 7         = %13.7e\n", Jeans_Coeff7                     );
      Aux_Message( stdout, "  coefficient 8         = %13.7e\n", Jeans_Coeff8                     );
      Aux_Message( stdout, "  initial phase shift   = %13.7e\n", Jeans_Phase0                     );
      Aux_Message( stdout, "  RhoBG                 = %13.7e\n", Jeans_RhoBG                      );
      Aux_Message( stdout, "  RhoBG fraction1       = %13.7e\n", Jeans_RhoBG_frac_1               );
      Aux_Message( stdout, "  RhoBG1                = %13.7e\n", Jeans_RhoBG*Jeans_RhoBG_frac_1   );
      Aux_Message( stdout, "  RhoBG fraction2       = %13.7e\n", Jeans_RhoBG_frac_2               );
      Aux_Message( stdout, "  RhoBG2                = %13.7e\n", Jeans_RhoBG*Jeans_RhoBG_frac_2   );
      if ( Jeans_AnalyticalForm == 1 ){
         Aux_Message( stdout, "  Analytical Form       = %d\n",     Jeans_AnalyticalForm              );
         Aux_Message( stdout, "             -> Single Exact\n"                                        );
         Aux_Message( stdout, "  real part 1 amplitude = %13.7e\n", 0.5*sqrt(Jeans_RhoBG*Jeans_RhoBG_frac_1)*Jeans_Single_deltaK(Jeans_Coeff0_1, y1));
         Aux_Message( stdout, "  imag part 1 amplitude = %13.7e\n", 0.5*sqrt(Jeans_RhoBG*Jeans_RhoBG_frac_1)*Jeans_Single_MddydeltaK(Jeans_Coeff0_1, y1));
         Aux_Message( stdout, "  real part 2 amplitude = %13.7e\n", 0.5*sqrt(Jeans_RhoBG*Jeans_RhoBG_frac_2)*Jeans_Single_deltaK(Jeans_Coeff0_2, y2));
         Aux_Message( stdout, "  imag part 2 amplitude = %13.7e\n", 0.5*sqrt(Jeans_RhoBG*Jeans_RhoBG_frac_2)*Jeans_Single_MddydeltaK(Jeans_Coeff0_2, y2));
         Aux_Message( stdout, "  real part 1 normalized amplitude = %13.7e\n", 0.5*Jeans_Single_deltaK(Jeans_Coeff0_1, y1));
         Aux_Message( stdout, "  imag part 1 normalized amplitude = %13.7e\n", 0.5*Jeans_Single_MddydeltaK(Jeans_Coeff0_1, y1));
         Aux_Message( stdout, "  real part 2 normalized amplitude = %13.7e\n", 0.5*Jeans_Single_deltaK(Jeans_Coeff0_2, y2));
         Aux_Message( stdout, "  imag part 2 normalized amplitude = %13.7e\n", 0.5*Jeans_Single_MddydeltaK(Jeans_Coeff0_2, y2));
      }
      if ( Jeans_AnalyticalForm == 2 ){
         Aux_Message( stdout, "  Analytical Form       = %d\n",     Jeans_AnalyticalForm              );
         Aux_Message( stdout, "  -> Two Large Y Approximation\n"                                      );
         Aux_Message( stdout, "  real part 1 amplitude = %13.7e\n", 0.5*sqrt(Jeans_RhoBG*Jeans_RhoBG_frac_1)*Jeans_Two_deltaK1_LargeY(Jeans_Coeff5, ELBDM_ETA1/ELBDM_ETA1, y));
         Aux_Message( stdout, "  imag part 1 amplitude = %13.7e\n", 0.5*sqrt(Jeans_RhoBG*Jeans_RhoBG_frac_1)*Jeans_Two_MddydeltaK1_LargeY(Jeans_Coeff5, ELBDM_ETA1/ELBDM_ETA1, y));
         Aux_Message( stdout, "  real part 2 amplitude = %13.7e\n", 0.5*sqrt(Jeans_RhoBG*Jeans_RhoBG_frac_2)*Jeans_Two_deltaK2_LargeY(Jeans_Coeff7, ELBDM_ETA1/ELBDM_ETA2, y));
         Aux_Message( stdout, "  imag part 2 amplitude = %13.7e\n", 0.5*sqrt(Jeans_RhoBG*Jeans_RhoBG_frac_2)*Jeans_Two_MddydeltaK2_LargeY(Jeans_Coeff7, ELBDM_ETA1/ELBDM_ETA2, y));
         Aux_Message( stdout, "  real part 1 normalized amplitude = %13.7e\n", 0.5*Jeans_Two_deltaK1_LargeY(Jeans_Coeff5, ELBDM_ETA1/ELBDM_ETA1, y));
         Aux_Message( stdout, "  imag part 1 normalized amplitude = %13.7e\n", 0.5*Jeans_Two_MddydeltaK1_LargeY(Jeans_Coeff5, ELBDM_ETA1/ELBDM_ETA1, y));
         Aux_Message( stdout, "  real part 2 normalized amplitude = %13.7e\n", 0.5*Jeans_Two_deltaK2_LargeY(Jeans_Coeff7, ELBDM_ETA1/ELBDM_ETA2, y));
         Aux_Message( stdout, "  imag part 2 normalized amplitude = %13.7e\n", 0.5*Jeans_Two_MddydeltaK2_LargeY(Jeans_Coeff7, ELBDM_ETA1/ELBDM_ETA2, y));
      }
      if ( Jeans_AnalyticalForm == 3 ){
         Aux_Message( stdout, "  Analytical Form       = %d\n",     Jeans_AnalyticalForm              );
         Aux_Message( stdout, "  -> Two Small Y Approximation\n"                                      );
         Aux_Message( stdout, "  real part 1 amplitude = %13.7e\n", 0.5*sqrt(Jeans_RhoBG*Jeans_RhoBG_frac_1)*Jeans_Two_deltaK1_SmallY(Jeans_Coeff1, Jeans_Coeff2, Jeans_RhoBG_frac_1, Jeans_RhoBG_frac_2, y));
         Aux_Message( stdout, "  imag part 1 amplitude = %13.7e\n", 0.5*sqrt(Jeans_RhoBG*Jeans_RhoBG_frac_1)*Jeans_Two_MddydeltaK1_SmallY(Jeans_Coeff1, ELBDM_ETA1/ELBDM_ETA1, y));
         Aux_Message( stdout, "  real part 2 amplitude = %13.7e\n", 0.5*sqrt(Jeans_RhoBG*Jeans_RhoBG_frac_2)*Jeans_Two_deltaK2_SmallY(Jeans_Coeff1, Jeans_Coeff2, Jeans_RhoBG_frac_1, Jeans_RhoBG_frac_2, y));
         Aux_Message( stdout, "  imag part 2 amplitude = %13.7e\n", 0.5*sqrt(Jeans_RhoBG*Jeans_RhoBG_frac_2)*Jeans_Two_MddydeltaK2_SmallY(Jeans_Coeff1, ELBDM_ETA2/ELBDM_ETA1, y));
         Aux_Message( stdout, "  real part 1 normalized amplitude = %13.7e\n", 0.5*Jeans_Two_deltaK1_SmallY(Jeans_Coeff1, Jeans_Coeff2, Jeans_RhoBG_frac_1, Jeans_RhoBG_frac_2, y));
         Aux_Message( stdout, "  imag part 1 normalized amplitude = %13.7e\n", 0.5*Jeans_Two_MddydeltaK1_SmallY(Jeans_Coeff1, ELBDM_ETA1/ELBDM_ETA1, y));
         Aux_Message( stdout, "  real part 2 normalized amplitude = %13.7e\n", 0.5*Jeans_Two_deltaK2_SmallY(Jeans_Coeff1, Jeans_Coeff2, Jeans_RhoBG_frac_1, Jeans_RhoBG_frac_2, y));
         Aux_Message( stdout, "  imag part 2 normalized amplitude = %13.7e\n", 0.5*Jeans_Two_MddydeltaK2_SmallY(Jeans_Coeff1, ELBDM_ETA2/ELBDM_ETA1, y));
      }
      if ( Jeans_AnalyticalForm == 4 ){
         Aux_Message( stdout, "  Analytical Form       = %d\n",     Jeans_AnalyticalForm              );
         Aux_Message( stdout, "    -> Single LargeY Approx \n"                                        );
         Aux_Message( stdout, "  real part 1 amplitude = %13.7e\n", 0.5*sqrt(Jeans_RhoBG*Jeans_RhoBG_frac_1)*Jeans_Single_deltaK_LargeY(Jeans_Coeff0_1, y1));
         Aux_Message( stdout, "  imag part 1 amplitude = %13.7e\n", 0.5*sqrt(Jeans_RhoBG*Jeans_RhoBG_frac_1)*Jeans_Single_MddydeltaK_LargeY(Jeans_Coeff0_1, y1));
         Aux_Message( stdout, "  real part 2 amplitude = %13.7e\n", 0.5*sqrt(Jeans_RhoBG*Jeans_RhoBG_frac_2)*Jeans_Single_deltaK_LargeY(Jeans_Coeff0_2, y2));
         Aux_Message( stdout, "  imag part 2 amplitude = %13.7e\n", 0.5*sqrt(Jeans_RhoBG*Jeans_RhoBG_frac_2)*Jeans_Single_MddydeltaK_LargeY(Jeans_Coeff0_2, y2));
         Aux_Message( stdout, "  real part 1 normalized amplitude = %13.7e\n", 0.5*Jeans_Single_deltaK_LargeY(Jeans_Coeff0_1, y1));
         Aux_Message( stdout, "  imag part 1 normalized amplitude = %13.7e\n", 0.5*Jeans_Single_MddydeltaK_LargeY(Jeans_Coeff0_1, y1));
         Aux_Message( stdout, "  real part 2 normalized amplitude = %13.7e\n", 0.5*Jeans_Single_deltaK_LargeY(Jeans_Coeff0_2, y2));
         Aux_Message( stdout, "  imag part 2 normalized amplitude = %13.7e\n", 0.5*Jeans_Single_MddydeltaK_LargeY(Jeans_Coeff0_2, y2));
      }
      if ( Jeans_AnalyticalForm == 5 ){
         Aux_Message( stdout, "  Analytical Form       = %d\n",     Jeans_AnalyticalForm              );
         Aux_Message( stdout, "     -> Single SmallY Approx\n"                                        );
         Aux_Message( stdout, "  real part 1 amplitude = %13.7e\n", 0.5*sqrt(Jeans_RhoBG*Jeans_RhoBG_frac_1)*Jeans_Single_deltaK_SmallY(Jeans_Coeff0_1, y1));
         Aux_Message( stdout, "  imag part 1 amplitude = %13.7e\n", 0.5*sqrt(Jeans_RhoBG*Jeans_RhoBG_frac_1)*Jeans_Single_MddydeltaK_SmallY(Jeans_Coeff0_1, y1));
         Aux_Message( stdout, "  real part 2 amplitude = %13.7e\n", 0.5*sqrt(Jeans_RhoBG*Jeans_RhoBG_frac_2)*Jeans_Single_deltaK_SmallY(Jeans_Coeff0_2, y2));
         Aux_Message( stdout, "  imag part 2 amplitude = %13.7e\n", 0.5*sqrt(Jeans_RhoBG*Jeans_RhoBG_frac_2)*Jeans_Single_MddydeltaK_SmallY(Jeans_Coeff0_2, y2));
         Aux_Message( stdout, "  real part 1 normalized amplitude = %13.7e\n", 0.5*Jeans_Single_deltaK_SmallY(Jeans_Coeff0_1, y1));
         Aux_Message( stdout, "  imag part 1 normalized amplitude = %13.7e\n", 0.5*Jeans_Single_MddydeltaK_SmallY(Jeans_Coeff0_1, y1));
         Aux_Message( stdout, "  real part 2 normalized amplitude = %13.7e\n", 0.5*Jeans_Single_deltaK_SmallY(Jeans_Coeff0_2, y2));
         Aux_Message( stdout, "  imag part 2 normalized amplitude = %13.7e\n", 0.5*Jeans_Single_MddydeltaK_SmallY(Jeans_Coeff0_2, y2));
      }
      Aux_Message( stdout, "  box size              = %13.7e\n", amr->BoxSize[0]                  );
      Aux_Message( stdout, "  wavelength            = %13.7e\n", Jeans_Wavelength                 );
      Aux_Message( stdout, "  wavenumber            = %13.7e\n", Jeans_WaveK                      );
      Aux_Message( stdout, "  ETA1                  = %13.7e\n", ELBDM_ETA1                       );
      Aux_Message( stdout, "  ETA2                  = %13.7e\n", ELBDM_ETA2                       );
      Aux_Message( stdout, "  initial a             = %13.7e\n", Time[0]                          );
      Aux_Message( stdout, "  end     a             = %13.7e\n", END_T                            );
      Aux_Message( stdout, "  initial y1            = %13.7e\n", y1                               );
      Aux_Message( stdout, "  end     y1            = %13.7e\n", y1_end                           );
      Aux_Message( stdout, "  initial y2            = %13.7e\n", y2                               );
      Aux_Message( stdout, "  end     y2            = %13.7e\n", y2_end                           );
      Aux_Message( stdout, "  critical y            = %13.7e\n", sqrt(6)                          );
      Aux_Message( stdout, "  critical wavenumber1  = %13.7e\n", Jeans_WaveKj_1                   );
      Aux_Message( stdout, "  critical wavenumber2  = %13.7e\n", Jeans_WaveKj_2                   );
      Aux_Message( stdout, "  Jeans stable1         = %s\n",     (Jeans_Stable_1)?"YES":"NO"      );
      Aux_Message( stdout, "  Jeans stable2         = %s\n",     (Jeans_Stable_2)?"YES":"NO"      );
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

   const double r       = 1.0/sqrt(3.0)*( x + y + z );
   const double Jeans_y = SQR(Jeans_WaveK)/ELBDM_ETA1*pow( Time, -0.5 );
   const double Jeans_y1= SQR(Jeans_WaveK)/ELBDM_ETA1*pow( Time, -0.5 );
   const double Jeans_y2= SQR(Jeans_WaveK)/ELBDM_ETA2*pow( Time, -0.5 );
   const double Phase   = Jeans_WaveK*r + Jeans_Phase0;

   if ( Jeans_AnalyticalForm == 1 ){
      // Exact Single
      fluid[REAL1] = sqrt(Jeans_RhoBG_frac_1 * Jeans_RhoBG) + 0.5*sqrt(Jeans_RhoBG*Jeans_RhoBG_frac_1)*Jeans_Single_deltaK(Jeans_Coeff0_1, Jeans_y1)*cos( Phase );
      fluid[IMAG1] =                                          0.5*sqrt(Jeans_RhoBG*Jeans_RhoBG_frac_1)*Jeans_Single_MddydeltaK(Jeans_Coeff0_1, Jeans_y1)*cos( Phase );
      fluid[DENS1] = SQR(fluid[REAL1]) + SQR(fluid[IMAG1]);

      fluid[REAL2] = sqrt(Jeans_RhoBG_frac_2 * Jeans_RhoBG) + 0.5*sqrt(Jeans_RhoBG*Jeans_RhoBG_frac_2)*Jeans_Single_deltaK(Jeans_Coeff0_2, Jeans_y2)*cos( Phase );
      fluid[IMAG2] =                                          0.5*sqrt(Jeans_RhoBG*Jeans_RhoBG_frac_2)*Jeans_Single_MddydeltaK(Jeans_Coeff0_2, Jeans_y2)*cos( Phase );
      fluid[DENS2] = SQR(fluid[REAL2]) + SQR(fluid[IMAG2]);
   }
   else if ( Jeans_AnalyticalForm == 2 ){
      // Two Large Y Approximation
      fluid[REAL1] = sqrt(Jeans_RhoBG_frac_1 * Jeans_RhoBG) + 0.5*sqrt(Jeans_RhoBG*Jeans_RhoBG_frac_1)*Jeans_Two_deltaK1_LargeY(Jeans_Coeff5, ELBDM_ETA1/ELBDM_ETA1, Jeans_y)*cos( Phase );
      fluid[IMAG1] =                                          0.5*sqrt(Jeans_RhoBG*Jeans_RhoBG_frac_1)*Jeans_Two_MddydeltaK1_LargeY(Jeans_Coeff5, ELBDM_ETA1/ELBDM_ETA1, Jeans_y)*cos( Phase );
      fluid[DENS1] = SQR(fluid[REAL1]) + SQR(fluid[IMAG1]);

      fluid[REAL2] = sqrt(Jeans_RhoBG_frac_2 * Jeans_RhoBG) + 0.5*sqrt(Jeans_RhoBG*Jeans_RhoBG_frac_2)*Jeans_Two_deltaK2_LargeY(Jeans_Coeff7, ELBDM_ETA1/ELBDM_ETA2, Jeans_y)*cos( Phase );
      fluid[IMAG2] =                                          0.5*sqrt(Jeans_RhoBG*Jeans_RhoBG_frac_2)*Jeans_Two_MddydeltaK2_LargeY(Jeans_Coeff7, ELBDM_ETA1/ELBDM_ETA2, Jeans_y)*cos( Phase );
      fluid[DENS2] = SQR(fluid[REAL2]) + SQR(fluid[IMAG2]);
   }
   else if ( Jeans_AnalyticalForm == 3 ){
     // Two Small Y Approximation
      fluid[REAL1] = sqrt(Jeans_RhoBG_frac_1 * Jeans_RhoBG) + 0.5*sqrt(Jeans_RhoBG*Jeans_RhoBG_frac_1)*Jeans_Two_deltaK1_SmallY(Jeans_Coeff1, Jeans_Coeff2, Jeans_RhoBG_frac_1, Jeans_RhoBG_frac_2, Jeans_y)*cos( Phase );
      fluid[IMAG1] =                                          0.5*sqrt(Jeans_RhoBG*Jeans_RhoBG_frac_1)*Jeans_Two_MddydeltaK1_SmallY(Jeans_Coeff1, ELBDM_ETA1/ELBDM_ETA1, Jeans_y)*cos( Phase );
      fluid[DENS1] = SQR(fluid[REAL1]) + SQR(fluid[IMAG1]);

      fluid[REAL2] = sqrt(Jeans_RhoBG_frac_2 * Jeans_RhoBG) + 0.5*sqrt(Jeans_RhoBG*Jeans_RhoBG_frac_2)*Jeans_Two_deltaK2_SmallY(Jeans_Coeff1, Jeans_Coeff2, Jeans_RhoBG_frac_1, Jeans_RhoBG_frac_2, Jeans_y)*cos( Phase ); // Be careful yhere should be y1.
      fluid[IMAG2] =                                          0.5*sqrt(Jeans_RhoBG*Jeans_RhoBG_frac_2)*Jeans_Two_MddydeltaK2_SmallY(Jeans_Coeff1, ELBDM_ETA2/ELBDM_ETA1, Jeans_y)*cos( Phase ); // Be careful yhere should be y1.
      fluid[DENS2] = SQR(fluid[REAL2]) + SQR(fluid[IMAG2]);
   }
   else if ( Jeans_AnalyticalForm == 4 ){
      // Single y>>1
      fluid[REAL1] = sqrt(Jeans_RhoBG_frac_1 * Jeans_RhoBG) + 0.5*sqrt(Jeans_RhoBG*Jeans_RhoBG_frac_1)*Jeans_Single_deltaK_LargeY(Jeans_Coeff0_1, Jeans_y1)*cos( Phase );
      fluid[IMAG1] =                                          0.5*sqrt(Jeans_RhoBG*Jeans_RhoBG_frac_1)*Jeans_Single_MddydeltaK_LargeY(Jeans_Coeff0_1, Jeans_y1)*cos( Phase );
      fluid[DENS1] = SQR(fluid[REAL1]) + SQR(fluid[IMAG1]);

      fluid[REAL2] = sqrt(Jeans_RhoBG_frac_2 * Jeans_RhoBG) + 0.5*sqrt(Jeans_RhoBG*Jeans_RhoBG_frac_2)*Jeans_Single_deltaK_LargeY(Jeans_Coeff0_2, Jeans_y2)*cos( Phase );
      fluid[IMAG2] =                                          0.5*sqrt(Jeans_RhoBG*Jeans_RhoBG_frac_2)*Jeans_Single_MddydeltaK_LargeY(Jeans_Coeff0_2, Jeans_y2)*cos( Phase );
      fluid[DENS2] = SQR(fluid[REAL2]) + SQR(fluid[IMAG2]);
   }
   else if ( Jeans_AnalyticalForm == 5 ){
      // Single y<<1
      fluid[REAL1] = sqrt(Jeans_RhoBG_frac_1 * Jeans_RhoBG) + 0.5*sqrt(Jeans_RhoBG*Jeans_RhoBG_frac_1)*Jeans_Single_deltaK_SmallY(Jeans_Coeff0_1, Jeans_y1)*cos( Phase );
      fluid[IMAG1] =                                          0.5*sqrt(Jeans_RhoBG*Jeans_RhoBG_frac_1)*Jeans_Single_MddydeltaK_SmallY(Jeans_Coeff0_1, Jeans_y1)*cos( Phase );
      fluid[DENS1] = SQR(fluid[REAL1]) + SQR(fluid[IMAG1]);

      fluid[REAL2] = sqrt(Jeans_RhoBG_frac_2 * Jeans_RhoBG) + 0.5*sqrt(Jeans_RhoBG*Jeans_RhoBG_frac_2)*Jeans_Single_deltaK_SmallY(Jeans_Coeff0_2, Jeans_y2)*cos( Phase );
      fluid[IMAG2] =                                          0.5*sqrt(Jeans_RhoBG*Jeans_RhoBG_frac_2)*Jeans_Single_MddydeltaK_SmallY(Jeans_Coeff0_2, Jeans_y2)*cos( Phase );
      fluid[DENS2] = SQR(fluid[REAL2]) + SQR(fluid[IMAG2]);
   }
   else{
      fluid[REAL1] = 0.0;
      fluid[IMAG1] = 0.0;
      fluid[DENS1] = SQR(fluid[REAL1]) + SQR(fluid[IMAG1]);

      fluid[REAL2] = 0.0;
      fluid[IMAG2] = 0.0;
      fluid[DENS2] = SQR(fluid[REAL2]) + SQR(fluid[IMAG2]);

   }

} // FUNCTION : SetGridIC



//-------------------------------------------------------------------------------------------------------
// Function    :  OutputError
// Description :  Output the L1 error
//
// Note        :  1. Invoke Output_L1Error()
//                2. Use SetGridIC() to provide the analytical solution at any given time
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void OutputError()
{

   const char Prefix[100]     = "JeansInstabilityComoving";
   const OptOutputPart_t Part = OUTPUT_DIAG;

   Output_L1Error( SetGridIC, Prefix, Part, NULL_REAL, NULL_REAL, NULL_REAL );

} // FUNCTION : OutputError
#endif // #if ( MODEL == ELBDM  &&  defined GRAVITY  &&  defined COMOVING )



//-------------------------------------------------------------------------------------------------------
// Function    :  Init_TestProb_ELBDM_JeansInstabilityComoving
// Description :  Test problem initializer
//
// Note        :  None
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Init_TestProb_ELBDM_JeansInstabilityComoving()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );


// validate the compilation flags and runtime parameters
   Validate();


#  if ( MODEL == ELBDM  &&  defined GRAVITY  &&  defined COMOVING )
// set the problem-specific runtime parameters
   SetParameter();


   Init_Function_User_Ptr   = SetGridIC;
   Flag_User_Ptr            = NULL;
   Mis_GetTimeStep_User_Ptr = NULL;
   BC_User_Ptr              = NULL;
   Flu_ResetByUser_Func_Ptr = NULL;
   Output_User_Ptr          = OutputError;
   Aux_Record_User_Ptr      = NULL;
   End_User_Ptr             = NULL;
   Init_ExternalAcc_Ptr     = NULL;
   Init_ExternalPot_Ptr     = NULL;
#  endif // #if ( MODEL == ELBDM  &&  defined GRAVITY  &&  defined COMOVING )


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : Init_TestProb_ELBDM_JeansInstabilityComoving
