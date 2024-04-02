#include "GAMER.h"

#ifdef MASSIVE_PARTICLES

#include "Par_EquilibriumIC.h"



Par_EquilibriumIC::Par_EquilibriumIC()
{
}

Par_EquilibriumIC::~Par_EquilibriumIC()
{
}



//-------------------------------------------------------------------------------------------------------
// Function    :  Read_Filenames
// Description :  Read in file names of physical parameters
//
// Note        :
//
// Parameter   :  filename_para : pointer of the file name
//
// Return      :
//-------------------------------------------------------------------------------------------------------
void Par_EquilibriumIC::Read_Filenames( const char *filename_para )
{

   vector <string> EMPTY;

   filenames.Cloud_Num = GetParams( filename_para, "ParEqmIC_Cloud_Num", 1, "int", EMPTY );

   GetParams( filename_para, "ParEqmIC_Params_Filenames", filenames.Cloud_Num, "string", filenames.Params_Filenames );


   Aux_Message( stdout, "Checking Params_Filenames\n" );

   fstream file;

   for (int k=0; k<filenames.Cloud_Num; k++)
   {

      const char * c = filenames.Params_Filenames[k].c_str();

      file.open( c, ios::in );

      if ( !file )
      {
         Aux_Message( stdout, "Test Problem parameter file %s cannot be found !!\n", filenames.Params_Filenames[k].c_str() );
         Aux_Error( ERROR_INFO, "Error in the input of Params_Filenames !!\n" );
      }

      file.close();
   }

} // FUNCTION : Read_Filenames



//-------------------------------------------------------------------------------------------------------
// Function    :  Load_Physical_Params
// Description :  Load the physical parameters from a file
//
// Note        :
//
// Parameter   :  filename_para : the file name loader's FP (file parameters) structure
//                cloud_idx     : index of the cloud
//                NPar_AllRank  : particle total number (including all clouds)
//
// Return      :
//-------------------------------------------------------------------------------------------------------
void Par_EquilibriumIC::Load_Physical_Params( const FP filename_para, const int cloud_idx, const long NPar_AllRank )
{

   params.Cloud_Center   = new double[3]; // central coordinates
   params.Cloud_BulkVel  = new double[3]; // bulk velocity

   Aux_Message( stdout, "Reading physical parameters input file: %s\n", filename_para.Params_Filenames[cloud_idx].c_str() );

   // (1) load the problem-specific runtime parameters
   params.Cloud_Num        = filename_para.Cloud_Num;
   params.Params_Filenames = filename_para.Params_Filenames[cloud_idx];

   const char* FileName  = filename_para.Params_Filenames[cloud_idx].c_str();
   ReadPara_t *ReadPara  = new ReadPara_t;
   double ratio;

   // (1-1) add parameters in the following format:
   // --> note that VARIABLE, DEFAULT, MIN, and MAX must have the same data type
   // --> some handy constants (e.g., Useless_bool, Eps_double, NoMin_int, ...) are defined in "include/ReadPara.h"
   // ********************************************************************************************************************************
   // ReadPara->Add( "KEY_IN_THE_FILE",         &VARIABLE,                          DEFAULT,       MIN,              MAX               );
   // ********************************************************************************************************************************
   ReadPara->Add( "Cloud_RSeed",                &params.Cloud_RSeed,                123,           0,                NoMax_int         );
   ReadPara->Add( "Cloud_Rho0",                 &params.Cloud_Rho0,                 1.0,           Eps_double,       NoMax_double      );
   ReadPara->Add( "Cloud_R0",                   &params.Cloud_R0,                   0.1,           Eps_double,       NoMax_double      );
   ReadPara->Add( "Cloud_MaxR",                 &params.Cloud_MaxR,                 0.375,         Eps_double,       NoMax_double      );
   ReadPara->Add( "Cloud_CenterX",              &params.Cloud_Center[0],            NoDef_double,  NoMin_double,     NoMax_double      );
   ReadPara->Add( "Cloud_CenterY",              &params.Cloud_Center[1],            NoDef_double,  NoMin_double,     NoMax_double      );
   ReadPara->Add( "Cloud_CenterZ",              &params.Cloud_Center[2],            NoDef_double,  NoMin_double,     NoMax_double      );
   ReadPara->Add( "Cloud_BulkVelX",             &params.Cloud_BulkVel[0],           0.0,           NoMin_double,     NoMax_double      );
   ReadPara->Add( "Cloud_BulkVelY",             &params.Cloud_BulkVel[1],           0.0,           NoMin_double,     NoMax_double      );
   ReadPara->Add( "Cloud_BulkVelZ",             &params.Cloud_BulkVel[2],           0.0,           NoMin_double,     NoMax_double      );
   ReadPara->Add( "Cloud_MassProfNBin",         &params.Cloud_MassProfNBin,         1000,          2,                NoMax_int         );
   ReadPara->Add( "Cloud_Par_Num_Ratio",        &ratio,                             0.,            0.,               1.0               );
   ReadPara->Add( "Cloud_Type",                  params.Cloud_Type,                 NoDef_str,     Useless_str,      Useless_str       );
   ReadPara->Add( "Density_Table_Name",          params.Density_Table_Name,         NoDef_str,     Useless_str,      Useless_str       );
   ReadPara->Add( "AddExtPot",                  &params.AddExtPot,                  0,             0,                1                 );
   ReadPara->Add( "ExtPot_Table_Name",           params.ExtPot_Table_Name,          NoDef_str,     Useless_str,      Useless_str       );
   ReadPara->Add( "Cloud_Einasto_Power_Factor", &params.Cloud_Einasto_Power_Factor, 1.0,           0.1,              10.0              );

   ReadPara->Read( FileName );
   delete ReadPara;

   // Convert Cloud_Par_Num_Ratio to Cloud_Par_Num
   params.Cloud_Par_Num = long(ratio*NPar_AllRank);

   // Check whether user forgot to fill in Cloud_Par_Num_Ratio
   if ( params.Cloud_Par_Num == 0 )
   {
      Aux_Error( ERROR_INFO, "Cloud_Par_Num_Ratio is 0! There is no particle in this cloud!!" );
   }

   // (1-2) set the default values
   for (int d=0; d<3; d++)
      if ( params.Cloud_Center[d] == NoDef_double )  params.Cloud_Center[d] = 0.5*amr->BoxSize[d];

   // (2) make a note
   if ( MPI_Rank == 0 )
   {
      Aux_Message( stdout, "=============================================================================\n" );
      Aux_Message( stdout, "  test problem ID                           = %d\n",     TESTPROB_ID );
      Aux_Message( stdout, "  random seed for setting particle position = %d\n",     params.Cloud_RSeed );
      Aux_Message( stdout, "  peak density                              = %13.7e\n", params.Cloud_Rho0 );
      Aux_Message( stdout, "  scale radius                              = %13.7e\n", params.Cloud_R0 );
      Aux_Message( stdout, "  maximum radius of particles               = %13.7e\n", params.Cloud_MaxR );

      for (int d=0; d<3; d++)
      {
      Aux_Message( stdout, "  central coordinate [%d]                    = %13.7e\n", d, params.Cloud_Center[d] );
      Aux_Message( stdout, "  bulk velocity [%d]                         = %13.7e\n", d, params.Cloud_BulkVel[d] );
      }

      if ( strcmp( params.Cloud_Type, "Table" ) != 0 )
      Aux_Message( stdout, "  number of radial bins in the mass profile = %d\n",     params.Cloud_MassProfNBin );

      Aux_Message( stdout, "  Cloud_Type                                = %s\n",     params.Cloud_Type );
      Aux_Message( stdout, "=============================================================================\n" );
   }//if ( MPI_Rank == 0 )

   // (3) Warn against small R0
   if ( params.Cloud_R0 < amr->dh[MAX_LEVEL] )
      Aux_Message( stdout, "WARNING : Characteristic length R0:%f is smaller than spatial resolution %f!\n", params.Cloud_R0, amr->dh[MAX_LEVEL] );

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Setting runtime parameters ... done\n" );

   // (4) Check Cloud_Type and table filenames
   // Checking Cloud_Type
   Aux_Message( stdout, "Checking Cloud_Type\n" );

   int flag = 0;

   if      ( strcmp( params.Cloud_Type, "Plummer"   ) == 0 )   flag = 1;
   else if ( strcmp( params.Cloud_Type, "NFW"       ) == 0 )   flag = 1;
   else if ( strcmp( params.Cloud_Type, "Burkert"   ) == 0 )   flag = 1;
   else if ( strcmp( params.Cloud_Type, "Jaffe"     ) == 0 )   flag = 1;
   else if ( strcmp( params.Cloud_Type, "Hernquist" ) == 0 )   flag = 1;
   else if ( strcmp( params.Cloud_Type, "Einasto"   ) == 0 )   flag = 1;
   else if ( strcmp( params.Cloud_Type, "Table"     ) == 0 )   flag = 1;

   if ( flag == 0 )
   {
      Aux_Message( stdout, "%s is not a Model Type\n", params.Cloud_Type );
      Aux_Error( ERROR_INFO, "Error in the input of Cloud_Type !!\n" );
   }

   // Checking Density_Table_Name
   Aux_Message( stdout, "Checking Density_Table_Name\n" );

   if ( strcmp( params.Cloud_Type, "Table" ) == 0 )
   {

      fstream file;
      file.open( params.Density_Table_Name, ios::in );

      if ( !file )
      {
         Aux_Message( stdout, "Density profile %s cannot be found !!\n", params.Density_Table_Name );
         Aux_Error( ERROR_INFO, "Error in the input of Density_Table_Name !!\n" );
      }

      file.close();
   }

   // Checking ExtPot_Table_Name
   Aux_Message( stdout, "Checking ExtPot_Table_Name\n" );

   if ( params.AddExtPot )
   {

      fstream file;
      file.open( params.ExtPot_Table_Name, ios::in );

      if ( !file )
      {
         Aux_Message( stdout, "External potential profile %s cannot be found !!\n", params.ExtPot_Table_Name );
         Aux_Error( ERROR_INFO, "Error in the input of ExtPot_Table_Name!!\n" );
      }

      file.close();
   }

} // FUNCTION : Load_Physical_Params



//-------------------------------------------------------------------------------------------------------
// Function    :  Init
// Description :  Initialize all necessary tables of physical parameters, including radius, mass, density, gravitational potential
//
// Note        :
//
// Parameter   :
//
// Return      :
//-------------------------------------------------------------------------------------------------------
void Par_EquilibriumIC::Init()
{

#  ifndef SUPPORT_GSL
   Aux_Error( ERROR_INFO, "Must enable SUPPORT_GSL for Par_EquilibriumIC !!\n" );
#  endif

   Table_r                 = NULL;
   Table_Enclosed_Mass     = NULL;
   Table_Density           = NULL;
   Table_dRho_dr           = NULL;
   Table_dRho_dx           = NULL;
   Table_Gravity_Field     = NULL;
   Table_Gravity_Potential = NULL;

   prob_dens               = NULL;
   int_prob_dens           = NULL;
   psi                     = NULL;

   //Set random seeds
   Random_Num_Gen = new RandomNumber_t( 1 );
   Random_Num_Gen->SetSeed( 0, params.Cloud_RSeed );

   //Initialize densities with Table
   if ( strcmp( params.Cloud_Type, "Table" ) == 0 )
   {

      int Tcol_r[1]   =  {0};
      int Tcol_rho[1] =  {1};
      int Row_r_Table;

      Aux_Message( stdout, "Loading Density Profile Table: %s\n", params.Density_Table_Name );

      Row_r_Table = Aux_LoadTable( Table_r, params.Density_Table_Name, 1, Tcol_r,true,true );

      int Row_Density_Table;
      Row_Density_Table = Aux_LoadTable( Table_Density, params.Density_Table_Name , 1, Tcol_rho,true,true );

      if ( Row_r_Table != Row_Density_Table )
         Aux_Error( ERROR_INFO, "Density row number is not equal to radius row number in the profile file !! Please check this file.\n" );

      params.Cloud_MassProfNBin = Row_r_Table;

      // Radii in the density table must be no less than Cloud_MaxR
      if ( Table_r[params.Cloud_MassProfNBin-1] < params.Cloud_MaxR )
      {
         Aux_Error( ERROR_INFO, "Maximum radius in your density table is smaller then Cloud_MaxR! Please check!\n" );
      }

      Table_Enclosed_Mass     = new double [params.Cloud_MassProfNBin];
      Table_dRho_dr           = new double [params.Cloud_MassProfNBin];
      Table_Gravity_Field     = new double [params.Cloud_MassProfNBin];
      Table_Gravity_Potential = new double [params.Cloud_MassProfNBin];
      Table_dRho_dx           = new double [params.Cloud_MassProfNBin];

      prob_dens               = new double [params.Cloud_MassProfNBin];
      int_prob_dens           = new double [params.Cloud_MassProfNBin];
      psi                     = new double [params.Cloud_MassProfNBin];

      Init_Mass_Table();
      Init_Pot_Table ();
      Add_Ext_Pot    ();
      Init_Prob_Dens ();

   }

   else
   {
      Table_r                 = new double [params.Cloud_MassProfNBin];
      Table_Enclosed_Mass     = new double [params.Cloud_MassProfNBin];
      Table_Density           = new double [params.Cloud_MassProfNBin];
      Table_dRho_dr           = new double [params.Cloud_MassProfNBin];
      Table_Gravity_Field     = new double [params.Cloud_MassProfNBin];
      Table_Gravity_Potential = new double [params.Cloud_MassProfNBin];
      Table_dRho_dx           = new double [params.Cloud_MassProfNBin];

      prob_dens               = new double [params.Cloud_MassProfNBin];
      int_prob_dens           = new double [params.Cloud_MassProfNBin];
      psi                     = new double [params.Cloud_MassProfNBin];

      Init_Mass     ();
      Init_Pot      ();
      Add_Ext_Pot   ();
      Init_Prob_Dens();
   }

} // FUNCTION : Init



//-------------------------------------------------------------------------------------------------------
// Function    :  Par_SetEquilibriumIC
// Description :  Set particle's initial conditions (IC) for a cloud that is in equilibrium state
//
// Note        :
//
// Parameter   :  Mass_AllRank : An array of all particles' masses
//                Pos_AllRank  : An array of all particles' position vectors
//                Vel_AllRank  : An array of all particles' velocity vectors
//                Par_Idx0     : Starting index of particles in this cloud
//
// Return      :  Mass_AllRank
//                Pos_AllRank
//                Vel_AllRank
//-------------------------------------------------------------------------------------------------------
void Par_EquilibriumIC::Par_SetEquilibriumIC( real *Mass_AllRank, real *Pos_AllRank[3], real *Vel_AllRank[3], const long Par_Idx0 )
{

   double *Table_MassProf_r = NULL;
   double *Table_MassProf_M = NULL;
   double  TotM, ParM, dr, RanM, RanR, EstM, ErrM, ErrM_Max=-1.0, RanVec[3];
   double  Vmax, RanV, RanProb, Prob;


   // determine the total enclosed mass within the maximum radius
   TotM = Set_Mass( params.Cloud_MaxR );
   ParM = TotM / (params.Cloud_Par_Num);

   // construct the mass profile table
   Table_MassProf_r = new double [params.Cloud_MassProfNBin];
   Table_MassProf_M = new double [params.Cloud_MassProfNBin];

   dr = params.Cloud_MaxR / (params.Cloud_MassProfNBin-1);

   for (int b=0; b<params.Cloud_MassProfNBin; b++)
   {
      Table_MassProf_r[b] = dr*b;
      Table_MassProf_M[b] = Set_Mass( Table_MassProf_r[b] );
   }

   // set particle attributes
   for (long p=Par_Idx0; p<Par_Idx0+params.Cloud_Par_Num; p++)
   {
      // mass
      Mass_AllRank[p] = ParM;

      // position
      // --> sample from the cumulative mass profile with linear interpolation
      RanM = Random_Num_Gen->GetValue( 0, 0.0, 1.0 )*TotM;
      RanR = Mis_InterpolateFromTable( params.Cloud_MassProfNBin, Table_MassProf_M, Table_MassProf_r, RanM );

      // record the maximum error
      EstM     = Set_Mass( RanR );
      ErrM     = fabs( (EstM-RanM)/RanM );
      ErrM_Max = fmax( ErrM, ErrM_Max );

      // randomly set the position vector with a given radius
      RanVec_FixRadius( RanR, RanVec );
      for (int d=0; d<3; d++)    Pos_AllRank[d][p] = RanVec[d] + params.Cloud_Center[d];

      // check periodicity
      for (int d=0; d<3; d++)
      {
         if ( OPT__BC_FLU[d*2] == BC_FLU_PERIODIC )
            Pos_AllRank[d][p] = FMOD( Pos_AllRank[d][p]+(real)amr->BoxSize[d], (real)amr->BoxSize[d] );
      }

      // velocity
      double a3 = RanR/params.Cloud_R0;

      RanV = Set_Velocity(a3);

      // randomly set the velocity vector with the given amplitude (RanV*Vmax)
      RanVec_FixRadius( RanV, RanVec );
      for (int d=0; d<3; d++)    Vel_AllRank[d][p] = RanVec[d] + params.Cloud_BulkVel[d];

   } // for (long p=0; p<NPar_AllRank; p++)

   Aux_Message( stdout, "   Total enclosed mass within MaxR  = %13.7e\n",  TotM );
   Aux_Message( stdout, "   Particle mass                    = %13.7e\n",  ParM );
   Aux_Message( stdout, "   Maximum mass interpolation error = %13.7e\n",  ErrM_Max );

   // free memory
   delete [] Table_MassProf_r;
   delete [] Table_MassProf_M;

} // FUNCTION : Par_SetEquilibriumIC



//Calculate Mass
//Different Model Type
//-------------------------------------------------------------------------------------------------------
// Function    :
// Description :
//
// Note        :
//
// Parameter   :
//
// Return      :
//-------------------------------------------------------------------------------------------------------
//Plummer
double mass_base_Plummer( double x, void* nothing )
{
   return pow( x, 2 )*pow( 1+x*x, -2.5 );
}



//-------------------------------------------------------------------------------------------------------
// Function    :
// Description :
//
// Note        :
//
// Parameter   :
//
// Return      :
//-------------------------------------------------------------------------------------------------------
//NFW
double mass_base_NFW( double x, void* nothing )
{
   return ( x/pow( 1+x, 2 ) );
}



//-------------------------------------------------------------------------------------------------------
// Function    :
// Description :
//
// Note        :
//
// Parameter   :
//
// Return      :
//-------------------------------------------------------------------------------------------------------
//Burkert
double mass_base_Burkert( double x, void* nothing )
{
   return pow( x, 2 )*( 1/( (1+x)*(1+x*x) ) );
}




//-------------------------------------------------------------------------------------------------------
// Function    :
// Description :
//
// Note        :
//
// Parameter   :
//
// Return      :
//-------------------------------------------------------------------------------------------------------
//Jaffe
double mass_base_Jaffe( double x, void* nothing )
{
   return (x/(1+x));
}



//-------------------------------------------------------------------------------------------------------
// Function    :
// Description :
//
// Note        :
//
// Parameter   :
//
// Return      :
//-------------------------------------------------------------------------------------------------------
//Hernquist
double mass_base_Hernquist( double x, void* nothing )
{
   return ( x/pow( 1+x, 3 ) );
}



//-------------------------------------------------------------------------------------------------------
// Function    :
// Description :
//
// Note        :
//
// Parameter   :
//
// Return      :
//-------------------------------------------------------------------------------------------------------
//Einasto
double mass_base_Einasto( double x, void *Einasto_Power_Factor )
{

   double a = *(double *) Einasto_Power_Factor;

   return pow( x, 2 ) *exp( -pow( x, a ) );
}



//-------------------------------------------------------------------------------------------------------
// Function    :  Set_Mass
// Description :  Calculate the enclosed mass of this cloud within radius r
//
// Note        :
//
// Parameter   :  r : radius
//
// Return      :  Enclosed mass of this cloud within radius r
//-------------------------------------------------------------------------------------------------------
double Par_EquilibriumIC::Set_Mass( double r )
{

   double x = r/params.Cloud_R0;

   if ( strcmp( params.Cloud_Type, "Table" ) == 0 )
   {

      if ( r >= Table_r[params.Cloud_MassProfNBin-1] )
         return Table_Enclosed_Mass[params.Cloud_MassProfNBin-1];

      else if ( r <= Table_r[0] )
         return Table_Enclosed_Mass[0];

      else
         return Mis_InterpolateFromTable( params.Cloud_MassProfNBin, Table_r, Table_Enclosed_Mass, r );
   }

   else
   {

      double result = NULL_REAL;
      double M0 = 4*M_PI*pow(params.Cloud_R0,3)*(params.Cloud_Rho0);

#     ifdef SUPPORT_GSL
      gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);

      double  error;
      gsl_function F;

      if      ( strcmp( params.Cloud_Type, "Plummer"   ) == 0 ) F.function = &mass_base_Plummer;
      else if ( strcmp( params.Cloud_Type, "NFW"       ) == 0 ) F.function = &mass_base_NFW;
      else if ( strcmp( params.Cloud_Type, "Burkert"   ) == 0 ) F.function = &mass_base_Burkert;
      else if ( strcmp( params.Cloud_Type, "Jaffe"     ) == 0 ) F.function = &mass_base_Jaffe;
      else if ( strcmp( params.Cloud_Type, "Hernquist" ) == 0 ) F.function = &mass_base_Hernquist;
      else if ( strcmp( params.Cloud_Type, "Einasto"   ) == 0 )
      {
         F.function = &mass_base_Einasto;
         F.params   = &params.Cloud_Einasto_Power_Factor;
      } // see reference here: https://www.gnu.org/software/gsl/doc/html/roots.html?highlight=gsl_function#providing-the-function-to-solve

      gsl_integration_qag( &F, 0, x, 0, 1e-7, 1000, 1, w, &result, &error );
      gsl_integration_workspace_free (w);
#     endif // #ifdef SUPPORT_GSL

      return result*M0;
   }

} // FUNCTION : SetMass



//-------------------------------------------------------------------------------------------------------
// Function    :  Set_Density
// Description :  Calculate the density of this cloud at radius r
//
// Note        :
//
// Parameter   :  r : radius
//
// Return      :  Density of this cloud at radius r
//-------------------------------------------------------------------------------------------------------
double Par_EquilibriumIC::Set_Density( double x )
{

   if ( strcmp( params.Cloud_Type, "Table" ) == 0 )
   {

      if ( x >= Table_r[params.Cloud_MassProfNBin-1] )
      {
         return Table_Density[params.Cloud_MassProfNBin-1];
      }

      return Mis_InterpolateFromTable( params.Cloud_MassProfNBin, Table_r, Table_Density, x );

   }

   else
   {

      double rho;
      double* nothing;

      if      ( strcmp( params.Cloud_Type, "Plummer"   ) == 0 ) rho = mass_base_Plummer( x, nothing );
      else if ( strcmp( params.Cloud_Type, "NFW"       ) == 0 ) rho = mass_base_NFW( x, nothing );
      else if ( strcmp( params.Cloud_Type, "Burkert"   ) == 0 ) rho = mass_base_Burkert( x, nothing );
      else if ( strcmp( params.Cloud_Type, "Jaffe"     ) == 0 ) rho = mass_base_Jaffe( x, nothing );
      else if ( strcmp( params.Cloud_Type, "Hernquist" ) == 0 ) rho = mass_base_Hernquist( x, nothing );
      else if ( strcmp( params.Cloud_Type, "Einasto"   ) == 0 ) rho = mass_base_Einasto( x, &params.Cloud_Einasto_Power_Factor );

      return rho*pow(x,-2);

   }

} // FUNCTION : Set_Density



//-------------------------------------------------------------------------------------------------------
// Function    :  Set_Velocity
// Description :  Set the velocity of a particle at radius r
//
// Note        :
//
// Parameter   :  r : radius
//
// Return      :  Particle velocity
//-------------------------------------------------------------------------------------------------------
double Par_EquilibriumIC::Set_Velocity( const double x )
{

   double index, sum = 0;
   double psi_per = -potential(x);

   for (int k=0; k<params.Cloud_MassProfNBin; k++)
   {

      if ( psi[k] > psi_per )
      {
         index = k-1;
         break;
      }

      sum += prob_dens[k] *pow( psi_per-psi[k], 0.5 ) *delta;
   }

   double sum_rad, sum_mes=0, par, psi_ass;
   int index_ass = 0;

   sum_rad = Random_Num_Gen->GetValue( 0, 0.0, 1.0 );

   sum_rad *= sum;

   for (int k=0; k<params.Cloud_MassProfNBin; k++)
   {

      if ( sum_mes > sum_rad )
      {
         index_ass = k-1;
         par = (sum_mes-sum_rad)/( prob_dens[index_ass] *pow( psi_per-psi[index_ass], 0.5 ) *delta );
         break;
      }

      sum_mes += prob_dens[k] *pow( psi_per-psi[k], 0.5 ) *delta;

      if ( k == params.Cloud_MassProfNBin-1 )
         index_ass = params.Cloud_MassProfNBin-1;
   }

   psi_ass = psi[index_ass] + delta*par;

   double kim = -2*(psi_ass-psi_per);

   if ( kim < 0.0 )
   {
      return 0;
   }

   double v = pow( kim, 0.5 );

   return v;

} // FUNCTION : Set_Velocity



//-------------------------------------------------------------------------------------------------------
// Function    :
// Description :
//
// Note        :
//
// Parameter   :
//
// Return      :
//-------------------------------------------------------------------------------------------------------
// Solve Eddington's equation
double Par_EquilibriumIC::potential( const double x )
{
   const double r = x*params.Cloud_R0;

   if ( r >= Table_r[params.Cloud_MassProfNBin-1] )
   {
      return Table_Gravity_Potential[params.Cloud_MassProfNBin-1]*Table_r[params.Cloud_MassProfNBin-1]/r;
   }

   if ( r <= Table_r[0] )
   {
      return Table_Gravity_Potential[0];
   }

   return Mis_InterpolateFromTable( params.Cloud_MassProfNBin, Table_r, Table_Gravity_Potential, r );

} // FUNCTION : potential



//-------------------------------------------------------------------------------------------------------
// Function    :
// Description :
//
// Note        :
//
// Parameter   :
//
// Return      :
//-------------------------------------------------------------------------------------------------------
double Par_EquilibriumIC::inverse_psi_to_index( const double psi )
{

   int max = params.Cloud_MassProfNBin-1;
   int min = 0;
   int mid = (max+min)/2;
   double mid_psi;

   while ( true )
   {

      if ( max == min+1 )
      {
         return max;
      }

      mid_psi = -Table_Gravity_Potential[mid] ;

      if ( psi > mid_psi )
      {
         max = mid;
         mid = (max+min)/2;
      }
      else
      {
         min = mid;
         mid = (max+min)/2;
      }
   }

} // FUNCTION : inverse_psi_to_index



//-------------------------------------------------------------------------------------------------------
// Function    :
// Description :
//
// Note        :
//
// Parameter   :
//
// Return      :
//-------------------------------------------------------------------------------------------------------
double Par_EquilibriumIC::integration_eng_base( const double eng )
{

   double min = eng_min;
   double max = eng;
   int    num = 1000;

   double dx = (max-min)/num;

   double result = 0;

   for (int i=0; i<num; i++)
   {

      double psi_l = min+i*dx;
      double psi_r = min+(i+1)*dx;

      int ind0 = inverse_psi_to_index(min+(i+0.5)*dx);

      if ( i == num-1 )
         result += -2*Table_dRho_dx[ind0]*( pow(eng-psi_l,0.5) );
      else
         result += -2*Table_dRho_dx[ind0]*( pow(eng-psi_l,0.5) - pow(eng-psi_r,0.5) );

   }

   return result;

} // FUNCTION : integration_eng_base



// Initialize physical parameter tables
//-------------------------------------------------------------------------------------------------------
// Function    :  Init_Mass
// Description :  Calculate the table of enclosed masses (vs. radius) of the cloud
//                by giving a known analytical density function of the cloud
//
// Note        :
//
// Parameter   :
//
// Return      :
//-------------------------------------------------------------------------------------------------------
void Par_EquilibriumIC::Init_Mass()
{

   double dr = params.Cloud_MaxR / (params.Cloud_MassProfNBin-1);

   //Radius & Mass
   for (int b=0; b<params.Cloud_MassProfNBin; b++)
   {
      Table_r[b] = dr*b;
      Table_Enclosed_Mass[b] = Set_Mass( Table_r[b] );
   }

   //Rho
   for (int b=1; b<params.Cloud_MassProfNBin; b++)
   {
      double x = dr*b/params.Cloud_R0;
      Table_Density[b] = Set_Density(x);
   }

   Table_Density[0] = Table_Density[1];

   //Rhodr
   Table_dRho_dr[0] = (Table_Density[1]-Table_Density[0])/dr;
   for (int b=1; b<params.Cloud_MassProfNBin-1; b++)
   {
      int num=3;

      if      ( b == 0 )
         Table_dRho_dr[b] = Slope_LinearRegression( Table_r, Table_Density, 0, num/2+1 );

      else if ( b == 1 )
         Table_dRho_dr[b] = Slope_LinearRegression( Table_r, Table_Density, 0, num/2+2 );

      else if ( b == params.Cloud_MassProfNBin-2 )
         Table_dRho_dr[b] = Slope_LinearRegression( Table_r, Table_Density, params.Cloud_MassProfNBin-num/2-1, num/2+1 );

      else
         Table_dRho_dr[b] = Slope_LinearRegression( Table_r, Table_Density, b-num/2, num+1 );

   }

   Table_dRho_dr[params.Cloud_MassProfNBin-1] = Table_dRho_dr[params.Cloud_MassProfNBin-2];

} // FUNCTION : Init_Mass



//-------------------------------------------------------------------------------------------------------
// Function    :  Init_Pot
// Description :  Calculate the table of potential (vs. radius) of the cloud
//                by giving a known analytical density function of the cloud
//
// Note        :
//
// Parameter   :
//
// Return      :
//-------------------------------------------------------------------------------------------------------
void Par_EquilibriumIC::Init_Pot()
{

   double dr = params.Cloud_MaxR / (params.Cloud_MassProfNBin-1);

   Table_Gravity_Field[0] =0;

   for (int b=1; b<params.Cloud_MassProfNBin; b++)
   {
      Table_Gravity_Field[b] = -NEWTON_G*Table_Enclosed_Mass[b]/pow( Table_r[b], 2 );
   }

   //Pot
   Table_Gravity_Potential[params.Cloud_MassProfNBin-1] = -NEWTON_G*Table_Enclosed_Mass[params.Cloud_MassProfNBin-1]/Table_r[params.Cloud_MassProfNBin-1];
   eng_min = -Table_Gravity_Potential[params.Cloud_MassProfNBin-1];

   for (int b=params.Cloud_MassProfNBin-2; b>0; b--)
   {
      Table_Gravity_Potential[b] = Table_Gravity_Potential[b+1] + Table_Gravity_Field[b]*dr;
   }

   Table_Gravity_Potential[0] = Table_Gravity_Potential[1];

   //derho_overdx
   for (int b=0; b<params.Cloud_MassProfNBin; b++)
   {
      Table_dRho_dx[b] = -Table_dRho_dr[b]/Table_Gravity_Field[b];
   }

} // FUNCTION : Init_Pot



// Initialization through loading a file of table
//-------------------------------------------------------------------------------------------------------
// Function    :  Init_Mass_Table
// Description :  Calculate the table of potential (vs. radius) of the cloud
//                by loading a file of table
//
// Note        :
//
// Parameter   :
//
// Return      :
//-------------------------------------------------------------------------------------------------------
void Par_EquilibriumIC::Init_Mass_Table()
{

   //Mass
   Table_Enclosed_Mass[0] = 0;
   double rho,dr,r;

   for (int b=1; b<params.Cloud_MassProfNBin; b++)
   {
      rho = (Table_Density[b] + Table_Density[b-1])/2;
      dr  = Table_r[b] - Table_r[b-1];
      r   = (Table_r[b] + Table_r[b-1])/2;

      Table_Enclosed_Mass[b] = Table_Enclosed_Mass[b-1] + 4*M_PI*pow(r,2) *rho * dr;
   }

   //Rhodr
   Table_dRho_dr[0] = (Table_Density[1]-Table_Density[0])/(Table_r[1]-Table_r[0]);

   for (int b=1; b<params.Cloud_MassProfNBin-1; b++)
   {
      int num = 3;
      if      ( b == 0 )
         Table_dRho_dr[b] = Slope_LinearRegression( Table_r, Table_Density, 0, num/2+1 );

      else if ( b == 1 )
         Table_dRho_dr[b] = Slope_LinearRegression( Table_r, Table_Density, 0, num/2+2 );

      else if ( b == params.Cloud_MassProfNBin-2 )
         Table_dRho_dr[b] = Slope_LinearRegression( Table_r, Table_Density, params.Cloud_MassProfNBin-num/2-1, num/2+1 );
      else
         Table_dRho_dr[b] = Slope_LinearRegression( Table_r, Table_Density, b-num/2, num+1 );

   }

   Table_dRho_dr[params.Cloud_MassProfNBin-1] = Table_dRho_dr[params.Cloud_MassProfNBin-2];

} // FUNCTION : Init_Mass_Table



//-------------------------------------------------------------------------------------------------------
// Function    :  Init_Pot_Table
// Description :  Calculate the table of potential (vs. radius) of the cloud
//                by loading a file of table
//
// Note        :
//
// Parameter   :
//
// Return      :
//-------------------------------------------------------------------------------------------------------
void Par_EquilibriumIC::Init_Pot_Table()
{

   Table_Gravity_Field[0] = 0;
   for (int b=1; b<params.Cloud_MassProfNBin; b++)
   {
      Table_Gravity_Field[b] = -NEWTON_G*Table_Enclosed_Mass[b]/pow( Table_r[b], 2 );
   }

   //Pot
   Table_Gravity_Potential[params.Cloud_MassProfNBin-1] = -NEWTON_G*Table_Enclosed_Mass[params.Cloud_MassProfNBin-1]/Table_r[params.Cloud_MassProfNBin-1];

   eng_min = -Table_Gravity_Potential[params.Cloud_MassProfNBin-1];

   for (int b=params.Cloud_MassProfNBin-2; b>0; b--)
   {
      double dr = Table_r[b+1] - Table_r[b];

      Table_Gravity_Potential[b] = Table_Gravity_Potential[b+1] + Table_Gravity_Field[b]*dr;
   }

   Table_Gravity_Potential[0] = Table_Gravity_Potential[1];

   //derho_overdx
   for (int b=0; b<params.Cloud_MassProfNBin; b++)
   {
      Table_dRho_dx[b] = -Table_dRho_dr[b]/Table_Gravity_Field[b];
   }

} // FUNCTION : Init_Pot_Table



//-------------------------------------------------------------------------------------------------------
// Function    :  Init_Prob_Dens
// Description :  Calculate the probability density function of particles' velocities
//
// Note        :
//
// Parameter   :
//
// Return      :
//-------------------------------------------------------------------------------------------------------
void Par_EquilibriumIC::Init_Prob_Dens()
{

   double min, max;
   min   = -Table_Gravity_Potential[params.Cloud_MassProfNBin-1];
   max   = -Table_Gravity_Potential[1];
   delta = (max-min)/params.Cloud_MassProfNBin;

   double eng = min;

   for (int k =0; k<params.Cloud_MassProfNBin; k++)
   {

      psi[k] = eng;

      int_prob_dens[k] = integration_eng_base(eng);

      eng += delta;

   }

   for (int k =0; k<params.Cloud_MassProfNBin; k++)
   {

      if      ( k == 0 )                           prob_dens[k] = Slope_LinearRegression( psi, int_prob_dens, k,   5 );
      else if ( k == 1 )                           prob_dens[k] = Slope_LinearRegression( psi, int_prob_dens, k-1, 5 );
      else if ( k == params.Cloud_MassProfNBin-2 ) prob_dens[k] = Slope_LinearRegression( psi, int_prob_dens, k-3, 5 );
      else if ( k == params.Cloud_MassProfNBin-1 ) prob_dens[k] = Slope_LinearRegression( psi, int_prob_dens, k-4, 5 );
      else                                         prob_dens[k] = Slope_LinearRegression( psi, int_prob_dens, k-2, 5 );

      if ( prob_dens[k] < 0 )                      prob_dens[k] = 0;

   }

   smooth_all( prob_dens, 0, params.Cloud_MassProfNBin );

} // FUNCTION : Init_Prob_Dens



//-------------------------------------------------------------------------------------------------------
// Function    :  Add_Ext_Pot
// Description :  Add exteranl potential through loading a file of table
//
// Note        :  Only will be activated if AddExtPot is turned on
//
// Parameter   :
//
// Return      :
//-------------------------------------------------------------------------------------------------------
void Par_EquilibriumIC::Add_Ext_Pot()
{

   if ( ! bool(params.AddExtPot) )  return;

   Aux_Message( stdout, "Loading External Potential Table...\n");

   int Tcol_r[1]={0};
   int Tcol_Pot[1]={1};

   double* Radius=NULL;
   double* Ext_Pot=NULL;

   int Row_r_Table;
   Row_r_Table = Aux_LoadTable( Radius, params.ExtPot_Table_Name, 1, Tcol_r, true, true );

   int Row_Ext_Pot_Table;
   Row_Ext_Pot_Table = Aux_LoadTable( Ext_Pot, params.ExtPot_Table_Name, 1, Tcol_Pot, true, true );

   Aux_Message( stdout, "Loading Ext_Pot Profile Table: %s\n", params.ExtPot_Table_Name );

   if ( Row_r_Table != Row_Ext_Pot_Table )
      Aux_Error( ERROR_INFO, "Ext_Pot row number is not equal to radius row number in the profile file !! Please check this file.\n" );

   if ( params.Cloud_MassProfNBin != Row_r_Table )
      Aux_Error( ERROR_INFO, "Cloud_MassProfNBin is not equal to the row number in profile file !!\n" );

   for (int i=0; i<params.Cloud_MassProfNBin; i++)
   {
      Table_Gravity_Potential[i] += Ext_Pot[i];
   }

} // FUNCTION : Add_Ext_Pot



// Auxiliary functions
//-------------------------------------------------------------------------------------------------------
// Function    :  GetParams
// Description :  Get the parameters from a file
//
// Note        :
//
// Parameter   :  filename   : file name of the input file
//                keyword    : the name of the parameter
//                para_num   : number of input parameters
//                containaer : container of return values
//
// Return      :  containaer
//
//-------------------------------------------------------------------------------------------------------
int Par_EquilibriumIC::GetParams( const char *filename, const char *keyword, const int para_num,
                                  const char *para_type, vector <string> &container )
{

   fstream file;
   file.open( filename, ios::in );

   string line;
   string para;

   if ( !file )
   {
      Aux_Error( ERROR_INFO, "Failed to open file : %s", filename );
   }
   else
   {
      do
      {
         getline( file, line );

         istringstream templine( line );

         string first_word;

         getline( templine, first_word, ' ' );

         if ( strcmp( first_word.c_str(), keyword ) == 0 )
         {
            for (int i=0; i<para_num; i++)
            {
               do
               {
                  getline( templine, para, ' ' );

               } while( para.length() == 0 );

               if ( strcmp( para_type, "string" ) == 0 )
                  container.push_back(para);
               else
                  return atoi( para.c_str() );
            }

            break;
         }

      } while( !file.eof() );
   }

   file.close();

   return 0;

} // FUNCTION : GetParams



//-------------------------------------------------------------------------------------------------------
// Function    :  RanVec_FixRadius
// Description :  Compute a random 3D vector with a fixed radius
//
// Note        :  Uniformly random sample in theta and phi does NOT give a uniformly random sample in 3D space
//                --> Uniformly random sample in a 3D sphere and then normalize all vectors to the given radius
//
// Parameter   :  r      : Input radius
//                RanVec : Array to store the random 3D vector
//
// Return      :  RanVec
//-------------------------------------------------------------------------------------------------------
void Par_EquilibriumIC::RanVec_FixRadius( const double r, double RanVec[] )
{

   double Norm, RanR2;

   do
   {
      RanR2 = 0.0;

      for (int d=0; d<3; d++)
      {
         RanVec[d]  = Random_Num_Gen->GetValue( 0, -1.0, +1.0 );
         RanR2     += SQR( RanVec[d] );
      }

   } while ( RanR2 > 1.0 );

   Norm = r / sqrt( RanR2 );

   for (int d=0; d<3; d++)    RanVec[d] *= Norm;

} // FUNCTION : RanVec_FixRadius



//-------------------------------------------------------------------------------------------------------
// Function    :
// Description :
//
// Note        :
//
// Parameter   :
//
// Return      :
//-------------------------------------------------------------------------------------------------------
void Par_EquilibriumIC::mask( double* x, int start, int fin )
{

   double standard = 3;

   for (int j=start; j<fin; j++)
   {

      bool flag = 0;

      for (int k=start; k<fin; k++)
      {
         if ( x[k] != 0 )
         {
            if ( fabs( x[j]/x[k] ) > standard ) flag = 1;
         }
      }

      if ( flag )
         x[j] = 0;
   }

} // FUNCTION : mask



//-------------------------------------------------------------------------------------------------------
// Function    :
// Description :
//
// Note        :
//
// Parameter   :
//
// Return      :
//-------------------------------------------------------------------------------------------------------
void Par_EquilibriumIC::add_num( double* x, int start, int fin )
{

   double sum = 0;
   int    num = 0;

   for (int j=start; j<fin; j++)
   {

      if ( x[j] != 0 )
      {
         sum += x[j];
         num++;
      }
   }

   if ( num != 0 )
   {
      double ave_x = sum/num;

      for (int j=start; j<fin; j++)
      {
         if ( x[j] == 0 ) x[j] = ave_x;
      }
   }

} // FUNCTION : add_num



//-------------------------------------------------------------------------------------------------------
// Function    :
// Description :
//
// Note        :
//
// Parameter   :
//
// Return      :
//-------------------------------------------------------------------------------------------------------
void Par_EquilibriumIC::smooth_all( double* x, int start, int fin )
{

   int num = 10;
   for (int k=start; k<fin-num+1; k++)
   {
      mask( x, k, k+num );
   }

   for (int k=start; k<fin-num+1; k++)
   {
      add_num( x, k, k+num );
   }

} // FUNCTION : smooth_all



//-------------------------------------------------------------------------------------------------------
// Function    : ArrayCovariance
// Description : Get the covariance between two arrays
//
// Note        : if x==y, then the covariance is the variance of x
//
// Parameter   : array_x     : array of x data
//               array_y     : array of y data
//               index_start : the first index in the array for the linear regression
//               n_elements  : number of elements for the linear regression
//
// Return      : covariance_xy
//-------------------------------------------------------------------------------------------------------
double Par_EquilibriumIC::ArrayCovariance( const double* array_x, const double* array_y,
                                           const int index_start, const int n_elements )
{
   const double normalized_factor = 1.0/n_elements;

   // average
   double average_x = 0.0;
   double average_y = 0.0;

   for (int i=index_start; i<index_start+n_elements; i++)
   {
      average_x += array_x[i];
      average_y += array_y[i];
   }

   average_x *= normalized_factor;
   average_y *= normalized_factor;


   // covariance
   double covariance_xy = 0.0;

   for (int i=index_start; i<index_start+n_elements; i++)
   {
      covariance_xy += (array_x[i]-average_x)*(array_y[i]-average_y);
   }

   covariance_xy *= normalized_factor;


   return covariance_xy;

} // FUNCTION : ArrayCovariance



//-------------------------------------------------------------------------------------------------------
// Function    : Slope_LinearRegression
// Description : Get the slope of y-x using linear regression
//
// Note        :
//
// Parameter   : array_x     : array of x data
//               array_y     : array of y data
//               index_start : the first index in the array for the linear regression
//               n_elements  : number of elements for the linear regression
//
// Return      : slope
//-------------------------------------------------------------------------------------------------------
double Par_EquilibriumIC::Slope_LinearRegression( const double* array_x, const double* array_y,
                                                  const int index_start, const int n_elements )
{

   const double variance_x    = ArrayCovariance( array_x, array_x, index_start, n_elements );
   const double covariance_xy = ArrayCovariance( array_x, array_y, index_start, n_elements );

   return covariance_xy/variance_x;

} // FUNCTION : Slope_LinearRegression



#endif // #ifdef PARTICLE
