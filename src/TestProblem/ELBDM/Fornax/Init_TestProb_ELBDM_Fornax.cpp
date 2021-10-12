#include "GAMER.h"
#include "TestProb.h"

#include"Particle_IC_Constructor.h"
#include"string"

// problem-specific global variables
// =======================================================================================
static int      Fornax_N;                               // total number of Fornaxs
static int      Fornax_RSeed;                           // random seed for setting Fornax_Center[]
                                                         //    (<0 --> hard coding each Fornax)
static double   Fornax_CoreRadiusAll;                   // core radius for all Fornaxs
                                                         //    (<=0.0 --> hard coding each Fornax)
static double   Fornax_EmptyRegion;                     // Fornax-free region from the boundary
                                                         //    (useful only when Fornax_RSeed>=0)
static char     Fornax_DensProf_Filename[MAX_STRING];   // filename of the reference Fornax density profile

static int      Fornax_DensProf_NBin;                   // number of radial bins of the Fornax density profile
static double  *Fornax_DensProf   = NULL;               // Fornax density profile [radius/density]
static double  *Fornax_CoreRadius = NULL;               // core radius of each Fornax
static double (*Fornax_Center)[3] = NULL;               // center coordinates of each Fornax
static double  *Fornax_ScaleL     = NULL;               // L/D: length/density scale factors of each Fornax
                                                         //      (defined as the ratio between the core radii/peak
                                                         //      density of the target and reference Fornax profiles)
static double  *Fornax_ScaleD     = NULL;
// =======================================================================================

//GC particle ic constructor
#  ifdef PARTICLE
void Par_Init_ByFunction_Fornax( const long NPar_ThisRank, const long NPar_AllRank,
                                  real *ParMass, real *ParPosX, real *ParPosY, real *ParPosZ,
                                  real *ParVelX, real *ParVelY, real *ParVelZ, real *ParTime,
                                  real *AllAttribute[PAR_NATT_TOTAL] );
extern Particle_IC_Constructor constructor_Fornax;//First Declared in Par_Init_ByFunction_Fornax
#  endif

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

//#  ifdef PARTICLE
//   Aux_Error( ERROR_INFO, "PARTICLE must be disabled !!\n" );
//#  endif

#  ifdef PARTICLE
   if ( OPT__INIT == INIT_BY_FUNCTION  &&  amr->Par->Init != PAR_INIT_BY_FUNCTION )
      Aux_Error( ERROR_INFO, "please set PAR_INIT = 1 (by FUNCTION) !!\n" );
#  endif

#  ifdef GRAVITY
   if ( OPT__BC_POT != BC_POT_ISOLATED )
      Aux_Error( ERROR_INFO, "must adopt isolated BC for gravity --> reset OPT__BC_POT !!\n" );
#  endif


// warnings
   if ( MPI_Rank == 0 )
   {
      if ( !OPT__INIT_RESTRICT )
         Aux_Message( stderr, "WARNING : it's recommended to enable OPT__INIT_RESTRICT !!\n" );
   }
// additional warnings
   if ( MPI_Rank == 0 )
   {
      for (int f=0; f<6; f++)
      if ( OPT__BC_FLU[f] == BC_FLU_PERIODIC )
         Aux_Message( stderr, "WARNING : periodic BC for fluid is not recommended for this test !!\n" );

#     ifdef GRAVITY
      if ( OPT__BC_POT == BC_POT_PERIODIC )
         Aux_Message( stderr, "WARNING : periodic BC for gravity is not recommended for this test !!\n" );
#     endif
   } // if ( MPI_Rank == 0 )


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
   ReadPara->Add( "Fornax_N",                 &Fornax_N,                 -1,             1,                NoMax_int         );
   ReadPara->Add( "Fornax_RSeed",             &Fornax_RSeed,              0,             NoMin_int,        NoMax_int         );
   ReadPara->Add( "Fornax_CoreRadiusAll",     &Fornax_CoreRadiusAll,      NoDef_double,  NoMin_double,     NoMax_double      );
   ReadPara->Add( "Fornax_EmptyRegion",       &Fornax_EmptyRegion,        0.0,           NoMin_double,     NoMax_double      );
   ReadPara->Add( "Fornax_DensProf_Filename",  Fornax_DensProf_Filename,  Useless_str,   Useless_str,      Useless_str       );

   ReadPara->Read( FileName );

   delete ReadPara;

// (1-2) set the default values

// (1-3) check the runtime parameters
   if ( Fornax_RSeed >= 0  &&  Fornax_EmptyRegion < 0.0 )
      Aux_Error( ERROR_INFO, "Fornax_EmptyRegion (%14.7e) < 0.0 !!\n", Fornax_EmptyRegion );

   if ( Fornax_CoreRadiusAll == NoDef_double )
      Aux_Error( ERROR_INFO, "Runtime parameter \"Fornax_CoreRadiusAll\" is not set !!\n" );


// (2) set the problem-specific derived parameters
// (2-1) allocate memory
   Fornax_CoreRadius = new double [Fornax_N];
   Fornax_Center     = new double [Fornax_N][3];
   Fornax_ScaleL     = new double [Fornax_N];
   Fornax_ScaleD     = new double [Fornax_N];

// (2-2) Fornax core radii
   if ( Fornax_CoreRadiusAll > 0.0 )
   {
      for (int t=0; t<Fornax_N; t++)  Fornax_CoreRadius[t] = Fornax_CoreRadiusAll;
   }

   else
   {
//    for Fornax_CoreRadiusAll <= 0.0, comment out the following line and hard code the core radius of each Fornax
      Aux_Error( ERROR_INFO, "for Fornax_CoreRadiusAll <= 0.0, please comment out this error check and hard code "
                             "the core radius of each Fornax !!\n" );
//    for (int t=0; t<Fornax_N; t++)  Fornax_CoreRadius[t] = XXX;
   }

// (2-3) Fornax centers
   if ( Fornax_RSeed >= 0 )
   {
      const double Coord_Min[3] = { Fornax_EmptyRegion, Fornax_EmptyRegion, Fornax_EmptyRegion };
      const double Coord_Max[3] = { amr->BoxSize[0] - Fornax_EmptyRegion,
                                    amr->BoxSize[1] - Fornax_EmptyRegion,
                                    amr->BoxSize[2] - Fornax_EmptyRegion };
      srand( Fornax_RSeed );

      for (int t=0; t<Fornax_N; t++)
      for (int d=0; d<3; d++)
         Fornax_Center[t][d] = ( (double)rand()/RAND_MAX )*(Coord_Max[d]-Coord_Min[d]) + Coord_Min[d];
   }

   else
   {
      if ( Fornax_N == 1 )
      {
         for (int d=0; d<3; d++)    Fornax_Center[0][d] = 0.5*amr->BoxSize[d];
      }

      else
      {
//       for Fornax_RSeed<0, comment out the following line and hard code the center of each Fornax
         Aux_Error( ERROR_INFO, "for Fornax_RSeed < 0 and Fornax_N > 1, please comment out this error check and hard code "
                                "the center of each Fornax !!\n" );

         /*
         for (int t=0; t<Fornax_N; t++)
         for (int d=0; d<3; d++)          Fornax_Center[t][d] = XXX;
         */
      }
   } // if ( Fornax_RSeed >= 0 ) ... else ...


// (3) load the reference Fornax density profile and evaluate the scale factors
   if ( OPT__INIT != INIT_BY_RESTART )
   {
//    load the reference profile
      const bool RowMajor_No  = false;    // load data into the column-major order
      const bool AllocMem_Yes = true;     // allocate memory for Fornax_DensProf
      const int  NCol         = 2;        // total number of columns to load
      const int  Col[NCol]    = {0, 1};   // target columns: (radius, density)

      Fornax_DensProf_NBin = Aux_LoadTable( Fornax_DensProf, Fornax_DensProf_Filename, NCol, Col, RowMajor_No, AllocMem_Yes );


//    get the core radius of the reference profile
      const double *RadiusRef = Fornax_DensProf + 0*Fornax_DensProf_NBin;
      const double *DensRef   = Fornax_DensProf + 1*Fornax_DensProf_NBin;
      const double  DensCore  = 0.5*DensRef[0];   // define core radius as the half-density radius

      double CoreRadiusRef = NULL_REAL;

      for (int b=1; b<Fornax_DensProf_NBin-1; b++)
      {
         if ( DensRef[b] >= DensCore  &&  DensRef[b+1] <= DensCore )
         {
            CoreRadiusRef = 0.5*( RadiusRef[b] + RadiusRef[b+1] );
            break;
         }
      }

      if ( CoreRadiusRef == NULL_REAL )
         Aux_Error( ERROR_INFO, "cannot determine the reference core radius !!\n" );


//    evaluate the scale factors of each Fornax
      for (int t=0; t<Fornax_N; t++)
      {
         Fornax_ScaleL[t] = Fornax_CoreRadius[t] / CoreRadiusRef;
         Fornax_ScaleD[t] = 1.0 / ( 4.0*M_PI*NEWTON_G*SQR(ELBDM_ETA)*POW4(Fornax_ScaleL[t]) );
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
      Aux_Message( stdout, "  total number of Fornaxs                  = %d\n",     Fornax_N                  );
      Aux_Message( stdout, "  random seed for setting the center coord. = %d\n",     Fornax_RSeed              );
      Aux_Message( stdout, "  size of the Fornax-free zone             = %13.7e\n", Fornax_EmptyRegion        );
      Aux_Message( stdout, "  density profile filename                  = %s\n",     Fornax_DensProf_Filename  );
      Aux_Message( stdout, "  number of bins of the density profile     = %d\n",     Fornax_DensProf_NBin      );
      Aux_Message( stdout, "\n" );
      Aux_Message( stdout, "  Fornax info:\n" );
      Aux_Message( stdout, "  %7s  %13s  %13s  %13s  %13s  %13s  %13s\n",
                   "ID", "CoreRadius", "ScaleL", "ScaleD", "Center_X", "Center_Y", "Center_Z" );
      for (int t=0; t<Fornax_N; t++)
      Aux_Message( stdout, "  %7d  %13.6e  %13.6e  %13.6e  %13.6e  %13.6e  %13.6e\n",
                   t, Fornax_CoreRadius[t], Fornax_ScaleL[t], Fornax_ScaleD[t],
                   Fornax_Center[t][0], Fornax_Center[t][1], Fornax_Center[t][2] );
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

   const double *Table_Radius  = Fornax_DensProf + 0*Fornax_DensProf_NBin;  // radius
   const double *Table_Density = Fornax_DensProf + 1*Fornax_DensProf_NBin;  // density

   double r_tar, r_ref, dens_ref;


// initialize density as zero since there may be multiple Fornaxs
   fluid[DENS] = 0.0;

// loop over all Fornaxs to get the total density
   for (int t=0; t<Fornax_N; t++)
   {
      r_tar = sqrt( SQR(x-Fornax_Center[t][0]) + SQR(y-Fornax_Center[t][1]) + SQR(z-Fornax_Center[t][2]) );

//    rescale radius (target radius --> reference radius)
      r_ref = r_tar / Fornax_ScaleL[t];

//    linear interpolation
      dens_ref = Mis_InterpolateFromTable( Fornax_DensProf_NBin, Table_Radius, Table_Density, r_ref );

      if ( dens_ref == NULL_REAL )
      {
         if      ( r_ref <  Table_Radius[0] )
            dens_ref = Table_Density[0];

         else if ( r_ref >= Table_Radius[Fornax_DensProf_NBin-1] )
            dens_ref = Table_Density[Fornax_DensProf_NBin-1];

         else
            Aux_Error( ERROR_INFO, "interpolation failed at radius %13.7e (min/max radius = %13.7e/%13.7e) !!\n",
                       r_ref, Table_Radius[0], Table_Radius[Fornax_DensProf_NBin-1] );
      }

//    rescale density (reference density --> target density) and add to the fluid array
      fluid[DENS] += dens_ref*Fornax_ScaleD[t];
   } // for (int t=0; t<Fornax_N; t++)


// set the real and imaginary parts
   fluid[REAL] = sqrt( fluid[DENS] );
   fluid[IMAG] = 0.0;                  // imaginary part is always zero --> no initial velocity

} // FUNCTION : SetGridIC



//-------------------------------------------------------------------------------------------------------
// Function    :  End_Fornax
// Description :  Free memory before terminating the program
//
// Note        :  1. Linked to the function pointer "End_User_Ptr" to replace "End_User()"
//
// Parameter   :  None
//-------------------------------------------------------------------------------------------------------
void End_Fornax()
{

   delete [] Fornax_DensProf;
   delete [] Fornax_CoreRadius;
   delete [] Fornax_Center;
   delete [] Fornax_ScaleL;
   delete [] Fornax_ScaleD;

} // FUNCTION : End_Fornax



//-------------------------------------------------------------------------------------------------------
// Function    :  BC_Fornax
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
void BC_Fornax( real fluid[], const double x, const double y, const double z, const double Time,
         const int lv, double AuxArray[] )
{

   fluid[REAL] = (real)0.0;
   fluid[IMAG] = (real)0.0;
   fluid[DENS] = (real)0.0;

} // FUNCTION : BC_Fornax
#endif // #if ( MODEL == ELBDM  &&  defined GRAVITY )



//-------------------------------------------------------------------------------------------------------
// Function    :  Init_TestProb_ELBDM_Fornax
// Description :  Test problem initializer
//
// Note        :  None
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Init_TestProb_ELBDM_Fornax()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );


// validate the compilation flags and runtime parameters
   Validate();


#  if ( MODEL == ELBDM  &&  defined GRAVITY )
// set the problem-specific runtime parameters
   SetParameter();

   //User Input for Particle IC Constructor
   int Fornax_num = 1;
   vector <string> TestProb_FileName;     //Test problem input parameter filenames
   vector <string> TypeName;                               //Type of Fornax 
   vector <string> Profile_FileName;
   TestProb_FileName.push_back("Input__TestProb_GC");
   TypeName.push_back("UNKNOWN");
   Profile_FileName.push_back("profile_GC.txt");

   
#  ifdef PARTICLE
   constructor_Fornax.construct_ic(Fornax_num,TestProb_FileName,TypeName,Profile_FileName);
#  endif


   Init_Function_User_Ptr = SetGridIC;
   BC_User_Ptr            = BC_Fornax;
   End_User_Ptr           = End_Fornax;
#  ifdef PARTICLE
   Par_Init_ByFunction_Ptr = Par_Init_ByFunction_Fornax;
#  endif

#  endif // #if ( MODEL == ELBDM  &&  defined GRAVITY )


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : Init_TestProb_ELBDM_Fornax
