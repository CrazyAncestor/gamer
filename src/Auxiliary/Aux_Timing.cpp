#include "Copyright.h"
#include "GAMER.h"

#ifdef TIMING

void Timing__EvolveLevel( const char FileName[], const double Time_LB_Main[][3] );
#ifdef TIMING_SOLVER
void Timing__Solver( const char FileName[] );
#endif


// global timing variables
// ----------------------------------------------------------
extern Timer_t *Timer_Main[7];
extern Timer_t *Timer_MPI[3];
extern Timer_t *Timer_Flu_Advance[NLEVEL];
extern Timer_t *Timer_Gra_Advance[NLEVEL];
extern Timer_t *Timer_Che_Advance[NLEVEL];
extern Timer_t *Timer_FixUp      [NLEVEL];
extern Timer_t *Timer_Flag       [NLEVEL];
extern Timer_t *Timer_Refine     [NLEVEL];
extern Timer_t *Timer_GetBuf     [NLEVEL][9];
extern Timer_t *Timer_Lv         [NLEVEL];
extern Timer_t *Timer_Par_Update [NLEVEL][3];
extern Timer_t *Timer_Par_2Sib   [NLEVEL];
extern Timer_t *Timer_Par_2Son   [NLEVEL];
extern Timer_t *Timer_Par_Collect[NLEVEL];
extern Timer_t *Timer_Par_MPI    [NLEVEL][6];

#ifdef TIMING_SOLVER
extern Timer_t *Timer_Pre         [NLEVEL][5];
extern Timer_t *Timer_Sol         [NLEVEL][5];
extern Timer_t *Timer_Clo         [NLEVEL][5];
extern Timer_t *Timer_Poi_PreRho  [NLEVEL];
extern Timer_t *Timer_Poi_PreFlu  [NLEVEL];
extern Timer_t *Timer_Poi_PrePot_C[NLEVEL];
extern Timer_t *Timer_Poi_PrePot_F[NLEVEL];
#endif

// accumulated timing results
static double Flu_Acc     [3] = { 0.0, 0.0, 0.0 };
static double Gra_Acc     [3] = { 0.0, 0.0, 0.0 };
static double Che_Acc     [3] = { 0.0, 0.0, 0.0 };
static double FixUp_Acc   [3] = { 0.0, 0.0, 0.0 };
static double Flag_Acc    [3] = { 0.0, 0.0, 0.0 };
static double Refine_Acc  [3] = { 0.0, 0.0, 0.0 };
static double MPI_Grid_Acc[3] = { 0.0, 0.0, 0.0 };
static double dt_Acc      [3] = { 0.0, 0.0, 0.0 };
static double Aux_Acc     [3] = { 0.0, 0.0, 0.0 };
static double Output_Acc  [3] = { 0.0, 0.0, 0.0 };
static double LB_Acc      [3] = { 0.0, 0.0, 0.0 };
static double Corr_Acc    [3] = { 0.0, 0.0, 0.0 };
static double Par_Acc     [3] = { 0.0, 0.0, 0.0 };
static double Sum_Acc     [3] = { 0.0, 0.0, 0.0 };
static double MPI_Par_Acc [3] = { 0.0, 0.0, 0.0 };




//-------------------------------------------------------------------------------------------------------
// Function    :  Aux_CreateTimer
// Description :  Create simulation timers by using the "Timer" structure
//-------------------------------------------------------------------------------------------------------
void Aux_CreateTimer()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "Aux_CreateTimer ... " );


   for (int t=0; t<7; t++)    Timer_Main[t] = new Timer_t( 1 );

   if ( OPT__TIMING_MPI )
   for (int t=0; t<3; t++)    Timer_MPI [t] = new Timer_t( 1 );

   for (int lv=0; lv<NLEVEL; lv++)
   {
#     ifdef INDIVIDUAL_TIMESTEP
      const int N = 1<<(1+lv);
#     else
      const int N = 1;
#     endif

      Timer_Flu_Advance[lv] = new Timer_t( N );
      Timer_Gra_Advance[lv] = new Timer_t( N );
      Timer_Che_Advance[lv] = new Timer_t( N );
      Timer_FixUp      [lv] = new Timer_t( N );
      Timer_Flag       [lv] = new Timer_t( N );
      Timer_Refine     [lv] = new Timer_t( N );
      Timer_Lv         [lv] = new Timer_t( N );
      for (int t=0; t<9; t++)    Timer_GetBuf    [lv][t] = new Timer_t( N );
      for (int t=0; t<3; t++)    Timer_Par_Update[lv][t] = new Timer_t( N );
      Timer_Par_2Sib   [lv] = new Timer_t( N );
      Timer_Par_2Son   [lv] = new Timer_t( N );
      Timer_Par_Collect[lv] = new Timer_t( N );
      for (int t=0; t<6; t++)    Timer_Par_MPI   [lv][t] = new Timer_t( N );

#     ifdef TIMING_SOLVER
      for (int v=0; v<5; v++)
      {
         Timer_Pre[lv][v]    = new Timer_t( 1 );
         Timer_Sol[lv][v]    = new Timer_t( 1 );
         Timer_Clo[lv][v]    = new Timer_t( 1 );
      }
      Timer_Poi_PreRho  [lv] = new Timer_t( 1 );
      Timer_Poi_PreFlu  [lv] = new Timer_t( 1 );
      Timer_Poi_PrePot_C[lv] = new Timer_t( 1 );
      Timer_Poi_PrePot_F[lv] = new Timer_t( 1 );
#     endif // #ifdef TIMING_SOLVER
   }


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "done\n" );

} // FUNCTION : Aux_CreateTimer



//-------------------------------------------------------------------------------------------------------
// Function    :  Aux_DeleteTimer
// Description :  Delete simulation timers by using the "Timer" structure
//-------------------------------------------------------------------------------------------------------
void Aux_DeleteTimer()
{

   for (int t=0; t<7; t++)    delete Timer_Main[t];

   if ( OPT__TIMING_MPI)
   for (int t=0; t<3; t++)    delete Timer_MPI [t];

   for (int lv=0; lv<NLEVEL; lv++)
   {
      delete Timer_Flu_Advance[lv];
      delete Timer_Gra_Advance[lv];
      delete Timer_Che_Advance[lv];
      delete Timer_FixUp      [lv];
      delete Timer_Flag       [lv];
      delete Timer_Refine     [lv];
      delete Timer_Lv         [lv];
      for (int t=0; t<9; t++)    delete Timer_GetBuf    [lv][t];
      for (int t=0; t<3; t++)    delete Timer_Par_Update[lv][t];
      delete Timer_Par_2Sib   [lv];
      delete Timer_Par_2Son   [lv];
      delete Timer_Par_Collect[lv];
      for (int t=0; t<6; t++)    delete Timer_Par_MPI   [lv][t];

#     ifdef TIMING_SOLVER
      for (int v=0; v<5; v++)
      {
         delete Timer_Pre      [lv][v];
         delete Timer_Sol      [lv][v];
         delete Timer_Clo      [lv][v];
      }
      delete Timer_Poi_PreRho  [lv];
      delete Timer_Poi_PreFlu  [lv];
      delete Timer_Poi_PrePot_C[lv];
      delete Timer_Poi_PrePot_F[lv];
#     endif // #ifdef TIMING_SOLVER
   }

} // FUNCTION : Aux_DeleteTimer



//-------------------------------------------------------------------------------------------------------
// Function    :  Aux_ResetTimer
// Description :  Reset simulation timers by using the "Timer" structure
//-------------------------------------------------------------------------------------------------------
void Aux_ResetTimer()
{

   for (int t=0; t<7; t++)    Timer_Main[t]->Reset();

   for (int lv=0; lv<NLEVEL; lv++)
   {
      Timer_Flu_Advance[lv]->Reset();
      Timer_Gra_Advance[lv]->Reset();
      Timer_Che_Advance[lv]->Reset();
      Timer_FixUp      [lv]->Reset();
      Timer_Flag       [lv]->Reset();
      Timer_Refine     [lv]->Reset();
      Timer_Lv         [lv]->Reset();
      for (int t=0; t<9; t++)    Timer_GetBuf    [lv][t]->Reset();
      for (int t=0; t<3; t++)    Timer_Par_Update[lv][t]->Reset();
      Timer_Par_2Sib   [lv]->Reset();
      Timer_Par_2Son   [lv]->Reset();
      Timer_Par_Collect[lv]->Reset();
      for (int t=0; t<6; t++)    Timer_Par_MPI   [lv][t]->Reset();

#     ifdef TIMING_SOLVER
      for (int v=0; v<5; v++)
      {
         Timer_Pre      [lv][v]->Reset();
         Timer_Sol      [lv][v]->Reset();
         Timer_Clo      [lv][v]->Reset();
      }
      Timer_Poi_PreRho  [lv]->Reset();
      Timer_Poi_PreFlu  [lv]->Reset();
      Timer_Poi_PrePot_C[lv]->Reset();
      Timer_Poi_PrePot_F[lv]->Reset();
#     endif // #ifdef TIMING_SOLVER
   }

} // FUNCTION : Aux_ResetTimer



//-------------------------------------------------------------------------------------------------------
// Function    :  Aux_Record_Timing
// Description :  Record the timing results (in second)
//
// Note        :  The option "TIMING_SOLVER" record the MAXIMUM values of all ranks
//-------------------------------------------------------------------------------------------------------
void Aux_Record_Timing()
{

   const char FileName[] = "Record__Timing";

   FILE *File = NULL;

   const char Comment_LB[][4] = { "Max", "Min", "Ave" };
   const int  NLB             = 7;
   double Time_LB_Main[NLB][3];     // Time[][0/1/2] = maximum/minimum/average

// only the root rank needs to output the timing results
   if ( MPI_Rank == 0 )
   {
//    check if file already exists
      static bool FirstTime = true;
      if ( FirstTime )
      {
         if ( Aux_CheckFileExist(FileName) )
            Aux_Message( stderr, "WARNING : file \"%s\" already exists !!\n", FileName );

         FirstTime = false;

//       output header
         File = fopen( FileName, "a" );

         fprintf( File, "# Flu_Adv    : evolve fluid variables\n" );
         fprintf( File, "# Gra_Adv    : calculate gravitational potential and acceleration\n" );
         fprintf( File, "# Che_Adv    : Grackle cooling and chemistry library\n" );
         fprintf( File, "# FixUp      : use fine-grid data to correct coarse-grid data\n" );
         fprintf( File, "# Flag       : check refinement criteria\n" );
         fprintf( File, "# Refine     : allocate/remove patches on a higher level\n" );
         fprintf( File, "# Buf_Rho    : MPI for exchanging density for the Poisson solver\n" );
         fprintf( File, "# Buf_Pot    : MPI for exchanging potential\n" );
         fprintf( File, "# Buf_Flu1   : MPI for exchanging fluid data after advance\n" );
         fprintf( File, "# Buf_Flu2   : MPI for exchanging fluid data after fix-up\n" );
         fprintf( File, "# Buf_Ref    : MPI for exchanging fluid data and potential after grid refinement\n" );
         fprintf( File, "# Buf_Flux   : MPI for exchanging fluxes across coarse-fine boundaries\n" );
         fprintf( File, "# Buf_Res    : MPI for exchanging restricted fluid data in the father buffer patches (LOAD_BALANCE only)\n" );
         fprintf( File, "# Buf_Che    : MPI for exchanging fluid data after calling Grackle\n" );
         fprintf( File, "# Par_KD     : kick-drift for particles\n" );
         fprintf( File, "# Par_K      : last kick for particles staying at the same level\n" );
         fprintf( File, "# Par_K-1    : last kick for particles entering coarser grids\n" );
         fprintf( File, "# Par_2Sib   : pass particles to sibling patches after drift\n" );
         fprintf( File, "# -MPI_Sib   : MPI for exchanging particles staying at the same level after drift (included in Par_2Sib)\n" );
         fprintf( File, "# -MPI_FaSib : MPI for exchanging particles entering coarser grids after drift (included in Par_2Sib)\n" );
         fprintf( File, "# Par_2Son   : pass particles to son patches after the last kick\n" );
         fprintf( File, "# -MPI       : MPI for exchanging particles after the last kick (included in Par_2Son)\n" );
         fprintf( File, "# Par_Coll   : collect particles from finer grids for the Poisson solver\n" );
         fprintf( File, "# -MPI_Real  : MPI for collecting particles from leaf patches in other ranks (included in Par_Coll)\n" );
         fprintf( File, "# -MPI_Sib   : MPI for collecting particles to sibling buffer patches (included in Par_Coll)\n" );
         fprintf( File, "# -MPI_FaSib : MPI for collecting particles to father-sibling buffer patches (included in Par_Coll)\n" );
         fprintf( File, "#--------------------------------------------------------------------------------------" );
         fprintf( File, "---------------------------------------\n\n" );


         fclose( File );
      } // if ( FirstTime )


      File = fopen( FileName, "a" );

      fprintf( File, "Time : %13.7e -> %13.7e,     Step : %8ld -> %8ld\n\n", Time[0]-dTime_Base, Time[0],
                                                                             Step-1, Step );
//    1. main loop
      fprintf( File, "Main Loop\n" );
      fprintf( File, "---------------------------------------------------------------------------------------" );
      fprintf( File, "---------------------------------------\n" );
      fprintf( File, "%3s%9s%17s%15s%13s%13s%15s%15s%15s\n",
               "", "Total", "TimeStep", "Integration", "Output", "Auxiliary", "LoadBalance", "CorrSync", "Sum" );
   } // if ( MPI_Rank == 0 )


   if ( OPT__TIMING_BALANCE )
   {
      double Send[NLB];
      double (*Recv)[NLB] = new double [MPI_NRank][NLB];

      for (int t=0; t<NLB; t++)  Send[t] = Timer_Main[t]->GetValue( 0 );

      MPI_Gather( Send, NLB, MPI_DOUBLE, Recv[0], NLB, MPI_DOUBLE, 0, MPI_COMM_WORLD );

      if ( MPI_Rank == 0 )
      {
         for (int t=0; t<NLB; t++)
         {
            Time_LB_Main[t][0] = __FLT_MIN__;
            Time_LB_Main[t][1] = __FLT_MAX__;
            Time_LB_Main[t][2] = 0.0;

            for (int r=0; r<MPI_NRank; r++)
            {
               Time_LB_Main[t][0]  = MAX( Time_LB_Main[t][0], Recv[r][t] );
               Time_LB_Main[t][1]  = MIN( Time_LB_Main[t][1], Recv[r][t] );
               Time_LB_Main[t][2] +=                          Recv[r][t];
            }

            Time_LB_Main[t][2] /= MPI_NRank;
         }

         for (int v=0; v<3; v++)
         fprintf( File, "%3s%9.4f%17.4f%15.4f%13.4f%13.4f%15.4f%15.4f%15.4f\n",
                  Comment_LB[v], Time_LB_Main[0][v], Time_LB_Main[1][v], Time_LB_Main[2][v], Time_LB_Main[3][v],
                  Time_LB_Main[4][v], Time_LB_Main[5][v], Time_LB_Main[6][v],
                  Time_LB_Main[1][v] + Time_LB_Main[2][v] + Time_LB_Main[3][v] +
                  Time_LB_Main[4][v] + Time_LB_Main[5][v] + Time_LB_Main[6][v] );

         fprintf( File, "\n\n" );

         fclose( File );
      } // if ( MPI_Rank == 0 )

      delete [] Recv;
   } // if ( OPT__TIMING_BALANCE )

   else
   {
      if ( MPI_Rank == 0 )
      {
         fprintf( File, "%3s%9.4f%17.4f%15.4f%13.4f%13.4f%15.4f%15.4f%15.4f\n", "",
                  Timer_Main[0]->GetValue( 0 ), Timer_Main[1]->GetValue( 0 ), Timer_Main[2]->GetValue( 0 ),
                  Timer_Main[3]->GetValue( 0 ), Timer_Main[4]->GetValue( 0 ), Timer_Main[5]->GetValue( 0 ),
                  Timer_Main[6]->GetValue( 0 ),
                  Timer_Main[1]->GetValue( 0 ) + Timer_Main[2]->GetValue( 0 ) + Timer_Main[3]->GetValue( 0 )
                + Timer_Main[4]->GetValue( 0 ) + Timer_Main[5]->GetValue( 0 ) + Timer_Main[6]->GetValue( 0 ) );

         fprintf( File, "\n\n" );

         fclose( File );
      }
   } // if ( OPT__TIMING_BALANCE ) ... else ...


// for timing the evolution of each AMR level
   Timing__EvolveLevel( FileName, Time_LB_Main );


// for timing GPU/CPU solvers
#  ifdef TIMING_SOLVER
   Timing__Solver( FileName );
#  endif


   if ( MPI_Rank == 0 )
   {
      FILE *File = fopen( FileName, "a" );
      fprintf( File, "=======================================================================================" );
      fprintf( File, "=======================================\n" );
      fprintf( File, "=======================================================================================" );
      fprintf( File, "=======================================\n\n\n\n" );
      fclose( File );
   }

} // FUNCTION : Aux_Record_Timing



//-------------------------------------------------------------------------------------------------------
// Function    :  Timing__EvolveLevel
// Description :  Record the timing results (in second) for different code sections in the function "EvolveLevel"
//
// Parameter   :  FileName     : Name of the output file
//                Time_LB_Main : Recorded elapsed time in the main function for the option "OPT__TIMING_BALANCE"
//-------------------------------------------------------------------------------------------------------
void Timing__EvolveLevel( const char FileName[], const double Time_LB_Main[][3] )
{

#  ifdef INDIVIDUAL_TIMESTEP
   const int NSubStep = 2;
#  else
   const int NSubStep = 1;
#  endif

   FILE *File = ( MPI_Rank == 0 ) ? fopen( FileName, "a" ) : NULL;

   double Total[NSubStep][NLEVEL], Flu_Advance[NSubStep][NLEVEL], Gra_Advance[NSubStep][NLEVEL], FixUp[NSubStep][NLEVEL];
   double Flag[NSubStep][NLEVEL], Refine[NSubStep][NLEVEL], GetBuf[NSubStep][NLEVEL][9], Sum[NSubStep][NLEVEL];
   double ParUpdate[NSubStep][NLEVEL][3], Par2Sib[NSubStep][NLEVEL], Par2Son[NSubStep][NLEVEL], ParCollect[NSubStep][NLEVEL];
   double ParMPI[NSubStep][NLEVEL][6], Che_Advance[NSubStep][NLEVEL];

   const char Comment_LB[][4] = { "Max", "Min", "Ave" };
   const int NLB = 27;
   double Time_LB[NSubStep][NLEVEL][NLB][3];    // [0/1/2] = maximum/minimum/average
   double Sum_LB[NSubStep][NLEVEL][3];


// 1. timing for each sub-step
   for (int SubStep=0; SubStep<NSubStep; SubStep++)
   {
      if ( MPI_Rank == 0 )
      {
         fprintf( File, "Integration sub-step %d\n", SubStep+1 );
         fprintf( File, "---------------------------------------------------------------------------------------" );
         fprintf( File, "---------------------------------------\n" );
         fprintf( File, "%3s%4s%8s %8s%8s%8s%8s%8s%8s%8s%8s%9s%9s%8s%9s%8s%8s%8s%8s%8s%9s%9s%11s%9s%8s%9s%10s%9s%11s%8s\n",
                  "", "Lv", "Total", "Flu_Adv", "Gra_Adv", "Che_Adv", "FixUp", "Flag", "Refine",
                  "Buf_Rho", "Buf_Pot", "Buf_Flu1", "Buf_Flu2", "Buf_Ref", "Buf_Flux", "Buf_Res", "Buf_Che",
                  "Par_KD", "Par_K", "Par_K-1", "Par_2Sib", "-MPI_Sib", "-MPI_FaSib", "Par_2Son", "-MPI",
                  "Par_Coll", "-MPI_Real", "-MPI_Sib", "-MPI_FaSib", "Sum" );
      }

      for (int lv=0; lv<NLEVEL; lv++)
      {
         Total      [SubStep][lv]    = 0.0;
         Flu_Advance[SubStep][lv]    = 0.0;
         Gra_Advance[SubStep][lv]    = 0.0;
         Che_Advance[SubStep][lv]    = 0.0;
         FixUp      [SubStep][lv]    = 0.0;
         Flag       [SubStep][lv]    = 0.0;
         Refine     [SubStep][lv]    = 0.0;
         GetBuf     [SubStep][lv][0] = 0.0;
         GetBuf     [SubStep][lv][1] = 0.0;
         GetBuf     [SubStep][lv][2] = 0.0;
         GetBuf     [SubStep][lv][3] = 0.0;
         GetBuf     [SubStep][lv][4] = 0.0;
         GetBuf     [SubStep][lv][5] = 0.0;
         GetBuf     [SubStep][lv][6] = 0.0;
         GetBuf     [SubStep][lv][7] = 0.0;
         GetBuf     [SubStep][lv][8] = 0.0;
         ParUpdate  [SubStep][lv][0] = 0.0;
         ParUpdate  [SubStep][lv][1] = 0.0;
         ParUpdate  [SubStep][lv][2] = 0.0;
         Par2Sib    [SubStep][lv]    = 0.0;
         Par2Son    [SubStep][lv]    = 0.0;
         ParCollect [SubStep][lv]    = 0.0;
         ParMPI     [SubStep][lv][0] = 0.0;
         ParMPI     [SubStep][lv][1] = 0.0;
         ParMPI     [SubStep][lv][2] = 0.0;
         ParMPI     [SubStep][lv][3] = 0.0;
         ParMPI     [SubStep][lv][4] = 0.0;
         ParMPI     [SubStep][lv][5] = 0.0;

#        ifdef INDIVIDUAL_TIMESTEP
         const int N = 1<<(1+lv);
#        else
         const int N = 1;
#        endif

         for (int t=SubStep; t<N; t+=2)
         {
            Total      [SubStep][lv]    += Timer_Lv         [lv]   ->GetValue( t );
            Flu_Advance[SubStep][lv]    += Timer_Flu_Advance[lv]   ->GetValue( t );
            Gra_Advance[SubStep][lv]    += Timer_Gra_Advance[lv]   ->GetValue( t );
            Che_Advance[SubStep][lv]    += Timer_Che_Advance[lv]   ->GetValue( t );
            FixUp      [SubStep][lv]    += Timer_FixUp      [lv]   ->GetValue( t );
            Flag       [SubStep][lv]    += Timer_Flag       [lv]   ->GetValue( t );
            Refine     [SubStep][lv]    += Timer_Refine     [lv]   ->GetValue( t );
            GetBuf     [SubStep][lv][0] += Timer_GetBuf     [lv][0]->GetValue( t );
            GetBuf     [SubStep][lv][1] += Timer_GetBuf     [lv][1]->GetValue( t );
            GetBuf     [SubStep][lv][2] += Timer_GetBuf     [lv][2]->GetValue( t );
            GetBuf     [SubStep][lv][3] += Timer_GetBuf     [lv][3]->GetValue( t );
            GetBuf     [SubStep][lv][4] += Timer_GetBuf     [lv][4]->GetValue( t );
            GetBuf     [SubStep][lv][5] += Timer_GetBuf     [lv][5]->GetValue( t );
            GetBuf     [SubStep][lv][6] += Timer_GetBuf     [lv][6]->GetValue( t );
            GetBuf     [SubStep][lv][7] += Timer_GetBuf     [lv][7]->GetValue( t );
            GetBuf     [SubStep][lv][8] += Timer_GetBuf     [lv][8]->GetValue( t );
            ParUpdate  [SubStep][lv][0] += Timer_Par_Update [lv][0]->GetValue( t );
            ParUpdate  [SubStep][lv][1] += Timer_Par_Update [lv][1]->GetValue( t );
            ParUpdate  [SubStep][lv][2] += Timer_Par_Update [lv][2]->GetValue( t );
            Par2Sib    [SubStep][lv]    += Timer_Par_2Sib   [lv]   ->GetValue( t );
            Par2Son    [SubStep][lv]    += Timer_Par_2Son   [lv]   ->GetValue( t );
            ParCollect [SubStep][lv]    += Timer_Par_Collect[lv]   ->GetValue( t );
            ParMPI     [SubStep][lv][0] += Timer_Par_MPI    [lv][0]->GetValue( t );
            ParMPI     [SubStep][lv][1] += Timer_Par_MPI    [lv][1]->GetValue( t );
            ParMPI     [SubStep][lv][2] += Timer_Par_MPI    [lv][2]->GetValue( t );
            ParMPI     [SubStep][lv][3] += Timer_Par_MPI    [lv][3]->GetValue( t );
            ParMPI     [SubStep][lv][4] += Timer_Par_MPI    [lv][4]->GetValue( t );
            ParMPI     [SubStep][lv][5] += Timer_Par_MPI    [lv][5]->GetValue( t );
         }

//       subtract the Par_CollectParticle2OneLevel time from the Gra_AdvanceDt time (only necessary for refinement levels)
         if ( lv > 0 )  Gra_Advance[SubStep][lv] -= ParCollect[SubStep][lv];

         if ( OPT__TIMING_BALANCE )
         {
            double Send[NLB];
            double (*Recv)[NLB] = new double [MPI_NRank][NLB];

            Send[ 0] = Total      [SubStep][lv];
            Send[ 1] = Flu_Advance[SubStep][lv];
            Send[ 2] = Gra_Advance[SubStep][lv];
            Send[ 3] = Che_Advance[SubStep][lv];
            Send[ 4] = FixUp      [SubStep][lv];
            Send[ 5] = Flag       [SubStep][lv];
            Send[ 6] = Refine     [SubStep][lv];
            Send[ 7] = GetBuf     [SubStep][lv][0];
            Send[ 8] = GetBuf     [SubStep][lv][1];
            Send[ 9] = GetBuf     [SubStep][lv][2];
            Send[10] = GetBuf     [SubStep][lv][3];
            Send[11] = GetBuf     [SubStep][lv][4] +
                       GetBuf     [SubStep][lv][5];
            Send[12] = GetBuf     [SubStep][lv][6];
            Send[13] = GetBuf     [SubStep][lv][7];
            Send[14] = GetBuf     [SubStep][lv][8];
            Send[15] = ParUpdate  [SubStep][lv][0];
            Send[16] = ParUpdate  [SubStep][lv][1];
            Send[17] = ParUpdate  [SubStep][lv][2];
            Send[18] = Par2Sib    [SubStep][lv];
            Send[19] = Par2Son    [SubStep][lv];
            Send[20] = ParCollect [SubStep][lv];
            Send[21] = ParMPI     [SubStep][lv][0];
            Send[22] = ParMPI     [SubStep][lv][1];
            Send[23] = ParMPI     [SubStep][lv][2];
            Send[24] = ParMPI     [SubStep][lv][3];
            Send[25] = ParMPI     [SubStep][lv][4];
            Send[26] = ParMPI     [SubStep][lv][5];

            MPI_Gather( Send, NLB, MPI_DOUBLE, Recv[0], NLB, MPI_DOUBLE, 0, MPI_COMM_WORLD );

            if ( MPI_Rank == 0 )
            {
               for (int t=0; t<NLB; t++)
               {
                  Time_LB[SubStep][lv][t][0] = __FLT_MIN__;
                  Time_LB[SubStep][lv][t][1] = __FLT_MAX__;
                  Time_LB[SubStep][lv][t][2] = 0.0;

                  for (int r=0; r<MPI_NRank; r++)
                  {
                     Time_LB[SubStep][lv][t][0]  = MAX( Time_LB[SubStep][lv][t][0], Recv[r][t] );
                     Time_LB[SubStep][lv][t][1]  = MIN( Time_LB[SubStep][lv][t][1], Recv[r][t] );
                     Time_LB[SubStep][lv][t][2] +=                                  Recv[r][t];
                  }

                  Time_LB[SubStep][lv][t][2] /= MPI_NRank;
               }

               for (int v=0; v<3; v++)
               {
                  Sum_LB[SubStep][lv][v] = 0.0;
                  for (int t=1; t<NLB; t++)     Sum_LB[SubStep][lv][v] += Time_LB[SubStep][lv][t][v];

//                don't add MPI time for particles in the sum
//                --> because they are already included in the corresponding host functions
                  for (int t=21; t<=26; t++)    Sum_LB[SubStep][lv][v] -= Time_LB[SubStep][lv][t][v];

                  fprintf( File, "%3s%4d%8.3f %8.3f%8.3f%8.3f%8.3f%8.3f%8.3f%8.3f%8.3f%9.3f%9.3f%8.3f%9.3f%8.3f%8.3f%8.3f%8.3f%8.3f%9.3f%9.3f%11.3f%9.3f%8.3f%9.3f%10.3f%9.3f%11.3f%8.3f\n",
                           Comment_LB[v], lv,
                           Time_LB[SubStep][lv][ 0][v], Time_LB[SubStep][lv][ 1][v], Time_LB[SubStep][lv][ 2][v],
                           Time_LB[SubStep][lv][ 3][v], Time_LB[SubStep][lv][ 4][v], Time_LB[SubStep][lv][ 5][v],
                           Time_LB[SubStep][lv][ 6][v], Time_LB[SubStep][lv][ 7][v], Time_LB[SubStep][lv][ 8][v],
                           Time_LB[SubStep][lv][ 9][v], Time_LB[SubStep][lv][10][v], Time_LB[SubStep][lv][11][v],
                           Time_LB[SubStep][lv][12][v], Time_LB[SubStep][lv][13][v], Time_LB[SubStep][lv][14][v],
                           Time_LB[SubStep][lv][15][v], Time_LB[SubStep][lv][16][v], Time_LB[SubStep][lv][19][v],
                           Time_LB[SubStep][lv][20][v], Time_LB[SubStep][lv][17][v], Time_LB[SubStep][lv][21][v],
                           Time_LB[SubStep][lv][18][v], Time_LB[SubStep][lv][22][v], Time_LB[SubStep][lv][23][v],
                           Time_LB[SubStep][lv][24][v], Time_LB[SubStep][lv][25][v], Time_LB[SubStep][lv][26][v],
                           Sum_LB [SubStep][lv][v] );
               }

               fprintf( File, "\n" );
            } // if ( MPI_Rank == 0 )

            delete [] Recv;
         } // if ( OPT__TIMING_BALANCE )

         else
         {
//          don't add MPI time for particles in the sum
//          --> because they are already included in the corresponding host functions
            Sum[SubStep][lv] = Flu_Advance[SubStep][lv] + Gra_Advance[SubStep][lv] + Che_Advance[SubStep][lv] +
                               FixUp[SubStep][lv] + Flag[SubStep][lv] + Refine[SubStep][lv] +
                               GetBuf[SubStep][lv][0] + GetBuf[SubStep][lv][1] + GetBuf[SubStep][lv][2] +
                               GetBuf[SubStep][lv][3] + GetBuf[SubStep][lv][4] + GetBuf[SubStep][lv][5] +
                               GetBuf[SubStep][lv][6] + GetBuf[SubStep][lv][7] + GetBuf[SubStep][lv][8] +
                               ParUpdate[SubStep][lv][0] + ParUpdate[SubStep][lv][1] + ParUpdate[SubStep][lv][2] +
                               Par2Sib[SubStep][lv] + Par2Son[SubStep][lv] + ParCollect[SubStep][lv];

            if ( MPI_Rank == 0 )
            fprintf( File, "%3s%4d%8.3f %8.3f%8.3f%8.3f%8.3f%8.3f%8.3f%8.3f%8.3f%9.3f%9.3f%8.3f%9.3f%8.3f%8.3f%8.3f%8.3f%8.3f%9.3f%9.3f%11.3f%9.3f%8.3f%9.3f%10.3f%9.3f%11.3f%8.3f\n",
                     "", lv, Total[SubStep][lv], Flu_Advance[SubStep][lv], Gra_Advance[SubStep][lv], Che_Advance[SubStep][lv],
                     FixUp[SubStep][lv], Flag[SubStep][lv], Refine[SubStep][lv],
                     GetBuf[SubStep][lv][0], GetBuf[SubStep][lv][1], GetBuf[SubStep][lv][2], GetBuf[SubStep][lv][3],
                     GetBuf[SubStep][lv][4]+GetBuf[SubStep][lv][5],  GetBuf[SubStep][lv][6], GetBuf[SubStep][lv][7],
                     GetBuf[SubStep][lv][8],
                     ParUpdate[SubStep][lv][0], ParUpdate[SubStep][lv][1], ParUpdate[SubStep][lv][2],
                     Par2Sib[SubStep][lv], ParMPI[SubStep][lv][0], ParMPI[SubStep][lv][1],
                     Par2Son[SubStep][lv], ParMPI[SubStep][lv][2],
                     ParCollect[SubStep][lv], ParMPI[SubStep][lv][3], ParMPI[SubStep][lv][4], ParMPI[SubStep][lv][5],
                     Sum[SubStep][lv] );
         } // if ( OPT__TIMING_BALANCE ) ... else ...

      } // for (int lv=0; lv<NLEVEL; lv++)


//    sum over all levels
      if ( MPI_Rank == 0 )
      {
         if ( OPT__TIMING_BALANCE )
         {
            for (int lv=1; lv<NLEVEL; lv++)
            {
               for (int t=0; t<NLB; t++)
               for (int v=0; v<3; v++)    Time_LB[SubStep][0][t][v] += Time_LB[SubStep][lv][t][v];

               for (int v=0; v<3; v++)    Sum_LB[SubStep][0][v] += Sum_LB[SubStep][lv][v];
            }

            for (int v=0; v<3; v++)
            fprintf( File, "%3s%4s%8.3f %8.3f%8.3f%8.3f%8.3f%8.3f%8.3f%8.3f%8.3f%9.3f%9.3f%8.3f%9.3f%8.3f%8.3f%8.3f%8.3f%8.3f%9.3f%9.3f%11.3f%9.3f%8.3f%9.3f%10.3f%9.3f%11.3f%8.3f\n",
                     Comment_LB[v], "Sum",
                     Time_LB[SubStep][0][ 0][v], Time_LB[SubStep][0][ 1][v], Time_LB[SubStep][0][ 2][v],
                     Time_LB[SubStep][0][ 3][v], Time_LB[SubStep][0][ 4][v], Time_LB[SubStep][0][ 5][v],
                     Time_LB[SubStep][0][ 6][v], Time_LB[SubStep][0][ 7][v], Time_LB[SubStep][0][ 8][v],
                     Time_LB[SubStep][0][ 9][v], Time_LB[SubStep][0][10][v], Time_LB[SubStep][0][11][v],
                     Time_LB[SubStep][0][12][v], Time_LB[SubStep][0][13][v], Time_LB[SubStep][0][14][v],
                     Time_LB[SubStep][0][15][v], Time_LB[SubStep][0][16][v], Time_LB[SubStep][0][19][v],
                     Time_LB[SubStep][0][20][v], Time_LB[SubStep][0][17][v], Time_LB[SubStep][0][21][v],
                     Time_LB[SubStep][0][18][v], Time_LB[SubStep][0][22][v], Time_LB[SubStep][0][23][v],
                     Time_LB[SubStep][0][24][v], Time_LB[SubStep][0][25][v], Time_LB[SubStep][0][26][v],
                     Sum_LB [SubStep][0][v] );

            fprintf( File, "\n\n" );
         } // if ( OPT__TIMING_BALANCE )

         else
         {
            for (int lv=1; lv<NLEVEL; lv++)
            {
               Total      [SubStep][0]    += Total      [SubStep][lv];
               Flu_Advance[SubStep][0]    += Flu_Advance[SubStep][lv];
               Gra_Advance[SubStep][0]    += Gra_Advance[SubStep][lv];
               Che_Advance[SubStep][0]    += Che_Advance[SubStep][lv];
               FixUp      [SubStep][0]    += FixUp      [SubStep][lv];
               Flag       [SubStep][0]    += Flag       [SubStep][lv];
               Refine     [SubStep][0]    += Refine     [SubStep][lv];
               GetBuf     [SubStep][0][0] += GetBuf     [SubStep][lv][0];
               GetBuf     [SubStep][0][1] += GetBuf     [SubStep][lv][1];
               GetBuf     [SubStep][0][2] += GetBuf     [SubStep][lv][2];
               GetBuf     [SubStep][0][3] += GetBuf     [SubStep][lv][3];
               GetBuf     [SubStep][0][4] += GetBuf     [SubStep][lv][4];
               GetBuf     [SubStep][0][5] += GetBuf     [SubStep][lv][5];
               GetBuf     [SubStep][0][6] += GetBuf     [SubStep][lv][6];
               GetBuf     [SubStep][0][7] += GetBuf     [SubStep][lv][7];
               GetBuf     [SubStep][0][8] += GetBuf     [SubStep][lv][8];
               ParUpdate  [SubStep][0][0] += ParUpdate  [SubStep][lv][0];
               ParUpdate  [SubStep][0][1] += ParUpdate  [SubStep][lv][1];
               ParUpdate  [SubStep][0][2] += ParUpdate  [SubStep][lv][2];
               Par2Sib    [SubStep][0]    += Par2Sib    [SubStep][lv];
               Par2Son    [SubStep][0]    += Par2Son    [SubStep][lv];
               ParCollect [SubStep][0]    += ParCollect [SubStep][lv];
               ParMPI     [SubStep][0][0] += ParMPI     [SubStep][lv][0];
               ParMPI     [SubStep][0][1] += ParMPI     [SubStep][lv][1];
               ParMPI     [SubStep][0][2] += ParMPI     [SubStep][lv][2];
               ParMPI     [SubStep][0][3] += ParMPI     [SubStep][lv][3];
               ParMPI     [SubStep][0][4] += ParMPI     [SubStep][lv][4];
               ParMPI     [SubStep][0][5] += ParMPI     [SubStep][lv][5];
               Sum        [SubStep][0]    += Sum        [SubStep][lv];
            }

            fprintf( File, "%3s%4s%8.3f %8.3f%8.3f%8.3f%8.3f%8.3f%8.3f%8.3f%8.3f%9.3f%9.3f%8.3f%9.3f%8.3f%8.3f%8.3f%8.3f%8.3f%9.3f%9.3f%11.3f%9.3f%8.3f%9.3f%10.3f%9.3f%11.3f%8.3f\n",
                     "", "Sum", Total[SubStep][0], Flu_Advance[SubStep][0], Gra_Advance[SubStep][0], Che_Advance[SubStep][0],
                     FixUp[SubStep][0], Flag[SubStep][0], Refine[SubStep][0],
                     GetBuf[SubStep][0][0], GetBuf[SubStep][0][1], GetBuf[SubStep][0][2], GetBuf[SubStep][0][3],
                     GetBuf[SubStep][0][4]+GetBuf[SubStep][0][5],  GetBuf[SubStep][0][6], GetBuf[SubStep][0][7],
                     GetBuf[SubStep][0][8],
                     ParUpdate[SubStep][0][0], ParUpdate[SubStep][0][1], ParUpdate[SubStep][0][2],
                     Par2Sib[SubStep][0], ParMPI[SubStep][0][0], ParMPI[SubStep][0][1],
                     Par2Son[SubStep][0], ParMPI[SubStep][0][2],
                     ParCollect[SubStep][0], ParMPI[SubStep][0][3], ParMPI[SubStep][0][4], ParMPI[SubStep][0][5],
                     Sum[SubStep][0] );

            fprintf( File, "\n" );
         } // if ( OPT__TIMING_BALANCE ) ... else ...
      } // if ( MPI_Rank == 0 )
   } // for (int SubStep=0; SubStep<NSubStep; SubStep++)


// 2. summary
   if ( MPI_Rank == 0 ) {
   if ( OPT__TIMING_BALANCE )
   {
//    _P : percentage; _IM : imbalance
      double Everything[3], MPI_Grid[3], dt[3], Aux[3], Corr[3], Output[3], LB[3], Par[3], MPI_Par[3];
      double Flu_P, Gra_P, Che_P, FixUp_P, Flag_P, Refine_P, Sum_P, MPI_Grid_P, dt_P, Aux_P, Corr_P, Output_P, LB_P, Par_P, MPI_Par_P;
      double Flu_IB, Gra_IB, Che_IB, FixUp_IB, Flag_IB, Refine_IB, Sum_IB, MPI_Grid_IB, dt_IB, Aux_IB, Corr_IB, Output_IB, LB_IB, Par_IB, MPI_Par_IB;

      fprintf( File, "Summary\n" );
      fprintf( File, "---------------------------------------------------------------------------------------" );
      fprintf( File, "---------------------------------------\n" );
      fprintf( File, "%3s%5s %11s%9s%9s%9s%9s%9s%9s%9s%9s%9s%9s%9s%9s%9s%12s\n",
               "", "", "Flu_Adv", "Gra_Adv", "Che_Adv", "FixUp", "Flag", "Refine", "MPI_Grid", "dt", "Output", "Aux", "LB",
               "CorrSync", "Par", "-MPI_Par", "Sum" );

      for (int v=0; v<3; v++)
      {
         Everything[v] = Time_LB_Main[0][v];
         dt        [v] = Time_LB_Main[1][v];
         Output    [v] = Time_LB_Main[3][v];
         Aux       [v] = Time_LB_Main[4][v];
         LB        [v] = Time_LB_Main[5][v];
         Corr      [v] = Time_LB_Main[6][v];

//       sum
         for (int t=1; t<NSubStep; t++)
         {
            for (int k=1; k<NLB; k++)  Time_LB[0][0][k][v] += Time_LB[t][0][k][v];

            Sum_LB[0][0][v] += Sum_LB[t][0][v];
         }

         MPI_Grid[v] = 0.0;
         for (int k=7; k<15; k++)   MPI_Grid[v] += Time_LB[0][0][k][v];

         Par[v] = 0.0;
         for (int k=15; k<21; k++)  Par[v] += Time_LB[0][0][k][v];

         MPI_Par[v] = 0.0;
         for (int k=21; k<27; k++)  MPI_Par[v] += Time_LB[0][0][k][v];

         Sum_LB[0][0][v] += dt[v] + Output[v] + Aux[v] + LB[v] + Corr[v];

//       2.1 time
         fprintf( File, "%3s%5s %11.4f%9.4f%9.4f%9.4f%9.4f%9.4f%9.4f%9.4f%9.4f%9.4f%9.4f%9.4f%9.4f%12.4f\n",
                  Comment_LB[v], "Time", Time_LB[0][0][1][v], Time_LB[0][0][2][v], Time_LB[0][0][3][v],
                  Time_LB[0][0][4][v], Time_LB[0][0][5][v], MPI_Grid[v], dt[v], Output[v], Aux[v], LB[v], Corr[v],
                  Par[v], MPI_Par[v], Sum_LB[0][0][v] );
      } // for (int v=0; v<3; v++)

      fprintf( File, "\n" );


//    2.2 "max" imbalance = (Max-Ave)/Ave
      Flu_IB      = 100.0*( Time_LB[0][0][1][0] - Time_LB[0][0][1][2] ) / ((Time_LB[0][0][1][2]==0.0)?1.0:Time_LB[0][0][1][2]);
      Gra_IB      = 100.0*( Time_LB[0][0][2][0] - Time_LB[0][0][2][2] ) / ((Time_LB[0][0][2][2]==0.0)?1.0:Time_LB[0][0][2][2]);
      Che_IB      = 100.0*( Time_LB[0][0][3][0] - Time_LB[0][0][3][2] ) / ((Time_LB[0][0][3][2]==0.0)?1.0:Time_LB[0][0][3][2]);
      FixUp_IB    = 100.0*( Time_LB[0][0][4][0] - Time_LB[0][0][4][2] ) / ((Time_LB[0][0][4][2]==0.0)?1.0:Time_LB[0][0][4][2]);
      Flag_IB     = 100.0*( Time_LB[0][0][5][0] - Time_LB[0][0][5][2] ) / ((Time_LB[0][0][5][2]==0.0)?1.0:Time_LB[0][0][5][2]);
      Refine_IB   = 100.0*( Time_LB[0][0][6][0] - Time_LB[0][0][6][2] ) / ((Time_LB[0][0][6][2]==0.0)?1.0:Time_LB[0][0][6][2]);
      MPI_Grid_IB = 100.0*( MPI_Grid        [0] - MPI_Grid        [2] ) / ((MPI_Grid        [2]==0.0)?1.0:MPI_Grid        [2]);
      dt_IB       = 100.0*( dt              [0] - dt              [2] ) / ((dt              [2]==0.0)?1.0:dt              [2]);
      Output_IB   = 100.0*( Output          [0] - Output          [2] ) / ((Output          [2]==0.0)?1.0:Output          [2]);
      Aux_IB      = 100.0*( Aux             [0] - Aux             [2] ) / ((Aux             [2]==0.0)?1.0:Aux             [2]);
      LB_IB       = 100.0*( LB              [0] - LB              [2] ) / ((LB              [2]==0.0)?1.0:LB              [2]);
      Corr_IB     = 100.0*( Corr            [0] - Corr            [2] ) / ((Corr            [2]==0.0)?1.0:Corr            [2]);
      Par_IB      = 100.0*( Par             [0] - Par             [2] ) / ((Par             [2]==0.0)?1.0:Par             [2]);
      MPI_Par_IB  = 100.0*( MPI_Par         [0] - MPI_Par         [2] ) / ((MPI_Par         [2]==0.0)?1.0:MPI_Par         [2]);
      Sum_IB      = 100.0*( Sum_LB[0][0]    [0] - Sum_LB[0][0]    [2] ) / ((Sum_LB[0][0]    [2]==0.0)?1.0:Sum_LB[0][0]    [2]);

      fprintf( File, "%9s%10.3f%%%8.3f%%%8.3f%%%8.3f%%%8.3f%%%8.3f%%%8.3f%%%8.3f%%%8.3f%%%8.3f%%%8.3f%%%8.3f%%%8.3f%%%8.3f%%%11.3f%%\n",
               "Imbalance", Flu_IB, Gra_IB, Che_IB, FixUp_IB, Flag_IB, Refine_IB, MPI_Grid_IB, dt_IB, Output_IB, Aux_IB, LB_IB,
               Corr_IB, Par_IB, MPI_Par_IB, Sum_IB );


//    2.3 "max" percentage
      Flu_P      = 100.0*Time_LB[0][0][1][0]/Everything[0];   // always divided by the maximum total time
      Gra_P      = 100.0*Time_LB[0][0][2][0]/Everything[0];
      Che_P      = 100.0*Time_LB[0][0][3][0]/Everything[0];
      FixUp_P    = 100.0*Time_LB[0][0][4][0]/Everything[0];
      Flag_P     = 100.0*Time_LB[0][0][5][0]/Everything[0];
      Refine_P   = 100.0*Time_LB[0][0][6][0]/Everything[0];
      MPI_Grid_P = 100.0*MPI_Grid        [0]/Everything[0];
      dt_P       = 100.0*dt              [0]/Everything[0];
      Output_P   = 100.0*Output          [0]/Everything[0];
      Aux_P      = 100.0*Aux             [0]/Everything[0];
      LB_P       = 100.0*LB              [0]/Everything[0];
      Corr_P     = 100.0*Corr            [0]/Everything[0];
      Par_P      = 100.0*Par             [0]/Everything[0];
      MPI_Par_P  = 100.0*MPI_Par         [0]/Everything[0];
      Sum_P      = 100.0*Sum_LB [0][0]   [0]/Everything[0];

      fprintf( File, "%3s%5s %10.3f%%%8.3f%%%8.3f%%%8.3f%%%8.3f%%%8.3f%%%8.3f%%%8.3f%%%8.3f%%%8.3f%%%8.3f%%%8.3f%%%8.3f%%%8.3f%%%11.3f%%\n",
               "Max", "Frac", Flu_P, Gra_P, Che_P, FixUp_P, Flag_P, Refine_P, MPI_Grid_P, dt_P, Output_P, Aux_P, LB_P,
               Corr_P, Par_P, MPI_Par_P, Sum_P );

      fprintf( File, "\n" );


//    2.4 record the accumulated timing results
      for (int v=0; v<3; v++)
      {
         Flu_Acc     [v] += Time_LB[0][0][1][v];
         Gra_Acc     [v] += Time_LB[0][0][2][v];
         Che_Acc     [v] += Time_LB[0][0][3][v];
         FixUp_Acc   [v] += Time_LB[0][0][4][v];
         Flag_Acc    [v] += Time_LB[0][0][5][v];
         Refine_Acc  [v] += Time_LB[0][0][6][v];
         MPI_Grid_Acc[v] += MPI_Grid        [v];
         dt_Acc      [v] += dt              [v];
         Output_Acc  [v] += Output          [v];
         Aux_Acc     [v] += Aux             [v];
         LB_Acc      [v] += LB              [v];
         Corr_Acc    [v] += Corr            [v];
         Par_Acc     [v] += Par             [v];
         Sum_Acc     [v] += Sum_LB[0][0]    [v];
         MPI_Par_Acc [v] += MPI_Par         [v];
      }
   } // if ( OPT__TIMING_BALANCE )

   else
   {
      double Everything, MPI_Grid, dt, Aux, Corr, Output, LB, Par, MPI_Par;
      double Flu_P, Gra_P, Che_P, FixUp_P, Flag_P, Refine_P, Sum_P, MPI_Grid_P, dt_P, Aux_P, Corr_P, Output_P, LB_P, Par_P, MPI_Par_P;

      Everything = Timer_Main[0]->GetValue( 0 );
      dt         = Timer_Main[1]->GetValue( 0 );
      Output     = Timer_Main[3]->GetValue( 0 );
      Aux        = Timer_Main[4]->GetValue( 0 );
      LB         = Timer_Main[5]->GetValue( 0 );
      Corr       = Timer_Main[6]->GetValue( 0 );

//    sum
      for (int t=1; t<NSubStep; t++)
      {
         Flu_Advance[0][0]    += Flu_Advance[t][0];
         Gra_Advance[0][0]    += Gra_Advance[t][0];
         Che_Advance[0][0]    += Che_Advance[t][0];
         FixUp      [0][0]    += FixUp      [t][0];
         Flag       [0][0]    += Flag       [t][0];
         Refine     [0][0]    += Refine     [t][0];
         GetBuf     [0][0][0] += GetBuf     [t][0][0];
         GetBuf     [0][0][1] += GetBuf     [t][0][1];
         GetBuf     [0][0][2] += GetBuf     [t][0][2];
         GetBuf     [0][0][3] += GetBuf     [t][0][3];
         GetBuf     [0][0][4] += GetBuf     [t][0][4];
         GetBuf     [0][0][5] += GetBuf     [t][0][5];
         GetBuf     [0][0][6] += GetBuf     [t][0][6];
         GetBuf     [0][0][7] += GetBuf     [t][0][7];
         GetBuf     [0][0][8] += GetBuf     [t][0][8];
         ParUpdate  [0][0][0] += ParUpdate  [t][0][0];
         ParUpdate  [0][0][1] += ParUpdate  [t][0][1];
         ParUpdate  [0][0][2] += ParUpdate  [t][0][2];
         Par2Sib    [0][0]    += Par2Sib    [t][0];
         Par2Son    [0][0]    += Par2Son    [t][0];
         ParCollect [0][0]    += ParCollect [t][0];
         ParMPI     [0][0][0] += ParMPI     [t][0][0];
         ParMPI     [0][0][1] += ParMPI     [t][0][1];
         ParMPI     [0][0][2] += ParMPI     [t][0][2];
         ParMPI     [0][0][3] += ParMPI     [t][0][3];
         ParMPI     [0][0][4] += ParMPI     [t][0][4];
         ParMPI     [0][0][5] += ParMPI     [t][0][5];
         Sum        [0][0]    += Sum        [t][0];
      } // for (int t=1; t<NSubStep; t++)

      MPI_Grid = 0.0;
      for (int x=0; x<=8; x++)   MPI_Grid += GetBuf[0][0][x];

      Par = ParUpdate[0][0][0] + ParUpdate[0][0][1] + ParUpdate[0][0][2] + Par2Sib[0][0] + Par2Son[0][0] + ParCollect[0][0];

      MPI_Par = 0.0;
      for (int x=0; x<=5; x++)   MPI_Par += ParMPI[0][0][x];

      Sum[0][0] += dt + Output + Aux + LB + Corr;

//    percentage
      Flu_P      = 100.0*Flu_Advance[0][0]/Everything;
      Gra_P      = 100.0*Gra_Advance[0][0]/Everything;
      Che_P      = 100.0*Che_Advance[0][0]/Everything;
      FixUp_P    = 100.0*FixUp      [0][0]/Everything;
      Flag_P     = 100.0*Flag       [0][0]/Everything;
      Refine_P   = 100.0*Refine     [0][0]/Everything;
      MPI_Grid_P = 100.0*MPI_Grid         /Everything;
      dt_P       = 100.0*dt               /Everything;
      Output_P   = 100.0*Output           /Everything;
      Aux_P      = 100.0*Aux              /Everything;
      LB_P       = 100.0*LB               /Everything;
      Corr_P     = 100.0*Corr             /Everything;
      Par_P      = 100.0*Par              /Everything;
      MPI_Par_P  = 100.0*MPI_Par          /Everything;
      Sum_P      = 100.0*Sum        [0][0]/Everything;

      fprintf( File, "\nSummary\n" );
      fprintf( File, "---------------------------------------------------------------------------------------" );
      fprintf( File, "---------------------------------------\n" );
      fprintf( File, "%3s%5s %11s%9s%9s%9s%9s%9s%9s%9s%9s%9s%9s%9s%9s%9s%12s\n",
               "", "", "Flu_Adv", "Gra_Adv", "Che_Adv", "FixUp", "Flag", "Refine", "MPI_Grid", "dt", "Output", "Aux", "LB",
               "CorrSync", "Par", "-MPI_Par", "Sum" );

//    2.1 time
      fprintf( File, "%3s%5s %11.4f%9.4f%9.4f%9.4f%9.4f%9.4f%9.4f%9.4f%9.4f%9.4f%9.4f%9.4f%9.4f%9.4f%12.4f\n",
               "", "Time", Flu_Advance[0][0], Gra_Advance[0][0], Che_Advance[0][0], FixUp[0][0], Flag[0][0], Refine[0][0], MPI_Grid,
               dt, Output, Aux, LB, Corr, Par, MPI_Par, Sum[0][0] );

//    2.2 percentage
      fprintf( File, "%3s%5s %10.3f%%%8.3f%%%8.3f%%%8.3f%%%8.3f%%%8.3f%%%8.3f%%%8.3f%%%8.3f%%%8.3f%%%8.3f%%%8.3f%%%8.3f%%%8.3f%%%11.3f%%\n",
               "", "Frac", Flu_P, Gra_P, Che_P, FixUp_P, Flag_P, Refine_P, MPI_Grid_P, dt_P, Output_P, Aux_P, LB_P,
               Corr_P, Par_P, MPI_Par_P, Sum_P );
      fprintf( File, "\n" );


//    2.3 record the accumulated timing results
      Flu_Acc     [0] += Flu_Advance[0][0];
      Gra_Acc     [0] += Gra_Advance[0][0];
      Che_Acc     [0] += Che_Advance[0][0];
      FixUp_Acc   [0] += FixUp      [0][0];
      Flag_Acc    [0] += Flag       [0][0];
      Refine_Acc  [0] += Refine     [0][0];
      MPI_Grid_Acc[0] += MPI_Grid;
      dt_Acc      [0] += dt;
      Output_Acc  [0] += Output;
      Aux_Acc     [0] += Aux;
      LB_Acc      [0] += LB;
      Corr_Acc    [0] += Corr;
      Par_Acc     [0] += Par;
      MPI_Par_Acc [0] += MPI_Par;
      Sum_Acc     [0] += Sum        [0][0];

   }} // if ( OPT__TIMING_BALANCE ) ... else ... if ( MPI_Rank == 0 )


   if ( MPI_Rank == 0 )    fclose( File );

} // FUNCTION : Timing__EvolveLevel



/* old shared time-step function
//-------------------------------------------------------------------------------------------------------
// Function    :  Timing__SharedTimestep
// Description :  Record the timing results (in second) for the shared time-step integration
//-------------------------------------------------------------------------------------------------------
void Timing__SharedTimestep( const char FileName[] )
{

   FILE *File = fopen( FileName, "a" );

   double Gra_Sum = 0.0;

   for (int lv=0; lv<NLEVEL; lv++)
   {
      Gra_Sum += Timer_Gra_Advance [lv]   ->GetValue( 0 );
      Gra_Sum += Timer_Gra_Restrict[lv]   ->GetValue( 0 );
      Gra_Sum += Timer_GetBuf      [lv][1]->GetValue( 0 );
      Gra_Sum += Timer_GetBuf      [lv][2]->GetValue( 0 );
   }

// a. overall
   fprintf( File, "Integration_SharedTimeStep\n" );
   fprintf( File, "-----------------------------------------------------------------------------------------" );
   fprintf( File, "------------------------\n" );
   fprintf( File, "%9s%15s%11s%11s%13s\n", "Total", "Fluid1", "Gravity", "Fluid2", "Sum" );
   fprintf( File, "%9.4f%15.4f%11.4f%11.4f%13.4f\n",
            Timer_Main[2]->GetValue( 0 ), Timer_Flu_Total[0]->GetValue( 0 ), Gra_Sum,
            Timer_Flu_Total[0]->GetValue( 1 ),
            Timer_Flu_Total[0]->GetValue( 0 ) + Timer_Flu_Total[0]->GetValue( 1 ) + Gra_Sum );
   fprintf( File, "\n" );


// b. hydro
   double Flu_Total;
   double Flu_Sum[NLEVEL];

   for (int t=0; t<2; t++)
   {
      fprintf( File, "Fluid %d\n", t+1 );
      fprintf( File, "---------------------------------------------------------------------------------------" );
      fprintf( File, "--------------------------\n" );
#     ifdef GRAVITY
      if ( t == 0 )
      fprintf( File, "%5s%13s%13s%9s%9s%9s%9s%11s\n",
               "Level", "Total", "Advance",  "FixUp", "Buf_Rho", "Flag", "Refine", "Sum" );
      else
      fprintf( File, "%5s%13s%13s%9s%9s%9s%9s%11s\n",
               "Level", "Total", "Advance",  "FixUp", "Buf_Flu", "Flag", "Refine", "Sum" );
#     else
      fprintf( File, "%5s%13s%13s%9s%9s%9s%9s%11s\n",
               "Level", "Total", "Advance",  "FixUp", "Buf_Flu", "Flag", "Refine", "Sum" );
#     endif

      for (int lv=0; lv<NLEVEL; lv++)
      {
         if ( lv != NLEVEL - 1 )
            Flu_Total = Timer_Flu_Total[lv]->GetValue( t ) - Timer_Flu_Total[lv+1]->GetValue( t );
         else
            Flu_Total = Timer_Flu_Total[lv]->GetValue( t );

         Flu_Sum[lv] =   Timer_Flu_Advance[lv]     ->GetValue( t ) + Timer_FixUp[lv]->GetValue( t )
                       + Timer_GetBuf     [lv][t*3]->GetValue( 0 ) + Timer_Flag [lv]->GetValue( t )
                       + Timer_Refine     [lv]     ->GetValue( t );

         fprintf( File, "%5d%13.4f%13.4f%9.4f%9.4f%9.4f%9.4f%11.4f\n",
                  lv, Flu_Total,
                  Timer_Flu_Advance[lv]     ->GetValue( t ),
                  Timer_FixUp      [lv]     ->GetValue( t ),
                  Timer_GetBuf     [lv][t*3]->GetValue( 0 ),
                  Timer_Flag       [lv]     ->GetValue( t ),
                  Timer_Refine     [lv]     ->GetValue( t ),
                  Flu_Sum          [lv] );
      }

      for (int lv=1; lv<NLEVEL; lv++)  Flu_Sum[0] += Flu_Sum[lv];
      fprintf( File, "%78.4f\n", Flu_Sum[0] );

      fprintf( File, "\n" );

   } // for (int t=0; t<2; t++)


// c. gravity
#  ifdef GRAVITY

   double PoiGra_Sum[NLEVEL];

   fprintf( File, "Poisson + Gravity\n" );
   fprintf( File, "-----------------------------------------------------------------------------------------" );
   fprintf( File, "------------------------\n" );
   fprintf( File, "%5s%13s%11s%11s%11s%13s\n",
            "Level", "Advance", "Restrict", "Buf_Pot", "Buf_Flu", "Sum" );

   for (int lv=0; lv<NLEVEL; lv++)
   {
      PoiGra_Sum[lv] =   Timer_Gra_Advance [lv]   ->GetValue( 0 ) + Timer_Gra_Restrict[lv]   ->GetValue( 0 )
                       + Timer_GetBuf      [lv][1]->GetValue( 0 ) + Timer_GetBuf      [lv][2]->GetValue( 0 );

      fprintf( File, "%5d%13.4f%11.4f%11.4f%11.4f%13.4f\n",
               lv,
               Timer_Gra_Advance [lv]   ->GetValue( 0 ), Timer_Gra_Restrict[lv]   ->GetValue( 0 ),
               Timer_GetBuf      [lv][1]->GetValue( 0 ), Timer_GetBuf      [lv][2]->GetValue( 0 ),
               PoiGra_Sum[lv] );
   }

   for (int lv=1; lv<NLEVEL; lv++)  PoiGra_Sum[0] += PoiGra_Sum[lv];
   fprintf( File, "%64.4f\n", PoiGra_Sum[0] );

   fprintf( File, "\n" );

#  endif // #ifdef GRAVITY

   fclose( File );

} // FUNCTION : Timing__SharedTimestep
*/



#ifdef TIMING_SOLVER
//-------------------------------------------------------------------------------------------------------
// Function    :  Timing__Solver
// Description :  Record the timing results (in second) for the option "TIMING_SOLVER"
//-------------------------------------------------------------------------------------------------------
void Timing__Solver( const char FileName[] )
{

// get the maximum values from all ranks
   int    ID;
   double Pre_loc[NLEVEL][5], Sol_loc[NLEVEL][5], Clo_loc[NLEVEL][5];
   double Pre_max[NLEVEL][5], Sol_max[NLEVEL][5], Clo_max[NLEVEL][5];
   double Poi_PreRho_loc[NLEVEL], Poi_PreFlu_loc[NLEVEL], Poi_PrePot_C_loc[NLEVEL], Poi_PrePot_F_loc[NLEVEL];
   double Poi_PreRho_max[NLEVEL], Poi_PreFlu_max[NLEVEL], Poi_PrePot_C_max[NLEVEL], Poi_PrePot_F_max[NLEVEL];

   for (int lv=0; lv<NLEVEL; lv++)
   {
      for (int v=0; v<5; v++)
      {
         Pre_loc[lv][v] = Timer_Pre[lv][v]->GetValue(0);
         Sol_loc[lv][v] = Timer_Sol[lv][v]->GetValue(0);
         Clo_loc[lv][v] = Timer_Clo[lv][v]->GetValue(0);
      }

      Poi_PreRho_loc  [lv] = Timer_Poi_PreRho  [lv]->GetValue(0);
      Poi_PreFlu_loc  [lv] = Timer_Poi_PreFlu  [lv]->GetValue(0);
      Poi_PrePot_C_loc[lv] = Timer_Poi_PrePot_C[lv]->GetValue(0);
      Poi_PrePot_F_loc[lv] = Timer_Poi_PrePot_F[lv]->GetValue(0);
   }

   MPI_Reduce( Pre_loc[0],       Pre_max[0],       NLEVEL*5, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD );
   MPI_Reduce( Sol_loc[0],       Sol_max[0],       NLEVEL*5, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD );
   MPI_Reduce( Clo_loc[0],       Clo_max[0],       NLEVEL*5, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD );

   MPI_Reduce( Poi_PreRho_loc,   Poi_PreRho_max,   NLEVEL,   MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD );
   MPI_Reduce( Poi_PreFlu_loc,   Poi_PreFlu_max,   NLEVEL,   MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD );
   MPI_Reduce( Poi_PrePot_C_loc, Poi_PrePot_C_max, NLEVEL,   MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD );
   MPI_Reduce( Poi_PrePot_F_loc, Poi_PrePot_F_max, NLEVEL,   MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD );


   if ( MPI_Rank == 0 )
   {
      FILE *File = fopen( FileName, "a" );

      fprintf( File, "\nGPU/CPU solvers\n" );
      fprintf( File, "---------------------------------------------------------------------------------------" );
      fprintf( File, "---------------------------------------\n" );
      fprintf( File, "%3s%10s%10s%10s%10s%10s%10s%10s%10s%10s%10s%10s%10s%10s\n",
               "Lv", "Flu_Pre", "Flu_Sol", "Flu_Clo", "Poi_Pre", "PreRho", "PreFlu", "PrePot_C", "Pre_Pot_F",
               "Poi_Sol", "Poi_Clo", "Che_Pre", "Che_Sol", "Che_Clo" );

      for (int lv=0; lv<NLEVEL; lv++)
      {
         if ( lv == 0 )    ID = 2;
         else              ID = 3;

         fprintf( File, "%3d%10.4f%10.4f%10.4f%10.4f%10.4f%10.4f%10.4f%10.4f%10.4f%10.4f%10.4f%10.4f%10.4f\n",
                  lv,
                  Pre_max[lv][ 0],                                             Sol_max[lv][ 0], Clo_max[lv][ 0],
                  Pre_max[lv][ID], Poi_PreRho_max  [lv], Poi_PreFlu_max  [lv],
                                   Poi_PrePot_C_max[lv], Poi_PrePot_F_max[lv], Sol_max[lv][ID], Clo_max[lv][ID],
                  Pre_max[lv][ 4],                                             Sol_max[lv][ 4], Clo_max[lv][ 4] );
      }

//    sum over all levels
      for (int lv=1; lv<NLEVEL; lv++)
      {
         Pre_max         [0][0] += Pre_max         [lv][0];
         Sol_max         [0][0] += Sol_max         [lv][0];
         Clo_max         [0][0] += Clo_max         [lv][0];
         Pre_max         [0][2] += Pre_max         [lv][3];
         Poi_PreRho_max  [0]    += Poi_PreRho_max  [lv];
         Poi_PreFlu_max  [0]    += Poi_PreFlu_max  [lv];
         Poi_PrePot_C_max[0]    += Poi_PrePot_C_max[lv];
         Poi_PrePot_F_max[0]    += Poi_PrePot_F_max[lv];
         Sol_max         [0][2] += Sol_max         [lv][3];
         Clo_max         [0][2] += Clo_max         [lv][3];
         Pre_max         [0][4] += Pre_max         [lv][4];
         Sol_max         [0][4] += Sol_max         [lv][4];
         Clo_max         [0][4] += Clo_max         [lv][4];
      }

      fprintf( File, "%3s%10.4f%10.4f%10.4f%10.4f%10.4f%10.4f%10.4f%10.4f%10.4f%10.4f%10.4f%10.4f%10.4f\n",
               "Sum",
               Pre_max[0][0],                                           Sol_max[0][0], Clo_max[0][0],
               Pre_max[0][2], Poi_PreRho_max  [0], Poi_PreFlu_max  [0],
                              Poi_PrePot_C_max[0], Poi_PrePot_F_max[0], Sol_max[0][2], Clo_max[0][2],
               Pre_max[0][4],                                           Sol_max[0][4], Clo_max[0][4] );
      fprintf( File, "\n" );

      fclose( File );
   } // if ( MPI_Rank == 0 )

} // FUNCTION : TimingSolver
#endif // TIMING_SOLVER



//-------------------------------------------------------------------------------------------------------
// Function    :  Aux_AccumulatedTiming
// Description :  Record the accumulated timing results (in second)
//
// Note        :  The timing results are accumulated in the function "Aux_Record_Timing"
//
// Parameter   :  TotalT : Total simulation time
//                InitT  : Initialization time
//                OtherT : Elapsed time in all other parts (Aux_Record_Performance, Aux_Record_Timing, Aux_ResetTimer)
//-------------------------------------------------------------------------------------------------------
void Aux_AccumulatedTiming( const double TotalT, double InitT, double OtherT )
{

   const char FileName[]      = "Record__Timing";
   const char Comment_LB[][4] = { "Max", "Min", "Ave" };
   const int  NNewTimer       = 2;

   double Flu_P, Gra_P, Che_P, FixUp_P, Flag_P, Refine_P, Sum_P, MPI_Grid_P, dt_P, Aux_P, Corr_P, Output_P, LB_P;
   double Par_P, MPI_Par_P, Init_P, Other_P;
   double Flu_IB, Gra_IB, Che_IB, FixUp_IB, Flag_IB, Refine_IB, Sum_IB, MPI_Grid_IB, dt_IB, Aux_IB, Corr_IB, Output_IB, LB_IB;
   double Par_IB, MPI_Par_IB, Init_IB, Other_IB;
   double NewTimer_Acc[NNewTimer][3], Send[NNewTimer];

   double (*Recv)[NNewTimer] = new double [MPI_NRank][NNewTimer];
   double *Init_Acc          = NewTimer_Acc[0];
   double *Other_Acc         = NewTimer_Acc[1];

   FILE *File = ( MPI_Rank == 0 ) ? fopen( FileName, "a" ) : NULL;

   if ( OPT__TIMING_BALANCE )
   {
      Send[0] = InitT;
      Send[1] = OtherT;

      MPI_Gather( Send, NNewTimer, MPI_DOUBLE, Recv[0], NNewTimer, MPI_DOUBLE, 0, MPI_COMM_WORLD );
   }

   if ( MPI_Rank != 0 )    return;


// collect the initialization time
   if ( OPT__TIMING_BALANCE )
   {
      for (int t=0; t<NNewTimer; t++)
      {
         NewTimer_Acc[t][0] = __FLT_MIN__;
         NewTimer_Acc[t][1] = __FLT_MAX__;
         NewTimer_Acc[t][2] = 0.0;

         for (int r=0; r<MPI_NRank; r++)
         {
            NewTimer_Acc[t][0]  = MAX( NewTimer_Acc[t][0], Recv[r][t] );
            NewTimer_Acc[t][1]  = MIN( NewTimer_Acc[t][1], Recv[r][t] );
            NewTimer_Acc[t][2] +=                          Recv[r][t];
         }

         NewTimer_Acc[t][2] /= MPI_NRank;
      }
   }

   else
   {
      NewTimer_Acc[0][0] = InitT;
      NewTimer_Acc[1][0] = OtherT;
   }

   for (int t=0; t<NNewTimer; t++)
   for (int v=0; v<3; v++)          Sum_Acc[v] += NewTimer_Acc[t][v];


// get the percentage
   Flu_P      = 100.0*Flu_Acc     [0]/TotalT;
   Gra_P      = 100.0*Gra_Acc     [0]/TotalT;
   Che_P      = 100.0*Che_Acc     [0]/TotalT;
   FixUp_P    = 100.0*FixUp_Acc   [0]/TotalT;
   Flag_P     = 100.0*Flag_Acc    [0]/TotalT;
   Refine_P   = 100.0*Refine_Acc  [0]/TotalT;
   MPI_Grid_P = 100.0*MPI_Grid_Acc[0]/TotalT;
   dt_P       = 100.0*dt_Acc      [0]/TotalT;
   Output_P   = 100.0*Output_Acc  [0]/TotalT;
   Aux_P      = 100.0*Aux_Acc     [0]/TotalT;
   LB_P       = 100.0*LB_Acc      [0]/TotalT;
   Corr_P     = 100.0*Corr_Acc    [0]/TotalT;
   Par_P      = 100.0*Par_Acc     [0]/TotalT;
   MPI_Par_P  = 100.0*MPI_Par_Acc [0]/TotalT;
   Init_P     = 100.0*Init_Acc    [0]/TotalT;
   Other_P    = 100.0*Other_Acc   [0]/TotalT;
   Sum_P      = 100.0*Sum_Acc     [0]/TotalT;


// record
   fprintf( File, "****************************************************************************************" );
   fprintf( File, "**************************************\n" );
   fprintf( File, "*************                                     Accumulated timing results            " );
   fprintf( File, "                         *************\n" );
   fprintf( File, "****************************************************************************************" );
   fprintf( File, "**************************************\n" );
   fprintf( File, "Total Number of Steps : %ld\n", Step );
   fprintf( File, "Final Physical Time   : %13.7e\n", Time[0] );
   fprintf( File, "Total Simulation Time : %lf sec\n", TotalT );
   fprintf( File, "Timing Diagnosis      :\n" );
   fprintf( File, "----------------------------------------------------------------------------------------" );
   fprintf( File, "--------------------------------------\n" );
   fprintf( File, "%3s%5s %9s%9s%9s%9s%9s%9s%9s%9s%9s%9s%9s%9s%9s%9s%9s%9s%9s\n",
            "", "", "Flu_Adv", "Gra_Adv", "Che_Adv", "FixUp", "Flag", "Refine", "MPI_Grid", "dt", "Output", "Aux", "LB",
            "CorrSync", "Par", "-MPI_Par", "Init", "Other", "Sum" );

   if ( OPT__TIMING_BALANCE )
   {
      for (int v=0; v<3; v++)
      fprintf( File, "%3s%5s %9.2f%9.2f%9.2f%9.2f%9.2f%9.2f%9.2f%9.2f%9.2f%9.2f%9.2f%9.2f%9.2f%9.2f%9.2f%9.2f%9.2f\n",
               Comment_LB[v], "Time", Flu_Acc[v], Gra_Acc[v], Che_Acc[v], FixUp_Acc[v], Flag_Acc[v], Refine_Acc[v], MPI_Grid_Acc[v],
               dt_Acc[v], Output_Acc[v], Aux_Acc[v], LB_Acc[v], Corr_Acc[v], Par_Acc[v], MPI_Par_Acc[v], Init_Acc[v], Other_Acc[v],
               Sum_Acc[v] );

//    "max" imbalance = (Max-Ave)/Ave
      Flu_IB      = 100.0*( Flu_Acc     [0] - Flu_Acc     [2] ) / ( (Flu_Acc     [2]==0.0) ? 1.0 : Flu_Acc     [2] );
      Gra_IB      = 100.0*( Gra_Acc     [0] - Gra_Acc     [2] ) / ( (Gra_Acc     [2]==0.0) ? 1.0 : Gra_Acc     [2] );
      Che_IB      = 100.0*( Che_Acc     [0] - Che_Acc     [2] ) / ( (Che_Acc     [2]==0.0) ? 1.0 : Che_Acc     [2] );
      FixUp_IB    = 100.0*( FixUp_Acc   [0] - FixUp_Acc   [2] ) / ( (FixUp_Acc   [2]==0.0) ? 1.0 : FixUp_Acc   [2] );
      Flag_IB     = 100.0*( Flag_Acc    [0] - Flag_Acc    [2] ) / ( (Flag_Acc    [2]==0.0) ? 1.0 : Flag_Acc    [2] );
      Refine_IB   = 100.0*( Refine_Acc  [0] - Refine_Acc  [2] ) / ( (Refine_Acc  [2]==0.0) ? 1.0 : Refine_Acc  [2] );
      MPI_Grid_IB = 100.0*( MPI_Grid_Acc[0] - MPI_Grid_Acc[2] ) / ( (MPI_Grid_Acc[2]==0.0) ? 1.0 : MPI_Grid_Acc[2] );
      dt_IB       = 100.0*( dt_Acc      [0] - dt_Acc      [2] ) / ( (dt_Acc      [2]==0.0) ? 1.0 : dt_Acc      [2] );
      Output_IB   = 100.0*( Output_Acc  [0] - Output_Acc  [2] ) / ( (Output_Acc  [2]==0.0) ? 1.0 : Output_Acc  [2] );
      Aux_IB      = 100.0*( Aux_Acc     [0] - Aux_Acc     [2] ) / ( (Aux_Acc     [2]==0.0) ? 1.0 : Aux_Acc     [2] );
      LB_IB       = 100.0*( LB_Acc      [0] - LB_Acc      [2] ) / ( (LB_Acc      [2]==0.0) ? 1.0 : LB_Acc      [2] );
      Corr_IB     = 100.0*( Corr_Acc    [0] - Corr_Acc    [2] ) / ( (Corr_Acc    [2]==0.0) ? 1.0 : Corr_Acc    [2] );
      Par_IB      = 100.0*( Par_Acc     [0] - Par_Acc     [2] ) / ( (Par_Acc     [2]==0.0) ? 1.0 : Par_Acc     [2] );
      MPI_Par_IB  = 100.0*( MPI_Par_Acc [0] - MPI_Par_Acc [2] ) / ( (MPI_Par_Acc [2]==0.0) ? 1.0 : MPI_Par_Acc [2] );
      Init_IB     = 100.0*( Init_Acc    [0] - Init_Acc    [2] ) / ( (Init_Acc    [2]==0.0) ? 1.0 : Init_Acc    [2] );
      Other_IB    = 100.0*( Other_Acc   [0] - Other_Acc   [2] ) / ( (Other_Acc   [2]==0.0) ? 1.0 : Other_Acc   [2] );
      Sum_IB      = 100.0*( Sum_Acc     [0] - Sum_Acc     [2] ) / ( (Sum_Acc     [2]==0.0) ? 1.0 : Sum_Acc     [2] );

      fprintf( File, "\n" );
      fprintf( File, "%9s%8.3f%%%8.3f%%%8.3f%%%8.3f%%%8.3f%%%8.3f%%%8.3f%%%8.3f%%%8.3f%%%8.3f%%%8.3f%%%8.3f%%%8.3f%%%8.3f%%%8.3f%%%8.3f%%%8.3f%%\n",
               "Imbalance", Flu_IB, Gra_IB, Che_IB, FixUp_IB, Flag_IB, Refine_IB, MPI_Grid_IB, dt_IB, Output_IB, Aux_IB, LB_IB,
               Corr_IB, Par_IB, MPI_Par_IB, Init_IB, Other_IB, Sum_IB );
      fprintf( File, "%3s%5s %8.3f%%%8.3f%%%8.3f%%%8.3f%%%8.3f%%%8.3f%%%8.3f%%%8.3f%%%8.3f%%%8.3f%%%8.3f%%%8.3f%%%8.3f%%%8.3f%%%8.3f%%%8.3f%%%8.3f%%\n",
               Comment_LB[0], "Frac", Flu_P, Gra_P, Che_P, FixUp_P, Flag_P, Refine_P, MPI_Grid_P, dt_P, Output_P, Aux_P, LB_P,
               Corr_P, Par_P, MPI_Par_P, Init_P, Other_P, Sum_P );
   } // if ( OPT__TIMING_BALANCE )

   else
   {
      fprintf( File, "%3s%5s %9.2f%9.2f%9.2f%9.2f%9.2f%9.2f%9.2f%9.2f%9.2f%9.2f%9.2f%9.2f%9.2f%9.2f%9.2f%9.2f%9.2f\n",
               "", "Time", Flu_Acc[0], Gra_Acc[0], Che_Acc[0], FixUp_Acc[0], Flag_Acc[0], Refine_Acc[0], MPI_Grid_Acc[0], dt_Acc[0],
               Output_Acc[0], Aux_Acc[0], LB_Acc[0], Corr_Acc[0], Par_Acc[0], MPI_Par_Acc[0], Init_Acc[0], Other_Acc[0], Sum_Acc[0] );
      fprintf( File, "%3s%5s %8.3f%%%8.3f%%%8.3f%%%8.3f%%%8.3f%%%8.3f%%%8.3f%%%8.3f%%%8.3f%%%8.3f%%%8.3f%%%8.3f%%%8.3f%%%8.3f%%%8.3f%%%8.3f%%%8.3f%%\n",
               "", "Frac", Flu_P, Gra_P, Che_P, FixUp_P, Flag_P, Refine_P, MPI_Grid_P, dt_P, Output_P, Aux_P, LB_P,
               Corr_P, Par_P, MPI_Par_P, Init_P, Other_P, Sum_P );
   } // if ( OPT__TIMING_BALANCE ) .. else ...

   fprintf( File, "****************************************************************************************" );
   fprintf( File, "**************************************\n\n\n" );
   fclose( File );


   delete [] Recv;

} // FUNCTION : Aux_AccumulatedTiming



#endif // #ifdef TIMING
