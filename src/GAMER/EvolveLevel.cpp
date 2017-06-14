#include "Copyright.h"
#include "GAMER.h"

#ifdef TIMING
extern Timer_t *Timer_Flu_Advance[NLEVEL];
extern Timer_t *Timer_Gra_Advance[NLEVEL];
extern Timer_t *Timer_FixUp      [NLEVEL];
extern Timer_t *Timer_Flag       [NLEVEL];
extern Timer_t *Timer_Refine     [NLEVEL];
extern Timer_t *Timer_Lv         [NLEVEL];
extern Timer_t *Timer_GetBuf     [NLEVEL][8];
extern Timer_t *Timer_Par_Update [NLEVEL][3];
extern Timer_t *Timer_Par_2Sib   [NLEVEL];
extern Timer_t *Timer_Par_2Son   [NLEVEL];
#endif




//-------------------------------------------------------------------------------------------------------
// Function    :  EvolveLevel
// Description :  Advance all physical attributes one time-step
//
// Note        :  1. Each step contains TWO sub-steps, and each of which will evolve the solution
//                   at lv by 0.5*dTime
//
// Parameter   :  lv    : Target refinement level
//                dTime : Time interval to update the physical time at this level
//-------------------------------------------------------------------------------------------------------
void EvolveLevel( const int lv, const double dTime )
{

#  ifdef GRAVITY
   const bool   SelfGravity   = ( OPT__GRAVITY_TYPE == GRAVITY_SELF  ||  OPT__GRAVITY_TYPE == GRAVITY_BOTH );
#  endif
#  ifdef INDIVIDUAL_TIMESTEP
   const int    NSubStep      = 2;
#  else // shared time-step
   const int    NSubStep      = 1;
#  endif
   const double dTime_SubStep = dTime/NSubStep;
#  ifdef PARTICLE
   const bool   StoreAcc_Yes      = true;
   const bool   StoreAcc_No       = false;
   const bool   UseStoredAcc_Yes  = true;
   const bool   UseStoredAcc_No   = false;
   const bool   TimingSendPar_Yes = true;
#  endif

   double dt_SubStep, TimeOld, TimeNew;


// sub-step loop for INDIVIDUAL_TIMESTEP
   for (int SubStep=0; SubStep<NSubStep; SubStep++ )
   {
#     ifdef TIMING
      MPI_Barrier( MPI_COMM_WORLD );
      Timer_Lv[lv]->Start();
#     endif

//    1. calculate the evolution time-step
// ===============================================================================================
      dt_SubStep = Mis_dTime2dt( Time[lv], dTime_SubStep );
      TimeOld    = Time[lv];
      TimeNew    = Time[lv] + dTime_SubStep;

//    synchronize the time array (avoid round-off errors)
      if (  lv > 0  &&  Mis_CompareRealValue( TimeNew, Time[lv-1], NULL, false )  )
      TimeNew    = Time[lv-1];
// ===============================================================================================


//    2. update particles (prediction for KDK) and exchange particles
// ===============================================================================================
#     ifdef PARTICLE
      if ( OPT__VERBOSE  &&  MPI_Rank == 0 )
         Aux_Message( stdout, "   Lv %2d: Par_UpdateParticle (predict) %5s... ", lv, "" );

#     ifdef STORE_PAR_ACC
      TIMING_FUNC(   Par_UpdateParticle( lv, TimeNew, TimeOld, PAR_UPSTEP_PRED,
                                         (amr->Par->Integ == PAR_INTEG_EULER) ? StoreAcc_Yes    : StoreAcc_No,
                                         (amr->Par->Integ == PAR_INTEG_EULER) ? UseStoredAcc_No : UseStoredAcc_Yes ),
                     Timer_Par_Update[lv][0],   true   );
#     else
      TIMING_FUNC(   Par_UpdateParticle( lv, TimeNew, TimeOld, PAR_UPSTEP_PRED,
                                         StoreAcc_No, UseStoredAcc_No ),
                     Timer_Par_Update[lv][0],   true   );
#     endif

      if ( OPT__VERBOSE  &&  MPI_Rank == 0 )    Aux_Message( stdout, "done\n" );

      if ( OPT__VERBOSE  &&  MPI_Rank == 0 )
         Aux_Message( stdout, "   Lv %2d: Par_PassParticle2Sibling %9s... ", lv, "" );

      TIMING_FUNC(   Par_PassParticle2Sibling( lv, TimingSendPar_Yes ),   Timer_Par_2Sib[lv],   true   );

      if ( OPT__VERBOSE  &&  MPI_Rank == 0 )    Aux_Message( stdout, "done\n" );
#     endif
// ===============================================================================================


//    3. fluid solver
// ===============================================================================================
      const int SaveSg_Flu = 1 - amr->FluSg[lv];

      if ( OPT__VERBOSE  &&  MPI_Rank == 0 )
         Aux_Message( stdout, "   Lv %2d: Flu_AdvanceDt, counter = %8ld ... ", lv, AdvanceCounter[lv] );

      if ( OPT__OVERLAP_MPI )
      {
//       enable OpenMP nested parallelism
#        ifdef OPENMP
         omp_set_nested( true );
#        endif

//       advance patches needed to be sent
         TIMING_FUNC(   Flu_AdvanceDt( lv, TimeNew, TimeOld, dt_SubStep, SaveSg_Flu, true, true ),
                        Timer_Flu_Advance[lv],   false   );

#        pragma omp parallel sections num_threads(2)
         {
#           pragma omp section
            {
//             transfer data simultaneously
#              ifdef GRAVITY
               if ( SelfGravity )
               TIMING_FUNC(   Buf_GetBufferData( lv, SaveSg_Flu, NULL_INT, DATA_GENERAL, _DENS,  Rho_ParaBuf,
                                                 USELB_YES ),
                              Timer_GetBuf[lv][0],   true   );
#              else
               TIMING_FUNC(   Buf_GetBufferData( lv, SaveSg_Flu, NULL_INT, DATA_GENERAL, _TOTAL, Flu_ParaBuf,
                                                 USELB_YES ),
                              Timer_GetBuf[lv][2],   true   );
#              endif
            }

#           pragma omp section
            {
//             advance patches not needed to be sent
               TIMING_FUNC(   Flu_AdvanceDt( lv, TimeNew, TimeOld, dt_SubStep, SaveSg_Flu, true, false ),
                              Timer_Flu_Advance[lv],   true   );
            }
         } // OpenMP parallel sections

//       disable OpenMP nested parallelism
#        ifdef OPENMP
         omp_set_nested( false );
#        endif
      } // if ( OPT__OVERLAP_MPI )

      else
      {
         TIMING_FUNC(   Flu_AdvanceDt( lv, TimeNew, TimeOld, dt_SubStep, SaveSg_Flu, false, false ),
                        Timer_Flu_Advance[lv],   true   );

#        ifdef GRAVITY
         if ( SelfGravity )
         TIMING_FUNC(   Buf_GetBufferData( lv, SaveSg_Flu, NULL_INT, DATA_GENERAL, _DENS,  Rho_ParaBuf,
                                           USELB_YES ),
                        Timer_GetBuf[lv][0],   true   );
#        else
         TIMING_FUNC(   Buf_GetBufferData( lv, SaveSg_Flu, NULL_INT, DATA_GENERAL, _TOTAL, Flu_ParaBuf,
                                           USELB_YES ),
                        Timer_GetBuf[lv][2],   true   );
#        endif
      } // if ( OPT__OVERLAP_MPI ) ... else ...

      amr->FluSg    [lv]             = SaveSg_Flu;
      amr->FluSgTime[lv][SaveSg_Flu] = TimeNew;

      if ( OPT__VERBOSE  &&  MPI_Rank == 0 )    Aux_Message( stdout, "done\n" );
// ===============================================================================================


//    4. Poisson + gravity solver
// ===============================================================================================
#     ifdef GRAVITY
      const int SaveSg_Pot = 1 - amr->PotSg[lv];

      if ( OPT__VERBOSE  &&  MPI_Rank == 0 )
         Aux_Message( stdout, "   Lv %2d: Gra_AdvanceDt, counter = %8ld ... ", lv, AdvanceCounter[lv] );

      if ( lv == 0 )
         Gra_AdvanceDt( lv, TimeNew, TimeOld, dt_SubStep, SaveSg_Flu, SaveSg_Pot, SelfGravity, true, false, false );

      else // lv > 0
      {
         if ( OPT__OVERLAP_MPI )
         {
//          enable OpenMP nested parallelism
#           ifdef OPENMP
            omp_set_nested( true );
#           endif

//          advance patches needed to be sent
            TIMING_FUNC(   Gra_AdvanceDt( lv, TimeNew, TimeOld, dt_SubStep, SaveSg_Flu, SaveSg_Pot,
                           SelfGravity, true, true, true ),
                           Timer_Gra_Advance[lv],   false   );

#           pragma omp parallel sections num_threads(2)
            {
#              pragma omp section
               {
//                transfer data simultaneously
                  if ( SelfGravity )
                  TIMING_FUNC(   Buf_GetBufferData( lv, NULL_INT, SaveSg_Pot, POT_FOR_POISSON, _POTE,
                                                    Pot_ParaBuf, USELB_YES ),
                                 Timer_GetBuf[lv][1],   true   );

                  TIMING_FUNC(   Buf_GetBufferData( lv, SaveSg_Flu, NULL_INT, DATA_GENERAL,    _TOTAL,
                                                    Flu_ParaBuf, USELB_YES ),
                                 Timer_GetBuf[lv][2],   true   );
               }

#              pragma omp section
               {
//                advance patches not needed to be sent
                  TIMING_FUNC(   Gra_AdvanceDt( lv, TimeNew, TimeOld, dt_SubStep, SaveSg_Flu, SaveSg_Pot,
                                 SelfGravity, true, true, false),
                                 Timer_Gra_Advance[lv],   true   );
               }
            } // OpenMP parallel sections

//          disable OpenMP nested parallelism
#           ifdef OPENMP
            omp_set_nested( false );
#           endif
         } // if ( OPT__OVERLAP_MPI )

         else
         {
            TIMING_FUNC(   Gra_AdvanceDt( lv, TimeNew, TimeOld, dt_SubStep, SaveSg_Flu, SaveSg_Pot,
                                          SelfGravity, true, false, false ),
                           Timer_Gra_Advance[lv],   true   );

            if ( SelfGravity )
            TIMING_FUNC(   Buf_GetBufferData( lv, NULL_INT, SaveSg_Pot, POT_FOR_POISSON, _POTE,  Pot_ParaBuf,
                                              USELB_YES ),
                           Timer_GetBuf[lv][1],   true   );

            TIMING_FUNC(   Buf_GetBufferData( lv, SaveSg_Flu, NULL_INT, DATA_GENERAL,    _TOTAL, Flu_ParaBuf,
                                              USELB_YES ),
                           Timer_GetBuf[lv][2],   true   );
         } // if ( OPT__OVERLAP_MPI ) ... else ...

//       note that the current implementation of external potential does NOT use PotSg/PotSgTime
         if ( SelfGravity )
         {
            amr->PotSg    [lv]             = SaveSg_Pot;
            amr->PotSgTime[lv][SaveSg_Pot] = TimeNew;
         }

      } // if ( lv == 0 ) ... else ...

      if ( OPT__VERBOSE  &&  MPI_Rank == 0 )    Aux_Message( stdout, "done\n" );
#     endif // #ifdef GRAVITY
// ===============================================================================================


//    5. correct particles velocity and send particles to lv+1
// ===============================================================================================
#     ifdef PARTICLE
      if ( amr->Par->Integ == PAR_INTEG_KDK )
      {
         if ( OPT__VERBOSE  &&  MPI_Rank == 0 )
            Aux_Message( stdout, "   Lv %2d: Par_UpdateParticle (correct Lv %2d)... ", lv, lv );

#        ifdef STORE_PAR_ACC
         TIMING_FUNC(   Par_UpdateParticle( lv, TimeNew, TimeOld, PAR_UPSTEP_CORR,
                                            StoreAcc_Yes, UseStoredAcc_No ),
                        Timer_Par_Update[lv][1],   true   );
#        else
         TIMING_FUNC(   Par_UpdateParticle( lv, TimeNew, TimeOld, PAR_UPSTEP_CORR,
                                            StoreAcc_No,  UseStoredAcc_No ),
                        Timer_Par_Update[lv][1],   true   );
#        endif

         if ( OPT__VERBOSE  &&  MPI_Rank == 0 )    Aux_Message( stdout, "done\n" );

         if ( lv > 0 )
         {
            if ( OPT__VERBOSE  &&  MPI_Rank == 0 )
               Aux_Message( stdout, "   Lv %2d: Par_UpdateParticle (correct Lv %2d)... ", lv, lv-1 );

//          apply velocity correction for particles just travelling from lv to lv-1
#           ifdef STORE_PAR_ACC
            TIMING_FUNC(   Par_UpdateParticle( lv-1, TimeNew, TimeOld, PAR_UPSTEP_CORR,
                                               StoreAcc_Yes, UseStoredAcc_No ),
                           Timer_Par_Update[lv][2],   true   );
#           else
            TIMING_FUNC(   Par_UpdateParticle( lv-1, TimeNew, TimeOld, PAR_UPSTEP_CORR,
                                               StoreAcc_No,  UseStoredAcc_No ),
                           Timer_Par_Update[lv][2],   true   );
#           endif

            if ( OPT__VERBOSE  &&  MPI_Rank == 0 )    Aux_Message( stdout, "done\n" );
         }
      }

      if ( OPT__VERBOSE  &&  MPI_Rank == 0 )
         Aux_Message( stdout, "   Lv %2d: Par_PassParticle2Son %12s ... ", lv, "" );

      TIMING_FUNC(   Par_PassParticle2Son_AllPatch( lv, TimingSendPar_Yes ),   Timer_Par_2Son[lv],   true   );

      if ( OPT__VERBOSE  &&  MPI_Rank == 0 )    Aux_Message( stdout, "done\n" );
#     endif
// ===============================================================================================


//    6. additional physics
// ===============================================================================================


// ===============================================================================================


      Time_Prev     [lv] = TimeOld;
      Time          [lv] = TimeNew;
      AdvanceCounter[lv] ++;

      if ( AdvanceCounter[lv] >= __LONG_MAX__ ) Aux_Message( stderr, "WARNING : AdvanceCounter overflow !!\n" );


      if ( lv != TOP_LEVEL  &&  NPatchTotal[lv+1] != 0 )
      {

//       7. enter the next refinement level
// ===============================================================================================
#        ifdef TIMING
         MPI_Barrier( MPI_COMM_WORLD );
         Timer_Lv[lv]->Stop( false );
#        endif

         EvolveLevel( lv+1, dTime_SubStep );

#        ifdef TIMING
         MPI_Barrier( MPI_COMM_WORLD );
         Timer_Lv[lv]->Start();
#        endif
// ===============================================================================================


//       8. correct the data at the current level with the data at the next finer level
// ===============================================================================================
         if ( OPT__VERBOSE  &&  MPI_Rank == 0 )    Aux_Message( stdout, "   Lv %2d: Flu_FixUp %24s... ", lv, "" );

         if ( OPT__FIXUP_FLUX )
         TIMING_FUNC(   Buf_GetBufferData( lv, NULL_INT, NULL_INT, COARSE_FINE_FLUX, _FLUX_TOTAL, NULL_INT,
                                           USELB_YES ),
                        Timer_GetBuf[lv][6],   true   );

         TIMING_FUNC(   Flu_FixUp( lv, dt_SubStep ),   Timer_FixUp[lv],   true   );

#        ifdef LOAD_BALANCE
         if ( OPT__FIXUP_RESTRICT )
         TIMING_FUNC(   LB_GetBufferData( lv, amr->FluSg[lv], NULL_INT, DATA_RESTRICT, _TOTAL, NULL_INT ),
                        Timer_GetBuf[lv][7],   true   );
#        endif

         if ( OPT__FIXUP_FLUX  ||  OPT__FIXUP_RESTRICT )
         TIMING_FUNC(   Buf_GetBufferData( lv, amr->FluSg[lv], NULL_INT, DATA_AFTER_FIXUP, _TOTAL, Flu_ParaBuf, USELB_YES  ),
                        Timer_GetBuf[lv][3],   true   );

         if ( OPT__VERBOSE  &&  MPI_Rank == 0 )    Aux_Message( stdout, "done\n" );
// ===============================================================================================

      } // if ( lv != TOP_LEVEL  &&  NPatchTotal[lv+1] != 0 )


//    9. flag the current level and create patches at the next finer level
// ===============================================================================================
      if ( lv != TOP_LEVEL  &&  AdvanceCounter[lv] % REGRID_COUNT == 0 )
      {
//       flag
         if ( OPT__VERBOSE  &&  MPI_Rank == 0 )    Aux_Message( stdout, "   Lv %2d: Flag %29s... ", lv, "" );

#        ifdef LOAD_BALANCE
         TIMING_FUNC(   Flag_Real( lv, USELB_YES ),       Timer_Flag[lv],   true    );
#        else
         TIMING_FUNC(   Flag_Real( lv, USELB_NO ),        Timer_Flag[lv],   false   );

         TIMING_FUNC(   MPI_ExchangeBoundaryFlag( lv ),   Timer_Flag[lv],   false   );

         TIMING_FUNC(   Flag_Buffer( lv ),                Timer_Flag[lv],   true    );
#        endif

         if ( OPT__VERBOSE  &&  MPI_Rank == 0 )    Aux_Message( stdout, "done\n" );


//       refine
         if ( OPT__VERBOSE  &&  MPI_Rank == 0 )    Aux_Message( stdout, "   Lv %2d: Refine %27s... ", lv, "" );

         TIMING_FUNC(   Refine( lv, USELB_YES ),   Timer_Refine[lv],   false   );

         Time          [lv+1]                     = Time[lv];
         amr->FluSgTime[lv+1][ amr->FluSg[lv+1] ] = Time[lv];
#        ifdef GRAVITY
//       note that the current implementation of external potential does NOT use PotSg/PotSgTime
         if ( SelfGravity )
         amr->PotSgTime[lv+1][ amr->PotSg[lv+1] ] = Time[lv];
#        endif

#        ifdef LOAD_BALANCE
         TIMING_FUNC(   Buf_GetBufferData( lv,   amr->FluSg[lv  ], NULL_INT, DATA_AFTER_REFINE, _TOTAL,  Flu_ParaBuf,
                                           USELB_YES ),
                        Timer_GetBuf[lv][4],    false   );
#        ifdef GRAVITY
         if ( SelfGravity )
         TIMING_FUNC(   Buf_GetBufferData( lv,   NULL_INT, amr->PotSg[lv  ], POT_AFTER_REFINE,  _POTE, Pot_ParaBuf,
                                           USELB_YES ),
                        Timer_GetBuf[lv][5],    false   );
#        endif
#        endif // #ifdef LOAD_BALANCE

         TIMING_FUNC(   Buf_GetBufferData( lv+1, amr->FluSg[lv+1], NULL_INT, DATA_AFTER_REFINE, _TOTAL,  Flu_ParaBuf,
                                           USELB_YES ),
                        Timer_GetBuf[lv][4],    true    );
#        ifdef GRAVITY
         if ( SelfGravity )
         TIMING_FUNC(   Buf_GetBufferData( lv+1, NULL_INT, amr->PotSg[lv+1], POT_AFTER_REFINE,  _POTE, Pot_ParaBuf,
                                           USELB_YES ),
                        Timer_GetBuf[lv][5],    true    );
#        endif

//       must call Poi_StorePotWithGhostZone AFTER collecting potential for buffer patches
#        ifdef STORE_POT_GHOST
         if ( SelfGravity )
         TIMING_FUNC(   Poi_StorePotWithGhostZone( lv+1, amr->PotSg[lv+1], false ),   Timer_Refine[lv],   false   );
#        endif

#        ifdef TIMING
         Timer_Refine[lv]->WorkingID ++;
#        endif

         if ( OPT__VERBOSE  &&  MPI_Rank == 0 )    Aux_Message( stdout, "done\n" );

         if ( OPT__PATCH_COUNT == 2 )  Aux_Record_PatchCount();
#        ifdef PARTICLE
         if ( OPT__PARTICLE_COUNT == 2 )  Par_Aux_Record_ParticleCount();
#        endif

      } // if ( lv != TOP_LEVEL  &&  AdvanceCounter[lv] % REGRID_COUNT == 0 )
// ===============================================================================================

#     ifdef TIMING
      else
      {
         Timer_Flag  [lv]   ->WorkingID ++;
         Timer_Refine[lv]   ->WorkingID ++;
         Timer_GetBuf[lv][4]->WorkingID ++;
         Timer_GetBuf[lv][5]->WorkingID ++;
      } // if ( lv != NLEVEL-1  &&  AdvanceCounter[lv] % REGRID_COUNT == 0 ) ... else ...
#     endif

#     ifdef TIMING
      MPI_Barrier( MPI_COMM_WORLD );
      Timer_Lv[lv]->Stop( true );
#     endif

   } // for (int SubStep=0; SubStep<NSubStep; SubStep++ )

} // FUNCTION : EvolveLevel

