#include "GAMER.h"




//-------------------------------------------------------------------------------------------------------
// Function    :  Flag_WLMDwarfGalaxy
// Description :  Flag cells for refinement for the WLM dwarf galaxy simulation
//
// Note        :  1. Linked to the function pointer "Flag_User_Ptr" by "Init_TestProb_Hydro_WLMDwarfGalaxy()"
//                   to replace "Flag_User()"
//                2. Please turn on the runtime option "OPT__FLAG_USER"
//
// Parameter   :  i,j,k       : Indices of the target element in the patch ptr[ amr->FluSg[lv] ][lv][PID]
//                lv          : Refinement level of the target patch
//                PID         : ID of the target patch
//                Threshold : User-provided threshold for the flag operation, which is loaded from the
//                            file "Input__Flag_User"
//
// Return      :  "true"  if the flag criteria are satisfied
//                "false" if the flag criteria are not satisfied
//-------------------------------------------------------------------------------------------------------
bool Flag_WLMDwarfGalaxy( const int i, const int j, const int k, const int lv, const int PID, const double *Threshold )
{

   bool Flag = false;

// 1. Cell Position
   const double dh     = amr->dh[lv];                                                  // grid size
   const double Pos[3] = { amr->patch[0][lv][PID]->EdgeL[0] + (i+0.5)*dh,              // x,y,z position
                           amr->patch[0][lv][PID]->EdgeL[1] + (j+0.5)*dh,
                           amr->patch[0][lv][PID]->EdgeL[2] + (k+0.5)*dh  };

// flag cells within the target region [Threshold ... BoxSize-Threshold]
   const double EdgeL = Threshold[0];
   const double EdgeR = amr->BoxSize[0]-Threshold[0];    // here we have assumed a cubic box

   Flag |=  (  Pos[0] >= EdgeL  &&  Pos[0] < EdgeR  &&
               Pos[1] >= EdgeL  &&  Pos[1] < EdgeR  &&
               Pos[2] >= EdgeL  &&  Pos[2] < EdgeR     );
   if ( Flag )    return Flag;


// 2. Hydro CFL
// fluid
   const real   Dens    = amr->patch[ amr->FluSg[lv] ][lv][PID]->fluid[DENS][k][j][i];
   const real   MomX    = amr->patch[ amr->FluSg[lv] ][lv][PID]->fluid[MOMX][k][j][i];
   const real   MomY    = amr->patch[ amr->FluSg[lv] ][lv][PID]->fluid[MOMY][k][j][i];
   const real   MomZ    = amr->patch[ amr->FluSg[lv] ][lv][PID]->fluid[MOMZ][k][j][i];
   const real   Engy    = amr->patch[ amr->FluSg[lv] ][lv][PID]->fluid[ENGY][k][j][i];

// pressure, assuming no contribution from the magnetic field and passive scalars
   const bool   CheckMinPres_Yes = true;
   const real   Emag    = NULL_REAL;
   const real  *Passive = NULL;

   const real   Pres    = Hydro_Con2Pres( Dens, MomX, MomY, MomZ, Engy, Passive,
                                          CheckMinPres_Yes, MIN_PRES, Emag,
                                          EoS_DensEint2Pres_CPUPtr, EoS_GuessHTilde_CPUPtr, EoS_HTilde2Temp_CPUPtr,
                                          EoS_AuxArray_Flt, EoS_AuxArray_Int, h_EoS_Table, NULL );

// sound speed squared
   const real   Cs2     = EoS_DensPres2CSqr_CPUPtr( Dens, Pres, Passive,
                                                    EoS_AuxArray_Flt, EoS_AuxArray_Int, h_EoS_Table );

// time-step from the CFL condition
   const real   CFLx    = SQRT( Cs2 ) + FABS( MomX/Dens );
   const real   CFLy    = SQRT( Cs2 ) + FABS( MomY/Dens );
   const real   CFLz    = SQRT( Cs2 ) + FABS( MomZ/Dens );

   const double dt_CFL  = DT__FLUID * dh / FMAX( CFLz, FMAX( CFLy, CFLx ) );


// flag the cell when its dt_CFL is smaller than the threshold
   Flag |= ( dt_CFL < Threshold[1] );


   return Flag;

} // FUNCTION : Flag_WLMDwarfGalaxy
