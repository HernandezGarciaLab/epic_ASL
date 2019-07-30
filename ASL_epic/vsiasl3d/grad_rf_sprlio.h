/* *************************************
 * grad_rf_sprlio.h		1/30/2008 	ghg
 * This structure is used to track the 
 * rf heating, SAR heating, Grad coil heating,
 * grad amplifier heating.
 * ********************************** */

/* only do this once in any given compilation.*/
#ifndef  grad_rf_sprlio__INCL
#define  grad_rf_sprlio_INCL

/*

RF_PULSE rfpulse[MAX_RFPULSE] = {
#include "rf_Prescan.h"
};
*/

RF_PULSE rfpulse[MAX_RFPULSE] = {

/*
 {
     (int *)&pw_vsitag1,   // pw 
     (float *)&a_vsitag1,  // amp 
     1.0,           // abswidth 
     0.0855,           // effwidth 
     0.117,           // area 
     0.2,              // dtycyc 
     0.2,              // maxpw 
     502,                // num 
     0.2925,         // maxb1 (from vel_sim_09.m) 
     0.0855,       // max_int_b1_sq 
     0.0855,        // max_rms_b1 
     180,             // nom_fa 
     &flip_vsitag1,        // act_fa 
     400,            // nom_pw (just the bit with the refocuser 180)
     2500,           // nom_bw 
     PSD_APS2_ON + PSD_MPS2_ON + PSD_SCAN_ON,      // activity 
     0,                // reference (not used) 
     0,                // isodelay 
     1.0,              // scale 
     (int*)&res_vsitag1,   // res (not used) 
     0,                // extgradflag (not used) 
     (int *)&wg_vsitag1    // waveform generator (e.g., TYPRHO1 or TYPRHO2) 
  }, 


  {  (int *)&pw_rf1,   // pw 
     (float *)&a_rf1,  // amp 
     0.1353,           // abswidth 
     0.0647,           // effwidth 
     0.5767,           // area 
     0.1753,              // dtycyc 
     0.1335,              // maxpw 
     502,                // num 
     0.147,         // maxb1 
     0.0226221,       // max_int_b1_sq 
     0.0868372,        // max_rms_b1 
     180,             // nom_fa 
     &flip_rf1,        // act_fa 
     3000,            // nom_pw 
     3000,           // nom_bw 
     PSD_APS2_ON + PSD_MPS2_ON + PSD_SCAN_ON,      // activity 
     0,                // reference (not used) 
     0,                // isodelay 
     1.0,              // scale 
     (int*)&res_rf1,   // res (not used) 
     0,                // extgradflag (not used) 
     (int *)&wg_rf1    // waveform generator (e.g., TYPRHO1 or TYPRHO2) 
  }, 
*/
  {  (int *)&pw_rfdummy,   // pw 
     (float *)&a_rfdummy,  // amp 
     0.1353,           // abswidth 
     0.0647,           // effwidth 
     0.5767,           // area 
     0.1753,              // dtycyc 
     0.1335,              // maxpw 
     502,                // num 
     //0.147,         // maxb1 ( this is the max for the unused ref dummy pulse) 
     0.200,         // maxb1 ( this is the max for the unused ref dummy pulse) 
     0.0226221,       // max_int_b1_sq 
     0.0868372,        // max_rms_b1 
     180,             // nom_fa 
     &flip_rfdummy,        // act_fa 
     3000,            // nom_pw 
     3000,           // nom_bw 
     PSD_APS2_ON + PSD_MPS2_ON + PSD_SCAN_ON,      // activity 
     0,                // reference (not used) 
     0,                // isodelay 
     1.0,              // scale 
     (int*)&res_rfdummy,   // res (not used) 
     0,                // extgradflag (not used) 
     (int *)&wg_rfdummy    // waveform generator (e.g., TYPRHO1 or TYPRHO2) 
  }, 
#include "rf_Prescan.h"



};

/* #define MAX_ENTRY_POINTS 15 */
float maxB1[MAX_ENTRY_POINTS], maxB1Seq;

GRAD_PULSE gradx[MAX_GRADX] = {
  {
  }
};

GRAD_PULSE grady[MAX_GRADY] = {
  {
  }
};

GRAD_PULSE gradz[MAX_GRADZ] = {
  {
  }
};

#endif  
