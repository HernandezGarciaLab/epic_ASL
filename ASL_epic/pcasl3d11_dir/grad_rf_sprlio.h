/* *************************************
 * grad_rf_sprlio.h		1/30/2008 	ghg
 * This structure is used to track the 
 * rf heating, SAR heating, Grad coil heating,
 * grad amplifier heating.
 * ********************************** */

/* only do this once in any given compilation.*/
#ifndef  grad_rf_sprlio__INCL
#define  grad_rf_sprlio_INCL

RF_PULSE rfpulse[MAX_RFPULSE] = {
#include "rf_Prescan.h"
};

#define MAX_ENTRY_POINTS 15 
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
