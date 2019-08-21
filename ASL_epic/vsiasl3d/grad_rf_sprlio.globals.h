/*@Start***********************************************************/
/* GEMSBG Include File
 * Copyright (C) 1995 The General Electric Company
 *
 *      Include File Name:  grad_rf_grass.globals
 *      Developer:              T. Hlaban        Original for 5.5
 *
 * $Source: grad_rf_grass.globals.h $
 * $Revision: 1.0 $  $Date: 4/18/95 15:39:04 $
 *      prescan.globals.h        10/1/95 ghg
 */

/*@Synopsis
  This has global #defines for ipg & host
*/

/*@Description

*/

/*@End*********************************************************/

/* only do this once in any given compilation.*/
#ifndef  grad_rf_globals_sprlio_INCL
#define  grad_rf_globals_sprlio_INCL

#define MAX_RFPULSE 25
#define MAX_GRADX 20
#define MAX_GRADY 20
#define MAX_GRADZ 20

#define RF1_SLOT 0
#define RF_FREE1 1

#define GX_FREE 0

#define GY_FREE 0

#define GZ_FREE 0

#include "rf_Prescan.globals.h"

#endif 
