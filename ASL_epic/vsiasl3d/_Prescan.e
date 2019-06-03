/*
 * GE Medical Systems
 * Copyright (C) 1993-1998 The General Electric Company
 * 
 * $Source: %M% $
 * $Revision: %I% $  $Date: %G% %U% $
 * 
 * Contains the epic code to support Prescan entry points for all types
 * of acquisition databases. The current entry points are cfl, mps1 & cfh.
 *
 * A brief description of the routines in this module:
 *
 *	PScvs: 
 *		Location - at the end of the host cv declarations
 *		Parameters -
 *
 *	        Description - define needed prescan CVS
 *
 *	PScvinit: 
 *		Location - In cvinit - at the end
 *		Parameters -
 *
 *	        Description - Load the RFPULSE structure
 *
 *	PScveval: 
 *		Location - In cveval - before minseqcoil calls
 *		Parameters -
 *
 *	        Description - Sets up gradient amplitudes
 *
 *	PSpredownload: 
 *		Location - In predownload - after initial RF peakB1 &
 *		                            orderslice & entrytabinit
 *		Parameters -
 *
 *	        Description - scales RF amplitudes, initializes timing
 *		              variables
 *
 *	PShost: 
 *		Location - After predownload - in host section
 *		Parameters -
 *
 *	        Description - has source code for host prescan routines
 *
 *	PSrspvar:
 *		Location - rspvar
 *	        Parameters -
 *
 *		Description - defines rspvars
 *
 *	PSpulsegen:
 *		Location - pulsegen
 *	        Parameters -
 *
 *		Description - Creates pulses for cfl, mps1 and cfh.
 *		 This needs the rsp variable cscfh_satindex
 *		 properly set in the calling pulse sequence 
 *		 if PSD_CFH_CHEMSAT is defined.
 *
 *	PScore:
 *	      	Location - rsp
 *		Parameters -
 *
 *		Description - executable code for cfl, mps1 and cfh.
 *
 *	PSeplist:
 *	      	Location - rsp
 *		Parameters -
 *
 *		Description - contains entry point ascii for simulation
 *
 *	PSipg: 
 *		Location - After core - in ipg section
 *		Parameters -
 *
 *	        Description - has source code for ipg prescan routines
 * 
 * Language : EPIC / ANSI C
 * Author   : Larry Ploetz & Al Li
 * Date     : 16/Mar/1993
 */
/* do not edit anything above this line */

/*
  Author  Date		Comments
------------------------------------------------------------------------------
  LP	   3/16/93  rewrite for simplification
  AL       4/02/93  Add autoshim entry point
  LP       5/12/93  SPR 16721 - 3d thickness fixed to 10
		            pw_sdixon2 initialized after cffield
  LP       6/10/93  added rspent initialization
  MHN      8/05/93  Made HSI compatible changes for optramp()
  MHN      8/11/93  Added amppwgrad for gykcfl/h - and changed
                    TRAPEZOID calls to TYPNDEF
5.5  MMS  11/22/93  Caught up 54 changes from HSI conversion
          CHANGES INCLUDED ARE:
             AL 8/11/93 SPR 17601 - make sure Autoshim don't        
                calculate sequence timing based on oblique plane    
             SS 8/11/93 SPR 17592 - mpsfov added and limited to     
                80 to avoid a_gxwmps1 overranging.                  
             AL 8/23/93 Change the number of concurrent             
               gradient in AutoShim to 1                           
             YS 9/08/93        Add MT CFH pulse                     
             YS 10/06/93 Change seq_type == TYP3D check to          
                opimode == PSD_3D or PSD_3DM check,                 
                so that all 3D prescans use 10mm slice.             
             LP 11/12/93 Changed cfh_tr for IR on change gscale_rf0 
5.5  MMS  12/8/93  HSI conversion continued, 
		SINSUSOID to TRAPEZOID2 and ampwencode to ampwencodet,
		pw_rampx to loggrd.xrt and pw_rampz to loggrd.zrt
5.5  MMS  12/10/93   chagned setiamp to setiampt for gylas and gylras
5.5  BES  01/04/94 Add #define amp_ and pw_killer to remove hard coded values
		   from amppwgrad calls.
    PL     01/13/94 amppwgrad mods; optramp fixed in FastTG section
                    trapezoidal phase encodes for autoshim fixed
                    initialize pulse widths in FastTG
    JM     01/21/94 more amppwgrad mods, deleted FTGdummy_area
    SS     02/05/94 Changed parameters in ACQUIREDATA call.
    JM     02/06/94 Argument changes for amppwgx1, amppwgz1, amppwlcrsh
                      and amppwencodet -> amppwtpe for AutoShim
    YS     03/14/94 Update with 5.4 spr fixes from 5.4.2. sccs 1.14.
    BES    03/15/94 Add new arg to ACQUIREDATA for xtr_pos parameter.
    PL     03/19/94 Add pscR1, pscR2, pscCF, and pscTG rsp vars
    PL     03/21/94 Add as_newgeo; obloptimize call line updated
                    New argument to loaddab calls
    JM     04/01/94 Added as_slquant declaration and changed
                    asslquant references to as_slquant.
    BES    04/22/94 Typo in fasttg on the loaddab call.
    PL     06/10/94 Reference asloggrd througout autoshim sections
                     (MRIge20350)
    BES    06/22/94 Change looping structure to include multi slab when
		    determining the prescan slice.  Change fasttg to
		    use PStloc instead of rpsinfo and to use 0 slice
		    offset in startseq.  This is consistent with other
	 	    entry points and avoids an error with 3d multi slab
		    and startseq (MRIge2023).
    BES    06/30/94 Define PSrot[9] as ipgexport in epic.h.  Need to do
		    this because we need to run scalerotmats on this
		    matrix on the host side and then setrotatearray
		    on the ipg side.  Initialize the matrix and scale
		    it in PSpredownload.  Setrotatearray in PSinit.
		    Remove temporary cvs PSrot0-9 and the PSrsprot. 
    LP     07/29/94 MRIge20781 - added minimum te times for CFH & PS1
                    sequencing to avoid .5T rf amp duty cycle problems
    LP     08/16/94 Change to base amp & pw selection on area & not
                    hard code them.
    BES    08/24/94 Reposition x and z killers in autoshim based on the
		    decay of readout instead of gy1ras.  With SR160 we
		    were seeing overlaps of the killers and readout.
    FHE    09/12/94 Fixed MRIge21243: Wrong argument in attenflag call.
7.0 YI     10/13/94 Modified for cerd and solaris 2.
    YH     10/24/94 Changed type short to int.
    YH     10/28/94 Cahnged F**_SLOT to filter*_slot.
    YI     11/15/94 Added changes of psdsource:55 94fw45.5 (5.5 BES  11/08/94.)
    YI     12/07/94 Change opthickPS to 3.5 mm if G-coil is switchable SGC.
                    Changed filter in Fast TG for cerd. Changed FTGopslthickz3
                    to 7mm if cfgradcoil equals 102 to avoid internal error that
                    the sum of time exceeds FTGtau1 with Vectra SGC coil. This
                    change is temorary because product VMX does not use
                    Vectra SGC coil.
    YI    03/06/95  Changed data type of asrottemp from short to int.
    YI    03/13/95  Added changes for new system safety.
    YI    03/23/95  Added changes of :55 95fw11.5b(to RJL 03/17/95.)
    YI    04/04/95  Added changes of :55 95fw14.1(to sccs 1.67 mgh 950402)
    YI    04/17/95  Changed data type of rsp variables for auto shim.
    YO    05/02/95  Changed CFH for separating the fat peak from the water peak
                    (Field strength <= 5000)
    YI    05/09/95  Changed elements in optramp() for pw_gzrf1mps1a
		    from loggrd.tz_xz to loggrd.tz
    NM    05/22/95  Added cffield support for 0.2T(porfile)
    YI    06/08/95  55 merged (King   05/26/95  MRIge23499)
		    Support 1:2 pulse width TG prescan.
    YI	  06/15/95  TG prescan will be done in axial plane if the system
		    is VMX.
   RJL    08/25/95  Merge with :55 95fw34.5 (QT 8/24/95)
   RJL    08/31/95  Merged with :55 (95fw35.3 MGH)
   RJL    10/04/95  Merged with :55(s95fw38.7) and :mrp(95fw39.5)
                    Change syncoff(&seqcore) to syncoff(&seqautshm) per LP.
		    It was incorrect.
   RJL    10/18/95  Made frequency offset dependant on entry point. Changes
                    to setfrequency and receive frequency. Added PSfreq_offset
		    array.
   RJL    10/28/95  Commented out YI changes (6/15/95) from PScveval() which caused
                    pw_rf2mps1 to be reinitialized after pulse stretched. This is
		    bad, and causes download failures for high patient weights.

   SVB    11/30/95  In sliceselz call changed `start time' to `RUP_GRD(start_time+pw_gzrf1cfha)'
                    -MRIge28625 fix.
                    The sliceselz call was setting the start time of the cfh pulse
                    to end of pw_gzrf1cfha instead of the beginning of pw_gzrf1cha.
   RJL    01/17/96  MRIge28581-Remove VMX use of axial TG if system is VMX. (6/15/95) Found
                    to be buggy and causes our system to always do axial prescan
		    also. Removed all code pertaining to this addition.
   RJL    02/26/96  MRIge30217, MRIge30365 - Change asfov to 320 for head coil per SRS
                    Autoshim Enhancements.

  JDM  08-Mar-96    MRIge30640
		    When cs is on and ir is off, need to account for 
		    pw_gzrf1cfha - already done in fix for MRIge28625

		    MRIge30641
		    Need check on sl thk. Otherwise causes dwnld failures 
		    w/0.8 G/cm. Assumes that gzrf1 gives the worst case min 
		    sl thk of all psc pulses. 

                    MRIge30642
                    Need check on fov.  Otherwise causes dwnld failures 
		    w/0.8 G/cm. 

                    MRIge30645
  		    Need to change posstart to eliminate dwnld failures with 
		    .8 G/cm.

   VB  28-Jun-96    MRIge33520: Fix suggested by Kathy Bahner
                    Incorrect number of baselines and disdaqs are played out during
		    autoshim prescan. This was resulting in incorrect baseline correction in 
		    subsequent scan. The fix is to change 
		    if ((view >= -asbaseline)&&(excitation > 0))
		    to
		    if ((view > -asbaseline)&&(excitation > 0))
		    See spr for more explanation.

  VB  21-Feb-97     Switched Inversion pulse to adiabatic pulse.
                    The inversion rf0cfh pulse definition and its pulsegen code
		    have been changed.

  VB  27-Mar-97     MRIge38829: Reset titime to opti if an invalid titime is
                    entered.   

sccs1.41 18-Oct-97  CMC   Fix for MRIge42359  CMC

  GFN 21-Oct-1997   Changed all references to ERMES_DEBUG to PSD_HW
                    when checking for Hardware compilations. The only
                    valid use of ERMES_DEBUG should be when setting
                    the use_ermes variable.
  MS  30-Oct-1997   Removed the RTrspvar inline. These variables are
                    declared in epic.h

sccs1.43 06-Jan-98  CMC   Internatioalization

sccs1.44 07-Jan-98  CMC   Syntax error

sccs1.48 17-Feb-98  BJM   MRIge43971: Modified handling of PSrot since
                          this is now a 2D array.

******************  Initial CV-1 version ******************

Version    Date       Person        Comment
-----------------------------------------------------------------
sccs1.48    12/12/97  JAP     Merged in changes from Lx-2.

sccs1.49    01/21/98  BJM     MRIge43968: modifed handling of PSrot in
                              PSpredownload since this variable is now a
                              matrix instead of an array.

sccs1.50    01/21/98  GFN     Changed cast in PSrot from (INT) to (int).
                              Removed (short) casts in setrotatearray()
                              and settriggerarray().
sccs1.53    06/08/98  AKG     Removed declaration of prescan rspvars pscR1,
                              pscR2, pscCF, pscTG from @rspvar 
                              and added them in epic.h
                              to ensure availability to all PSDs.
                              Fix for spr MRIge46008

	    09/02/98  SGT    MRIge45995 - Select the slice most close to the
 			     center of the prescribed region, instead of the
			     slice most close to the iso-center.

            10/14/98  GFN    MRIge47485 - Moved the @pg PSipg section right
                             after the @pg PSpulsegen section to keep all the
                             @pg sections together before the @rsp sections.
                             Converted to ANSI.

         13-Oct-1998  RJF    Fixing n32 clearcase build no image problem.
                             Was a coding error in PSfilter(). This has to
                             be a void function which makes use of the CV
                             from the PSD as an index to the array of filter
                             spec structures.

         28-Oct-1998  JAP    MRIge48530 - Fixed typo in 1 T cffield check.

         07-Dec-1998  PRA    MRIge47735 - Faster prescan cv no longer read
                             from file instead, it is set using the cv
                             cffastprescan.

         11-Dec-1998  GFN    MRIge49699 - Removed the set of pos_start and
                             start_time in FTGpulsegen(). pos_start is a CV
                             and should not be modified in IPG code. The
                             start_time variable is local and unused in
                             that function.

         07-Jan-1999  JAH    MRIge50128 - moved PSrsptrigger declaration
                             from inside PSinit() to top of PSCore

                             MRIge50129 - Added SPECTROSCOPY flag so
                             spectroscopy specific portions of prescan
                             could be used instead of the standard
                             autoshim and cfh routines. This was
                             required because necessary rsp functions in 
                             Prescan.e were moved to the PScore section
                             from the PSipg section. Spectroscopy
                             sequences did not inline PScore because they
                             used their own versions of cfh and autoshim.
                             In order to pick up the other portion of
                             PScore, the spectroscopy PSDs now inline PScore
                             leaving out the spectroscopy specific portions.

         25-Mar-1999  RAK    Changed to use update time macros instead of
                             hardcoded numbers.

         02-Apr-1999  AMJ    MRIge52384 - Cvs te_as and tr_as are updated only
                             for fastprescan.

         16-Apr-1999  PRA    MRIge52773 - Initialized the filter slots
                             for Autoshim and FTG in PSfilter section .

         06-May-1999  AMJ    MRIge42428 - Added maxte calculation for mps1. 
                             Mismatch was that ps1_tr was set to 2s and if max allowed 
                             te which is 2s is selected there would be download failure
                             due to increase in seqlenght beyond ps1_tr.

         12-Dec-2000  RJF    LxMGD conversion. This version is not backward compatible with IPG
                             MRIge63434 caputures this. To see what changed, see MRIge63432.

         22-Sep-1999  RAK    MRIge55994 - Added the SDL header files to allow
                             tools PSDs to compile.

         13-Oct-1999  AF     MRIge56170 - Moved optramp for pw_gzrf3ftg function from
	                     after the area_g2ftg function to before the area_g2ftg 
			     function.
         24-Apr-2000  MM     Merge from MFO.

         14-Jul-2000  BSA    MRIge60975 Compatibility with 16us WARP for nMR.

         25-Sep-2000  SK     YMSmr02591 - change PSD_OBL_OPTIMAL to PSD_OBL_RESTRICT

         05-Dec-2000  BWL    MRIge63262 - fixed gradient amplitude calculation for
                             opspf=1.

         19-Apr-2001  RAK    MRIge65870 - Added PScveval1( te ) to eliminate Prescan
                             using the opte to set the TE of the Prescan entry points.
                             The PScveval1 input argument is the TE that is used for 
                             the Prescan entry points. PScveval1() is called
                             by PScveval() with an input argument of exist(opte).

         05-Jun-2001  RAK    MRIge66667 - GE ssfse protocol download fails.

         06-Jun-2001  RAK    MRIge66767 - oppseq logic problems in Prescan.e and Inversion.e for ssfse

         28-Jan-2002  DCZ    The routine optramp for gzrf2ftga needs to be 
                             executed before the amplitude and pulse width
                             of the gzrf2ftga pulse are used for the 
                             calculation of area_g1ftg.

         10-Jun-2002  JAH    MRIge75690 - change the receive coil updates
                     /RJF    to the entry point table to a function call
                             that updates the previous entries plus the
                             new MGD USSM multi-coil driver control fields.

         13-Jun-2002  DCZ    MRIge76008: Merge the following SPR fix from 91.merge.
                             MRIge72999, MRIge74381 - TG & CFH Prescan entry points 
                             were having flat lines during LongTE prescriptions. SCAN_TR 
                             entry point alone needs to have the eff_te
                             or opte dependency. Removed PScveval1() & retained
                             the pre-9.0 PScveval(). Removed the dependencies 
                             of ps_te/opte in the MPS1/APS1 & CFH entry points.

         31-Oct-2002  SP     MRIge78369:  Change asxres and asyres to 128x64 when
                             fast prescan is enabled for 1.5T.  We found much better
                             fat suppression with this higher autoshim resolution.

         13-Sep-2001  GE     flip_rf0cfh stretch modification

         18-Nov-2001  GE     modifications of APS1/MPS1 for 3T body coil
                             (aps1_mod)

         03-May-2002  SK     MRIge74630 - Changing filter for CFH to 1kHz bandwidth

         8-June-2002  LS     MRIge73651 - Changing ftgtr from 2s to 4s to fix fast
                             tg failure problem.

         14-June-2002  LS    MRIge75192 - wrong cfh_acq_length calcuation cuased
                             prescan failure. cfh is using 1KHz filter, so the 
                             cfh_acq_length should be 2x tdaq of a 2K filter.

         22-Aug-2002  NDG    Merge VH3 with MGD for Corona 3T.

         06-Dec-2002  AMR    MRIge79721: Merge from 9.1 to MGD for MGD2.

         01/03/2003   AP     Added bilateral multislab prescan functionality 

         04/17/2003   AP     MRIge82455 - Initialize psc_vol_index to 0 for multi
                             volume prescan BBA feature

         08/01/2003   HD     Added Check for Gradient Overlapping issues in FTG
                             when derating B1

         02/24/2004     LS   MRIge90312 - adding IR pulse to do fat sat for cfh.
                                1. PSir is on for all cases 
                                2. cfh_ti is 120ms, cfh_tr 1.5sec, cfh_te 50ms
                                3. prescan will shift cf by 220Hz if water sat is ON      
                                4. added PSslice_inc for prescan cfh/ftg slice index.

         02/24/2004     LS   MRIge92242 - IR pulse for CFH is not stretched proporly, causes download failuer.

         05/17/2004     ZL   MRIge93862 - Entrypoint table needs to be updated for all
                                          phased array coils for Autoshim and FTG.

         05/19/2004     LS   MRIge93237 - start time error when mt and cs sat ON at the same time .

         05/22/2004     LS   MRIge94017 - CFH failure caused by cssat cfh index mistake.
         
         6/25/2004     LS   MRIhc01181 - Reduce cfh TE from 50ms to 30ms for 3T. 

         8/19/2004     AP   MRIhc02472 - Use the location of the scan volume to set the 
                            center frequency for one prescan volume case.
         9/14/2004      LS  MRIhc03365 - replace cfhtun with showfp.

         10/6/2004      LS  MRIhc03942 - VERSE request shim. set pidoshim to on for
                                VERSE type scan (opvrg).
12.0	 20-Oct-2004  HH     Added changes for MRIhc03410 - Fix axial orientation and location
                             for 3-Plane for CFl, CFH, FTG and APS1

14.0     02-Mar-2005  LS   MRIhc06378 - make psd_grd_wait, psd_rf_wait part of the 
                                                         ta_180 calculation. 
14.0	 08-Apr-2005  LS   MRIhc06679,06683 - Fix ps1, ftg readout axis and offset.

14.0	 06-May-2005  HH   MRIhc07062 Fix asfov, mpsfov and FTGfov to FOV_MAX for knee, wrist coils

14.0     14-Jun-2005  SVR  MRIhc07934: SWIFT changes phase-1.

14.0     27-May-2005  ARI  Add a new entry point for receiver noise - called RCVN
                           all new variable and function names include word "RCVN"
                           therefore search "RCVN" for changes. 

14.0     01-Jul-2005 SXZ/RBA MRIhc08321:Localized Volume Excitation with non-selective IR PRESS
                            sequence for CFH Entry Point. Henceforth this will be the new CFH 
                            Entry Point for all Pulse Sequences except Spectroscopy. 
                            Spectroscopy will be addressed later. Shim Volume is a must for the
                            new localized CFH to be run, else it will run the STIR (12.0) CFH.
                            PRESSCFH_NONE  will be used when the Shim Vol is not used.
                            PRESSCFH_SLICE will be used when the Shim Vol is used and the PSD 
                                           is a 2D PSD.
                            PRESSCFH_SLAB  will be used when the Shim Vol is used and the PSD is
                                           a 3D Slab.
                            PRESSCFH_SHIMVOL will be used when the Shim Vol is used and the
                                           (1) PSD is either a 2D/3D and does not intersect 
                                           the Shim Volume, This is determined by the geometry 
                                           calculations.
                                           (2) The Prescription is a Radial prescription (opcoax)
                                           is 0. This is determined by the opimode and opcoax.

14.0    11-Jul-2005  SVR   Changes for SWIFT prescan.

14.0    15-Jul-2005  HKC   MRIhc08595: Introduced echo1ptcfh for echo1 CFH 
                           filter output points; corrected two epic_error 
                           message typo.

14.0    02-Aug-2005  SVR   SWIFT rot matrix changes.

14.0    19-Aug-2005  SVR   SWIFT slice trigger changes.

14.0    18-Nov-2005  ARI   Set max RCVN number of points to 1K for 32ch fix

14.0    20-Dec-2005  CRM   MRIhc09608: For SWIFT, always setup prescan volume
                           and coil when switching entry points.  Also switch
                           coils during MPS1 since prescan uses surface coil
                           receive.

14.0    07-Feb-2006  CRM   MRIhc13033: Add SwiFT support for single connector 
                           PV Coil using scn provided biasEnable and 
                           switchSelects.

14.0    16-Feb-2006  RBA   MRIhc11621: Initialzied Variables for CFH

14.0    28-Feb-2006  HKC   MRIhc08633: Set epprexres for CFH entry point to 
                           echo1ptcfh.

14.0    01-Mar-2006  ARI   MRIhc13698: Set max RCVN bandwidth to 62.5kHz

HDe     26-Apr-2006  YI    YMSmr09211: Fixed scaling problem on rotation matrix
                           in PRESSCFH_SHIMVOL mode.

14.5    22-Jun-2006  Teja  coilInfo related changes 

22.0    31-Dec-2009 VSN/VAK MRIhc46886: SV to DV Apps Sync Up

22.0    11-Feb-2010 LS     MRIhc47602: fix for noise cal STD issue. Grad crushers
                           are added before noise cal acquisition to kill unwanted
                           MR signal from previous prescan entry point(s). 

22.0    25-Feb-2010 AKR/KN MRIhc48093: Added support for HD coil switching.

SV      19-Jan-2010  MHi   GEHmr03577: Subtract specir_delay to avoid 
                           download failure in case of SPECIAL with long TI.
                           For compatibility with normal ChemSat, this is
                           distinguished by PSD_CFH_CHEMSAT_SPECIAL,
                           which is newly defined in ChemSatSpecIR.e.

22.0    29-Apr-2010  AE    MRIhc49539: check current nucleus against coil DB
                           to prevent switching of nucleus.

23.0    12-Jan-2011  AE  MRIhc54366: introduced CV CFLxres to control number 
                         of points in CFL entry point. Set defaults for
                         CFLxres and echo1bwcfl to accommodate B0 drift.
                         Doubled CFL RF1 bandwidth to reduce B0-drift 
                         induced slice offsets.

23.0    01-Apr-2011  TAC MRIhc54029: Added support for Burst Mode.

23.0    06-Sep-2011  LS  HCSDM00096222: Enable dual shim localized CFH mode (slab).

HD23.0  28-Jul-2011  NB  HCSDM00089407: rcvn_xres to be changed for HD23 program

HD23.0  22-Aug-2011  AE  HCSDM00084103: removed function ReadFtgTr(). (23.0: HCSDM00125966)

23.0    07-Feb-2012  LS  HCSDM00113265: Enable local rec coil for Unilateral breast prescan TG.
 
24.0    25-Jul-2012  LS  HCSDM00143679: Add derating factor for Prescan

24.0    24-May-2012  AE  HCSDM00136895: added express TG.

24.0    26-Sep-2012  AE  HCSDM00155770: static defect fixes.

24.0    04-Oct-2012  AE  HCSDM00161809: calculate minimum FTGtau1 value.

24.0    28-Sep-2012  AE  HCSDM00157626: avoid stretching of XTG Fermi
                         pulses by scaling down their flip angle, as needed.

24.0    08-Oct-2012  AE  HCSDM00163151: fgre 3plane auto prescan failed
                         BECAUSE CFH was using Sag plane instead of Ax plane.
                         PSslice_num was set differently in cveval and
                         predownload. Moved the PSslice_num assignment to
                         PScvinit().

24.0    07-Dec-2012  HS  HCSDM00174876: correct amplitude setting of rf3xtg in
                         XTGpredownload(), because calcPulseParams() in previous
                         efgre3d.e appeared below XTGpredownload() and ia_rf3xtg
                         was overwritten in calcPulseParams().

24.0    29-Aug-2013  LS/HS  HCSDM00232515: Download failed on Silent MRA

24.0    16-Sep-2013  LS  HCSDM00235363: change XTGfov to opfov if aps1_mod = 0 to
                        reduce signal from outside ROI, which could lead to incorrect
                        TG when suboptimal coil mode is selected.

-------------------------------------------------------------------------------
*/

@global PSglobal
/*********************************************************************
 *                    PRESCAN.E GLOBAL SECTION                       *
 *                            PSglobal                               *
 *                                                                   *
 * Common code shared between the Host and Tgt PSD processes.  This  *
 * section contains all the #define's, global variables and function *
 * declarations (prototypes).                                        *
 *********************************************************************/
#include "Prescan.h"

#include <stdio.h>
#include <sysDep.h>
#include <sysDepSupport.h>

#define amp_killer 0.4
#define pw_killer 3.6ms
#define FA_FERMI_BLS 630

/* defines for pimrsaps CVs from op_prescan.h */
#define MRSAPS_OFF 0
#define MRSAPS_CFL 1
#define MRSAPS_TG 2
#define MRSAPS_VTG 2
#define MRSAPS_CFH 3
#define MRSAPS_TR 4
#define MRSAPS_FSEPS 9
#define MRSAPS_APA 10
#define MRSAPS_RCVN 12
#define MRSAPS_AWS 101
#define MRSAPS_AVS 102
#define MRSAPS_AS  103
#define MRSAPS_FTG 104
#define MRSAPS_XTG 116

/* defines for cfh_ti */
#define CFHTI_1HT 120000
#define CFHTI_3T  190000
#define CFHTE_1HT 50000
#define CFHTE_3T  30000

/* defines for rcvn_filter */
#define RCVN_MIN_BW  4.0
#define RCVN_MAX_BW  62.5
#define RCVN_MIN_TR  250ms

@host PShostVars
/* Structure definitions for prescan filters*/
FILTER_INFO echo1as_filt;
FILTER_INFO echo1mps1_filt;
FILTER_INFO echo1ftg_filt;
FILTER_INFO echo1xtg_filt;
FILTER_INFO echo1cfl;
FILTER_INFO echo1cfh;
FILTER_INFO echo1rcvn;

/* defines for pimrsaps CVs from op_prescan.h */
#ifndef PSC_MAX_CONTROL_ARRAY
#define PSC_MAX_CONTROL_ARRAY 15
#endif

#define COILLOG_LOG_MAXSIZE 262144 /* quarter-Meg */

int* pimrs[PSC_MAX_CONTROL_ARRAY];

/* YMSmr09211 04/26/2006 YI */
SCAN_INFO cfh_info[MAX_PSC_VQUANT];
SCAN_INFO ps1scan_info[PRESCAN_ROT_MAX];

@cv PScvs
/*********************************************************************
 *                     PRESCAN.E HOST SECTION                        *
 *                             PScvs                                 *
 *                                                                   *
 * Standard C variables of _limited_ types common for both the Host  *
 * and Tgt PSD processes. Declare here all the simple types, e.g,    *
 * int, float, and C structures containing the min and max values,   *
 * and ID description, etc.                                          *
 *                                                                   *
 * NOTE FOR Lx:                                                      *
 * Since the architectures between the Host and the Tgt sides are    *
 * different, the memory alignment for certain types varies. Hence,  *
 * the following types are "forbidden": short, char, and double.     *
 *********************************************************************/
float PSsr_derate_factor = 1.0 with {1.0, 50.0, 1.0, VIS, "PSC SR derating factor",};
float PSassr_derate_factor = 1.0 with {1.0, 50.0, 1.0, VIS, "AutoShim SR derating factor",};

/* ezi (GE) */
float mpsfov = 100 with {FOV_MIN,,,VIS, "mpsfov",} ;

int fastprescan = 0 with {0, 1, 0, VIS, "Fast prescan on/off",};
int pre_slice = 0 with {0,,0,INVIS, "prescan slice number",};
int PSslice_num;
float xmtaddAPS1, xmtaddCFL, xmtaddCFH, xmtaddFTG, xmtadd, xmtaddRCVN;
float ps1scale, cflscale, cfhscale, ftgscale;
float extraScale;  /* for external PSD use */
int PSdebugstate = 0 with {0, 1, 0, VIS, "Debug flag for Prescan",};
int PSfield_strength = B0_5000 with {0,,, VIS, "Prescan Field Strength",};
int PScs_sat = 1 with {0, 1, 1, VIS, "Prescan Chem-SAT flag",};
int PSir = 1 with {0, 1, 1, VIS, "Prescan IR flag",};
int PSmt = 1 with {0, 1, 1, VIS, "Prescan MT flag",};
int ps1_rxcoil = 0 with {0, 1, 0, VIS, "TG PS1 Coil: 1=Rxed coil, 0=default",};

/* vmx 06/14/95 YI */
int tg_1_2_pw = 1 with {0,1,1,VIS,"1:2 pw TG prescan flag",};
int tg_axial = 1 with {0,1,1,VIS,"0:user plane 2:axial",};
float coeff_pw_tg = 1.0;
float fov_lim_mps = 350.0 with {30,450,350,VIS,"FOV limt for MPS",};
/* end vmx */

int TGspf = 0 with {0, 1, 0, VIS, "TG swap phase/freq. flag", };

float flip_rf2cfh;
float flip_rf3cfh; /* For presscfh MRIhc08321 */
int ps1_tr=2s;
int cfl_tr=398ms;
int cfh_tr=398ms;
int rcvn_tr=398ms;

float cfh_ec_position = (16.0/256.0) with {0.0, 1.0, (16.0/256.0), VIS, "Position of the echo center",};
				/* vmx 05/02/95 YO */
int cfl_dda = 4 with {0, 4, 4, VIS, "Num. disdaqs in cfl",};
int cfl_nex = 2 with {1, 2, 2, VIS, "Num. nex in cfl",};
int cfh_dda = 4 with {0, 4, 4, VIS, "Num. disdaqs in cfh",};
int cfh_nex = 2 with {1, 2, 2, VIS, "Num. nex in cfh",};
int rcvn_dda = 0 with {0, 4, 4, VIS, "Num. disdaqs in RCVN",};
int rcvn_nex = 1 with {1, 2, 2, VIS, "Num. nex in RCVN",};

/* presscfh  cvs -  For  MRIhc08321 */ 
/* Definitions in epic.h
 * PRESSCFH_SLICE 1      
 * PRESSCFH_SLAB 2                                              
 * PRESSCFH_SHIMVOL 3
 * PRESSCFH_NONE 4
 * */
int presscfh_override = 0  with {0, PRESSCFH_NONE,0, VIS, "PSD Level Control for Overriding the default CFH selected by opimode",};
int presscfh = PRESSCFH_NONE  with {1,PRESSCFH_NONE,PRESSCFH_NONE, VIS, "PSD Level Control for CFH",};
int presscfh_ctrl = PRESSCFH_NONE  with {1,PRESSCFH_NONE ,PRESSCFH_NONE, VIS, "Geometry Level control for CFH",};
int presscfh_outrange = 0;
int presscfh_cgate = 0;
int presscfh_debug = 0 with {0,DEBUG_DEV,0,VIS,"IR PRESS CFH debugging information",};
int presscfh_wait_rf12 = 0;
int presscfh_minte = 20000;
float presscfh_fov = 0.0;
float presscfh_fov_ratio = 1.0;
float presscfh_pfov_ratio = 1.0;
float presscfh_slab_ratio = 1.0;
float presscfh_pfov = 0.0; 
float presscfh_slthick = 10.0;
float presscfh_ir_slthick = 10.0;
int presscfh_ir_noselect = 1;
/* SXZ: change from 0.5 to 0.3 to make sure presscfh is more probably
 * used for PRESSCFH_SLICE mode */
float presscfh_minfov_ratio = 0.3; /* SXZ: change from 0.5 to 0.3 */ 

/* steam_flag */
int cfh_steam_flag = 0;
int steam_pg_gap = 8;

float area_gykcfl;
float area_gykcfh;
float area_xtgzkiller;
float area_xtgykiller;
int PSoff90=80us;
int dummy_pw;
int min180te;

float PStloc;
float PSrloc;
float PSphasoff;
int PStrigger;

/* begin aps1_mod changes (GE) */
float PStloc_mod;
float PSrloc_mod;
float thickPS_mod;
/* end aps1_mod changes (GE) */

float asx_killer_area = 840.0; /* based on 5.4:  .7 * ( 800 + 400 ) */
float asz_killer_area = 840.0; /* based on 5.4:  .7 * ( 800 + 400 ) */
float cfhir_killer_area = 4086.0; /* based on 5.4:  .9 * ( 4000 + 540 ) */
float ps_crusher_area = 714.0; /* based on 5.4:  .7 * ( 600 + 420 ) */
float cfh_crusher_area = 4000.0; /* MRIhc57311: increased crusher area to kill signal from outside slice */
float target_cfh_crusher;
float target_cfh_crusher2;  /* For presscfh MRIhc08321 */ 

int cfh_newmode = 1;
float cfh_rf2freq = 0 with {,,,VIS,"",};
float cfh_rf3freq = 0 with {,,,VIS,"",}; /* For presscfh MRIhc08321 */
float cfh_rf1freq = 0 with {,,,VIS,"",};
float cfh_fov = 0 with {,,,VIS,"",};
int cfh_ti = CFHTI_1HT;
int eff_cfh_te = CFHTE_1HT;

/**** FastTG CVs   *******/
float FTGslthk = 20 with {,,,VIS,"",};

float FTGopslthickz1=80 with {MINTHICK,80,80,VISONLY,"Slice thickness in mm.",};
float FTGopslthickz2=80 with {MINTHICK,80,80,VISONLY,"Slice thickness in mm.",};
float FTGopslthickz3=20 with {MINTHICK,80,20,VISONLY,"Slice thickness in mm.",};
int   ftgtr = 2s with {TR_MIN,TR_MAX,2s,VISONLY,"Fast TG time of repetition",};
float FTGfov = 480.0 with {FOV_MIN,,,VISONLY,"FastTG fov",};
float FTGau  = 4 with {,,4,VISONLY,"Tau scale factor"};
float FTGtecho = 4 with {,,4,VISONLY,"",};
int FTGtau1   = 8192us with {0,64ms,8192us,VISONLY,"Theta1 to Theta2 time (center to center)"};
int FTGtau2   = 32768us with {0,64ms,32768us,VISONLY, "Theta1 to Theta3 time (center to center)"};
int FTGacq1   = 0 with {0,1,0,VISONLY, "1=window one active, 0=disabled",};
int FTGacq2   = 1 with {0,1,1,VISONLY, "1=window two active, 0=disabled",};
int epi_ir_on = 0 with {0,1,1,VIS, "1=IR EPI",};	/* ypd */
int ssfse_ir_on = 0 with {0,1,1,VIS, "1=IR EPI",};	/* MRIge66767 */
int ftg_dda = 0 with {0,16,0,INVIS, "Num. disdaqs in fast TG",};

float FTGecho1bw = 3.90625 with {2,32,4,VISONLY, "Fast TG Echo1 filter bw. in KHz",};
int FTGtestpulse  = 0 with {0,1,0,VISONLY,"Test pulse for gradient moment tests."
,};
int FTGxres    = 256 with {16,512,256,VISONLY, "X(frequency) resolution",};
float FTGxmtadd;
int pw_gxw2ftgleft = 4096; /* HCSDM00161809: time of 2nd readout window to S1 echo */

/**** eXpressTG CVs   *******/
int   xtgtr = 200ms with {TR_MIN,TR_MAX,200ms,VISONLY,"eXpress TG time of repetition",};
int XTGtau1   = 8192us with {0,64ms,8192us,VISONLY,"Theta1 to Theta2 time (center to center)"};
float XTGfov = 480.0 with {FOV_MIN,,,VISONLY,"eXpress TG fov",};
int pw_bsrf = 4ms;
int xtg_volRecCoil = 0 with {0,1,0,VIS, "XTG coil 0=Rx'ed Coil; 1=vol. rec. coil",};
int xtg_offres_freq = 2000; /* 2kHz off-resonance */
float XTGecho1bw = 15.625 with {2,32,15.625,VISONLY, "eXpress TG Echo1 filter bw. in KHz",};
int XTGxres    = 256 with {16,512,256,VISONLY, "X(frequency) resolution",};
float xmtaddXTG, xtgscale;
int xtg_dda = 0 with {0,16,0,INVIS, "Num. disdaqs in express TG",};
int XTGacq1   = 0 with {0,1,0,VISONLY, "1=window one active, 0=disabled",};
float XTGopslthick = 10 with {,,,VIS,"",};


/**** CFL/CFH CVs   *******/
int CFLxres = 256 with {16,4096,256,VISONLY, "CFL X(frequency) resolution",}; /* MRIhc54366 */
 
float echo1bwcfl = 2.016129 with {,,,INVIS, "Echo1 CFL filter bw. in KHz",};
float echo1bwcfh = 0.50 with {,,,INVIS, "Echo1 CFH filter bw. in KHz",};

int echo1ptcfh = 256 with {,,,INVIS, "Echo1 CFH filter output points",};  /* MRIhc08595 */

float echo1bwrcvn = 15.625 with {,, 15.625, INVIS, "Echo1 RCVN filter BW in KHz",};
int   rcvn_xres   = 4096 with {128, 4096, 4096, VIS, "X(frequency) resolution for RCVN filter",};
int   rcvn_loops  = 10;  /* number of iterations RCVN entry point will loop */

/**** AutoShim CVs *******/

float echo1bwas = 15.625 with {,,,INVIS, "Echo1 auto-shim filter bw. in KHz",};

int off90as  = 80us with {,,80,INVIS, "Comp factor for real 90",};
int td0as  = 4 with {0,,1,INVIS, "Init deadtime",};
int t_exaas = 0 with {0,,0,INVIS,"time from start of 90 to mid 90"};
int time_ssias = 400us with {0,,400us,INVIS, "time from eos to ssi in intern trig",};
int tleadas = 25us with {0,,25us,INVIS, "Init deadtime",};

int te_as;
int tr_as;
int as_dda = 4 with {0,4,4,INVIS, "Num. disdaqs in autoshim",};

/**** Receive Gain CVs ****/

int rgfeature_enable  = PSD_OFF with {PSD_OFF, PSD_ON, PSD_OFF, INVIS, "Enable RG Reduction Feature (0:No, 1:Yes)",};

/********** Required for mrsaps/opt prescan ****************/
/********** pimrsapsflg must be on *************************/
float aslenap = 200 with {0, , 200, VIS, "AP length of shim voxel (mm)", };
float aslenrl = 200 with {0, , 200, VIS, "RL length of shim voxel (mm)", };
float aslensi = 200 with {0, , 200, VIS, "SI length of shim voxel (mm)", };

float aslocap = 0 with { , , 0, VIS, "AP location of shim voxel (mm)", };
float aslocrl = 0 with { , , 0, VIS, "RL location of shim voxel (mm)", };
float aslocsi = 0 with { , , 0, VIS, "SI location of shim voxel (mm)", };
/***********************************************************/

/* temp crusher amplitudes */
float area_gxwas;                  /* readout pulse area */
float area_gz1as;
float area_readrampas;             /* area of left readout ramp */
int avail_pwgx1as;                 /* avail time for gx1as pulse */
int avail_pwgy1as;                 /* avail time for gy1as pulse */
int avail_pwgz1as;               /* avail time for gz1as pulse */
int bw_rf1as;                    /* bandwidth of rf pulses */

/* filter info for 1st, 2nd echo */
float flip_pctas=1.0;                 /* flip angle % for rf scaling */

int dix_timeas;              /* dixon delay for even excitations */
float xmtaddas,xmtlogas;     /* rf attenuation */
int ps1obl_debug = 0 with {0,1,0,INVIS,
                        "On(=1) to print messages for obloptimize",};
int asobl_debug = 0 with {0,1,0,INVIS,
                        "On(=1) to print messages for obloptimize",};
int ps1_newgeo = 1;

int as_newgeo = 1;
int pw_gy1as_tot;
int endview_iampas;
float endview_scaleas;

/* YMSmr09211  04/26/2006 YI */
int cfh_newgeo = 1; 
int cfhobl_debug = 0 with {0,1,0,INVIS,
                        "On(=1) to print messages for obloptimize",};

float deltf = 1.0 with {,,1.0,VIS,"Frequency shift",};

int IRinCFH = 0 with {0, 0, 0, INVIS, "YMS IR in CFH flag", };
int cfh_each = 0 with {0, 0, 0, INVIS, "YMS CFH per slice flag", };
int cfh_slquant = 0 with {0, 0, 0, INVIS, "YMS CFH slice quantity", };

int noswitch_slab_psc = 0 with{ 0,1,0,VIS,"No slab switch for psc(sotf)",};
int noswitch_coil_psc = 0 with{ 0,1,0,VIS,"No coil switch for psc(sotf)",};
int PStest_slab = 1 with{ 1,2,1,INVIS,"Testing slab for psc(sotf)",};
/******************** Communication cv's (Prescan)**********/
int pimrsapsflg = 0 with {0, 1, 0, VIS, "flag for MRS AutoPrescan", };
int pimrsaps1 = MRSAPS_CFL with  {0, 116, 1, VIS, "MRS AutoPrescan step 1: "
                                 "0=Stop, 1=CFL, 2=APS1, 3=CFH, 4=TR,"
                                 "9=PC, 10=APA, 12=RCVN, 101=AWS, 103=SHIM, "
                                 "104=FTG, 116=XTG", };
int pimrsaps2 = MRSAPS_FTG with  {0, 116, 2, VIS, "MRS AutoPrescan step 2: "
                                 "0=Stop, 1=CFL, 2=APS1, 3=CFH, 4=TR, "
                                 "9=PC, 10=APA, 12=RCVN, 101=AWS, 103=SHIM, "
                                 "104=FTG, 116=XTG", };
int pimrsaps3 = MRSAPS_AS with   {0, 116, 3, VIS, "MRS AutoPrescan step 3: "
                                 "0=Stop, 1=CFL, 2=APS1, 3=CFH, 4=TR, "
                                 "9=PC, 10=APA, 12=RCVN, 101=AWS, 103=SHIM, "
                                 "104=FTG, 116=XTG", };
int pimrsaps4 = MRSAPS_TR with   {0, 116, 4, VIS, "MRS AutoPrescan step 4: "
                                 "0=Stop, 1=CFL, 2=APS1, 3=CFH, 4=TR, "
                                 "9=PC, 10=APA, 12=RCVN, 101=AWS, 103=SHIM, "
                                 "104=FTG, 116=XTG", };
int pimrsaps5 = MRSAPS_RCVN with  {0, 116, 1, VIS, "MRS AutoPrescan step 5: "
                                 "0=Stop, 1=CFL, 2=APS1, 3=CFH, 4=TR, "
                                 "9=PC, 10=APA, 12=RCVN, 101=AWS, 103=SHIM, "
                                 "104=FTG, 116=XTG", };
int pimrsaps6 = MRSAPS_CFH with {0, 116, 3, VIS, "MRS AutoPrescan step 6: "
                                 "0=Stop, 1=CFL, 2=APS1, 3=CFH, 4=TR, "
                                 "9=PC, 10=APA, 12=RCVN, 101=AWS, 103=SHIM, "
                                 "104=FTG, 116=XTG", };
int pimrsaps7 = MRSAPS_OFF with  {0, 116, 3, VIS, "MRS AutoPrescan step 7: "
                                 "0=Stop, 1=CFL, 2=APS1, 3=CFH, 4=TR, "
                                 "9=PC, 10=APA, 12=RCVN, 101=AWS, 103=SHIM, "
                                 "104=FTG, 116=XTG", };
int pimrsaps8 = MRSAPS_OFF with  {0, 116, 101, VIS, "MRS AutoPrescan step 8: "
                                 "0=Stop, 1=CFL, 2=APS1, 3=CFH, 4=TR, "
                                 "9=PC, 10=APA, 12=RCVN, 101=AWS, 103=SHIM, "
                                 "104=FTG, 116=XTG", };
int pimrsaps9 = MRSAPS_OFF with  {0, 116, 0, VIS, "MRS AutoPrescan step 9: "
                                 "0=Stop, 1=CFL, 2=APS1, 3=CFH, 4=TR, "
                                 "9=PC, 10=APA, 12=RCVN, 101=AWS, 103=SHIM, "
                                 "104=FTG, 116=XTG", };
int pimrsaps10 = MRSAPS_OFF with  {0, 116, 0, VIS, "MRS AutoPrescan step 10: "
                                  "0=Stop, 1=CFL, 2=APS1, 3=CFH, 4=TR, "
                                  "9=PC, 10=APA, 12=RCVN, 101=AWS, 103=SHIM, "
                                  "104=FTG, 116=XTG", };
int pimrsaps11 = MRSAPS_OFF with  {0, 116, 0, VIS, "MRS AutoPrescan step 11: "
                                  "0=Stop, 1=CFL, 2=APS1, 3=CFH, 4=TR, "
                                  "9=PC, 10=APA, 12=RCVN, 101=AWS, 103=SHIM, "
                                  "104=FTG, 116=XTG", };
int pimrsaps12 = MRSAPS_OFF with  {0, 116, 0, VIS, "MRS AutoPrescan step 12: "
                                  "0=Stop, 1=CFL, 2=APS1, 3=CFH, 4=TR, "
                                  "9=PC, 10=APA, 12=RCVN, 101=AWS, 103=SHIM, "
                                  "104=FTG, 116=XTG", };
int pimrsaps13 = MRSAPS_OFF with  {0, 116, 0, VIS, "MRS AutoPrescan step 13: "
                                  "0=Stop, 1=CFL, 2=APS1, 3=CFH, 4=TR, "
                                  "9=PC, 10=APA, 12=RCVN, 101=AWS, 103=SHIM, "
                                  "104=FTG, 116=XTG", };
int pimrsaps14 = MRSAPS_OFF with  {0, 116, 0, VIS, "MRS AutoPrescan step 14: "
                                  "0=Stop, 1=CFL, 2=APS1, 3=CFH, 4=TR, "
                                  "9=PC, 10=APA, 12=RCVN, 101=AWS, 103=SHIM, "
                                  "104=FTG, 116=XTG", };
int pimrsaps15 = MRSAPS_OFF with  {0, 116, 0, VIS, "MRS AutoPrescan step 15: "
                                  "0=Stop, 1=CFL, 2=APS1, 3=CFH, 4=TR, "
                                  "9=PC, 10=APA, 12=RCVN, 101=AWS, 103=SHIM, "
                                  "104=FTG, 116=XTG", };

/* MRIhc15304: CVs related to coil switching
 * By changing the hub index for the coil and sending
 * it through SSP packet (RFHUB ssp packet)
 */
int pw_contrfhubsel = 4 with {0, , 4, INVIS, "Width of the change hub index packet",};
int delay_rfhubsel = 20;
int pw_contrfsel = 4 with {0, , 4, INVIS, "Width of the modify receiver port packet",};
int csw_tr = 0 with {0, ,0, VIS,"seq length for receiver coil switch core",};
int csw_wait_sethubindeximm = 250ms with
    {0, , 250ms, VIS, "Additional time for coil switch when calling sethubindeximm",};
int csw_wait_setrcvportimm = 100ms with
    {0, , 100ms, VIS, "Additional time for coil switch when calling setrcvportimm ",};
int csw_wait_before = 10ms with {0, , 10ms,INVIS,"Delay for coil switching startup",};
int csw_time_ssi = 50us with
    {0, , 50ms, VIS, "time from eos to ssi in intern trig for coil switch",};

/* MRIhc47602/MRIhc47515/GEHmr03545 : Killer gradient in the RCVN sequence */
float area_gxkrcvn = 10000;
float area_gykrcvn = 10000;
float area_gzkrcvn = 10000;
int pre_rcvn_tr = 20ms with {0, , 0, VIS, "Pre sequence before RCVN",};
int rcvn_flag = 2 with {0, 2, 1, VIS, "0: OFF; 1: crusher; 2: delay b4 RCVN",};

/***********************************************************/

@host PSipgexport
/*********************************************************************
 *                  PRESCAN.E IPGEXPORT SECTION                      *
 *                                                                   *
 * Standard C variables of _any_ type common for both the Host and   *
 * Tgt PSD processes. Declare here all the complex type, e.g.,       *
 * structures, arrays, files, etc.                                   *
 *                                                                   *
 * NOTE FOR Lx:                                                      *
 * Since the architectures between the Host and the Tgt sides are    *
 * different, the memory alignment for certain types varies. Hence,  *
 * the following types are "forbidden": short, char, and double.     *
 *********************************************************************/
int PSfreq_offset[ENTRY_POINT_MAX];
int cfl_tdaq;
int cfh_tdaq;
int rcvn_tdaq;
long rsp_PSrot[MAX_PSC_VQUANT] [9];

/* For presscfh MRIhc08321 */
PSC_INFO presscfh_info[MAX_PSC_VQUANT]={ {0,0,0,{0},0,0,0} };

/* YMSmr09211  04/26/2006 YI */
LOG_GRAD  cflloggrd = {0,0,0,0,0,0,{0,0,0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 };
LOG_GRAD  ps1loggrd = {0,0,0,0,0,0,{0,0,0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 };
LOG_GRAD  cfhloggrd = {0,0,0,0,0,0,{0,0,0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 };
LOG_GRAD  rcvnloggrd = {0,0,0,0,0,0,{0,0,0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 };


WF_PROCESSOR read_axis = XGRAD;
WF_PROCESSOR killer_axis = YGRAD;

@host PScvinit
    { /* Start of code inlined from Prescan.e PScvinit */
        /*********************************************************************
         *                     PRESCAN.E HOST SECTION                        *
         *                           PScvinit                                *
         *                                                                   *
         * Write here the code unique to the Host PSD process. All code to   *
         * be executed in cvinit() must be written here.                     *
         *********************************************************************/
        
        cvdef(mpsfov, cfsystemmaxfov);
        cvdef(FTGfov, cfsystemmaxfov);
        
        FTGfov = cfsystemmaxfov;

        PScvinit();
        FTGcvinit();
        XTGcvinit();
        AScvinit();
        RGcvinit();
     } /* End of code inlined from Prescan.e PScvinit */

@host PScveval
    { /* Start of code inlined from Prescan.e PScveval */
        /*********************************************************************
         *                     PRESCAN.E HOST SECTION                        *
         *                           PScveval                                *
         *                                                                   *
         * Write here the code unique to the Host PSD process. All code to   *
         * be executed in cveval() must be written here.                     *
         *********************************************************************/

        /* MRIhc49539: check current nucleus against coil DB: */
        if(specnuc != coilInfo[0].rxNucleus)
        {
            epic_error(use_ermes, "%s is incompatible with %s.", EM_PSD_INCOMPATIBLE, EE_ARGS(2), STRING_ARG, "This PSD", STRING_ARG, "the selected coil");
            return FAILURE; 
        }
         
        TGspf = ( (0 == getAps1Mod()) && opspf );
    	read_axis = TGspf ? YGRAD : XGRAD;
        if(read_axis == XGRAD)
        {
            killer_axis = YGRAD;
        }
        else
        {
            killer_axis = XGRAD;
        }

        if (psddebugcode)
        {
            psd_dump_coil_info();
        }

        if (FAILURE==PScveval())
        {
            return FAILURE;
        }
        if (FAILURE==FTGcveval())
        {
            return FAILURE;
        }
        if (FAILURE==XTGcveval())
        {
            return FAILURE;
        }
        if (FAILURE==AScveval())
        {
            return FAILURE;
        }
        if (FAILURE==RGcveval())
        {
            return FAILURE;
        }

    } /* End of code inlined from Prescan.e PScveval */


@host PSpredownload
{ /* Start of code inlined from Prescan.e PSpredownload */
    /*********************************************************************
     *                     PRESCAN.E HOST SECTION                        *
     *                         PSpredownload                             *
     *                                                                   *
     * Write here the code unique to the Host PSD process. All code to   *
     * be executed in predownload() must be written here.                *
     *********************************************************************/
    if (FAILURE==PSpredownload())
    {
        return FAILURE;
    }
    if (FAILURE==FTGpredownload())
    {
        return FAILURE;
    }
    if (FAILURE==XTGpredownload())
    {
        return FAILURE;
    }
    if (FAILURE==ASpredownload())
    {
        return FAILURE;
    }
} /* End of code inlined from Prescan.e PSpredownload */

@host PSfilter
{ /* Start of code inlined from Prescan.e PSfilter */
    /*********************************************************************
     *                     PRESCAN.E HOST SECTION                        *
     *                            PSfilter                               *
     *                                                                   *
     * Write here the code unique to the Host PSD process.               *
     *********************************************************************/
    /* PS filter must be called with no arguments. num_filter_slot 
       is a CV which comes from individual PSDs, which get incremented 
       in the PSD so that Prescan filter generation takes place for 
       the next slot in psd_filter_spec. - RJF 13/Oct/1998 */
    /* vmx 10/13/94 YI */
    PSfilter();
    /* end vmx */

    if(psddebugcode)
    {
        dump_runtime_filter_info(psd_filt_spec);
    }
} /* End of code inlined from Prescan.e PSfilter */

@host PShost
/*********************************************************************
 *                     PRESCAN.E HOST SECTION                        *
 *                             PShost                                *
 *                                                                   *
 * Write here the code unique to the Host PSD process.               *
 *********************************************************************/
#include <string.h>
#include "sar_pm.h" 
#include "epic_usercv.h"

/*
 *  set_presscfh_mode
 *  
 *  Type: Private Function
 *  
 *  Description: Set the modes for presscfh
 * 			PRESSCFH_SLICE 1      
 * 			PRESSCFH_SLAB 2                                              
 * 			PRESSCFH_SHIMVOL 3
 * 			PRESSCFH_NONE 4
 * 
 * */

STATUS
set_presscfh_mode (void)
{
    if(exist(oppscvquant) >= 1)
    {
        if (presscfh_override > 0) 
        {
            cvoverride(presscfh, presscfh_override, PSD_FIX_ON, PSD_EXIST_ON);
        }
        else
        {
            switch(exist(opimode)) 
            {
                case PSD_2D:
                case PSD_CINE: 
                    if( (opplane == PSD_3PLANE) || ((exist(oprealtime) == PSD_ON)) )
                    {   /* For 3-Plane or Realtime mode use SHIMVOL Mode*/
                        cvoverride(presscfh, PRESSCFH_SHIMVOL, PSD_FIX_ON, PSD_EXIST_ON);
                    }
                    else 
                    {
                        if( exist(opcoax) != 0 ) 
                        {   /* NON-RADIAL MODE */
                            cvoverride(presscfh, PRESSCFH_SLICE, PSD_FIX_ON, PSD_EXIST_ON);
                            if(exist(opassetcal) && existcv(opassetcal) 
                               && ( (1==exist(opasset)) || (exist(opasset) == ASSET_REG_CAL) ) && existcv(opasset)) 
                            {
                                cvoverride(presscfh, PRESSCFH_SHIMVOL, PSD_FIX_ON, PSD_EXIST_ON);
                            }
                        } 
                        else if(exist(opcoax) == 0) 
                        {   /* RADIAL MODE */
                            cvoverride(presscfh, PRESSCFH_SHIMVOL, PSD_FIX_ON, PSD_EXIST_ON);
                        } 
                        else 
                        {
                            cvoverride(presscfh, PRESSCFH_NONE, PSD_FIX_ON, PSD_EXIST_ON);	
                        }
                    }
                    break;
                case PSD_3D:
                case PSD_3DM:
                case PSD_ANGIO:
                    /* If Shim Volume is placed */
                    if( exist(opcoax) != 0 )      /* Non-Radial/Oblique MODE */
                    {
                        if( exist(opvquant) > 1 )     /* Multiple Slabs */
                        {
                            cvoverride(presscfh, PRESSCFH_SHIMVOL, PSD_FIX_ON, PSD_EXIST_ON);	
                        }
                        else     /* Single Slab */
                        {
                            cvoverride(presscfh, PRESSCFH_SLAB, PSD_FIX_ON, PSD_EXIST_ON);
                        }
                    }
                    else     /* Radial Case or Orthogonal Slabs */
                    {
                        cvoverride(presscfh, PRESSCFH_SHIMVOL, PSD_FIX_ON, PSD_EXIST_ON);
                    }
                    break;
                default: 
                    cvoverride(presscfh, PRESSCFH_NONE, PSD_FIX_ON, PSD_EXIST_ON);
                    break;
            }
        } 
    }
    else
    {
        cvoverride(presscfh, PRESSCFH_NONE, PSD_FIX_ON, PSD_EXIST_ON);
    }
    
    if(presscfh_debug) 
    {
        printf("\n The presscfh is %d (1 - SLICE, 2 - SLAB, 3 - SHIMVOL 4- NONE)\n",presscfh);
        printf("\n The presscfh_ctrl is %d (1 - SLICE, 2 - SLAB, 3 - SHIMVOL 4- NONE), \n",presscfh_ctrl);
        printf("\n The presscfh_override is %d (1 - SLICE, 2 - SLAB, 3 - SHIMVOL 4- NONE)\n",presscfh_override);
        fflush(stdout);
    }
    return SUCCESS;
}

/*
 *  sr_derate
 *  
 *  Type: Private Function
 *  
 *  Description:
 *    update the ramp time for loggrad by *sr_derate_factor  
 *    @param[out] lgrad, logical gradient characteristics
 *    @param[in]  sc_derate_factor, SR derating factor ( >=1.0 ) 
 *
 */

STATUS
sr_derate (LOG_GRAD *lgrad, const float sr_derate_factor)
{
    STATUS status = SUCCESS;

    if(sr_derate_factor < 1.0)
    {
        status = FAILURE;
    }
    else
    {
        lgrad->xrt = lgrad->xrt*sr_derate_factor; 
        lgrad->yrt = lgrad->yrt*sr_derate_factor; 
        lgrad->zrt = lgrad->zrt*sr_derate_factor; 
        lgrad->xft = lgrad->xft*sr_derate_factor;
        lgrad->yft = lgrad->yft*sr_derate_factor;
        lgrad->zft = lgrad->zft*sr_derate_factor;

        status = SUCCESS;
    }

    return status;
}

/*
 *  PS1cvinit
 *  
 *  Type: Private Function
 *  
 *  Description:
 *  
 */
STATUS
PS1cvinit( void )
{
    INT index;

    /* initialize pulse widths */
    pw_rf1mps1 = rfpulse[RF1_APS1_SLOT].nom_pw;
    pw_rf2mps1 = rfpulse[RF2_APS1_SLOT].nom_pw;

    /* initialize resolutions */
    res_rf1mps1 = 0;
    res_rf2mps1 = 0;

    /* initialize amplitudes */
    a_rf1mps1 = 0.5;
    a_rf2mps1 = 1.0;

    /* initialize flip angles */
    flip_rf1mps1 = 90;
    flip_rf2mps1 = 180;

    /* initialize sinc cycles */
    cyc_rf1mps1 = 1;
    cyc_rf2mps1 = 1;

    /* initialize gscale values */
    gscale_rf1mps1 = 0.90909;
    gscale_rf2mps1 = 0.4545;

    /* begin aps1_mod changes (GE) */
    if ( 1 == getAps1Mod() )
    {
        for( index = 0; index < 9; index++)
        {
            ps1scan_info[0].oprot[index] = 0.0;
        }

        if (1==getAps1ModPlane()) /* Axial, read=x */
        {
            ps1scan_info[0].oprot[0] = ps1scan_info[0].oprot[4] = ps1scan_info[0].oprot[8] = 1.0;
        }
        if (2==getAps1ModPlane()) /* Sagittal, read=z */
        {
            ps1scan_info[0].oprot[2] = ps1scan_info[0].oprot[4] = ps1scan_info[0].oprot[6] = 1.0;
        }
        if (3==getAps1ModPlane()) /* Coronal, read=z */
        {
            ps1scan_info[0].oprot[1] = ps1scan_info[0].oprot[5] = ps1scan_info[0].oprot[6] = 1.0;
        }
        if (4==getAps1ModPlane()) /* Axial, read=y */
        {
            ps1scan_info[0].oprot[1] = ps1scan_info[0].oprot[3] = ps1scan_info[0].oprot[8] = 1.0;
        }
        if (5==getAps1ModPlane()) /* Sagittal, read=y */
        {
            ps1scan_info[0].oprot[2] = ps1scan_info[0].oprot[3] = ps1scan_info[0].oprot[7] = 1.0;
        }
        if (6==getAps1ModPlane()) /* Coronal, read=x */
        {
            ps1scan_info[0].oprot[0] = ps1scan_info[0].oprot[5] = ps1scan_info[0].oprot[7] = 1.0;
        }
    }
    else
    {
        for (index = 0; index < 9; index++)
        {
            ps1scan_info[0].oprot[index] = scan_info[PSslice_num].oprot[index]; 
        }
    }
    /* end aps1_mod changes (GE) */

    ps1_newgeo = 1;
    if (obloptimize(&ps1loggrd, &phygrd, ps1scan_info, 1, PSD_OBL,
                    0, obl_method, ps1obl_debug, &ps1_newgeo,
                    cfsrmode) == FAILURE)
    {
        epic_error(use_ermes, "%s failed in PS1cvinit.", EM_PSD_SUPPORT_FAILURE,
                   EE_ARGS(1), STRING_ARG, "obloptimize"); 
        return FAILURE;
    }

    /* derate SR for quiet PSC */ 
    sr_derate(&ps1loggrd, PSsr_derate_factor);

    return SUCCESS;
}

/*
 *  CFLcvinit
 *  
 *  Type: Private Function
 *  
 *  Description:
 *  
 */
STATUS
CFLcvinit( void )
{
    /* initialize pulse widths */
    pw_rf1cfl = rfpulse[RF1_CFL_SLOT].nom_pw;

    /* initialize resolutions */
    res_rf1cfl =  0;

    /* initialize amplitudes */
    a_rf1cfl = 0.5;

    /* initialize flip angles */
    flip_rf1cfl = 90;

    /* initialize sinc cycles */
    cyc_rf1cfl = 2; /* MRIhc54366: increased from 1 */

    /* initialize gscale values */
    gscale_rf1cfl = 0.90909;

    return SUCCESS;
}


/*
 *  RCVNcvinit
 *  
 *  Type: Private Function
 *  
 *  Description:
 *  
 */
STATUS
RCVNcvinit( void )
{
    if (CFG_VAL_RECEIVER_RRF == cfreceivertype)
    {
        rcvn_xres = 1024;
    }
    else
    {
        rcvn_xres = 4096;
    }
    return SUCCESS;
}


/*
 *  CFHcvinit
 *  
 *  Type: Private Function
 *  
 *  Description:
 *  
 */
STATUS
CFHcvinit( void )
{
    /* initialize pulse widths */
    pw_rf0cfh = rfpulse[RF0_CFH_SLOT].nom_pw;
    pw_rf1cfh = rfpulse[RF1_CFH_SLOT].nom_pw;
    pw_rf2cfh = rfpulse[RF2_CFH_SLOT].nom_pw;
    pw_rf3cfh = rfpulse[RF3_CFH_SLOT].nom_pw; /* For presscfh MRIhc08321 */

    /* initialize resolutions */
    res_rf0cfh =  RES_SH_ADIABATIC; /* Adiabatic pulse */
    res_rf1cfh =  0;
    res_rf2cfh =  0;
    res_rf3cfh =  0; /* For presscfh MRIhc08321 */

    /* initialize amplitudes */
    a_rf0cfh = 0.61;
    a_rf1cfh = 0.5;
    a_rf2cfh = 1.0;
    a_rf3cfh = 1.0; /* For presscfh MRIhc08321 */

    /* initialize flip angles */
    flip_rf0cfh = 180;
    flip_rf1cfh = 90;
    flip_rf2cfh = 180;
    flip_rf3cfh = 180; /* For presscfh MRIhc08321 */

    /* initialize sinc cycles */
    cyc_rf0cfh = 2; /* Adiabatic pulse */
    cyc_rf1cfh = 1;
    cyc_rf2cfh = 1;
    cyc_rf3cfh = 1; /* For presscfh MRIhc08321 */
    a_gyrf3cfh = 0.0;/*For presscfh MRIhc11621 */
    a_gxrf2cfh = 0.0;/*For  presscfh MRIhc11621 */
    a_gzrf1cfh = 0.0;/*For presscfh MRIhc11621 */
    a_gzrf0cfh = 0.0;/*For presscfh MRIhc11621 */


    /* initialize gscale values */
    gscale_rf1cfh = 0.90909;

#ifdef PSD_CFH_CHEMSAT
    rfpulse[RFCSSAT_CFH_SLOT].pw = &pw_rfcssatcfh;
    rfpulse[RFCSSAT_CFH_SLOT].amp = &a_rfcssatcfh;
    rfpulse[RFCSSAT_CFH_SLOT].act_fa = &flip_rfcssatcfh;
#endif
%ifdef PSD_CFH_MT
    rfpulse[RFMT_CFH_SLOT].pw = &pw_rfmtcfh;
    rfpulse[RFMT_CFH_SLOT].amp = &a_rfmtcfh;
    rfpulse[RFMT_CFH_SLOT].act_fa = &flip_rfmtcfh;
%endif /* PSD_CFH_MT */

    if( presscfh_ctrl != PRESSCFH_NONE && cfh_steam_flag == PSD_ON )
    {
        setuprfpulse(RF1_CFH_SLOT, &pw_rf1cfh, &a_rf1cfh, SAR_ABS_SINC3, SAR_PSINC3,
                     SAR_ASINC3, SAR_DTYCYC_SINC3, SAR_MAXPW_SINC3, 1,
                     MAX_B1_SINC3_90, MAX_INT_B1_SQ_SINC3_90,
                     MAX_RMS_B1_SINC3_90, 90.0, &flip_rf1cfh, 3200.0,
                     3750, PSD_CFH_ON, 0,
                     0, 0, &res_rf1cfh, 0, &wg_rf1cfh, rfpulse);

        setuprfpulse(RF2_CFH_SLOT, &pw_rf2cfh, &a_rf2cfh, SAR_ABS_SINC3, SAR_PSINC3,
                     SAR_ASINC3, SAR_DTYCYC_SINC3, SAR_MAXPW_SINC3, 1,
                     MAX_B1_SINC3_90, MAX_INT_B1_SQ_SINC3_90,
                     MAX_RMS_B1_SINC3_90, 90.0, &flip_rf2cfh, 3200.0,
                     3750, PSD_CFH_ON, 0,
                     0, 0, &res_rf2cfh, 0, &wg_rf2cfh, rfpulse);

        setuprfpulse(RF3_CFH_SLOT, &pw_rf3cfh, &a_rf3cfh, SAR_ABS_SINC3, SAR_PSINC3,
                     SAR_ASINC3, SAR_DTYCYC_SINC3, SAR_MAXPW_SINC3, 1,
                     MAX_B1_SINC3_90, MAX_INT_B1_SQ_SINC3_90,
                     MAX_RMS_B1_SINC3_90, 90.0, &flip_rf3cfh, 3200.0,
                     3750, PSD_CFH_ON, 0,
                     0, 0, &res_rf3cfh, 0, &wg_rf3cfh, rfpulse);

        a_rf0cfh = 1;
        a_rf1cfh = 0.5464; /* 0.5/0.61*(60/90) */
        a_rf2cfh = 0.5464;
        a_rf3cfh = 0.5464; 

        flip_rf1cfh = 60;
        flip_rf2cfh = 60;
        flip_rf3cfh = 60; 

        cyc_rf1cfh = 3;
        cyc_rf2cfh = 3;
        cyc_rf3cfh = 3; 
    }
    else
    {
        setuprfpulse(RF1_CFH_SLOT, &pw_rf1cfh, &a_rf1cfh, SAR_ABS_SINC1, SAR_PSINC1,
                     SAR_ASINC1, SAR_DTYCYC_SINC1, SAR_MAXPW_SINC1, 1,
                     MAX_B1_SINC1_90, MAX_INT_B1_SQ_SINC1_90,
                     MAX_RMS_B1_SINC1_90, 90.0, &flip_rf1cfh, 3200.0,
                     1250, PSD_CFH_ON, 0,
                     0, 0, &res_rf1cfh, 0, &wg_rf1cfh, rfpulse);

        setuprfpulse(RF2_CFH_SLOT, &pw_rf2cfh, &a_rf2cfh, SAR_ABS_SINC1, SAR_PSINC1,
                     SAR_ASINC1, SAR_DTYCYC_SINC1, SAR_MAXPW_SINC1, 1,
                     MAX_B1_SINC1_90, MAX_INT_B1_SQ_SINC1_90,
                     MAX_RMS_B1_SINC1_90, 90.0, &flip_rf2cfh, 3200.0,
                     1250, PSD_CFH_ON, 0,
                     0, 0, &res_rf2cfh, 0, &wg_rf2cfh, rfpulse);

        setuprfpulse(RF3_CFH_SLOT, &pw_rf3cfh, &a_rf3cfh, SAR_ABS_SINC1, SAR_PSINC1,
                     SAR_ASINC1, SAR_DTYCYC_SINC1, SAR_MAXPW_SINC1, 1,
                     MAX_B1_SINC1_90, MAX_INT_B1_SQ_SINC1_90,
                     MAX_RMS_B1_SINC1_90, 90.0, &flip_rf3cfh, 3200.0,
                     1250, PSD_CFH_ON, 0,
                     0, 0, &res_rf3cfh, 0, &wg_rf3cfh, rfpulse);

        a_rf0cfh = 0.61;
        a_rf1cfh = 0.5;
        a_rf2cfh = 1.0;
        a_rf3cfh = 1.0; 

        flip_rf1cfh = 90;
        flip_rf2cfh = 180;
        flip_rf3cfh = 180; 

        cyc_rf1cfh = 1;
        cyc_rf2cfh = 1;
        cyc_rf3cfh = 1; 
    }

    return SUCCESS;
}



/*
 *  PScvinit
 *  
 *  Type: Public Function
 *  
 *  Description:
 *  
 */
STATUS
PScvinit( void )
{

    if( (exist(opimode) == PSD_SPECTRO) || (PSD_ON == pimrsapsflg) )
    {
        fastprescan = 0;
    }
    else
    {
        fastprescan = cffastprescan;
    }

    /* Set the modes for presscfh */
    set_presscfh_mode();
    
    /* Fix to set CFH,CFL,APS1, FTG slice location to be axial mid slice to ensure
       the presence of signal when sample position is off iso center such as wrist scout scans */
    if (opplane == PSD_3PLANE)
    {
        PSslice_num = (int)(opaxial_slice/2);
    }
    else
    {  /* Start of code moved up from PSpredownload section below by SL */
        FLOAT minloc;
        FLOAT maxloc;         /* For MRIge45995 */
        FLOAT midloc;         /* For MRIge45995 */
        FLOAT minimum;        /* For MRIge45995 */
        FLOAT loc;
        INT index;

        /* search through the scan table to find the location nearest
           isocenter. Save the location information for prescan since
           graphic sat and different imaging techniques can alter the
           content and format of the rspinfo structure
        */

        minloc = MAXFLOAT;  /* For MRIge45995 */
        maxloc = -MAXFLOAT; /* For MRIge45995 */

        /* First to find the center of the prescribed region, MRIge45995 */
        for (index=0;index<opslquant*opvquant;index++)
        {
            loc = scan_info[index].optloc;
            if (loc < minloc)
            {
                minloc = loc;
            }
            if (loc > maxloc)
            {
                maxloc = loc;
            }
        }
        midloc = 0.5*(minloc + maxloc);

        /* Find the slice most close to the center of the prescribed region, MRIge45995 */
        minimum = MAXFLOAT;
        for (index=0;index<opslquant*opvquant;index++)
        {
            loc = scan_info[index].optloc;
            if (fabs(loc-midloc) < minimum)
            {
                minimum = fabs(loc-midloc);
                PSslice_num = index;
            }
        }
    } /* End of code moved up from PSpredownload section below */

    /* Comment from /vobs/scan/MrScan/SrxControl/SrxGeos.m on setting opcoax */
    /* slightly modified */

    /* This next section determines the value of the opcoax CV  */
    /* The rules are as follows:  */
    /*   */
    /* coaxial      N   Y   Y  */
    /* offcenter    -   N   Y  */
    /* ======================  */
    /* opcoax       0   1   2  */
    /*  */
    /* The decision of NOT coaxial is determined by the following three  */
    /* conditions: (graphic rx imaging option selected) AND (oblique   */
    /* prescription) AND ( (more than one group has been prescribed) OR (Number
     * of radial slices > 1 ) ).  */
    /*  */
    /* The decision of offcenter is made based on whether or not   */
    /* all of the slices have 0 offset in their phase and frequency  */
    /* directions.  If even one is not 0, then the prescription is  */
    /* said to be offcenter.  */

    /* Use the below in the PSD Code if the HDMR2 changes for MRIhc08321
     * need to be overridden. RBA for MRIhc08321.
     */
    
    
    PS1cvinit();
    CFLcvinit();
    CFHcvinit();
    RCVNcvinit();

    /* initialize field strength and PS variables */ /* vmx 05/02/95 YO */
    PSfield_strength = (int) cffield;
    if(PSfield_strength <= B0_5000)
    {
        PScs_sat = 1;
        PSmt = 0;
    }
    else
    {
        if (aspir_flag || (PSD_ON == exist(opspecir)))
        {   /* turn OFF ASPIR/SPECIAL during CFH */
            PScs_sat = 0;
        } else {
            PScs_sat = 1;
        }

%ifdef PSD_SPSPFATSAT
        if (use_spsp_fatsat)
        {   /* turn OFF SPSP FATSAT during CFH */
            PScs_sat = 0;
        } 
%endif /* PSD_SPSPFATSAT */
        PSmt = 1;
    }

    pimrs[0] = &pimrsaps1;
    pimrs[1] = &pimrsaps2;
    pimrs[2] = &pimrsaps3;
    pimrs[3] = &pimrsaps4;
    pimrs[4] = &pimrsaps5;
    pimrs[5] = &pimrsaps6;
    pimrs[6] = &pimrsaps7;
    pimrs[7] = &pimrsaps8;
    pimrs[8] = &pimrsaps9;
    pimrs[9] = &pimrsaps10;
    pimrs[10] = &pimrsaps11;
    pimrs[11] = &pimrsaps12;
    pimrs[12] = &pimrsaps13;
    pimrs[13] = &pimrsaps14;
    pimrs[14] = &pimrsaps15;

    
    /* MRIhc15304: we will keep asfov as cv and fill it with the value
     * from coil. Some application (Spectro related) need to decide the
     * asfov based on the application (overriding the value decided by
     * coil). */ 
    asfov = coilInfo[0].autoshimFov;

    /* 
     * Set the wait time for sethubindeximm.  A 100 ms delay is
     * sufficent to apply the settings on the driver module.  On systems
     * with an RRF receiver, a longer delay (250 ms) is needed for the
     * new receiver channel map to be loaded into the DRF. 
     */
    if( cfcoilswitchmethod & COIL_SWITCH_RSP_SETHUBINDEXIMM )
    {
        if( CFG_VAL_RECEIVER_RRF == cfreceivertype )
        {
            csw_wait_sethubindeximm  = 250ms;
        }
        else
        {
            csw_wait_sethubindeximm  = 100ms;
        }
    }

    return SUCCESS;

}   /* end PScvinit() */


/*
 *  FTGcvinit
 *  
 *  Type: Public Function
 *  
 *  Description:
 *  
 */
STATUS
FTGcvinit( void )
{
    a_rf1ftg    = 0.5;
    a_rf2ftg    = 1.0;
    a_rf3ftg    = 1.0;
    pw_rf1ftg   = (int)rfpulse[RF1_FTG_SLOT].nom_pw;
    pw_rf2ftg   = (int)rfpulse[RF2_FTG_SLOT].nom_pw;
    pw_rf3ftg   = (int)rfpulse[RF3_FTG_SLOT].nom_pw;
    cyc_rf1ftg  = 1;
    res_rf1ftg  = 800;
    cyc_rf2ftg  = 1;
    res_rf2ftg  = 800;
    cyc_rf3ftg  = 1;
    res_rf3ftg  = 800;
    flip_rf1ftg =  90.0;
    flip_rf2ftg = 180.0;
    flip_rf3ftg = 180.0;

    return SUCCESS;
} /* end FTGcvinit() */


/*
 *  XTGcvinit
 *  
 *  Type: Public Function
 *  
 *  Description:
 *  
 */
STATUS
XTGcvinit( void )
{    
    double coil_maxb1, fermi_b1; /* HCSDM00157626 */
    a_rf1xtg    = 0.5;
    a_rf2xtg    = 1.0;
    pw_rf1xtg   = (int)rfpulse[RF1_XTG_SLOT].nom_pw;
    pw_rf2xtg   = (int)rfpulse[RF2_XTG_SLOT].nom_pw;
    cyc_rf1xtg  = 1;
    res_rf1xtg  = 800;
    cyc_rf2xtg  = 1;
    res_rf2xtg  = 800;
    flip_rf1xtg =  90.0;
    flip_rf2xtg = 180.0;

    a_rf3xtg = -1.0; 
    a_rf4xtg = 1.0;
    res_rf3xtg = RES_FERMI_BLS;
    res_rf4xtg = res_rf3xtg;
    pw_rf3xtg   = pw_bsrf;
    pw_rf4xtg   = pw_rf3xtg;
    
    flip_rf3xtg = (FA_FERMI_BLS*pw_rf3xtg)/SAR_FERMI_BLS_NOM_PW;  /* scale flip to keep its maxB1 constant */
    /* 630 degree flip, 4ms pw: max B1 = 0.071 */

    /* HCSDM00157626: scale rf3xtg and rf4xtg flip angle if coil B1 limit
       would otherwise force pulses to stretch */
    coil_maxb1 = txCoilInfo[getTxIndex(coilInfo[0])].maxB1Peak;
    fermi_b1 = 100*FA_FERMI_BLS/NOM_FA_RFMT*SAR_MAXB1_FERMI_BLS;

    if(fermi_b1 > coil_maxb1)
    {
        flip_rf3xtg = (float)((int)(flip_rf3xtg*coil_maxb1/fermi_b1));
    }

    flip_rf4xtg = flip_rf3xtg;       
    
    rfpulse[RF3_XTG_SLOT].num = 1;  /* New TG */
    rfpulse[RF3_XTG_SLOT].activity = rfpulse[RF1_XTG_SLOT].activity;
    rfpulse[RF4_XTG_SLOT].num = 1;  /* New TG */
    rfpulse[RF4_XTG_SLOT].activity = rfpulse[RF1_XTG_SLOT].activity;    
        
    xtgtr = 200ms;
    
    xtg_dda = 2;
    XTGacq1 = PSD_ON;
        
    return SUCCESS;
} /* end XTGcvinit() */

/*
 *  AScvinit
 *  
 *  Type: Public Function
 *  
 *  Description:
 *  
 */
STATUS
AScvinit( void )
{
    INT asplane;
    INT ascell;

    td0as = GRAD_UPDATE_TIME;

    /* Create a pseudo scan_info table for obloptimize to work with.
       The rotation matrices for the planes are what is important. */
    for( asplane = 0; asplane < 3; asplane++ )
    {
        for( ascell = 0; ascell < 9; ascell++ )
        {
            asscan_info[asplane].oprot[ascell] = 0.0;
        }
    }

    /* Axial */
    asscan_info[0].oprot[0] = asscan_info[0].oprot[4] = asscan_info[0].oprot[8]
        = 1.0;
    /* Sagittal */
    asscan_info[1].oprot[2] = asscan_info[1].oprot[4] = asscan_info[1].oprot[6]
        = 1.0;
    /* Coronal */
    asscan_info[2].oprot[1] = asscan_info[2].oprot[5] = asscan_info[2].oprot[6]
        = 1.0;
  
    as_newgeo = 1;
    if (FAILURE==obloptimize(&asloggrd, &phygrd, asscan_info, 3, PSD_OBL,
                             0, PSD_OBL_RESTRICT, asobl_debug, &as_newgeo,
                             cfsrmode))
    {
        epic_error(use_ermes, "%s failed in AScvinit.", EM_PSD_SUPPORT_FAILURE,
                   EE_ARGS(1), STRING_ARG, "obloptimize"); 
        return FAILURE;
    }

    /* derate SR for quiet PSC */ 
    sr_derate(&asloggrd, PSassr_derate_factor);

    /* X Killer CVs */
    if (FAILURE==amppwgrad(asx_killer_area, asloggrd.tx_xz, 0.0, 0.0, asloggrd.xrt,
                           MIN_PLATEAU_TIME, &a_gxkas, &pw_gxkasa,
                           &pw_gxkas, &pw_gxkasd ))
    {
        epic_error(use_ermes, "%s failed in AScvinit.",
                   EM_PSD_SUPPORT_FAILURE,EE_ARGS(1),
                   STRING_ARG,"amppwgrad:gxkas"); 
        return FAILURE;
    }

    /* Z Killer CVs */
    if (FAILURE==amppwgrad(asz_killer_area, asloggrd.tz_xz, 0.0, 0.0, asloggrd.zrt,
                           MIN_PLATEAU_TIME, &a_gzkas, &pw_gzkasa,
                           &pw_gzkas, &pw_gzkasd ))
    {
        epic_error(use_ermes, "%s failed in AScvinit.",
                   EM_PSD_SUPPORT_FAILURE,EE_ARGS(1),
                   STRING_ARG,"amppwgrad:gzkas"); 
        return FAILURE;
    }

    /* rf1 cvs  */
    a_rf1as = 1.0;
    pw_rf1as = rfpulse[RF1_AUTOSHIM].nom_pw;
    gscale_rf1as = .90909;
    cyc_rf1as = 1;
    res_rf1as = 0; /* initialized to zero for system safety check in cveval */

    /* gzrf1 cvs */
    pw_gzrf1as = pw_rf1as;
    flip_rf1as = asflip;

    /*******************/
    /* Starting point  */
    /*******************/
    tleadas  = RUP_GRD(24us);
    bw_rf1as = 4 * cyc_rf1as / ((float)pw_rf1as / (float)1.0s);
    t_exaas  = pw_gzrf1asa + pw_rf1as / 2;

    return SUCCESS;
} /* end AScvinit() */


/*
 *  AScveval
 *  
 *  Type: Public Function
 *  
 *  Description:
 *  
 */
STATUS
AScveval( void )
{
    /* REDFLAG : Making the same as 1.0T and 1.5T (3T/4T) */
    if ((fastprescan == 1) &&
        ((cffield == B0_10000) || (cffield == B0_15000) ||
         (cffield == B0_30000) || (cffield == B0_40000)))
    {
        as_dda = 0;
        echo1bwas = 62.5;
        asbaseline = 0;
        asxres = 128;
        asyres =  64;
        /* asres=128 and asyres=64 for all T fields */
    }
    else
    {
        as_dda = 4;
        echo1bwas = 15.625;
        asxres = 256;
        asyres = 128;
        asbaseline = 8;
        te_as = 9000;
        if (cffield == B0_15000)
        {
            tr_as = 25000;
        }
        else if (cffield == B0_10000)
        {
            tr_as = 30000;
        }
        else if (cffield == B0_2000) /* profile 05/22/95 NM */
        {
            tr_as = 40000; /* Profile 09/29/95 NM */
        }
        else if (cffield == B0_3500) /* MFO,Hino, Feb/02/00 MM */
        {
            tr_as = 40000;
        }
        else if (cffield == B0_5000)
        {
            tr_as = 35000;
        }
        else if (cffield == B0_40000)
        {
            /* REDFLAG : Using same value as 1.5T. */
            tr_as = 25000;
            DEBUG_4_0(SD_PSD_SUPPORT,__FILE__,__LINE__);
        }
        else if (cffield == B0_30000)
        {
            /* REDFLAG : Using same value as 1.5T. */
            tr_as = 25000;
            DEBUG_3_0(SD_PSD_SUPPORT,__FILE__,__LINE__);
        }
        else if (cffield == B0_7000)
        {
            tr_as = 32000;
            DEBUG_0_7(SD_PSD_SUPPORT,__FILE__,__LINE__);
        }
        else 
        {
            SDL_PrintFStrengthWarning(SD_PSD_SUPPORT,cffield,__FILE__,__LINE__);
        }
    }

    /* MRIge21914 - moved deltf to @cv */

    /* Call the SDL function to compute fat-water separation. */
    deltf = SDL_GetChemicalShift( cffield );
    /* In phase delta TE will be used in MFO. Hino, Feb/02/00 MM */

    /* dixon time shift.  put it on grad boundary. */
    dix_timeas = RUP_GRD((int)(1.0s / deltf));
    pw_sdixon2 = GRAD_UPDATE_TIME + dix_timeas;

    /********************************************************/
    /*   Z Board                                            */
    /*   Slice Selection                                    */
    /********************************************************/
    if (FAILURE==ampslice(&a_gzrf1as, bw_rf1as, asslthick, gscale_rf1as, TYPDEF))
    {
        epic_error(use_ermes, "%s failed", EM_PSD_SUPPORT_FAILURE,
                   EE_ARGS(1), STRING_ARG, "ampslice for gzrf1as.");

        return FAILURE;
    }

    if (FAILURE==optramp(&pw_gzrf1asa, a_gzrf1as, asloggrd.tz_xyz, asloggrd.zrt, TYPDEF))
    {
        epic_error(use_ermes, "%s failed", EM_PSD_SUPPORT_FAILURE,
                   EE_ARGS(1), STRING_ARG, "optramp for gzrf1asa.");
        return FAILURE;
    }
    pw_gzrf1asd = pw_gzrf1asa;

    /******************************************/
    /*   Calc area needed for z rephaser      */
    /******************************************/
    area_gz1as = (off90as + pw_rf1as/2.0 +pw_gzrf1asd/2.0)*a_gzrf1as;

    /* availible time for rephaser */
    avail_pwgz1as = 1s;

    if ( FAILURE==amppwgz1(&a_gz1as, &pw_gz1as, &pw_gz1asa, &pw_gz1asd,
                           area_gz1as, avail_pwgz1as, MIN_PLATEAU_TIME,
                           asloggrd.zrt, asloggrd.tz_xyz) )
    {
        epic_error(use_ermes, "%s failed", EM_PSD_SUPPORT_FAILURE,
                   EE_ARGS(1), STRING_ARG, "amppwgz1 for gz1as.");
        return FAILURE;
    }

    /*****************************************************************/
    /*  X Board - readout and dephaser                               */
    /*****************************************************************/

    if ( FAILURE==calcfilter( &echo1as_filt, echo1bwas, asxres, OVERWRITE_NONE) ) 
    {
        epic_error(use_ermes, "%s failed.", EM_PSD_SUPPORT_FAILURE,
                   EE_ARGS(1), STRING_ARG, "calcfilter for echo1as_filt");

        return FAILURE;
    }

    if ( FAILURE==ampfov(&a_gxwas, echo1as_filt.bw, asfov) )
    {
        epic_error(use_ermes, "%s failed", EM_PSD_SUPPORT_FAILURE,
                   EE_ARGS(1), STRING_ARG, "ampfov for gxwas.");
        return FAILURE;
    }

    if (FAILURE==optramp(&pw_gxwasa, a_gxwas, asloggrd.tx_xyz, asloggrd.xrt, TYPDEF)) 
    {
        epic_error(use_ermes, "%s failed", EM_PSD_SUPPORT_FAILURE,
                   EE_ARGS(1), STRING_ARG, "optramp for pw_gxwasa.");
        return FAILURE;
    }

    pw_gxwasd = pw_gxwasa;
    pw_gxwas = echo1as_filt.tdaq;

    avail_pwgx1as = 1s;

    area_readrampas = 0.5*pw_gxwasa*a_gxwas;
    area_gxwas = pw_gxwas*a_gxwas;

    if ( FAILURE==amppwgx1(&a_gx1as, &pw_gx1as, &pw_gx1asa, &pw_gx1asd,
                           (int)TYPGRAD, area_gxwas, area_readrampas,
                           avail_pwgx1as, 1.0, MIN_PLATEAU_TIME,
                           asloggrd.xrt, asloggrd.tx_xyz) )
    {
        epic_error(use_ermes, "%s failed", EM_PSD_SUPPORT_FAILURE,
                   EE_ARGS(1), STRING_ARG, "amppwgx1 for gx1as.");
        return FAILURE;
    }

    /***************************************************************/
    /*   Y Board - Dephaser and Killer                             */
    /*                                                             */
    /*   Calculate Y Phase encode amp and pw.                      */
    /*                                                             */
    /***************************************************************/

    /* find min and max time for gy1as */
    avail_pwgy1as = 1s;

    /* Scale the waveform amps for the phase encodes 
     * so each phase instruction jump is an integer step */
    if ( FAILURE==endview((int)(asyres), &endview_iampas) )
    {
        epic_error(use_ermes, supfailfmt, EM_PSD_SUPPORT_FAILURE,
                   EE_ARGS(1), STRING_ARG, "endview:autoshim");
        return FAILURE;
    } 
  
    endview_scaleas = (float)max_pg_iamp / (float)endview_iampas;

    if ( FAILURE==amppwtpe(&a_gy1asa, &a_gy1asb, &pw_gy1as, &pw_gy1asa, &pw_gy1asd,
                           asloggrd.ty_xyz/endview_scaleas,asloggrd.yrt,
                           (0.5 * (FLOAT)(asyres-1))/(asfov * 0.1) * 1.0e6 / GAM) ) 
    {
        epic_error(use_ermes, supfailfmt, EM_PSD_SUPPORT_FAILURE,
                   EE_ARGS(1), STRING_ARG, "amppwtpe:autoshim");
        return FAILURE;
    }

    /* phase rewinders */
    a_gy1ras = a_gy1as;
    a_gy1rasa = a_gy1asa;
    a_gy1rasb = a_gy1asb;
    pw_gy1ras = pw_gy1as;
    pw_gy1rasa = pw_gy1asa;
    pw_gy1rasd = pw_gy1asd;

    if(1==fastprescan) 
    {

        te_as = RUP_GRD(pw_rf1as/2 + off90as + pw_gzrf1asd + pw_gz1asa + pw_gz1as 
                        + pw_gz1asd + pw_gy1asa + pw_gy1as + pw_gy1asd + pw_gx1asa 
                        + pw_gx1as + pw_gx1asd + pw_gxwasa + pw_gxwas/2);

        tr_as = RUP_GRD(te_as + dix_timeas + pw_gzrf1as/2 + pw_gzrf1asa + td0as 
                        + tleadas - rfupa + pw_gxwas/2 + pw_gxwasd + pw_gzkasa 
                        + pw_gzkas + pw_gzkasd + 1ms);
    }

    return SUCCESS;
}   /* end AScveval() */


/*
 *  PS1cveval
 *  
 *  Type: Private Function
 *  
 *  Description:
 *  
 */
STATUS
PS1cveval( FLOAT *opthickPS )
{
    INT bw_rf1mps1, bw_rf2mps1;    /* band widths of rf pulses */
    FLOAT area_pulse;
    FLOAT area_readrampmps1;
    FLOAT area_gxwmps1;
    FLOAT av_temp_float = 0;
    int ps1_xrt;
    float ps1_tx;
    float ps1_tx_xz;
    float ps1_tz_xz;

    /* check for breast L/R coil */
    if( (strstr(coilInfo[0].coilName, "R_BREAST") != NULL) ||
        (strstr(coilInfo[0].coilName, "L_BREAST") != NULL) ||
        (strstr(coilInfo[0].coilName, "BreastL") != NULL) ||
        (strstr(coilInfo[0].coilName, "BreastR") != NULL) ||
        (strstr(coilInfo[0].coilName, "RtBreast") != NULL) ||
        (strstr(coilInfo[0].coilName, "LtBreast") != NULL) ||
        (strstr(coilInfo[0].coilName, "breast L") != NULL) ||
        (strstr(coilInfo[0].coilName, "breast R") != NULL) ||
        (strstr(coilInfo[0].coilName, "BrstL") != NULL) ||
        (strstr(coilInfo[0].coilName, "BrstR") != NULL) ||
        (strstr(coilInfo[0].coilName, "BREASTPA R") != NULL) ||
        (strstr(coilInfo[0].coilName, "BREASTPA L") != NULL) )
    {
        ps1_rxcoil = PSD_ON;  /* Flag for R or L breast coil for TG */
    }
    else
    {
        ps1_rxcoil = PSD_OFF;
    }

    ps1_xrt = (TGspf ? ps1loggrd.yrt : ps1loggrd.xrt);
    ps1_tx = (TGspf ? ps1loggrd.ty : ps1loggrd.tx);
    ps1_tx_xz = (TGspf ? ps1loggrd.ty_yz : ps1loggrd.tx_xz);
    ps1_tz_xz = (TGspf ? ps1loggrd.tz_yz : ps1loggrd.tz_xz);

    /* Z slice select for 90 pulse */
    pw_gzrf1mps1 = pw_rf1mps1;
    bw_rf1mps1 = rfpulse[RF1_APS1_SLOT].nom_bw*rfpulse[RF1_APS1_SLOT].nom_pw/(float)pw_rf1mps1;

    /* MRIge30641 */
    /* Need check on sl thk.  Otherwise causes dwnld failures w/0.8 G/cm. */
    /* Assumes that gzrf1 gives the worst case min sl thk of all psc pulses. */ 
    if (FAILURE==ampslice(&av_temp_float, bw_rf1mps1, ps1loggrd.tz, gscale_rf1mps1, TYPDEF))
    {
        epic_error(use_ermes, "%s failed", EM_PSD_SUPPORT_FAILURE,
                   EE_ARGS(1), STRING_ARG, "ampslice");
        return FAILURE;
    }

    av_temp_float = ceil(av_temp_float*10.0)/10.0;
    if (av_temp_float > *opthickPS)
    {
        *opthickPS = av_temp_float;
    }

    if (0 == getAps1Mod()) {
        cvoverride(thickPS_mod, *opthickPS, PSD_FIX_ON, PSD_EXIST_ON);
    } else {
        float fov = FMax(2, getAps1ModSlThick(), av_temp_float);
        cvoverride(thickPS_mod, fov, PSD_FIX_ON, PSD_EXIST_ON);
    }

    if (FAILURE==ampslice(&a_gzrf1mps1, bw_rf1mps1,thickPS_mod,gscale_rf1mps1,TYPDEF))
    {
        epic_error(use_ermes, "%s failed", EM_PSD_SUPPORT_FAILURE,
                   EE_ARGS(1), STRING_ARG, "ampslice for gzrf1mps1.");
        return FAILURE;
    }
    /* end aps1_mod changes (GE) */

    /* slice selection ramp */
    if (optramp(&pw_gzrf1mps1a, a_gzrf1mps1, ps1loggrd.tz, ps1loggrd.zrt,
                TYPDEF)==FAILURE) /* vmx 5/9/95 YI */ 
    {
        epic_error(use_ermes, "%s failed", EM_PSD_SUPPORT_FAILURE,
                   EE_ARGS(1), STRING_ARG, "optramp for gzrf1mps1a.");
        return FAILURE;
    }
    pw_gzrf1mps1d = pw_gzrf1mps1a;

    /* Z gradient refocus */
    /* available time not calculated, defaulted to 10ms */
    area_pulse = a_gzrf1mps1*(pw_gzrf1mps1/2 + PSoff90 + pw_gzrf1mps1d/2);
    if (amppwgz1(&a_gz1mps1,&pw_gz1mps1,&pw_gz1mps1a,&pw_gz1mps1d,area_pulse,
                 (int)(1s),MIN_PLATEAU_TIME,ps1loggrd.zrt,ps1_tz_xz) == FAILURE)
    {
        epic_error(use_ermes, "%s failed", EM_PSD_SUPPORT_FAILURE,
                   EE_ARGS(1), STRING_ARG, "amppwgz1 for gz1mps1.");
        return FAILURE;
    }

    /* Z gradient crushers for 180 pulse */
    /* Left crusher. Denoted by the "l" after the "2"  in "gzrf2lmps1" */
    if (amppwgrad(ps_crusher_area, ps1_tz_xz, 0.0, 0.0, ps1loggrd.zrt,
                  MIN_PLATEAU_TIME, &a_gzrf2lmps1, &pw_gzrf2lmps1a,
                  &pw_gzrf2lmps1, &pw_gzrf2lmps1d) == FAILURE)
    {
        epic_error(use_ermes, "%s failed in PScveval.", EM_PSD_SUPPORT_FAILURE,
                   EE_ARGS(1), STRING_ARG, "amppwgrad:gzrf2lmps1"); 
        return FAILURE;
    }
  
    /* Right crusher. Denoted by the "r" after the "2" in "gzrf2rmps1"*/
    /* This is identical to left crusher */
    pw_gzrf2rmps1 = pw_gzrf2lmps1;
    a_gzrf2rmps1 =  a_gzrf2lmps1;

    /* right crusher ramps */
    pw_gzrf2rmps1a = pw_gzrf2lmps1a;
    pw_gzrf2rmps1d = pw_gzrf2lmps1d;

    /* Z slice select for 180 pulse */
    pw_gzrf2mps1 = pw_rf2mps1;
    bw_rf2mps1 = rfpulse[RF2_APS1_SLOT].nom_bw*rfpulse[RF2_APS1_SLOT].nom_pw/(float)pw_rf2mps1;

    /* begin aps1_mod changes (GE) */
    if (FAILURE==ampslice(&a_gzrf2mps1, bw_rf2mps1, thickPS_mod, gscale_rf2mps1, TYPDEF))
    {
        epic_error(use_ermes, "%s failed", EM_PSD_SUPPORT_FAILURE,
                   EE_ARGS(1), STRING_ARG, "ampslice for gzrf2mps1.");
        return FAILURE;
    }
    /* end aps1_mod changes (GE) */

    /* match ramps so gradient can be bridged in pulsegen */
    pw_gzrf2mps1a = pw_gzrf2lmps1d;
    pw_gzrf2mps1d = pw_gzrf2rmps1a;

    /* readout gradient */

    /* begin aps1_mod changes (GE) */
    if( (1 == getAps1Mod()) && (PSD_ON == ps1_rxcoil) )
    {
        mpsfov = getAps1ModFov();
    }
    else
    {
        mpsfov = cfsystemmaxfov;
    }

    if (FAILURE==calcfilter( &echo1mps1_filt, 15.625, 256, OVERWRITE_NONE))
    {
        epic_error(use_ermes, "%s failed.", EM_PSD_SUPPORT_FAILURE,
                   EE_ARGS(1), STRING_ARG, "calcfilter for echo1mps1_filt");
        return FAILURE;
    }

    /* MRIge30642 */
    /* Need check on fov.  Otherwise causes dwnld failures w/0.8 G/cm. */ 
    if (ampfov(&av_temp_float, echo1mps1_filt.bw, ps1_tx) == FAILURE)
    {
        epic_error(use_ermes, "%s failed", EM_PSD_SUPPORT_FAILURE,
                   EE_ARGS(1), STRING_ARG, "ampfov");
        return FAILURE;
    }

    av_temp_float = ceil(av_temp_float / 10.0) * 10.0;
    if( av_temp_float > mpsfov )
    {
        mpsfov = av_temp_float;
    }

    if (ampfov(&a_gxwmps1, echo1mps1_filt.bw, mpsfov) == FAILURE)
    {
        epic_error(use_ermes, "%s failed", EM_PSD_SUPPORT_FAILURE,
                   EE_ARGS(1), STRING_ARG, "ampfov for gxwmps1.");
        return FAILURE;
    }

    /* attack and decay ramps */
    if (optramp(&pw_gxwmps1a, a_gxwmps1, ps1_tx, ps1_xrt,
                TYPDEF)==FAILURE)
    {
        epic_error(use_ermes, "%s failed", EM_PSD_SUPPORT_FAILURE,
                   EE_ARGS(1), STRING_ARG, "optramp for gxwmps1a.");
        return FAILURE;
    }
  
    pw_gxwmps1d = pw_gxwmps1a;

    pw_gxwmps1 = echo1mps1_filt.tdaq;
  
    /* dephaser */
    area_gxwmps1 = a_gxwmps1*(pw_gxwmps1);
    area_readrampmps1 = 0.5*pw_gxwmps1a*a_gxwmps1;

    if (amppwgx1(&a_gx1mps1, &pw_gx1mps1, &pw_gx1mps1a ,&pw_gx1mps1d, TYPSPIN,
                 area_gxwmps1, (float)area_readrampmps1, 
                 (int)1s, 1.0, MIN_PLATEAU_TIME, ps1_xrt, ps1_tx_xz) == FAILURE)
    {
        epic_error(use_ermes, "%s failed", EM_PSD_SUPPORT_FAILURE,
                   EE_ARGS(1), STRING_ARG, "amppwgx1 for gx1mps1.");
        return FAILURE;
    }

    /* Y gradient is not used in MPS1 */

    return SUCCESS;
}

/*
 *  CFLcveval
 *  
 *  Type: Private Function
 *  
 *  Description:
 *  
 */
STATUS
CFLcveval( FLOAT opthickPS )
{
    LONG bw_rf1cfl;
    FLOAT area_gz1cfl;

    cflloggrd = loggrd; /* same as imaging loggrd */
    /* derate SR for quiet PSC */ 
    sr_derate(&cflloggrd, PSsr_derate_factor);

    if ((fastprescan == 1) &&
        ((cffield == B0_10000) || (cffield == B0_15000) || 
         (cffield == B0_30000) || (cffield == B0_40000)))
    {
        cfl_dda = 2;      /* BJM MRIge80347: Changed from 0 -> 2 for MGD, coil switch problem */
        cfl_nex = 1;
    }
    else
    {
        cfl_dda = 4;
        cfl_nex = 2;
    }

    if(cffield == B0_2000) 
    {
        echo1bwcfl = 10.41666;
    } 
    else if(cffield == B0_15000)
    {
        /* MRIhc54366: accommodate B0 drift using larger receive bandwidth 
           with matched spectral resolution */
        echo1bwcfl = 7.8125;         
        CFLxres = 1024;
    }
    else
    {
        echo1bwcfl = 15.625;         
        CFLxres = 1024;
    }

    /* MRIhc54366: changed hard coded number of output points to CFLxres */
    if ( FAILURE==calcfilter( &echo1cfl, echo1bwcfl, CFLxres, OVERWRITE_NONE) )
    {
        epic_error(use_ermes, "%s failed.", EM_PSD_SUPPORT_FAILURE,
                   EE_ARGS(1), STRING_ARG, "calcfilter for echo1cfl_filt");
        return FAILURE;
    }

    /* MRIhc54366: dynamic CFL excitation pulse selection: */
    if(cyc_rf1cfl == 2)
    {
        rfpulse[RF1_CFL_SLOT].abswidth = SAR_ABS_SINC2;
        rfpulse[RF1_CFL_SLOT].effwidth = SAR_PSINC2;
        rfpulse[RF1_CFL_SLOT].area = SAR_ASINC2;
        rfpulse[RF1_CFL_SLOT].dtycyc = SAR_DTYCYC_SINC2;
        rfpulse[RF1_CFL_SLOT].maxpw = SAR_MAXPW_SINC2;
        rfpulse[RF1_CFL_SLOT].max_b1 = SAR_MAXB1_SINC2_90;
        rfpulse[RF1_CFL_SLOT].max_int_b1_sq = SAR_MAX_INT_B1_SQ_SINC2_90;
        rfpulse[RF1_CFL_SLOT].max_rms_b1 = SAR_MAX_RMS_B1_SINC2_90;
        rfpulse[RF1_CFL_SLOT].nom_bw = 2500;
    }
    else
    {
        rfpulse[RF1_CFL_SLOT].abswidth = SAR_ABS_SINC1;
        rfpulse[RF1_CFL_SLOT].effwidth = SAR_PSINC1;
        rfpulse[RF1_CFL_SLOT].area = SAR_ASINC1;
        rfpulse[RF1_CFL_SLOT].dtycyc = SAR_DTYCYC_SINC1;
        rfpulse[RF1_CFL_SLOT].maxpw = SAR_MAXPW_SINC1;
        rfpulse[RF1_CFL_SLOT].max_b1 = MAX_B1_SINC1_90;
        rfpulse[RF1_CFL_SLOT].max_int_b1_sq = MAX_INT_B1_SQ_SINC1_90;
        rfpulse[RF1_CFL_SLOT].max_rms_b1 = MAX_RMS_B1_SINC1_90;
        rfpulse[RF1_CFL_SLOT].nom_bw = 1250;
    }
    /* MRIhc54366: END pulse selection. */

    /* CFL acq duration needed for attenuator setting */
    cfl_tdaq = echo1cfl.tdaq;

    pw_gzrf1cfl = pw_rf1cfl;
    bw_rf1cfl = rfpulse[RF1_CFL_SLOT].nom_bw*rfpulse[RF1_CFL_SLOT].nom_pw/(float)pw_rf1cfl;

    /* MRIhc54366: new lower limit */
    opthickPS = (exist(opslthick) < 5.0) ? 5.0 : exist(opslthick);

    if ( FAILURE==ampslice(&a_gzrf1cfl, bw_rf1cfl, opthickPS, gscale_rf1cfl, TYPDEF) ) 
    {
        epic_error(use_ermes, "%s failed for gzrf1cfl.",
                   EM_PSD_SUPPORT_FAILURE, EE_ARGS(1), STRING_ARG, "ampslice");
        return FAILURE;
    }

    if ( FAILURE==optramp(&pw_gzrf1cfla, a_gzrf1cfl, cflloggrd.tz, cflloggrd.zrt, TYPDEF) )  
    {
        epic_error(use_ermes, "%s failed for gzrf1cfl.", 
                   EM_PSD_SUPPORT_FAILURE, EE_ARGS(1), STRING_ARG, "optramp");
        return FAILURE;
    }

    pw_gzrf1cfld = pw_gzrf1cfla;

    /* Find Params for refocusing pulse */
    area_gz1cfl =  a_gzrf1cfl *0.5* ( pw_gzrf1cfl + pw_gzrf1cfld);
    if ( FAILURE==amppwgz1(&a_gz1cfl, &pw_gz1cfl, &pw_gz1cfla, &pw_gz1cfld, 
                           area_gz1cfl, (INT)1s, MIN_PLATEAU_TIME,
                           cflloggrd.zrt, cflloggrd.tz) ) 
    {
        epic_error(use_ermes, "%s failed in cfl.", EM_PSD_SUPPORT_FAILURE,
                   EE_ARGS(1), STRING_ARG, "amppwgz1 for gz1cfl");
        return FAILURE;
    }

    /* Find Params for killer pulse */
    area_gykcfl = amp_killer*pw_killer;
    if ( FAILURE==amppwgrad(area_gykcfl, cflloggrd.ty, 0.0, 0.0, cflloggrd.yrt,
                            MIN_PLATEAU_TIME, &a_gykcfl, &pw_gykcfla,
                            &pw_gykcfl, &pw_gykcfld) ) 
    {
        epic_error(use_ermes, "%s failed in cfl.", EM_PSD_SUPPORT_FAILURE,
                   EE_ARGS(1), STRING_ARG, "amppwgrad:gykcfl");
        return FAILURE;
    }

    return SUCCESS;
}


/*
 *  RCVNcveval
 *  
 *  Type: Private Function
 *  
 *  Description:
 *  
 */
STATUS
RCVNcveval( void )
{
    rcvnloggrd = loggrd; /* same as imaging loggrd */
    /* derate SR for quiet PSC */ 
    sr_derate(&rcvnloggrd, PSsr_derate_factor);

    /* MRIhc47602/MRIhc47515/GEHmr03545 : Killer gradient before Receiver noise sequence */

    if ( FAILURE==amppwgrad( area_gxkrcvn, rcvnloggrd.tx_xyz, 0.0, 0.0, rcvnloggrd.xrt,
                             MIN_PLATEAU_TIME, &a_gxkrcvn, &pw_gxkrcvna,
                             &pw_gxkrcvn, &pw_gxkrcvnd ) )
    {
        epic_error(use_ermes, "Support routine %s failed", EM_PSD_SUPPORT_FAILURE,
                   EE_ARGS(1), STRING_ARG, "amppwgrad:gxkrcvn");
        return FAILURE;
    }

    if ( FAILURE==amppwgrad( area_gykrcvn, rcvnloggrd.ty_xyz, 0.0, 0.0, rcvnloggrd.yrt,
                             MIN_PLATEAU_TIME, &a_gykrcvn, &pw_gykrcvna,
                             &pw_gykrcvn, &pw_gykrcvnd ) )
    {
        epic_error(use_ermes, "Support routine %s failed", EM_PSD_SUPPORT_FAILURE,
                   EE_ARGS(1), STRING_ARG, "amppwgrad:gykrcvn");
        return FAILURE;
    }

    if ( FAILURE==amppwgrad( area_gzkrcvn, rcvnloggrd.tz_xyz, 0.0, 0.0, rcvnloggrd.zrt,
                             MIN_PLATEAU_TIME, &a_gzkrcvn, &pw_gzkrcvna,
                             &pw_gzkrcvn, &pw_gzkrcvnd ) )
    {
        epic_error(use_ermes, "Support routine %s failed", EM_PSD_SUPPORT_FAILURE,
                   EE_ARGS(1), STRING_ARG, "amppwgrad:gzkrcvn");
        return FAILURE;
    }

    if ( rcvn_flag == 1 )
    {
        pre_rcvn_tr = IMax(2, 20ms, RUP_GRD(pw_gxkrcvna+pw_gxkrcvn+pw_gxkrcvnd+1ms));
    }
    else if(rcvn_flag == 2)  /* extra delay before rcvn */
    {
        pre_rcvn_tr = 1000ms;
    }

    if (existcv(oprbw))
    {
        echo1bwrcvn = exist(oprbw);
    }

    echo1bwrcvn = FMax(2, echo1bwrcvn, (float) RCVN_MIN_BW);

    echo1bwrcvn = FMin(2, echo1bwrcvn, (float) RCVN_MAX_BW);

    if ( FAILURE==calcfilter( &echo1rcvn, echo1bwrcvn, rcvn_xres, OVERWRITE_NONE) ) 
    {
        epic_error(use_ermes, "%s failed.", EM_PSD_SUPPORT_FAILURE,
                   EE_ARGS(1), STRING_ARG, "calcfilter for echo1rcvn");
        return FAILURE;
    }

    /* RCVN acq duration needed for attenuator setting */
    rcvn_tdaq = echo1rcvn.tdaq;

    /* Add 20ms dead time to prevent runtime errors */
    rcvn_tr = IMax(2, rcvn_tdaq + 20ms, (int) RCVN_MIN_TR);

    /* To make sure RCVN loops enough to acquire 4K points */
    rcvn_loops = IMax(2, 2 * (int)(4096 / rcvn_xres), 1);

    return SUCCESS;
}


/*
 *  CFHfilter
 *  
 *  Type: Private Function
 *  
 *  Description: Separate function for CFH for inclusion in 
 *               in Spectroscopy volume localized CFH
 *  
 */
STATUS
CFHfilter( void )
{
    /* MRIhc08595: Replaced cfh filter output points, hard coded 256, with echo1ptcfh CV. */
    echo1ptcfh = 256;

    if(cffield <= B0_5000) 
    {
        echo1bwcfh = 0.25;
    }
    else if (cffield >= B0_30000) 
    {
        echo1bwcfh = 1.0;
    }
    else 
    {
        echo1bwcfh = 0.50;
    }

    if ( FAILURE==calcfilter( &echo1cfh, echo1bwcfh, echo1ptcfh, OVERWRITE_NONE) ) 
    {
        epic_error(use_ermes, "%s failed.", EM_PSD_SUPPORT_FAILURE,
                   EE_ARGS(1), STRING_ARG, "calcfilter for echo1cfh_filt");
        return FAILURE;
    }

    cfh_tdaq = echo1cfh.tdaq;

    return SUCCESS;
}

/*
 *  CFHcveval
 *  
 *  Type: Private Function
 *  
 *  Description:
 *  
 */
STATUS
CFHcveval( FLOAT opthickPS )
{
    INT bw_rf0cfh;
    INT bw_rf1cfh;

    GRAD_PULSE psd_cfhrightcrush;
    GRAD_PULSE psd_cfhleftcrush;

    FLOAT area_gz1cfh;
    a_gyrf3cfh = 0.0;/*For MRIhc11621 */
    /* For presscfh MRIhc08321 */
    if( presscfh == PRESSCFH_SLICE && PSfield_strength > B0_5000 && 
        cfh_newmode && exist(oppscvquant)>= 1 ) 
    {
        cvoverride(presscfh_ctrl, PRESSCFH_SLICE, PSD_FIX_ON, PSD_EXIST_ON);
        presscfh_slthick = opthickPS;
        presscfh_fov_ratio = 0.8;
        presscfh_pfov_ratio = 0.8;
        presscfh_slab_ratio = 0.5;
    } 
    else if(presscfh == PRESSCFH_SLAB && PSfield_strength > B0_5000 &&
            cfh_newmode && exist(oppscvquant)>= 1 ) 
    {
        cvoverride(presscfh_ctrl, PRESSCFH_SLAB, PSD_FIX_ON, PSD_EXIST_ON);
        presscfh_fov_ratio = 0.8;
        presscfh_pfov_ratio = 0.8;
        presscfh_slab_ratio = 0.75;
    }
    else if( presscfh == PRESSCFH_SHIMVOL && PSfield_strength > B0_5000 &&
             cfh_newmode && exist(oppscvquant) >= 1 ) 
    {
        cvoverride(presscfh_ctrl, PRESSCFH_SHIMVOL, PSD_FIX_ON, PSD_EXIST_ON); 
        presscfh_fov_ratio = 0.5;
        presscfh_pfov_ratio = 0.5;
        presscfh_slab_ratio = 0.5;
    }
    else
    { 
        cvoverride(presscfh_ctrl, PRESSCFH_NONE, PSD_FIX_ON, PSD_EXIST_ON); 
    }
    
    if(presscfh_debug) 
    {
        printf("\n CFHcveval : presscfh = %d,presscfh_ctrl = %d,presscfh_override = %d\n",presscfh,presscfh_ctrl,presscfh_override);
        fflush(stdout);
    }

    /* REDFLAG : Making the same as 1.0T and 1.5T (3T/4T) */
    if ((fastprescan == 1) &&
        ((cffield == B0_10000) || (cffield == B0_15000) || 
         (cffield == B0_30000) || (cffield == B0_40000)))
    {
        cfh_dda = 0;
        cfh_nex = 1;
    }
    else
    {
        cfh_dda = 4;
        cfh_nex = 2;
    }

    if(cffield == B0_15000) 
    {
        cfh_ti = CFHTI_1HT;
        eff_cfh_te = CFHTE_1HT;
    }
    else if(cffield >= B0_30000) 
    {
        cfh_ti = CFHTI_3T;
        eff_cfh_te = CFHTE_3T;
    }

    CFHfilter();

    /* Initialize some grad structures 
       so we can use the psdsupport routine amppwlcrsh */
    psd_cfhleftcrush.attack = &pw_gzrf2lcfha;
    psd_cfhleftcrush.decay = &pw_gzrf2lcfhd;
    psd_cfhleftcrush.pw = &pw_gzrf2lcfh;
    psd_cfhleftcrush.amp = &a_gzrf2lcfh;
 
    psd_cfhrightcrush.attack = &pw_gzrf2rcfha;
    psd_cfhrightcrush.decay = &pw_gzrf2rcfhd;
    psd_cfhrightcrush.pw = &pw_gzrf2rcfh;
    psd_cfhrightcrush.amp = &a_gzrf2rcfh;

    if( presscfh_ctrl != PRESSCFH_NONE && cfh_steam_flag == PSD_ON )
    {
        psd_cfhrightcrush.attack = &pw_gzrf3rcfha;
        psd_cfhrightcrush.decay = &pw_gzrf3rcfhd;
        psd_cfhrightcrush.pw = &pw_gzrf3rcfh;
        psd_cfhrightcrush.amp = &a_gzrf3rcfh;
    }

    if( presscfh_ctrl != PRESSCFH_NONE && cfh_steam_flag == PSD_ON )
    {
        setuprfpulse(RF1_CFH_SLOT, &pw_rf1cfh, &a_rf1cfh, SAR_ABS_SINC3, SAR_PSINC3,
                     SAR_ASINC3, SAR_DTYCYC_SINC3, SAR_MAXPW_SINC3, 1,
                     MAX_B1_SINC3_90, MAX_INT_B1_SQ_SINC3_90,
                     MAX_RMS_B1_SINC3_90, 90.0, &flip_rf1cfh, 3200.0,
                     3750, PSD_CFH_ON, 0,
                     0, 0, &res_rf1cfh, 0, &wg_rf1cfh, rfpulse);

        setuprfpulse(RF2_CFH_SLOT, &pw_rf2cfh, &a_rf2cfh, SAR_ABS_SINC3, SAR_PSINC3,
                     SAR_ASINC3, SAR_DTYCYC_SINC3, SAR_MAXPW_SINC3, 1,
                     MAX_B1_SINC3_90, MAX_INT_B1_SQ_SINC3_90,
                     MAX_RMS_B1_SINC3_90, 90.0, &flip_rf2cfh, 3200.0,
                     3750, PSD_CFH_ON, 0,
                     0, 0, &res_rf2cfh, 0, &wg_rf2cfh, rfpulse);

        setuprfpulse(RF3_CFH_SLOT, &pw_rf3cfh, &a_rf3cfh, SAR_ABS_SINC3, SAR_PSINC3,
                     SAR_ASINC3, SAR_DTYCYC_SINC3, SAR_MAXPW_SINC3, 1,
                     MAX_B1_SINC3_90, MAX_INT_B1_SQ_SINC3_90,
                     MAX_RMS_B1_SINC3_90, 90.0, &flip_rf3cfh, 3200.0,
                     3750, PSD_CFH_ON, 0,
                     0, 0, &res_rf3cfh, 0, &wg_rf3cfh, rfpulse);

        a_rf0cfh = 1;
        a_rf1cfh = 0.5464;
        a_rf2cfh = 0.5464;
        a_rf3cfh = 0.5464; 

        flip_rf1cfh = 60;
        flip_rf2cfh = 60;
        flip_rf3cfh = 60; 

        cyc_rf1cfh = 3;
        cyc_rf2cfh = 3;
        cyc_rf3cfh = 3; 
    }
    else
    {
        setuprfpulse(RF1_CFH_SLOT, &pw_rf1cfh, &a_rf1cfh, SAR_ABS_SINC1, SAR_PSINC1,
                     SAR_ASINC1, SAR_DTYCYC_SINC1, SAR_MAXPW_SINC1, 1,
                     MAX_B1_SINC1_90, MAX_INT_B1_SQ_SINC1_90,
                     MAX_RMS_B1_SINC1_90, 90.0, &flip_rf1cfh, 3200.0,
                     1250, PSD_CFH_ON, 0,
                     0, 0, &res_rf1cfh, 0, &wg_rf1cfh, rfpulse);

        setuprfpulse(RF2_CFH_SLOT, &pw_rf2cfh, &a_rf2cfh, SAR_ABS_SINC1, SAR_PSINC1,
                     SAR_ASINC1, SAR_DTYCYC_SINC1, SAR_MAXPW_SINC1, 1,
                     MAX_B1_SINC1_90, MAX_INT_B1_SQ_SINC1_90,
                     MAX_RMS_B1_SINC1_90, 90.0, &flip_rf2cfh, 3200.0,
                     1250, PSD_CFH_ON, 0,
                     0, 0, &res_rf2cfh, 0, &wg_rf2cfh, rfpulse);

        setuprfpulse(RF3_CFH_SLOT, &pw_rf3cfh, &a_rf3cfh, SAR_ABS_SINC1, SAR_PSINC1,
                     SAR_ASINC1, SAR_DTYCYC_SINC1, SAR_MAXPW_SINC1, 1,
                     MAX_B1_SINC1_90, MAX_INT_B1_SQ_SINC1_90,
                     MAX_RMS_B1_SINC1_90, 90.0, &flip_rf3cfh, 3200.0,
                     1250, PSD_CFH_ON, 0,
                     0, 0, &res_rf3cfh, 0, &wg_rf3cfh, rfpulse);

        a_rf0cfh = 0.61;
        a_rf1cfh = 0.5;
        a_rf2cfh = 1.0;
        a_rf3cfh = 1.0; 

        flip_rf1cfh = 90;
        flip_rf2cfh = 180;
        flip_rf3cfh = 180; 

        cyc_rf1cfh = 1;
        cyc_rf2cfh = 1;
        cyc_rf3cfh = 1; 

    }

    if( presscfh_ctrl == PRESSCFH_SLICE || presscfh_ctrl == PRESSCFH_SLAB ) 
    {
        /* it is assumed: 
         * (a) cubicle local shim volume 
         */

        FLOAT av;
        FLOAT dxv, dyv, dzv, dxs, dys, dzs, dxs0, dys0, dzs0; 
        FLOAT Dz, al2, al, dxl, dyl, dzl;
        FLOAT rs[9], rv[9], rstrv[9];
        INT ii, vidx;
       
        for( vidx = 0; vidx < exist(oppscvquant); vidx ++ )
        {
            av = (psc_info[vidx].oppsclenx)/2.0;
            dxv = psc_info[vidx].oppscrloc;
            dyv = psc_info[vidx].oppscphasoff;
            dzv = psc_info[vidx].oppsctloc;

            dxs = scan_info[PSslice_num].oprloc;
            dys = scan_info[PSslice_num].opphasoff;
            dzs = scan_info[PSslice_num].optloc;

            for( ii = 0; ii < 9; ii++ )
            {
                rs[ii] = scan_info[PSslice_num].oprot[ii];
                rv[ii] = psc_info[vidx].oppscrot[ii];
            }

            for( ii = 0; ii < 9; ii++ )
            {
                int ir, ic;
                ir = ii / 3;
                ic = ii % 3;
                rstrv[ii] = rs[3*0+ir]*rv[3*0+ic] 
                    + rs[3*1+ir]*rv[3*1+ic] 
                    + rs[3*2+ir]*rv[3*2+ic]; 
            }

            if( presscfh_debug )
            {
                printf("rot for shim volume %d\n", vidx);
                printf("%8.2f %8.2f %8.2f\n", rv[0], rv[3], rv[6]);
                printf("%8.2f %8.2f %8.2f\n", rv[1], rv[4], rv[7]);
                printf("%8.2f %8.2f %8.2f\n", rv[2], rv[5], rv[8]);
                printf("rot for slice\n");
                printf("%8.2f %8.2f %8.2f\n", rs[0], rs[3], rs[6]);
                printf("%8.2f %8.2f %8.2f\n", rs[1], rs[4], rs[7]);
                printf("%8.2f %8.2f %8.2f\n", rs[2], rs[5], rs[8]);
            }

            dxs0 = rstrv[0]*dxv + rstrv[1]*dyv + rstrv[2]*dzv;
            dys0 = rstrv[3]*dxv + rstrv[4]*dyv + rstrv[5]*dzv;
            dzs0 = rstrv[6]*dxv + rstrv[7]*dyv + rstrv[8]*dzv;

            if( presscfh == PRESSCFH_SLICE ) 
            {
                presscfh_slthick = opthickPS;
                if( abs(dzs0 - dzs) >= av ) 
                {
                    presscfh_outrange = 1;
                } else {
                    presscfh_outrange = 0;
                }
            } 
            else 
            {   /* presscfh == PRESSCFH_SLAB */
                FLOAT dzss, dzse, dzs0s, dzs0e;
                FLOAT dz1, dz2;
                dzss = scan_info[0].optloc;
                dzse = scan_info[opslquant*opvquant-1].optloc;
                dzs0s = dzs0 - av/2;
                dzs0e = dzs0 + av/2;

                /* find dz1 and dz2 */
                if( dzs0e - dzs0s < 0 ) 
                {
                    FLOAT temp;
                    temp = dzs0e;
                    dzs0e = dzs0s;
                    dzs0s = temp;
                }
                if( dzse - dzss < 0 ) 
                {
                    FLOAT temp;
                    temp = dzse;
                    dzse = dzss;
                    dzss = temp;
                }

                dz1 = dzs0s;
                dz2 = dzs0e;

                if( dzss > dz1 ) 
                { 
                    dz1 = dzss;
                }
                if( dzse < dz2 ) 
                {
                    dz2 = dzse;
                }

                if( dzss >= dz2 || dzse <= dz1 ) 
                {
                    presscfh_outrange = 1;
                }
                else
                {
                    dzs = dz1 + (dz2-dz1)/2;   
                    presscfh_slthick = presscfh_slab_ratio*(dz2-dz1);
                    if( presscfh_slthick < opthickPS ) 
                    {
                        presscfh_slthick = opthickPS;
                    }
                    presscfh_outrange = 0;
                }         
            }

            if( presscfh_outrange == 0 ) 
            {
                Dz = dzs0 - dzs;
                al2 = av*av - Dz*Dz;
                if( al2 > 0 ) 
                {
                    al = sqrt(al2);
                } else {
                    al = 0;
                }
            }

            if( al >= presscfh_minfov_ratio*av ) 
            {
                dxl = dxs0;
                dyl = dys0;
                dzl = dzs;

                presscfh_fov = presscfh_fov_ratio*al*2.0;
                presscfh_pfov = presscfh_pfov_ratio*al*2.0;

                presscfh_info[vidx].oppsctloc = dzl;
                presscfh_info[vidx].oppscrloc = dxl;
                presscfh_info[vidx].oppscphasoff = dyl;
                for( ii = 0; ii < 9; ii++ ) 
                {
                    presscfh_info[vidx].oppscrot[ii] = scan_info[PSslice_num].oprot[ii]; 
                }
                presscfh_info[vidx].oppsclenx = (INT)(presscfh_fov_ratio*al*2.0);
                presscfh_info[vidx].oppscleny = (INT)(presscfh_pfov_ratio*al*2.0);
                presscfh_info[vidx].oppsclenz = (INT)presscfh_slthick;
            } 
            else 
            {
                presscfh_outrange = 1;
            }

            if( presscfh_debug ) 
            {
                printf("av, dxv, dyv, dzv: %8.2f, %8.2f, %8.2f, %8.2f\n", av, dxv, dyv, dzv);    
                printf("dxs, dys, dzs: %8.2f, %8.2f, %8.2f\n", dxs, dys, dzs);    
                printf("al, dxl, dyl, dzl: %8.2f, %8.2f, %8.2f, %8.2f\n", al, dxl, dyl, dzl);    
                printf("presscfh_slthick: %8.2f\n", presscfh_slthick);    
                printf("outrange: %d\n", presscfh_outrange);    
            }

            if( presscfh_outrange == 1 ) 
            {
                cvoverride(presscfh_ctrl, PRESSCFH_NONE, PSD_FIX_ON, PSD_EXIST_ON); 
                if(presscfh_debug) 
                {
                    printf("\nCFHcveval : THIS TURNS OUT TO BE PRESSCFH_NONE,but was initially %d\n",presscfh);
                    fflush(stdout);
                }
                break;
            }
        }
    } 

    if (presscfh_ctrl == PRESSCFH_SHIMVOL ) 
    {
        INT vidx = 0;
        for( vidx = 0; vidx < exist(oppscvquant); vidx ++ )
        {
            presscfh_fov = presscfh_fov_ratio * psc_info[vidx].oppsclenx;
            presscfh_pfov = presscfh_pfov_ratio * psc_info[vidx].oppscleny;
            presscfh_slthick = presscfh_slab_ratio * psc_info[vidx].oppsclenz;

            presscfh_info[vidx].oppsclenx = (INT)(presscfh_fov);
            presscfh_info[vidx].oppscleny = (INT)(presscfh_pfov);
            presscfh_info[vidx].oppsclenz = (INT)(presscfh_slthick);

            if( presscfh_debug ) 
            {
                printf(": SHIM VOLUME %d presscfh=2; presscfh_fov=%8.2f\n", vidx, presscfh_fov);
                fflush(stdout);
            }
        }
    }

    if( presscfh_ctrl != PRESSCFH_NONE ) 
    {
        INT vidx = 0;
        for( vidx = 0; vidx < exist(oppscvquant); vidx ++ )
        {
            if( presscfh_fov < presscfh_info[vidx].oppsclenx )
            {
                presscfh_fov = presscfh_info[vidx].oppsclenx;
            }
            if( presscfh_pfov < presscfh_info[vidx].oppscleny )
            {
                presscfh_pfov = presscfh_info[vidx].oppscleny;
            }
            if( presscfh_slthick < presscfh_info[vidx].oppsclenz )
            {
                presscfh_slthick = presscfh_info[vidx].oppsclenz;
            }
        }

        presscfh_ir_slthick = presscfh_slthick;

        rfpulse[RF3_CFH_SLOT].num = 1;
        rfpulse[RF3_CFH_SLOT].activity = PSD_CFH_ON;

    } else {
        rfpulse[RF3_CFH_SLOT].num = 0;
        rfpulse[RF3_CFH_SLOT].activity = PSD_PULSE_OFF;

    }

    /* If inversion recovery image mode, Calcs for the Inversion pulse */
    if (PSD_ON == PSir)
    {
        res_rf0cfh = RES_SH_ADIABATIC;  /* Adiabatic pulse */
        /* MRIge90312 -- use 1.5sec TR for IR cfh */
        cfh_tr = 1500ms;
        gscale_rf0cfh = 0.87; /* Changed from .65 to .87 for adiabatic pulse */

        pw_gzrf0cfh   = pw_rf0cfh;

        /* Y Killer CVs */ /* YMSmr09211  04/26/2006 YI */
        if(amppwgrad(cfhir_killer_area, cfhloggrd.ty_yz, 0.0, 0.0, cfhloggrd.yrt,
                     MIN_PLATEAU_TIME, &a_gyrf0kcfh, &pw_gyrf0kcfha, &pw_gyrf0kcfh, &pw_gyrf0kcfhd) == FAILURE) 
        {
            epic_error(use_ermes, "%s failed in PScveval.", EM_PSD_SUPPORT_FAILURE,
                       EE_ARGS(1), STRING_ARG, "amppwgrad:gyrf0kcfh"); 
            return FAILURE;
        }

        rfpulse[RF0_CFH_SLOT].num = 1;
        rfpulse[RF0_CFH_SLOT].activity = PSD_CFH_ON;
        bw_rf0cfh = 5.12*cyc_rf0cfh/((FLOAT)pw_rf0cfh/(FLOAT)1.0s); /* adiabatic pulse */
 
        if(ampslice(&a_gzrf0cfh, bw_rf0cfh, ((presscfh_ctrl == PRESSCFH_NONE) ? opthickPS : presscfh_ir_slthick),
                    gscale_rf0cfh, TYPDEF) == FAILURE) 
        {
            epic_error(use_ermes, "%s failed for gzrf0cfh.",EM_PSD_SUPPORT_FAILURE, 
                       EE_ARGS(1), STRING_ARG, "ampslice");
            return FAILURE;
        }
        /* Non Selective IR */
        if( (presscfh_ctrl != PRESSCFH_NONE) && presscfh_ir_noselect ) 
        {
            a_gzrf0cfh = 0;
        }
        /* YMSmr09211  04/26/2006 YI */
        if(optramp(&pw_gzrf0cfha,a_gzrf0cfh, cfhloggrd.tz, cfhloggrd.zrt, TYPDEF)==FAILURE) 
        {
            epic_error(use_ermes, "%s failed", EM_PSD_SUPPORT_FAILURE, EE_ARGS(1),
                       STRING_ARG, "optramp for gzrf0cfha.");
            return FAILURE;
        }
        pw_gzrf0cfhd = pw_gzrf0cfha;
    } else {
        cfh_tr = 398ms;
        rfpulse[RF0_CFH_SLOT].num = 0;
        rfpulse[RF0_CFH_SLOT].activity = PSD_PULSE_OFF;
    }

    /* Calculations for the 90 pulse */
    pw_gzrf1cfh = pw_rf1cfh;
    bw_rf1cfh = rfpulse[RF1_CFH_SLOT].nom_bw*rfpulse[RF1_CFH_SLOT].nom_pw/(float)pw_rf1cfh;
       
    if (ampslice(&a_gzrf1cfh, bw_rf1cfh, ( (presscfh_ctrl == PRESSCFH_NONE) ? opthickPS : presscfh_slthick),
                 gscale_rf1cfh, TYPDEF) == FAILURE) 
    {
        epic_error(use_ermes, "%s failed for gzrf1cfh.", EM_PSD_SUPPORT_FAILURE,
                   EE_ARGS(1), STRING_ARG, "ampslice");
        return FAILURE;
    }

    /* YMSmr09211  04/26/2006 YI */
    if (optramp(&pw_gzrf1cfha, a_gzrf1cfh, cfhloggrd.tz, cfhloggrd.zrt, TYPDEF) == FAILURE) 
    {
        epic_error(use_ermes, "%s failed for gzrf1cfh.", 
                   EM_PSD_SUPPORT_FAILURE, EE_ARGS(1), STRING_ARG, "optramp");
        return FAILURE;
    }

    pw_gzrf1cfhd = pw_gzrf1cfha;

    /* Calculations for the 180 pulse */
    /* Have to rely on the PSD to keep these up to date - can't tell from
       here whether or not the PSD is a mempoid. */
    pw_gxrf2cfh = pw_rf2cfh;
    pw_gyrf2cfh = pw_rf2cfh;
    pw_gyrf3cfh = pw_rf3cfh; /* for presscfh */

    /* Find amplitudes for the FOV selective traps */
    if (opcoax==1) 
    {
        /* If coaxial through isocenter */
        cfh_fov = FMin(2, 40.0, opfov/10.0);
    } else {
        /* Otherwise open it up */
        cfh_fov = 40.0;
    }

    if(cfh_newmode) 
    {   /*override with new mode*/
        float cfh_new_fov = opspf ? (opphasefov*opfov/10.0) : (opfov/10.0);

        cfh_fov = FMin(2, 40.0, cfh_new_fov);
        /* For non-coaxials (multi angle), center cfh_rf2freq and up the
           excitation region to 40 cm to cover lots of ground */
        cfh_fov = (opcoax == 0) ? 40.0 : FMax(2,(float)FOV_MIN,cfh_fov);
    }

    if( presscfh_ctrl != PRESSCFH_NONE ) 
    {
        /* X FOV Selective */ /* YMSmr09211  04/26/2006 YI */
        a_gxrf2cfh = 4*cyc_rf2cfh/(GAM*(float)pw_rf2cfh/(float)(1.0s)*presscfh_fov/10.0);
        if (optramp(&pw_gxrf2cfha, a_gxrf2cfh, cfhloggrd.tx, cfhloggrd.xrt, TYPDEF) == FAILURE) 
        {
            epic_error(use_ermes, "%s failed for gxrf2cfh.", EM_PSD_SUPPORT_FAILURE, EE_ARGS(1),
                       STRING_ARG, "optramp for gxrf2cfh"); 
            return FAILURE;
        }
        pw_gxrf2cfha = IMax(2,(int)(2*GRAD_UPDATE_TIME), pw_gxrf2cfha);
        pw_gxrf2cfhd = pw_gxrf2cfha;
        target_cfh_crusher = cfhloggrd.tz_xz; /* YMSmr09211  04/26/2006 YI */

        /* Y FOV Selective */ /* YMSmr09211  04/26/2006 YI */
        a_gyrf3cfh = 4*cyc_rf3cfh/(GAM*(float)pw_rf3cfh/(float)(1.0s)*presscfh_pfov/10.0);
        if (optramp(&pw_gyrf3cfha, a_gyrf3cfh, cfhloggrd.ty, cfhloggrd.yrt, TYPDEF) == FAILURE) 
        {
            epic_error(use_ermes, "%s failed for gyrf3cfh.", EM_PSD_SUPPORT_FAILURE, EE_ARGS(1),
                       STRING_ARG, "optramp for gyrf3cfh");
            return FAILURE;
        }
        pw_gyrf3cfha = IMax(2,(int)(2*GRAD_UPDATE_TIME), pw_gyrf3cfha);
        pw_gyrf3cfhd = pw_gyrf3cfha;
        target_cfh_crusher2 = cfhloggrd.tz_yz; /* YMSmr09211  04/26/2006 YI */

        if( cfh_steam_flag != PSD_ON )
        {
            /* Z CRUSHER CVs */ /* YMSmr09211  04/26/2006 YI */
            if (amppwgrad(cfh_crusher_area, target_cfh_crusher, 0.0, 0.0, cfhloggrd.zrt,
                          MIN_PLATEAU_TIME, &a_gzrf2rcfh, &pw_gzrf2rcfha,
                          &pw_gzrf2rcfh, &pw_gzrf2rcfhd) == FAILURE) 
            {
                epic_error(use_ermes, "%s failed in PScveval.",
                           EM_PSD_SUPPORT_FAILURE, EE_ARGS(1),
                           STRING_ARG, "amppwgrad:gzrf2rcfh"); 
                return FAILURE;
            }
            pw_gzrf2lcfha = pw_gzrf2rcfha;
            pw_gzrf2lcfhd = pw_gzrf2rcfhd;
        }
        else
        {
            if (amppwgrad(cfh_crusher_area, target_cfh_crusher2, 0.0, 0.0, cfhloggrd.zrt,
                          MIN_PLATEAU_TIME, &a_gzrf3rcfh, &pw_gzrf3rcfha,
                          &pw_gzrf3rcfh, &pw_gzrf3rcfhd) == FAILURE) 
            {
                epic_error(use_ermes, "%s failed in PScveval.",
                           EM_PSD_SUPPORT_FAILURE, EE_ARGS(1),
                           STRING_ARG, "amppwgrad:gzrf3rcfh"); 
                return FAILURE;
            }
            pw_gzrf2lcfha = pw_gzrf3rcfha;
            pw_gzrf2lcfhd = pw_gzrf3rcfhd;
        }

        if( cfh_steam_flag != PSD_ON )
        {
            /* Z CRUSHER CVs */ /* YMSmr09211  04/26/2006 YI */
            if (amppwgrad(cfh_crusher_area, target_cfh_crusher2, 0.0, 0.0, cfhloggrd.zrt,
                          MIN_PLATEAU_TIME, &a_gzrf3rcfh, &pw_gzrf3rcfha,
                          &pw_gzrf3rcfh, &pw_gzrf3rcfhd) == FAILURE) 
            {
                epic_error(use_ermes, "%s failed in PScveval.", EM_PSD_SUPPORT_FAILURE, EE_ARGS(1),
                           STRING_ARG, "amppwgrad:gzrf3rcfh"); 
                return FAILURE;
            }

            /* Z CRUSHER CVs */ /* YMSmr09211  04/26/2006 YI */
            if (amppwgrad(cfh_crusher_area, target_cfh_crusher2, 0.0, 0.0, cfhloggrd.zrt,
                          MIN_PLATEAU_TIME, &a_gzrf3lcfh, &pw_gzrf3lcfha,
                          &pw_gzrf3lcfh, &pw_gzrf3lcfhd) == FAILURE) 
            {
                epic_error(use_ermes, "%s failed in PScveval.", EM_PSD_SUPPORT_FAILURE, EE_ARGS(1),
                           STRING_ARG, "amppwgrad:gzrf3lrcfhd"); 
                return FAILURE;
            }
        }
        else
        {
            if (amppwgrad(cfh_crusher_area, target_cfh_crusher, 0.0, 0.0, cfhloggrd.zrt,
                          MIN_PLATEAU_TIME, &a_gzrf2rcfh, &pw_gzrf2rcfha,
                          &pw_gzrf2rcfh, &pw_gzrf2rcfhd) == FAILURE) 
            {
                epic_error(use_ermes, "%s failed in PScveval.", EM_PSD_SUPPORT_FAILURE, EE_ARGS(1),
                           STRING_ARG, "amppwgrad:gzrf2rcfh"); 
                return FAILURE;
            }

            if (amppwgrad(cfh_crusher_area, target_cfh_crusher2, 0.0, 0.0, cfhloggrd.zrt,
                          MIN_PLATEAU_TIME, &a_gzrf3lcfh, &pw_gzrf3lcfha,
                          &pw_gzrf3lcfh, &pw_gzrf3lcfhd) == FAILURE) 
            {
                epic_error(use_ermes, "%s failed in PScveval.", EM_PSD_SUPPORT_FAILURE, EE_ARGS(1),
                           STRING_ARG, "amppwgrad:gzrf3lrcfhd"); 
                return FAILURE;
            }
        }
    } else {
        if (opspf == 0) 
        {
            /* X FOV Selective */ /* YMSmr09211  04/26/2006 YI */
            a_gxrf2cfh = 4*cyc_rf2cfh/(GAM*(float)pw_rf2cfh/(float)(1.0s)*cfh_fov);
            if (optramp(&pw_gxrf2cfha, a_gxrf2cfh, cfhloggrd.tx, cfhloggrd.xrt, TYPDEF) == FAILURE) 
            {
                epic_error(use_ermes, "%s failed for gxrf2cfh.", 
                           EM_PSD_SUPPORT_FAILURE, EE_ARGS(1),
                           STRING_ARG, "optramp for gxrf2cfh");
                return FAILURE;
            }
            pw_gxrf2cfha = IMax(2,(int)(2*GRAD_UPDATE_TIME), pw_gxrf2cfha);
            pw_gxrf2cfhd = pw_gxrf2cfha;
            target_cfh_crusher = cfhloggrd.tz_xz; /* YMSmr09211  04/26/2006 YI */
        } 
        else 
        {
            /* Y FOV Selective */ /* YMSmr09211  04/26/2006 YI */
            a_gyrf2cfh = 4*cyc_rf2cfh/(GAM*(float)pw_rf2cfh/(float)(1.0s)*cfh_fov);
            if (optramp(&pw_gyrf2cfha, a_gyrf2cfh, cfhloggrd.ty, cfhloggrd.yrt, TYPDEF) == FAILURE) 
            {
                epic_error(use_ermes, "%s failed for gyrf2cfh.", EM_PSD_SUPPORT_FAILURE, EE_ARGS(1),
                           STRING_ARG, "optramp for gyrf2cfh");
                return FAILURE;
            }
            pw_gyrf2cfha = IMax(2,(int)(2*GRAD_UPDATE_TIME), pw_gyrf2cfha);
            pw_gyrf2cfhd = pw_gyrf2cfha;
            target_cfh_crusher = cfhloggrd.tz_yz; /* YMSmr09211  04/26/2006 YI */
        }

        /* Z CRUSHER CVs */ /* YMSmr09211  04/26/2006 YI */
        if (amppwgrad(cfh_crusher_area, target_cfh_crusher, 0.0, 0.0, cfhloggrd.zrt,
                      MIN_PLATEAU_TIME, &a_gzrf2rcfh, &pw_gzrf2rcfha,
                      &pw_gzrf2rcfh, &pw_gzrf2rcfhd) == FAILURE) 
        {
            epic_error(use_ermes, "%s failed in PScveval.",
                       EM_PSD_SUPPORT_FAILURE, EE_ARGS(1),
                       STRING_ARG, "amppwgrad:gzrf2rcfh"); 
            return FAILURE;
        }

        pw_gzrf2lcfha = pw_gzrf2rcfha;
        pw_gzrf2lcfhd = pw_gzrf2rcfhd;
    }

    area_gz1cfh = (PSoff90 + pw_rf1cfh/2.0 + pw_gzrf1cfhd/2.0)*a_gzrf1cfh;
    /* YMSmr09211  04/26/2006 YI */
    if (amppwlcrsh(&psd_cfhleftcrush, &psd_cfhrightcrush,
                   area_gz1cfh, (float)0, cfhloggrd.tz_xz, 
                   MIN_PLATEAU_TIME, cfhloggrd.zrt, &dummy_pw) == FAILURE) 
    {
        epic_error(use_ermes, "%s failed in cfh.", EM_PSD_SUPPORT_FAILURE,
                   EE_ARGS(1), STRING_ARG, "amppwlcrsh for gzrf2lcfh");
        return FAILURE;
    }

    if( presscfh_ctrl != PRESSCFH_NONE && cfh_steam_flag == PSD_ON )
    {
        FLOAT area_g1 = (pw_gxrf2cfh + pw_gxrf2cfha)*a_gxrf2cfh/2.0; 
        if (amppwgrad(area_g1, cfhloggrd.tx, 0.0, 0.0, cfhloggrd.xrt,
                      MIN_PLATEAU_TIME, &a_gx1cfh, &pw_gx1cfha,
                      &pw_gx1cfh, &pw_gx1cfhd) == FAILURE) 
        {
            epic_error(use_ermes, "%s failed in PScveval.",
                       EM_PSD_SUPPORT_FAILURE, EE_ARGS(1),
                       STRING_ARG, "amppwgrad:gx1cfh"); 
            return FAILURE;
        }

        area_g1 = (pw_gyrf3cfh + pw_gyrf3cfhd)*a_gyrf3cfh/2.0; 
        if (amppwgrad(area_g1, cfhloggrd.ty, 0.0, 0.0, cfhloggrd.yrt,
                      MIN_PLATEAU_TIME, &a_gy1cfh, &pw_gy1cfha,
                      &pw_gy1cfh, &pw_gy1cfhd) == FAILURE) 
        {
            epic_error(use_ermes, "%s failed in PScveval.",
                       EM_PSD_SUPPORT_FAILURE, EE_ARGS(1),
                       STRING_ARG, "amppwgrad:gy1cfh"); 
            return FAILURE;
        }
    }

    /* Find Params for killer pulse */ /* YMSmr09211  04/26/2006 YI */
    area_gykcfh = amp_killer*pw_killer;
    if (amppwgrad(area_gykcfh, cfhloggrd.ty, 0.0, 0.0, cfhloggrd.yrt,
                  MIN_PLATEAU_TIME, &a_gykcfh, &pw_gykcfha,
                  &pw_gykcfh, &pw_gykcfhd) == FAILURE) 
    {
        epic_error(use_ermes, "%s failed in cfh.", EM_PSD_SUPPORT_FAILURE,
                   EE_ARGS(1), STRING_ARG, "amppwgrad:gykcfh");
        return FAILURE;
    }

    return SUCCESS;
}

/*
 *  PScveval
 *  
 *  Type: Public Function
 *  
 *  Description:
 *  
 */
STATUS
PScveval( void )
{
    STATUS ps_status;
    FLOAT opthickPS;
    INT i,j; /* YMSmr09211  04/26/2006 YI */
    INT num_cfhlocs;

    /*********************************************************************
     * Generic SECTION
     *********************************************************************/
   
    /* Set the modes for presscfh */
    set_presscfh_mode();

    /* YMSmr09211  04/26/2006 YI */
    if( presscfh_ctrl == PRESSCFH_SHIMVOL )
    {
        num_cfhlocs = exist(oppscvquant);
        for (i=0; i< num_cfhlocs; i++)
        {
            for (j=0; j<9; j++)
            {
                cfh_info[i].oprot[j] = psc_info[i].oppscrot[j];
            }
        }
    }
    else
    {
        num_cfhlocs = IMax(2, 1, exist(oppscvquant));
        for (i=0; i< num_cfhlocs; i++)
        {
            for (j=0; j<9; j++)
            {
                cfh_info[i].oprot[j] = scan_info[PSslice_num].oprot[j];
            }
        }
    }

    cfh_newgeo = 1;
    if (obloptimize(&cfhloggrd, &phygrd, cfh_info, num_cfhlocs,
                    PSD_OBL, 0, obl_method, cfhobl_debug, &cfh_newgeo, cfsrmode)==FAILURE)
    {
        epic_error(use_ermes,"obloptimize failed in PScveval()",
                   EM_PSD_FUNCTION_FAILURE,2,STRING_ARG,"obloptimize",STRING_ARG,"PScveval()");
        return FAILURE; 

    }

    /* derate SR for quiet PSC */ 
    sr_derate(&cfhloggrd, PSsr_derate_factor);

    if (opimode == PSD_3D)  
    {
        opthickPS = 10.0;
    } 
    else 
    {
        /* vmx 12/07/94 YI */
        if( (cfcoilshld == 1) && (pigradcoil == 1) )
        {
            opthickPS = (exist(opslthick) < 3.5) ? 3.5 : exist(opslthick);
        }
        else
        {
            opthickPS = (exist(opslthick) < 3.0) ? 3.0 : exist(opslthick);
        }
        /* end vmx */
    }

    /***********************************************************************
     * MPS1/APS1 SECTION
     ***********************************************************************/

    if ( SUCCESS != (ps_status = PS1cveval( &opthickPS )) ) {
        epic_error(use_ermes, "%s failed.", EM_PSD_SUPPORT_FAILURE,
                   EE_ARGS(1), STRING_ARG, "PS1cveval");
        return ps_status;
    }

    /***********************************************************************
     * CFL SECTION
     ***********************************************************************/

    if ( SUCCESS != (ps_status = CFLcveval( opthickPS )) ) {
        epic_error(use_ermes, "%s failed.", EM_PSD_SUPPORT_FAILURE,
                   EE_ARGS(1), STRING_ARG, "CFLcveval");
        return ps_status;
    }

    /***********************************************************************
     * RCVN SECTION
     ***********************************************************************/

    if ( SUCCESS != (ps_status = RCVNcveval( )) ) {
        epic_error(use_ermes, "%s failed.", EM_PSD_SUPPORT_FAILURE,
                   EE_ARGS(1), STRING_ARG, "RCVNcveval");
        return ps_status;
    }

    /***********************************************************************
     * CFH SECTION
     ***********************************************************************/

    if ( SUCCESS != (ps_status = CFHcveval( opthickPS )) ) {
        epic_error(use_ermes, "%s failed.", EM_PSD_SUPPORT_FAILURE,
                   EE_ARGS(1), STRING_ARG, "CFHcveval");
        return ps_status;
    }

    /***********************************************************************
     * Programmable Prescan SECTION
     ***********************************************************************/

    /* Check for valid Optprescan entry points */
    if(pimrsapsflg) {
        int ipi;
        char strmrs[16] = "pimrsaps1";

        for (ipi = 0; ipi < PSC_MAX_CONTROL_ARRAY; ipi++) 
        {
            switch(*pimrs[ipi]) 
            {
            case MRSAPS_OFF: case MRSAPS_CFL: case MRSAPS_TG:
            case MRSAPS_CFH: case MRSAPS_TR: case MRSAPS_FSEPS:
            case MRSAPS_APA: case MRSAPS_AWS: case MRSAPS_AVS:
            case MRSAPS_AS: case MRSAPS_FTG: case MRSAPS_RCVN:
            case MRSAPS_XTG:
                break;
            default:
                sprintf(strmrs, "pimrsaps%d", ipi + 1);
                epic_error(use_ermes, "%s is out of range.",
                           EM_PSD_CV_OUT_OF_RANGE, EE_ARGS(1),
                           STRING_ARG, strmrs);
                return FAILURE;
            }
        }
    } /* end pimrsapsflg check */

    return SUCCESS;
}   /* end PScveval() */


/*
 *  FTGcveval
 *  
 *  Type: Public Function
 *  
 *  Description:
 *  
 */
STATUS
FTGcveval( void )
{
    INT bw_rf1ftg;
    INT bw_rf2ftg;
    INT bw_rf3ftg;
    FLOAT area_g1ftg;
    FLOAT area_g2ftg;
    FLOAT area_g2bftg;
    FLOAT area_g3ftg;
    FLOAT area_postgxw1ftg;
    FLOAT area_gxw1ftg;
    FLOAT area_gx1ftg;
    FLOAT area_gx2ftg;
    INT PosReadoutWindow;  /* Readout window location            */
    int ftg_xrt;
    float ftg_tx;
    float ftg_tx_xz;
    float ftg_tz_xz;
    float ftg_tx_xyz;
    int FTGtau1check1, FTGtau1check2, FTGtau1check3,
        FTGtau1check4, FTGtau1check5, FTGtau1check6, temp_FTGtau1;

    /* HCSDM00161809: for calculation of minimum FTGtau1 */
    FTGtau1 = 8192us;
    FTGtau1check1   = FTGtau1;
    FTGtau1check2   = FTGtau1;
    FTGtau1check3   = FTGtau1;
    FTGtau1check4   = FTGtau1;
    FTGtau1check5   = FTGtau1;
    FTGtau1check6   = FTGtau1;
    temp_FTGtau1    = FTGtau1;
    
    ftg_xrt = (TGspf ? ps1loggrd.yrt : ps1loggrd.xrt);
    ftg_tx = (TGspf ? ps1loggrd.ty : ps1loggrd.tx);
    ftg_tx_xz = (TGspf ? ps1loggrd.ty_yz : ps1loggrd.tx_xz);
    ftg_tz_xz = (TGspf ? ps1loggrd.tz_yz : ps1loggrd.tz_xz);
    ftg_tx_xyz = (TGspf ? ps1loggrd.ty_xyz : ps1loggrd.tx_xyz);

    FTGfov = cfsystemmaxfov;

    FTGopslthickz1 = 4*FTGslthk;
    FTGopslthickz2 = 4*FTGslthk;
    FTGopslthickz3 = FTGslthk;
    gscale_rf1ftg  = 0.90909;
    gscale_rf2ftg  = gscale_rf1ftg;
    gscale_rf3ftg  = gscale_rf1ftg;
    pw_gzrf1ftg = pw_rf1ftg;
    pw_gzrf2ftg = pw_rf2ftg;
    pw_gzrf3ftg = pw_rf3ftg;

    /* 12/07/94 YI temporary change to avoid internal error with Vectra SGC coil */
    if( cfgradcoil == GCOIL_VECTRA )
    {
        FTGopslthickz3 = 7;
    }
    
    bw_rf1ftg = 4 * cyc_rf1ftg / ((FLOAT)pw_rf1ftg / (FLOAT)1.0s);

    if( !existcv(FTGtau2) )
    {
        FTGtau2 = exist(FTGtau1)*exist(FTGau);
    }

    if( ampslice(&a_gzrf1ftg, bw_rf1ftg, FTGopslthickz1, gscale_rf1ftg, TYPDEF)
        == FAILURE )
    {
        epic_error(use_ermes, "%s failed", EM_PSD_SUPPORT_FAILURE,
                   EE_ARGS(1), STRING_ARG, "ampslice for gzrf1ftg.");
        return FAILURE;
    }
    
    bw_rf2ftg = 4 * cyc_rf2ftg / ((FLOAT)pw_rf2ftg / (FLOAT)1.0s);
 
    if( ampslice(&a_gzrf2ftg, bw_rf2ftg, FTGopslthickz2, gscale_rf2ftg, TYPDEF)
        == FAILURE )
    {
        epic_error(use_ermes, "%s failed", EM_PSD_SUPPORT_FAILURE,
                   EE_ARGS(1), STRING_ARG, "ampslice for gzrf2ftg.");
        return FAILURE;
    }
    
 
    bw_rf3ftg = 4 * cyc_rf3ftg / ((FLOAT)pw_rf3ftg / (FLOAT)1.0s);
 
    if( ampslice(&a_gzrf3ftg, bw_rf3ftg, FTGopslthickz3, gscale_rf3ftg, TYPDEF)
        == FAILURE ) {
        epic_error(use_ermes, "%s failed", EM_PSD_SUPPORT_FAILURE,
                   EE_ARGS(1), STRING_ARG, "ampslice for gzrf3ftg.");
        return FAILURE;
    }
    if( optramp(&pw_gzrf1ftga, a_gzrf1ftg, ps1loggrd.tz, ps1loggrd.zrt, TYPDEF)
        == FAILURE ) {
        epic_error(use_ermes, "%s failed", EM_PSD_SUPPORT_FAILURE,
                   EE_ARGS(1), STRING_ARG, "optramp for gzrf1ftga.");
        return FAILURE;
    }    
    
    pw_gzrf1ftgd = pw_gzrf1ftga;

    if (optramp(&pw_gzrf2ftga, a_gzrf2ftg, ps1loggrd.tz, ps1loggrd.zrt, TYPDEF)
        == FAILURE) {
        epic_error(use_ermes, "%s failed", EM_PSD_SUPPORT_FAILURE,
                   EE_ARGS(1), STRING_ARG, "optramp for gzrf2ftga.");
        return FAILURE;
    }    
    
    pw_gzrf2ftgd = pw_gzrf2ftga;
 
 
    /* Find Params for first zgrad refocusing pulse */
    area_g1ftg = ( a_gzrf1ftg *.5* ( pw_gzrf1ftg + pw_gzrf1ftgd)
                   + a_gzrf2ftg *.5 * (pw_gzrf2ftga + pw_gzrf2ftg) );
    if( amppwgz1(&a_gz1ftg, &pw_gz1ftg, &pw_gz1ftga, &pw_gz1ftgd,
                 area_g1ftg, (INT)1s, MIN_PLATEAU_TIME,
                 ps1loggrd.zrt, ps1loggrd.tz) == FAILURE ) {
        epic_error(use_ermes, "%s failed", EM_PSD_SUPPORT_FAILURE,
                   EE_ARGS(1), STRING_ARG, "amppwgz1 for gz1ftg.");
        return FAILURE;
    }
    

    /* Find Params for refocusing pulse gz2ftg */
    area_g2ftg =  a_gzrf2ftg *.5* ( pw_gzrf2ftg + pw_gzrf2ftgd);
    if( amppwgz1(&a_gz2ftg, &pw_gz2ftg, &pw_gz2ftga, &pw_gz2ftgd,
                 area_g2ftg, (INT)1s, MIN_PLATEAU_TIME,
                 ps1loggrd.zrt, ftg_tz_xz) == FAILURE ) {
	epic_error(use_ermes, "%s failed", EM_PSD_SUPPORT_FAILURE,
                   EE_ARGS(1), STRING_ARG, "amppwgz1 for gz2ftg.");
	return FAILURE;
    }
 
    /* Find Params for refocusing pulse gz2btg */

    /* MRIge56170  AF  10/13/99 */
    if( optramp(&pw_gzrf3ftga, a_gzrf3ftg, ps1loggrd.tz, ps1loggrd.zrt, TYPDEF)
        == FAILURE ) {
        epic_error(use_ermes, "%s failed", EM_PSD_SUPPORT_FAILURE,
                   EE_ARGS(1), STRING_ARG, "optramp for gzrf3ftga.");
        return FAILURE;
    }

    area_g2bftg =  a_gzrf3ftg * .5 *(pw_gzrf3ftga + pw_gzrf3ftg);
    if( amppwgz1(&a_gz2bftg, &pw_gz2bftg, &pw_gz2bftga, &pw_gz2bftgd,
                 area_g2bftg, (INT)1s, MIN_PLATEAU_TIME,
                 ps1loggrd.zrt, ps1loggrd.tz) == FAILURE ) {
        epic_error(use_ermes, "%s failed", EM_PSD_SUPPORT_FAILURE,
                   EE_ARGS(1), STRING_ARG, ":amppwgz1 for gz2bftg.");
        return FAILURE;
    }
    
    pw_gzrf3ftgd = pw_gzrf3ftga;
 
    /* Find Params for refocusing pulse */
    area_g3ftg =  a_gzrf3ftg *.5* ( pw_gzrf3ftg + pw_gzrf3ftgd);
    if( amppwgz1(&a_gz3ftg, &pw_gz3ftg, &pw_gz3ftga, &pw_gz3ftgd,
                 area_g3ftg, (INT)1s, MIN_PLATEAU_TIME,
                 ps1loggrd.zrt, ftg_tz_xz) == FAILURE )
    {

        epic_error(use_ermes, "%s failed", EM_PSD_SUPPORT_FAILURE,
                   EE_ARGS(1), STRING_ARG, "amppwgz1 for gz3ftg.");
        return FAILURE;
    }
    
    if( calcfilter( &echo1ftg_filt, FTGecho1bw, FTGxres, OVERWRITE_NONE)
        == FAILURE ) {
        epic_error(use_ermes, "%s failed.", EM_PSD_SUPPORT_FAILURE,
                   EE_ARGS(1),STRING_ARG,"calcfilter for echo1ftg_filt");
        return FAILURE;
    }

    if( ampfov(&a_gxw1ftg, echo1ftg_filt.bw, FTGfov) == FAILURE )
    {
        epic_error(use_ermes, "%s failed", EM_PSD_SUPPORT_FAILURE,
                   EE_ARGS(1), STRING_ARG, "ampfov for gxw1ftg.");
        return FAILURE;
 
    }
 
    pw_gxw1ftg = RUP_GRD((int)(echo1ftg_filt.tdaq)/4);

    area_gxw1ftg = 2.0*a_gxw1ftg*(float)(pw_gxw1ftg);
 
    if( amppwgx1(&a_gx1ftg, &pw_gx1ftg, &pw_gx1ftga, &pw_gx1ftgd,
                 TYPSPIN, area_gxw1ftg, .0,
                 (int)1s, 1.0, MIN_PLATEAU_TIME,
                 ftg_xrt, ftg_tx_xz) == FAILURE )
    {
        epic_error(use_ermes, "%s failed", EM_PSD_SUPPORT_FAILURE,
                   EE_ARGS(1), STRING_ARG, "amppwgx1 for gx1ftg.");
        return FAILURE;
    }
    
    if( amppwgx1(&a_gx2test, &pw_gx2test, &pw_gx2testa, &pw_gx2testd,
                 TYPSPIN, area_gxw1ftg, .0,
                 (int)1s, 1.0, MIN_PLATEAU_TIME,
                 ftg_xrt, ftg_tx_xz) == FAILURE )
    {
        epic_error(use_ermes, "%s failed", EM_PSD_SUPPORT_FAILURE,
                   EE_ARGS(1), STRING_ARG, "amppwgx1 for gx1ftg.");

        return FAILURE;
    }
    
 
    area_gxw1ftg = a_gxw1ftg*(float)(pw_gxw1ftg);
    if( amppwgx1(&a_gx1bftg, &pw_gx1bftg, &pw_gx1bftga, &pw_gx1bftgd,
                 TYPSPIN, area_gxw1ftg, .0,
                 (int)1s, 1.0, MIN_PLATEAU_TIME,
                 ftg_xrt, ftg_tx_xyz) == FAILURE )
    {
        epic_error(use_ermes, "%s call failed", EM_PSD_SUPPORT_FAILURE,
                   EE_ARGS(1), STRING_ARG, "amppwgx1:gx1btg");
    }

    if( optramp(&pw_gxw1ftga, a_gxw1ftg, ftg_tx,
                ftg_xrt, TYPDEF) == FAILURE )
    {
        epic_error(use_ermes, "%s failed", EM_PSD_SUPPORT_FAILURE,
                   EE_ARGS(1), STRING_ARG, "optramp for gxw1ftga.");
        return FAILURE;
    }    
    
    pw_gxw1ftgd = pw_gxw1ftga;
 
    area_postgxw1ftg = ( 3.0 * a_gxw1ftg * (float)(pw_gxw1ftg)
                         - 0.5 * a_gxw1ftg * (pw_gxw1ftga + pw_gxw1ftgd) );
    if( amppwgx1(&a_postgxw1ftg, &pw_postgxw1ftg, &pw_postgxw1ftga,
                 &pw_postgxw1ftgd, TYPSPIN, area_postgxw1ftg, .0,
                 (int)1s, 1.0, MIN_PLATEAU_TIME,
                 ftg_xrt, ftg_tx) == FAILURE )
    {
        epic_error(use_ermes, "%s failed", EM_PSD_SUPPORT_FAILURE,
                   EE_ARGS(1), STRING_ARG, "amppwgx1 for postgxw1ftg.");
        return FAILURE;
    }
    
    if( ampfov(&a_gxw2ftg, echo1ftg_filt.bw, FTGfov) == FAILURE )
    {
        epic_error(use_ermes, "%s failed", EM_PSD_SUPPORT_FAILURE,
                   EE_ARGS(1), STRING_ARG, "ampfov for gxw2ftg.");
        return FAILURE;
    }
    
    pw_gxw2ftg = RUP_GRD((int)(4.0*(float)pw_gxw1ftg));
    if( optramp(&pw_gxw2ftga, a_gxw2ftg, ftg_tx,
                ftg_xrt, TYPDEF) == FAILURE )
    { 
        epic_error(use_ermes, "%s failed", EM_PSD_SUPPORT_FAILURE,
                   EE_ARGS(1), STRING_ARG, "optramp for gxw2ftga.");
 
        return FAILURE;
    }
    
    pw_gxw2ftgd = pw_gxw2ftga;
 
    /* We want S1 to refocus tau1 ms after rf3.  The gx2tg refocusing
       pulse accomplishes this.  gxw2tg pulse can't start before the
       end of the gz3tg refocussing pulse (why bother collecting
       corrupt data? ). */
    area_gx1ftg = a_gx1ftg * (0.5 * (float)(pw_gx1ftga + pw_gx1ftgd) + pw_gx1ftg);
 
    /* Starting readout window position relative to center of rf3 */
    PosReadoutWindow = ( pw_rf3ftg/2 + pw_gzrf3ftgd + pw_gz3ftga
                         + pw_gz3ftg + pw_gz3ftgd );
 
    pw_gxw2ftgleft = RUP_GRD(echo1ftg_filt.tdaq/8);   /* 1/8th of readout window to left of S1 */

    area_gx2ftg = area_gx1ftg - a_gxw2ftg * (float)(pw_gxw2ftga / 2 + pw_gxw2ftgleft);
 
    if( amppwgrad(area_gx2ftg, ftg_tx_xz,
                  0.0, 0.0, ftg_xrt,
                  MIN_PLATEAU_TIME, &a_gx2ftg, &pw_gx2ftga,
                  &pw_gx2ftg, &pw_gx2ftgd) == FAILURE ) {
        epic_error(use_ermes, "%s failed in fasttg.", EM_PSD_SUPPORT_FAILURE,
                   EE_ARGS(1), STRING_ARG, "amppwgrad:gx2ftg");
        return FAILURE;
    }
   
    if( (pw_rf3ftg/2+pw_gx2ftga+pw_gx2ftg+pw_gx2ftgd+pw_gxw2ftga+pw_gxw2ftgleft
         > FTGtau1) || (pw_rf3ftg/2+pw_gzrf3ftgd+pw_gz3ftga+pw_gz3ftg+pw_gz3ftgd+pw_gxw2ftga+pw_gxw2ftgleft > FTGtau1) )
    {
        /* return FAILURE; */
    }
    
    /* HCSDM00161809: calculate minimum FTGtau1 */
    FTGtau1check1 = RUP_GRD(pw_rf1ftg/2+pw_gx1ftga+pw_gx1ftg+pw_gx1ftgd+pw_rf2ftg/2);
    FTGtau1check2 = RUP_GRD(pw_rf1ftg/2+pw_gzrf1ftgd+pw_gz1ftga+pw_gz1ftg+pw_gz1ftgd+pw_gzrf2ftga+pw_rf2ftg/2);
    FTGtau1check3 = RUP_GRD(pw_rf2ftg/2+pw_gx1bftga+pw_gx1bftg+pw_gx1bftgd+pw_gxw1ftga+pw_gxw1ftg/2);
    FTGtau1check4 = RUP_GRD(pw_rf2ftg/2+pw_gzrf2ftgd+pw_gz2ftga+pw_gz2ftg+pw_gz2ftgd+pw_gxw1ftga+pw_gxw1ftg/2);
    FTGtau1check5 = RUP_GRD(pw_rf3ftg/2+pw_gx2ftga+pw_gx2ftg+pw_gx2ftgd+pw_gxw2ftga+pw_gxw2ftgleft);
    FTGtau1check6 = RUP_GRD(pw_rf3ftg/2+pw_gzrf3ftgd+pw_gz3ftga+pw_gz3ftg+pw_gz3ftgd+pw_gxw2ftga+pw_gxw2ftgleft);

    temp_FTGtau1 = IMax(7, FTGtau1, FTGtau1check1, FTGtau1check2, FTGtau1check3, FTGtau1check4, FTGtau1check5, FTGtau1check6);

    cvoverride(FTGtau1, temp_FTGtau1, PSD_FIX_OFF, PSD_EXIST_ON);

    if( !existcv(FTGtau2) )
    {
        FTGtau2 = exist(FTGtau1)*exist(FTGau); 
    }

    flip_rf1ftg =  90.0;
    flip_rf2ftg = 180.0;
    flip_rf3ftg = 180.0;
    
    strcpy(entry_point_table[L_FTG].epname, "fasttg");
 
    return SUCCESS;
}   /* end FTGcveval() */


/*
 *  XTGcveval
 *  
 *  Type: Public Function
 *  
 *  Description:
 *  
 */
STATUS
XTGcveval( void )
{
    INT bw_rf1xtg;
    INT bw_rf2xtg;
    FLOAT area_g1xtg;
    FLOAT area_gxw1xtg;
    INT ta_180, tb_180;
    FLOAT area_gxw1rampxtg;
    INT xtg_xrt, xtg_yrt;
    FLOAT xtg_tx;
    FLOAT xtg_tx_xz;
    FLOAT xtg_tz_xz;
    FLOAT xtg_tx_xyz, xtg_ty_xyz, xtg_tz_xyz;
    FLOAT xtg_temp_float;
    INT temp_max_pw;

    XTGfov = FMin(2, 480.0, cfsystemmaxfov);

    if(0 == getAps1Mod())
    {
        XTGfov = FMin(2, opfov, XTGfov);
        XTGfov = FMax(2, XTGfov, 50.0);
    }

    xtg_xrt = (TGspf ? ps1loggrd.yrt : ps1loggrd.xrt);
    xtg_yrt = (TGspf ? ps1loggrd.xrt : ps1loggrd.yrt);
    xtg_tx = (TGspf ? ps1loggrd.ty : ps1loggrd.tx);
    xtg_tx_xz = (TGspf ? ps1loggrd.ty_yz : ps1loggrd.tx_xz);
    xtg_tz_xz = (TGspf ? ps1loggrd.tz_yz : ps1loggrd.tz_xz);
    xtg_tz_xyz = ps1loggrd.tz_xyz;
    xtg_tx_xyz = (TGspf ? ps1loggrd.ty_xyz : ps1loggrd.tx_xyz);
    xtg_ty_xyz = (TGspf ? ps1loggrd.tx_xyz : ps1loggrd.ty_xyz);

    XTGopslthick = 10.0;

    gscale_rf1xtg  = 0.90909;
    gscale_rf2xtg  = gscale_rf1xtg;
    pw_gzrf1xtg = pw_rf1xtg;
    pw_gzrf2xtg = pw_rf2xtg;

    XTGecho1bw = 15.625;

    bw_rf1xtg = 4 * cyc_rf1xtg / ((FLOAT)pw_rf1xtg / (FLOAT)1.0s);

    if( ampslice(&a_gzrf1xtg, bw_rf1xtg, XTGopslthick, gscale_rf1xtg, TYPDEF)
        == FAILURE )
    {
        epic_error(use_ermes, "%s failed", EM_PSD_SUPPORT_FAILURE,
                   EE_ARGS(1), STRING_ARG, "ampslice for gzrf1xtg.");
        return FAILURE;
    }

    bw_rf2xtg = 4 * cyc_rf2xtg / ((FLOAT)pw_rf2xtg / (FLOAT)1.0s);

    if( ampslice(&a_gzrf2xtg, bw_rf2xtg, XTGopslthick, gscale_rf2xtg, TYPDEF)
        == FAILURE )
    {
        epic_error(use_ermes, "%s failed", EM_PSD_SUPPORT_FAILURE,
                   EE_ARGS(1), STRING_ARG, "ampslice for gzrf2xtg.");
        return FAILURE;
    }

    if( optramp(&pw_gzrf1xtga, a_gzrf1xtg, ps1loggrd.tz, ps1loggrd.zrt, TYPDEF)
        == FAILURE ) {
        epic_error(use_ermes, "%s failed", EM_PSD_SUPPORT_FAILURE,
                   EE_ARGS(1), STRING_ARG, "optramp for gzrf1xtga.");
        return FAILURE;
    }    

    pw_gzrf1xtgd = pw_gzrf1xtga;
    
    if (optramp(&pw_gzrf2xtga, a_gzrf2xtg, ps1loggrd.tz, ps1loggrd.zrt, TYPDEF)
        == FAILURE) {
        epic_error(use_ermes, "%s failed", EM_PSD_SUPPORT_FAILURE,
                   EE_ARGS(1), STRING_ARG, "optramp for gzrf2xtga.");
        return FAILURE;
    }    

    pw_gzrf2xtgd = pw_gzrf2xtga;

    /* Z CRUSHER CVs */ 
    area_g1xtg = (PSoff90 + pw_gzrf1xtg/2.0 + pw_gzrf1xtgd/2.0)*a_gzrf1xtg;
    area_xtgzkiller = amp_killer*pw_killer+area_g1xtg;
    if (amppwgrad(area_xtgzkiller, xtg_tz_xyz, 0.0, 0.0, ps1loggrd.zrt,
                  MIN_PLATEAU_TIME, &a_gz2xtg, &pw_gz2xtga,
                  &pw_gz2xtg, &pw_gz2xtgd) == FAILURE) {
        epic_error(use_ermes, "%s failed in XTGcveval.",
                   EM_PSD_SUPPORT_FAILURE, EE_ARGS(1),
                   STRING_ARG, "amppwgrad: gz2xtg"); 
        return FAILURE;
    }

    if (amppwgrad((area_xtgzkiller-area_g1xtg), xtg_tz_xyz, 0.0, 0.0, ps1loggrd.zrt,
                  MIN_PLATEAU_TIME, &a_gz1xtg, &pw_gz1xtga,
                  &pw_gz1xtg, &pw_gz1xtgd) == FAILURE) {
        epic_error(use_ermes, "%s failed in XTGcveval.",
                   EM_PSD_SUPPORT_FAILURE, EE_ARGS(1),
                   STRING_ARG, "amppwgrad: gz1xtg"); 
        return FAILURE;
    }

    /* use 16kHz, 256 xres */
    if( calcfilter( &echo1xtg_filt, XTGecho1bw, XTGxres, OVERWRITE_NONE)
        == FAILURE ) {
        epic_error(use_ermes, "%s failed.", EM_PSD_SUPPORT_FAILURE,
                   EE_ARGS(1),STRING_ARG,"calcfilter for echo1xtg_filt");
        return FAILURE;
    }
    
    if (ampfov(&xtg_temp_float, echo1xtg_filt.bw, xtg_tx) == FAILURE)
    {
        epic_error(use_ermes, "%s failed", EM_PSD_SUPPORT_FAILURE,
                   EE_ARGS(1), STRING_ARG, "xtg ampfov");
        return FAILURE;
    }

    xtg_temp_float = ceil(xtg_temp_float / 10.0) * 10.0;
    if( xtg_temp_float > XTGfov )
    {
        XTGfov = xtg_temp_float;
    }

    if( ampfov(&a_gxw1xtg, echo1xtg_filt.bw, XTGfov) == FAILURE )
    {
        epic_error(use_ermes, "%s failed", EM_PSD_SUPPORT_FAILURE,
                   EE_ARGS(1), STRING_ARG, "ampfov for gxw1xtg.");
        return FAILURE;        
    }

    if( optramp(&pw_gxw1xtga, a_gxw1xtg, xtg_tx,
                xtg_xrt, TYPDEF) == FAILURE )
    {
        epic_error(use_ermes, "%s failed", EM_PSD_SUPPORT_FAILURE,
                   EE_ARGS(1), STRING_ARG, "optramp for gxw1xtga.");
        return FAILURE;
    }    

    pw_gxw1xtgd = pw_gxw1xtga;
    
    pw_gxw1xtg = RUP_GRD((int)(echo1xtg_filt.tdaq));
    
    area_gxw1xtg = a_gxw1xtg*(float)(pw_gxw1xtg);
    area_gxw1rampxtg = 0.5*a_gxw1xtg*(float)(pw_gxw1xtga);
    
    if( amppwgx1(&a_gx1bxtg, &pw_gx1bxtg, &pw_gx1bxtga, &pw_gx1bxtgd,
                 TYPSPIN, area_gxw1xtg, area_gxw1rampxtg,
                 10000, 1.0, MIN_PLATEAU_TIME, xtg_xrt, xtg_tx_xyz) == FAILURE )
    {
        epic_error(use_ermes, "%s failed", EM_PSD_SUPPORT_FAILURE,
                   EE_ARGS(1), STRING_ARG, "amppwgx1 for gx1bxtg.");
        return FAILURE;
    }
    a_gx1bxtg = -a_gx1bxtg;

    /* Find Params for killer pulse */
    area_xtgykiller = amp_killer*pw_killer;
    if (amppwgrad(area_xtgykiller, xtg_ty_xyz, 0.0, 0.0, xtg_yrt,
                  MIN_PLATEAU_TIME, &a_gykxtgl, &pw_gykxtgla,
                  &pw_gykxtgl, &pw_gykxtgld) == FAILURE) {
        epic_error(use_ermes, "%s failed in xtg.", EM_PSD_SUPPORT_FAILURE,
                   EE_ARGS(1), STRING_ARG, "amppwgrad:gykxtgl");
        return FAILURE;
    }
    
    if (amppwgrad(area_xtgykiller, xtg_ty_xyz, 0.0, 0.0, xtg_yrt,
                  MIN_PLATEAU_TIME, &a_gykxtgr, &pw_gykxtgra,
                  &pw_gykxtgr, &pw_gykxtgrd) == FAILURE) {
        epic_error(use_ermes, "%s failed in xtg.", EM_PSD_SUPPORT_FAILURE,
                   EE_ARGS(1), STRING_ARG, "amppwgrad:gykxtgr");
        return FAILURE;
    }
    flip_rf1xtg =  90.0;
    flip_rf2xtg = 180.0;
    
    temp_max_pw = IMax(2,pw_gzrf1xtgd, pw_gykxtgla+pw_gykxtgl+pw_gykxtgld);
    ta_180 = pw_gzrf1xtg/2+pw_gz1xtga+pw_gz1xtg+pw_gz1xtgd+temp_max_pw+pw_rf3xtg+rfupd+pw_gzrf2xtga+pw_gzrf2xtg/2;
    
    temp_max_pw = IMax(2,(pw_gx1bxtga+pw_gx1bxtg+pw_gx1bxtgd+pw_gxw1xtga), pw_gykxtgra+pw_gykxtgr+pw_gykxtgrd);
    tb_180 = pw_gzrf2xtg/2+pw_gzrf2xtgd+pw_gz2xtga+pw_gz2xtg+pw_gz2xtgd+temp_max_pw+pw_rf4xtg+rfupd+pw_gxw1xtg/2;
    XTGtau1 = RUP_GRD(IMax(2, ta_180, tb_180)); 
    
    strcpy(entry_point_table[L_XTG].epname, "expresstg");
    
    return SUCCESS;
}   /* end XTGcveval() */


/*
 *  RGcvinit
 *  
 *  Type: Public Function
 *  
 *  Description:
 *  
 */
STATUS
RGcvinit( void )
{
    return SUCCESS;
}


/*
 *  RGcveval
 *
 *  Type: Public Function
 *
 *  Description:
 *
 */
STATUS
RGcveval( void )
{
    if (PSD_ON == rgfeature_enable) 
    { 
        cvmod( opuser39, RG_CAL_MODE_MEASURED, RG_CAL_MODE_HIGH_FIXED, RG_CAL_MODE_HIGH_FIXED, 
               "Receiver Gain (0:Prescan Measured, 1:Predefined)", 0, "" );

        if ( (PSD_2D == exist(opimode)) &&
             ((PSD_SE == exist(oppseq)) || (PSD_IR == exist(oppseq))) &&
             (exist(opptsize) == 4) &&
             (exist(opslthick) <= 10) ) 
        {
            opuser39 = _opuser39.defval;
            activate_reserved_usercv(39);

            if( existcv(opuser39) && (exist(opuser39) != _opuser39.minval) &&
                (exist(opuser39) != _opuser39.maxval) )
            {
                epic_error(use_ermes, "%s must be 0 or 1", EM_PSD_CV_0_OR_1,
                           EE_ARGS(1), STRING_ARG, "UserCV39");

                return FAILURE;
            }

            oprgcalmode = exist(opuser39);
        }
        else 
        {
            opuser39 = 0;
            deactivate_reserved_usercv(39);

            oprgcalmode = RG_CAL_MODE_MEASURED;
        }
    } 
    else 
    {
        opuser39 = 0;
        deactivate_reserved_usercv(39);

        oprgcalmode = RG_CAL_MODE_MEASURED;
    }

    return SUCCESS;
}   /* end RGcveval() */


/*
 *  PSfilter
 *  
 *  Type: Public Function
 *  
 *  Description:
 *  
 */
STATUS
PSfilter( void )
{
    if (setfilter( &echo1cfl,PRESCAN) == FAILURE) {
        epic_error(use_ermes, "%s failed", EM_PSD_SUPPORT_FAILURE,
                   EE_ARGS(1), STRING_ARG, "setfilter");
        return FAILURE;
    }
    filter_cfl_fid = echo1cfl.fslot;

    if (setfilter( &echo1rcvn,PRESCAN) == FAILURE) {
        epic_error(use_ermes, "%s failed", EM_PSD_SUPPORT_FAILURE,
                   EE_ARGS(1), STRING_ARG, "setfilter");
        return FAILURE;
    }
    filter_rcvn_fid = echo1rcvn.fslot;

    if (setfilter(&echo1cfh, PRESCAN) == FAILURE) {
        epic_error(use_ermes, "%s failed", EM_PSD_SUPPORT_FAILURE,
                   EE_ARGS(1), STRING_ARG, "setfilter");
        return FAILURE;
    }
    filter_cfh_fid = echo1cfh.fslot;


    if (setfilter(&echo1mps1_filt, PRESCAN) == FAILURE) {
        epic_error(use_ermes, "%s failed", EM_PSD_SUPPORT_FAILURE,
                   EE_ARGS(1), STRING_ARG, "setfilter");
        return FAILURE;
    }
    filter_echo1mps1 = echo1mps1_filt.fslot;

    if (setfilter(&echo1ftg_filt, PRESCAN) == FAILURE) {
        epic_error(use_ermes, "%s failed", EM_PSD_SUPPORT_FAILURE,
                   EE_ARGS(1), STRING_ARG, "setfilter");
        return FAILURE;
    }
    filter_echo1ftg  = echo1ftg_filt.fslot; /* 11/24/94 YI */
    filter_echo2ftg =  filter_echo1ftg;

    if (setfilter(&echo1xtg_filt, PRESCAN) == FAILURE) {
        epic_error(use_ermes, "%s failed", EM_PSD_SUPPORT_FAILURE,
                   EE_ARGS(1), STRING_ARG, "setfilter");
        return FAILURE;
    }
    filter_echo1xtg  = echo1xtg_filt.fslot;

    if (setfilter(&echo1as_filt, PRESCAN)  == FAILURE) {
        epic_error(use_ermes, "%s failed", EM_PSD_SUPPORT_FAILURE,
                   EE_ARGS(1), STRING_ARG, "setfilter");
        return FAILURE;
    }
    filter_echo1as = echo1as_filt.fslot;

    return SUCCESS;
}   /* end PSfilter() */


/*
 *  PS1predownload
 *  
 *  Type: Private Function
 *  
 *  Description:
 *  
 */
STATUS
PS1predownload( void )
{
    /* Set xmtaddAPS1 according to maximum B1 and rescale for powermon,
       adding additional (audio) scaling if xmtaddAPS1 is too big.
       Add in coilatten, too. */
    xmtaddAPS1 = -200*log10(maxB1[L_APS1]/maxB1Seq) + getCoilAtten();
    if (xmtaddAPS1 > cfdbmax)
    {
        ps1scale = (float) pow(10.0, (cfdbmax - xmtaddAPS1)/200.0);
        xmtaddAPS1 = cfdbmax;
    }
    else
    {
        ps1scale = 1.0;
    }
  
    if( (B0_30000 == cffield) && (PSD_XRMW_COIL == cfgcoiltype) )
    {
        /* MRIhc57081: Limit TG to coil peak B1 on MR750w */
        calcTGLimit(&tgcap, &tgwindow, maxB1Seq, txCoilInfo[getTxIndex(coilInfo[0])]);
    }
    else
    {
        /* Otherwise use defaults */
        tgcap = _tgcap.defval;
        tgwindow = _tgwindow.defval;
    }

    if (setScale(L_APS1, RF_FREE, rfpulse, maxB1[L_APS1], 
                 ps1scale) == FAILURE)
    {
        epic_error(use_ermes, "%s failed.",
                   EM_PSD_SUPPORT_FAILURE, EE_ARGS(1),
                   STRING_ARG, "setScale ps1");
        return FAILURE;
    }

    ia_rf1mps1 = max_pg_iamp*(*rfpulse[RF1_APS1_SLOT].amp);
    ia_rf2mps1 = max_pg_iamp*(*rfpulse[RF2_APS1_SLOT].amp);

    entry_point_table[L_APS1].epxmtadd = (short) rint((double)xmtaddAPS1);
    /* APS1 & MPS1 */
    strcpy(entry_point_table[L_APS1].epname,"aps1");
    entry_point_table[L_APS1].epfilter=(n8)filter_echo1mps1;

    /* MRIge75651 */
    if( powermon( &entry_point_table[L_APS1],
                  L_APS1,
                  (int)RF_FREE,
                  rfpulse,
                  ps1_tr ) == FAILURE )
    {
        epic_error(use_ermes, "%s failed.", EM_PSD_SUPPORT_FAILURE,
                   EE_ARGS(1), STRING_ARG, "powermon ps1");
        return FAILURE;
    }

    {
        double pmfsrf[NUM_POWERMON_PORT];
        double pmfspw[NUM_POWERMON_PORT];
        double pmfsdc[NUM_POWERMON_PORT];

        /* Check for multiple Tx coils. This routine is optimized only
         * for single Tx coil. */ 
        if (1 != getNumTxCoils(coilInfo, opncoils))
        {
            epic_error( use_ermes,
                        "%s does not support more than one transmit coil",
                        EM_PSD_MULTI_TX_NOT_SUPPORTED, EE_ARGS(1), STRING_ARG,
                        "Prescan");
            return FAILURE;
        }
 
        if( cfpwrmontyp >= PMTYP_UPM ) 
        {
            if( SUCCESS != getActivePowerMonPeakLimits( pmfsrf, pmfspw,
                                                        pmfsdc, (int)cffield))
            {
                epic_error( use_ermes, "Support routine %s failed.",
                            EM_PSD_SUPPORT_FAILURE, EE_ARGS(1),
                            STRING_ARG, "getActivePowerMonPeakLimits" );
                return FAILURE;
            }
        }
        else
        {
            int i = 0;
            for( i = 0; i < NUM_POWERMON_PORT; ++i )
            {
                pmfsrf[i] = PMFULL;
                pmfspw[i] = PMFULL;
                pmfsdc[i] = PMFULL;
            }
        }

        /* for APS1, set amp values to full scale */
        /* Use the first Tx coil */
        switch(activePowerMonChannel(txCoilInfo[getTxIndex(coilInfo[0])]))
        {
            case PMCH_HEAD:
                entry_point_table[L_APS1].epamph = pmfsrf[PMCH_HEAD];
                break;
            case PMCH_BODY:
                entry_point_table[L_APS1].epampb = pmfsrf[PMCH_BODY];
                break;
            case PMCH_SPECTRO:
                entry_point_table[L_APS1].epamps = pmfsrf[PMCH_SPECTRO];
                break;
        }
    }

#ifdef BURST_MODE_SUPPORT
    /* Call minseq() to calculate Burst Mode model inputs.
       Disable writing corner points files for prescan entry points. */
    int gradHeatFile_save = gradHeatFile;
    gradHeatFile = FALSE;
    if ( FAILURE == minseq( &min_seqgrad,
                            gradx, GX_FREE,
                            grady, GY_FREE,
                            gradz, GZ_FREE,
                            &loggrd, idx_seqmps1, tsamp, tmin,
                            use_ermes, seg_debug ) ) {
        epic_error( use_ermes, "%s failed.",
                    EM_PSD_ROUTINE_FAILURE, EE_ARGS(1), STRING_ARG, "minseq" );
        return FAILURE;
    }
    gradHeatFile = gradHeatFile_save;
#endif

    /* If aps1_mod set to 1 and NOT breast R or L coil, use volRec coil to set TG */
    if ( (1 == getAps1Mod()) && (PSD_OFF == ps1_rxcoil) )
    {
        if (volRecCoilInfo[0].hubIndex != coilInfo[0].hubIndex)
        {
            UpdateEntryTabRecCoil(&entry_point_table[L_APS1], 0);
        }
    }

    /* copy APS1 to MPS1 */
    entry_point_table[L_MPS1] = entry_point_table[L_APS1];

    strcpy(entry_point_table[L_MPS1].epname, "mps1");

    /* This is usually equal to the scan entry point.
       Make sure it is continuous for manual prescan */
    entry_point_table[L_MPS1].eppmtable.pmContinuousUpdate = 1;
    entry_point_table[L_MPS2].eppmtable.pmContinuousUpdate = 1;

    return SUCCESS;
}

/*
 *  CFLpredownload
 *  
 *  Type: Private Function
 *  
 *  Description:
 *  
 */
STATUS
CFLpredownload( void )
{
    /* Sample time for cfl */
    pitsp1 = echo1cfl.tdaq/echo1cfl.outputs;

    xmtaddCFL = -200*log10(maxB1[L_CFL]/maxB1Seq) + getCoilAtten();
    if (xmtaddCFL > cfdbmax)
    {
        cflscale = (float) pow(10.0, (cfdbmax - xmtaddCFL)/200.0);
        xmtaddCFL = cfdbmax;
    }
    else
    {
        cflscale = 1.0;
    }

    if (setScale(L_CFL, RF_FREE, rfpulse, maxB1[L_CFL], 
                 cflscale) == FAILURE)
    {
        epic_error(use_ermes, "%s failed.", EM_PSD_SUPPORT_FAILURE,
                   EE_ARGS(1), STRING_ARG, "setScale cfl");
        return FAILURE;
    }

    ia_rf1cfl = max_pg_iamp*(*rfpulse[RF1_CFL_SLOT].amp);

    entry_point_table[L_CFL].epxmtadd = (short) rint((double)xmtaddCFL);

    entry_point_table[L_CFL].epprexres = (s16)CFLxres; /* MRIhc54366 */

    strcpy(entry_point_table[L_CFL].epname,"cfl");
    entry_point_table[L_CFL].epfilter=(n8)filter_cfl_fid;
    
    /* MRIge75651 */
    if( powermon( &entry_point_table[L_CFL],
                  L_CFL,
                  (int)RF_FREE,
                  rfpulse,
                  cfl_tr ) == FAILURE )
    {
        epic_error(use_ermes, "%s failed.", EM_PSD_SUPPORT_FAILURE,
                   EE_ARGS(1), STRING_ARG, "powermon cfl");
        return FAILURE;
    }

#ifdef BURST_MODE_SUPPORT
    /* Call minseq() to calculate Burst Mode model inputs.
       Disable writing corner points files for prescan entry points. */
    int gradHeatFile_save = gradHeatFile;
    gradHeatFile = FALSE;
    if ( FAILURE == minseq( &min_seqgrad,
                            gradx, GX_FREE,
                            grady, GY_FREE,
                            gradz, GZ_FREE,
                            &loggrd, idx_seqcfl, tsamp, tmin,
                            use_ermes, seg_debug ) ) {
        epic_error( use_ermes, "%s failed.",
                    EM_PSD_ROUTINE_FAILURE, EE_ARGS(1), STRING_ARG, "minseq" );
        return FAILURE;
    }
    gradHeatFile = gradHeatFile_save;
#endif

    return SUCCESS;
}

/*
 *  RCVNpredownload
 *  
 *  Type: Private Function
 *  
 *  Description:
 *  
 */
STATUS
RCVNpredownload( void )
{
    entry_point_table[L_RCVN] = entry_point_table[L_MPS2];

    xmtaddRCVN = xmtaddCFL;
    entry_point_table[L_RCVN].epxmtadd = (short) rint((double)xmtaddRCVN);

    strcpy(entry_point_table[L_RCVN].epname,"rcvn");
    
    entry_point_table[L_RCVN].epfilter  = (n8)filter_rcvn_fid;
    entry_point_table[L_RCVN].epprexres = rcvn_xres;

#ifdef BURST_MODE_SUPPORT
    /* Call minseq() to calculate Burst Mode model inputs.
       Disable writing corner points files for prescan entry points. */
    int gradHeatFile_save = gradHeatFile;
    gradHeatFile = FALSE;
    if ( FAILURE == minseq( &min_seqgrad,
                            gradx, GX_FREE,
                            grady, GY_FREE,
                            gradz, GZ_FREE,
                            &loggrd, idx_seqrcvn, tsamp, tmin,
                            use_ermes, seg_debug ) ) {
        epic_error( use_ermes, "%s failed.",
                    EM_PSD_ROUTINE_FAILURE, EE_ARGS(1), STRING_ARG, "minseq" );
        return FAILURE;
    }
    gradHeatFile = gradHeatFile_save;
#endif

    return SUCCESS;
}

/*
 *  CFHpredownload
 *  
 *  Type: Private Function
 *  
 *  Description:
 *  
 */
STATUS
CFHpredownload( void )
{
    /* Sample time for cfh */
    pitsp2 = echo1cfh.tdaq/echo1cfh.outputs;

    xmtaddCFH = -200*log10(maxB1[L_CFH]/maxB1Seq) + getCoilAtten();
    if (xmtaddCFH > cfdbmax)
    {
        cfhscale = (float) pow(10.0, (cfdbmax - xmtaddCFH)/200.0);
        xmtaddCFH = cfdbmax;
    }
    else
    {
        cfhscale = 1.0;
    }

    if (setScale(L_CFH, RF_FREE, rfpulse, maxB1[L_CFH], 
                 cfhscale) == FAILURE)
    {
        epic_error(use_ermes, "%s failed.", EM_PSD_SUPPORT_FAILURE,
                   EE_ARGS(1), STRING_ARG, "setScale cfh");
        return FAILURE;
    }

    ia_rf1cfh = max_pg_iamp*(*rfpulse[RF1_CFH_SLOT].amp);
    ia_rf2cfh = max_pg_iamp*(*rfpulse[RF2_CFH_SLOT].amp);
#ifdef PSD_CFH_CHEMSAT
    if ((cs_sat == PSD_ON) && PScs_sat)
    {
        ia_rfcssatcfh = max_pg_iamp*(*rfpulse[RFCSSAT_CFH_SLOT].amp);
    }
#endif
%ifdef PSD_CFH_MT
    if ((opmt == PSD_ON) && PSmt)
    {
        ia_rfmtcfh = max_pg_iamp*(*rfpulse[RFMT_CFH_SLOT].amp);
    }
%endif /* PSD_CFH_MT */
    if (PSD_ON == PSir)
    {
        ia_rf0cfh = max_pg_iamp*(*rfpulse[RF0_CFH_SLOT].amp);
    }

    cfh_rf1freq = ( GAM * a_gzrf1cfh * PStloc
                    / (10 * TARDIS_FREQ_RES) );
    /* factor 10 is because rsptloc is in mm */

    if( ( (opcoax != 0) && cfh_newmode ) || (presscfh != PRESSCFH_NONE) ) {
        if( presscfh_ctrl == PRESSCFH_SLICE || presscfh_ctrl == PRESSCFH_SLAB ) {
	    cfh_rf2freq = GAM * presscfh_info[0].oppscrloc * a_gxrf2cfh / (10 * TARDIS_FREQ_RES) ;
            cfh_rf3freq =  GAM * presscfh_info[0].oppscphasoff * a_gyrf3cfh / (10 * TARDIS_FREQ_RES) ;
	} else if(presscfh_ctrl == PRESSCFH_SHIMVOL){
	    cfh_rf2freq = GAM * psc_info[0].oppscrloc * a_gxrf2cfh / (10 * TARDIS_FREQ_RES) ;
            cfh_rf3freq =  GAM * psc_info[0].oppscphasoff * a_gyrf3cfh / (10 * TARDIS_FREQ_RES) ;
            if((presscfh != presscfh_ctrl) && presscfh_debug) {
                printf("\n  presscfh %d changes to presscfh_ctrl %d \n",presscfh, presscfh_ctrl);
                fflush(stdout);
            }
        } else {
            cfh_rf2freq = ( GAM * (opspf ? PSphasoff * a_gyrf2cfh : PSrloc * a_gxrf2cfh )
                        / (10 * TARDIS_FREQ_RES) );
        /* factor 10 is because rloc/phasoff is in mm */
    }
    } else {
        /* For non-coaxials (multi angle), center cfh_rf2freq and up the
           excitation region to 40 cm to cover lots of ground 
        */
        cfh_rf2freq = 0;
    }

    entry_point_table[L_CFH].epxmtadd = (short) rint((double)xmtaddCFH);
    strcpy(entry_point_table[L_CFH].epname,"cfh");
    entry_point_table[L_CFH].epfilter=(n8)filter_cfh_fid;
    entry_point_table[L_CFH].epprexres = (s16)echo1ptcfh; /* MRIhc08633 */
    
    /* MRIge75651 */
    if( powermon( &entry_point_table[L_CFH],
                  L_CFH,
                  (int)RF_FREE,
                  rfpulse,
                  cfh_tr) == FAILURE )
    {
        epic_error(use_ermes, "%s failed.", EM_PSD_SUPPORT_FAILURE,
                   EE_ARGS(1), STRING_ARG, "powermon cfh");
        return FAILURE;
    }

#ifdef BURST_MODE_SUPPORT
    /* Call minseq() to calculate Burst Mode model inputs.
       Disable writing corner points files for prescan entry points. */
    int gradHeatFile_save = gradHeatFile;
    gradHeatFile = FALSE;
    if ( FAILURE == minseq( &min_seqgrad,
                            gradx, GX_FREE,
                            grady, GY_FREE,
                            gradz, GZ_FREE,
                            &loggrd, idx_seqcfh, tsamp, tmin,
                            use_ermes, seg_debug ) ) {
        epic_error( use_ermes, "%s failed.",
                    EM_PSD_ROUTINE_FAILURE, EE_ARGS(1), STRING_ARG, "minseq" );
        return FAILURE;
    }
    gradHeatFile = gradHeatFile_save;
#endif

    switch (getTxCoilType())
    /* Only 0.5/0.2T will use. */ /* vmx 07/27/95 YO */
    {
        case TX_COIL_LOCAL:
            cfh_ec_position = (10.0 / 256.0);
            break;
        default:
            cfh_ec_position = (16.0 / 256.0);
            break;
    }

    return SUCCESS;
}

/*
 *  PSpredownload
 *  
 *  Type: Public Function
 *  
 *  Description:
 *  
 */
STATUS
PSpredownload( void )
{
    STATUS ps_status;
    INT i;
    INT j;
    INT index, vidx;
    
    /* Check for multiple Tx coils. This routine is optimized only for
     * single Tx coil. */
    if (1 != getNumTxCoils(coilInfo, opncoils))
    {
        epic_error( use_ermes,
                    "%s does not support more than one transmit coil",
                    EM_PSD_MULTI_TX_NOT_SUPPORTED, EE_ARGS(1), STRING_ARG,
                    "prescan");
        return FAILURE;
    }
    
    /***********************************************************************
     * Generic SECTION
     ***********************************************************************/

    if(PSD_ON == exist(oprtcgate)) {
        phys_record_flag = PSD_ON; /* flag for rt data recording */
    } else {
        phys_record_flag = PSD_OFF;
    }

    /* go through entry point table and set frequency offset based on receiver */
    for( i = 0; i < ENTRY_POINT_MAX; i++ )
    {
        PSfreq_offset[i] = cfreceiveroffsetfreq;
    }

    pw_omegarf0cfh = pw_rf0cfh; /* adiabatic pulse */

    PSslice_ind = PSslice_num;  /* MRIge90312 -- for smart prescan */

    PStloc = scan_info[PSslice_num].optloc;
    PSrloc = scan_info[PSslice_num].oprloc;
    PSphasoff = scan_info[PSslice_num].opphasoff;

    /* begin aps1_mod changes (GE) */
    if ( 1 == getAps1Mod() )
    {
        cvunlock(PStloc_mod);
        cvunlock(PSrloc_mod);
        PStloc_mod = getAps1ModPsTloc();
        PSrloc_mod = getAps1ModPsRloc();
    }
    else  /* use imaging locs */
    { 
        PStloc_mod = PStloc;
        PSrloc_mod = opspf ? PSphasoff : PSrloc;    
    }
    /* end aps1_mod changes (GE) */

   /* Create rsp_psc_info table based on psc_info */
    for ( i=0; i < exist(oppscvquant); i ++) {
        rsp_psc_info[i].rsppsctloc = psc_info[i].oppsctloc;
        rsp_psc_info[i].rsppscrloc = psc_info[i].oppscrloc;
        rsp_psc_info[i].rsppscphasoff = psc_info[i].oppscphasoff;

        rsp_psc_info[i].rsppsclenx = psc_info[i].oppsclenx;
        rsp_psc_info[i].rsppscleny = psc_info[i].oppscleny;
        rsp_psc_info[i].rsppsclenz = psc_info[i].oppsclenz; 
    }

    /* Check the rotation matrix for rsp_psc_info */
    for (i=0; i< exist(oppscvquant); i++) {
        for (j=0; j<9; j++) {
            rsp_psc_info[i].rsppscrot[j] = hostToRspRotMat(psc_info[0].oppscrot[j]);
        }
    } 

    /* fill in the prescan rotation array for the prescan slice.
       PSrot is an ipgexport defined in epic.h  */

    /*
     * MRIge43971 BJM: loop over 2D PSrot array to be consistent with other
     *                 rotation matrices and since scalerotmats() expects a 2D
     *                 argument.
     */
    for (index = 0; index < 9; index++)
    {
        PSrot[0][index] = hostToRspRotMat(scan_info[PSslice_num].oprot[index]);
        PSrot_mod[0][index] = hostToRspRotMat(ps1scan_info[0].oprot[index]); 

        rsp_PSrot[0][index] = hostToRspRotMat(cfh_info[0].oprot[index]);

        /* set up rot for CFH */
        for( vidx = 0; vidx < oppscvquant; vidx++ )
        {
            rsp_PSrot[vidx][index] = hostToRspRotMat(cfh_info[vidx].oprot[index]);

        }
    }

    /* Scale Rot matrix for CFH */
    if(scalerotmats(rsp_PSrot, &cfhloggrd, &phygrd, IMax(2,1,exist(oppscvquant)), obl_debug) == FAILURE) /* YMSmr09211  04/26/2006 YI */
    {
        epic_error(use_ermes,"System configuration data integrity violation detected in PSD.\nPlease try again or restart the system.", 
                   EM_PSD_PSDCRUCIAL_CONFIG_FAILURE,EE_ARGS(0));
        return FAILURE;
    }
    
    /* Scale Rot matrix for CFL */
    if(scalerotmats(PSrot, &loggrd, &phygrd, 1, obl_debug) == FAILURE)
    {
        epic_error(use_ermes,"System configuration data integrity violation detected in PSD.\nPlease try again or restart the system.", 
                   EM_PSD_PSDCRUCIAL_CONFIG_FAILURE,EE_ARGS(0));
        return FAILURE;
    }

    /* Scale Rot matrix for ps1 & FTG */
    if(scalerotmats(PSrot_mod, &ps1loggrd, &phygrd, 1, obl_debug) == FAILURE) 
    {
        epic_error(use_ermes,"System configuration data integrity violation detected in PSD.\nPlease try again or restart the system.", 
                   EM_PSD_PSDCRUCIAL_CONFIG_FAILURE,EE_ARGS(0));
        return FAILURE;
    }

    PStrigger = TRIG_LINE;

    /* For Prescan: Inform 'Auto' Prescan about prescan parameters 	*/
    pitr = 2s;	        /* 1st pass prescan TR 	*/
    pichop = 0;		/* No chop		*/

    /* find minimum rfamp te time based on duty cycle */
    min180te = RUP_GRD((int)((float)(pw_rf1mps1 + cfrfminblank)/
                             ((TX_COIL_LOCAL == getTxCoilType()) ? 
                              cfrfmdch : cfrfmdcb)))*2;

    /***********************************************************************
     * MPS1/APS1 SECTION
     ***********************************************************************/

    if ( SUCCESS != (ps_status = PS1predownload()) ) {
        epic_error(use_ermes, "%s failed.", EM_PSD_SUPPORT_FAILURE,
                   EE_ARGS(1), STRING_ARG, "PS1predownload");
        return ps_status;
    }


    /***********************************************************************
     * CFL SECTION
     ***********************************************************************/

    if ( SUCCESS != (ps_status = CFLpredownload()) ) {
        epic_error(use_ermes, "%s failed.", EM_PSD_SUPPORT_FAILURE,
                   EE_ARGS(1), STRING_ARG, "CFLpredownload");  /* MRIhc08595 */
        return ps_status;
    }

    /***********************************************************************
     * CFH SECTION
     ***********************************************************************/

    if ( SUCCESS != (ps_status = CFHpredownload()) ) {
        epic_error(use_ermes, "%s failed.", EM_PSD_SUPPORT_FAILURE,
                   EE_ARGS(1), STRING_ARG, "CFHpredownload");  /* MRIhc08595 */
        return ps_status;
    }

    /***********************************************************************
     * RCVN SECTION
     ***********************************************************************/

    if ( SUCCESS != (ps_status = RCVNpredownload()) ) {
        epic_error(use_ermes, "%s failed.", EM_PSD_SUPPORT_FAILURE,
                   EE_ARGS(1), STRING_ARG, "RCVNpredownload");
        return ps_status;
    }

    /* xmtaddRef is used for calculating TG value to be stored in smart prescan DB wrt a
     * reference maxB1 used in TG entry point */
    xmtaddRef = -200*log10(cfmaxb1ref/100.0/maxB1Seq) + getCoilAtten();

    /* HCSDM00184619 : Tools PSD dont inline PSpreDwonload.
     * This should move to a routine which  psdIF can handle. */

    /* Copy coilInfo, volRecCoilInfo, txCoilInfo to target side */
    memcpy(coilInfo_tgt, coilInfo, sizeof(coilInfo));
    memcpy(volRecCoilInfo_tgt, volRecCoilInfo, sizeof(volRecCoilInfo));
    memcpy(txCoilInfo_tgt, txCoilInfo, sizeof(txCoilInfo));
    memcpy(cttEntry_tgt, cttEntry, sizeof(cttEntry)); 
    chksum_rampdir_tgt = chksum_rampdir;
    cframpdir_tgt = cframpdir;

    return SUCCESS;
}   /* end PSpredownload() */


/*
 *  FTGpredownload
 *  
 *  Type: Public Function
 *  
 *  Description:
 *  
 */
STATUS
FTGpredownload( void )
{
    
    /* Set xmtaddFTG according to maximum B1 and rescale for powermon,
       adding additional (audio) scaling if xmtaddFTG is too big.
       We are assuming that the pulse shapes used in CFH are the
       same as in scan. */
    xmtaddFTG = -200*log10(maxB1[L_FTG]/maxB1Seq) + getCoilAtten();
    if (xmtaddFTG > cfdbmax)
    {
        ftgscale = (float) pow(10.0, (cfdbmax - xmtaddFTG)/200.0);
        xmtaddFTG = cfdbmax;
    }
    else
    {
        ftgscale = 1.0;
    }

    if (setScale(L_FTG,RF_FREE,rfpulse,maxB1[L_FTG],ftgscale) == FAILURE)
    {
        epic_error(use_ermes, "%s failed.", EM_PSD_SUPPORT_FAILURE,
                   EE_ARGS(1), STRING_ARG, "setScale ftg");
        return FAILURE;
    }
    
    /* Set the amplitude scale factors. */
    ia_rf1ftg = max_pg_iamp*(*rfpulse[RF1_FTG_SLOT].amp);
    ia_rf2ftg = max_pg_iamp*(*rfpulse[RF2_FTG_SLOT].amp);
    ia_rf3ftg = max_pg_iamp*(*rfpulse[RF3_FTG_SLOT].amp);
    
    entry_point_table[L_FTG].epxmtadd = (short)rint((double)xmtaddFTG);

    /* MRIge75651 */
    if( powermon( &entry_point_table[L_FTG],
                  L_FTG,
                  (int)RF_FREE,
                  rfpulse,
                  (int)exist(ftgtr) ) == FAILURE )
    {
        epic_error(use_ermes, "%s failed.", EM_PSD_SUPPORT_FAILURE,
                   EE_ARGS(1), STRING_ARG, "powermon ftg");
        return FAILURE;
    }
    xmtadd = xmtaddFTG;

    entry_point_table[L_FTG].epfilter = (n8)filter_echo1ftg; /* 11/24/94 YI */
    entry_point_table[L_FTG].epprexres = 256;
    
    FTGxmtadd = entry_point_table[L_APS1].epxmtadd-entry_point_table[L_FTG].epxmtadd;

#ifdef BURST_MODE_SUPPORT
    /* Call minseq() to calculate Burst Mode model inputs.
       Disable writing corner points files for prescan entry points. */
    int gradHeatFile_save = gradHeatFile;
    gradHeatFile = FALSE;
    if ( FAILURE == minseq( &min_seqgrad,
                            gradx, GX_FREE,
                            grady, GY_FREE,
                            gradz, GZ_FREE,
                            &loggrd, idx_seqftg, tsamp, tmin,
                            use_ermes, seg_debug ) ) {
        epic_error( use_ermes, "%s failed.",
                    EM_PSD_ROUTINE_FAILURE, EE_ARGS(1), STRING_ARG, "minseq" );
        return FAILURE;
    }
    gradHeatFile = gradHeatFile_save;
#endif

    /* Use volRec coil for FTG */
    if (volRecCoilInfo[0].hubIndex != coilInfo[0].hubIndex)
    {
        UpdateEntryTabRecCoil(&entry_point_table[L_FTG], 0);
    } 

    return SUCCESS;
}   /* end FTGpredownload() */


/*
 *  XTGpredownload
 *  
 *  Type: Public Function
 *  
 *  Description:
 *  
 */
STATUS
XTGpredownload( void )
{
    
    /* Set xmtaddXTG according to maximum B1 and rescale for powermon,
       adding additional (audio) scaling if xmtaddXTG is too big. */
    xmtaddXTG = -200*log10(maxB1[L_XTG]/maxB1Seq) + getCoilAtten();
    if (xmtaddXTG > cfdbmax)
    {
        xtgscale = (float) pow(10.0, (cfdbmax - xmtaddXTG)/200.0);
        xmtaddXTG = cfdbmax;
    }
    else
    {
        xtgscale = 1.0;
    }

    if (setScale(L_XTG,RF_FREE,rfpulse,maxB1[L_XTG],xtgscale) == FAILURE)
    {
        epic_error(use_ermes, "%s failed.", EM_PSD_SUPPORT_FAILURE,
                   EE_ARGS(1), STRING_ARG, "setScale xtg");
        return FAILURE;
    }
    
    /* Set the amplitude scale factors. */
    ia_rf1xtg = max_pg_iamp*(*rfpulse[RF1_XTG_SLOT].amp);
    ia_rf2xtg = max_pg_iamp*(*rfpulse[RF2_XTG_SLOT].amp);
    ia_rf4xtg = max_pg_iamp*(*rfpulse[RF4_XTG_SLOT].amp);
    ia_rf3xtg = -ia_rf4xtg;
    a_rf3xtg  = -a_rf4xtg;
    
    entry_point_table[L_XTG].epxmtadd = (short)rint((double)xmtaddXTG);

    /* MRIge75651 */
    if( powermon( &entry_point_table[L_XTG],
                  L_XTG,
                  (int)RF_FREE,
                  rfpulse,
                  (int)exist(xtgtr) ) == FAILURE )
    {
        epic_error(use_ermes, "%s failed.", EM_PSD_SUPPORT_FAILURE,
                   EE_ARGS(1), STRING_ARG, "powermon xtg");
        return FAILURE;
    }
    xmtadd = xmtaddXTG;

    entry_point_table[L_XTG].epfilter = (n8)filter_echo1xtg; /* 11/24/94 YI */
    entry_point_table[L_XTG].epprexres = 256;

#ifdef BURST_MODE_SUPPORT
    /* Call minseq() to calculate Burst Mode model inputs.
       Disable writing corner points files for prescan entry points. */
    int gradHeatFile_save = gradHeatFile;
    gradHeatFile = FALSE;
    if ( FAILURE == minseq( &min_seqgrad,
                            gradx, GX_FREE,
                            grady, GY_FREE,
                            gradz, GZ_FREE,
                            &loggrd, idx_seqxtg, tsamp, tmin,
                            use_ermes, seg_debug ) ) {
        epic_error( use_ermes, "%s failed.",
                    EM_PSD_ROUTINE_FAILURE, EE_ARGS(1), STRING_ARG, "minseq" );
        return FAILURE;
    }
    gradHeatFile = gradHeatFile_save;
#endif

    if(PSD_ON == xtg_volRecCoil)
    {
        /* Use volRec coil for XTG */
        if (volRecCoilInfo[0].hubIndex != coilInfo[0].hubIndex)
        {
            UpdateEntryTabRecCoil(&entry_point_table[L_XTG], 0);
        } 
    }

    return SUCCESS;
}   /* end XTGpredownload() */


/*
 *  ASpredownload
 *  
 *  Type: Public Function
 *  
 *  Description:
 *  
 */
STATUS
ASpredownload( void )
{
    FLOAT asscale;
 
    /******************************************************************/
    /* Set xmtaddas according to maximum B1 and rescale for powermon, */
    /* adding additional (audio) scaling if xmtaddas is too big.      */
    /* We are assuming that the pulse shapes used in CFH are the      */
    /* same as in scan.                                               */
    /******************************************************************/

    strcpy(entry_point_table[L_AUTOSHIM].epname, "autoshim");

    xmtaddas = -200*log10(maxB1[L_AUTOSHIM]/maxB1Seq) + getCoilAtten();
    if (xmtaddas > cfdbmax) 
    {
        asscale = (float) pow(10.0, (cfdbmax - xmtaddas)/200.0);
        xmtaddas = cfdbmax;
    } 
    else
    {
        asscale = 1.0;
    }

    if (setScale(L_AUTOSHIM, RF_FREE, rfpulse, maxB1[L_AUTOSHIM],
                 asscale) == FAILURE) 
    {
        epic_error(use_ermes, "%s failed.", EM_PSD_SUPPORT_FAILURE,
                   EE_ARGS(1), STRING_ARG, "setScale autoshim");
        return FAILURE;
    }

    /* Set the amplitude scale factors. */
    ia_rf1as = max_pg_iamp*(*rfpulse[RF1_AUTOSHIM].amp);

    entry_point_table[L_AUTOSHIM].epxmtadd=(short)rint((double)xmtaddas);
    
    /* MRIge75651 */
    if( powermon( &entry_point_table[L_AUTOSHIM],
                  L_AUTOSHIM,
                  (int)RF_FREE,
                  rfpulse,
                  (int)tr_as ) == FAILURE )
    {
        epic_error(use_ermes, "%s failed.", EM_PSD_SUPPORT_FAILURE,
                   EE_ARGS(1), STRING_ARG, "powermon autoshim");
        return FAILURE;
    }

#ifdef BURST_MODE_SUPPORT
    /* Call minseq() to calculate Burst Mode model inputs.
       Disable writing corner points files for prescan entry points. */
    int gradHeatFile_save = gradHeatFile;
    gradHeatFile = FALSE;
    if ( FAILURE == minseq( &min_seqgrad,
                            gradx, GX_FREE,
                            grady, GY_FREE,
                            gradz, GZ_FREE,
                            &loggrd, idx_seqaushim, tsamp, tmin,
                            use_ermes, seg_debug ) ) {
        epic_error( use_ermes, "%s failed.",
                    EM_PSD_ROUTINE_FAILURE, EE_ARGS(1), STRING_ARG, "minseq" );
        return FAILURE;
    }
    gradHeatFile = gradHeatFile_save;
#endif

    entry_point_table[L_AUTOSHIM].epfilter = (n8)filter_echo1as;
    entry_point_table[L_AUTOSHIM].epprexres = asxres;
    
    /* Use volRec coil for autoshim */
    if (volRecCoilInfo[0].hubIndex != coilInfo[0].hubIndex)
    {
        UpdateEntryTabRecCoil(&entry_point_table[L_AUTOSHIM], 0);
    } 

    pidoshim = PSD_OFF;  /* MRIge92917 - initialize to OFF */
    if( pimrsapsflg == PSD_ON || exist(opepi) == PSD_ON 
        || exist(opspiral) == PSD_ON || oppseq == PSD_SSFP
        || exist(opvrg) == PSD_ON )
    {
        pidoshim = PSD_ON;
    }  /* Auto Shim is required for MRS, EPI, Spiral, VERSE and fiesta */

    if( exist(opfat) || exist(opfatcl)
        || exist(opspecir) || exist(opwater) )
    {
        pidoshim = PSD_ON;
    }  /* Auto Shim is required for Chem Sat scans */

    return SUCCESS;
}   /* end ASpredownload() */

/* CoilSwitchGetTR (MRIhc15304)
 * 
 * Description:
 *   This function returns the TR for the coilSwitch SSP sequence based on
 *   the setRcvPortFlag as passed to CoilSwitchSetCoil.
 *
 *  Parameters:
 *  (I: for input parameters, O: for output parameters)
 *  (O) return: TR in usec
 *  (I) setRcvPortFlag - flag to indicate that setrcvportimm will be
 *        executed when switching coils.  See CoilSwitchSetCoil
 */
int
CoilSwitchGetTR(const int setRcvPortFlag)
{
    int wait_rspimm = 0;
    
    /* When asynchronous RSP function calls are needed, the wait time is
     * extended to accomodate the worst case RSP time */
    if( COIL_SWITCH_RSP_SETHUBINDEXIMM & cfcoilswitchmethod )
    {
        wait_rspimm = csw_wait_sethubindeximm;
    }

    if( setRcvPortFlag )
    {
        wait_rspimm = IMax(2, wait_rspimm, csw_wait_setrcvportimm);
    }

    return csw_tr + wait_rspimm;
}

@pg PSpulsegen
/*********************************************************************
 *                      PRESCAN.E PG SECTION                         *
 *                          PSpulsegen                               *
 *                                                                   *
 * Write here the functional code that loads hardware sequencer      *
 * memory with data that will allow it to play out the sequence.     *
 * These functions call pulse generation macros previously defined   *
 * with @pulsedef, and must return SUCCESS or FAILURE.               *
 *********************************************************************/
PSpulsegen();
FTGpulsegen();
XTGpulsegen();
ASpulsegen();

@pg PSipg
/*********************************************************************
 *                      PRESCAN.E PG SECTION                         *
 *                             PSipg                                 *
 *                                                                   *
 * Write here the functional code that loads hardware sequencer      *
 * memory with data that will allow it to play out the sequence.     *
 * These functions call pulse generation macros previously defined   *
 * with @pulsedef, and must return SUCCESS or FAILURE.               *
 *********************************************************************/

/*
 *  PS1pulsegen
 *  
 *  Type: Private Function
 *  
 *  Description:
 *  
 */
STATUS
PS1pulsegen( INT posstart )
{
    INT postemp;
    INT ta_180, tb_180, te_180;
    INT end_encode;

    /***********************************************************************
     * MPS1/APS1 SECTION
     ***********************************************************************/


    SLICESELZ2( rf1mps1, RF1_APS1_SLOT, posstart, pw_rf1mps1, thk_rf1mps1,
                flip_rf1mps1, cyc_rf1mps1, TYPNDEF, ps1loggrd );

    /* Refocus on Z gradient */
    postemp = RUP_GRD(pend(&gzrf1mps1d,"gzrf1mps1d",0));
    TRAPEZOID(ZGRAD, gz1mps1, postemp+pw_gz1mps1a, 0, TYPNDEF, ps1loggrd);
    end_encode = pend(&gz1mps1d,"gz1mps1d",0);    

    /* read dephaser */
    postemp = RUP_GRD(pend(&gzrf1mps1, "gzrf1mps1", 0)+rfupd+pw_gx1mps1a);
    TRAPEZOID(read_axis, gx1mps1, postemp, 0, TYPNDEF, ps1loggrd);
    if(pend(&gx1mps1d,"gx1mps1d",0)>end_encode)
    {
        end_encode = pend(&gx1mps1d,"gx1mps1d",0);
    }

    /****** figure out minimum te from z grad, compare needed 
      time before and after 180 *********/
    tb_180  =  end_encode - ( RUP_GRD(1ms) + pw_rf1mps1/2 - PSoff90)
        + pw_gzrf2lmps1a + pw_gzrf2lmps1 + pw_gzrf2lmps1d + pw_rf2mps1/2;

    ta_180  = pw_rf2mps1/2 + pw_gzrf2rmps1a + pw_gzrf2rmps1 + pw_gzrf2rmps1d
        + pw_gxwmps1/2 - psd_rf_wait + psd_grd_wait + DABSETUP;

    te_180 = RUP_GRD(2*(IMax(3, ta_180, tb_180, min180te/2)));

    postemp = RUP_GRD( (RUP_GRD(1ms) + pw_rf1mps1/2 - PSoff90)
                       + (te_180/2) - pw_rf2mps1/2 );

    SLICESELZ2( rf2mps1, RF2_APS1_SLOT, postemp, pw_rf2mps1, thk_rf2mps1,
                flip_rf2mps1, cyc_rf2mps1, TYPNDEF, ps1loggrd );

    /* crushers */
    postemp = pbeg(&gzrf2mps1,"gzrf2mps1",0) - pw_gzrf2lmps1 - pw_gzrf2lmps1d;
    TRAPEZOID(ZGRADB, gzrf2lmps1, postemp, 0, TYPNDEF, ps1loggrd );
  
    TRAPEZOID( ZGRADB, gzrf2rmps1, pendall(&gzrf2mps1,0), 0, TYPNDEF, ps1loggrd );
  
    postemp = RUP_GRD(pmid(&gzrf2mps1,"gzrf2mps1",0)+ (te_180/2) - pw_gxwmps1/2);
    TRAPEZOID(read_axis, gxwmps1, postemp, 0, TYPNDEF, ps1loggrd);    

    ACQUIREDATA( echo1mps1, pbeg(&gxwmps1,"gxwmps1",0)+psd_grd_wait, , ,DABNORM);

    ATTENUATOR(attenuator_keymps1, pend(&gxwmps1, "gxwmps1",0));

    SEQLENGTH(seqmps1, ps1_tr, seqmps1);

    return SUCCESS;
}

/*
 *  CFLpulsegen
 *  
 *  Type: Private Function
 *  
 *  Description:
 *  
 */
STATUS
CFLpulsegen( INT posstart )
{
    INT postemp;
    INT tsamp_delay_cfl;

    tsamp_delay_cfl = RUP_GRD(1ms);

    /*  90 slice sel pulse  */
    SLICESELZ2(rf1cfl, RF1_CFL_SLOT, posstart, pw_rf1cfl,
               thk_rf1cfl, flip_rf1cfl, cyc_rf1cfl, TYPNDEF, cflloggrd);
  
    /* Refocusing Pulse */
    TRAPEZOID(ZGRAD, gz1cfl, pbeg(&gzrf1cfl,"gzrf1cfl",0) + pw_gzrf1cfl
              + pw_gzrf1cfld + pw_gz1cfla, 0, TYPNDEF, cflloggrd);
  
    /* Data Acquisiton with 2K filter */
    ACQUIREDATA(cfl_fid, pendall(&gz1cfl,0) + tsamp_delay_cfl, , ,DABNORM);
    /* Assert the ESSP flag on the rf1cfl pulse */
    attenflagon(&rf1cfl, 0);

    postemp = RUP_GRD(pendall(&gz1cfl,0) + tsamp_delay_cfl + cfl_tdaq + pw_gykcfla);

    ATTENUATOR(cfl_attenkey, postemp);
    TRAPEZOID(YGRAD, gykcfl, postemp, 0, TYPNDEF, cflloggrd);
  
    SEQLENGTH(seqcfl, cfl_tr, seqcfl);

    return SUCCESS;
}

/*
 *  RCVNpulsegen
 *  
 *  Type: Private Function
 *  
 *  Description:
 *  
 */
STATUS
RCVNpulsegen( INT posstart )
{
    INT postemp, postemp2;
    INT tsamp_delay_rcvn;
    INT unblank_delay_rcvn;

    /* MRIhc47602/MRIhc47515/GEHmr03545 */
    if ( rcvn_flag == 1 )
    {
        TRAPEZOID(XGRAD, gxkrcvn, RUP_GRD(posstart + pw_gxkrcvna), 0, TYPNDEF, rcvnloggrd);
        TRAPEZOID(YGRAD, gykrcvn, RUP_GRD(posstart + pw_gykrcvna), 0, TYPNDEF, rcvnloggrd);
        TRAPEZOID(ZGRAD, gzkrcvn, RUP_GRD(posstart + pw_gzkrcvna), 0, TYPNDEF, rcvnloggrd);

    }
    else if ( rcvn_flag == 2 ) /* extra delay before rcvn */
    {
        WAIT(ZGRAD,rcvn_wait,RUP_GRD(posstart),GRAD_UPDATE_TIME);
    }
    SEQLENGTH(pre_rcvn, RUP_GRD(pre_rcvn_tr+posstart), pre_rcvn);

    unblank_delay_rcvn = RUP_GRD(1ms);
    tsamp_delay_rcvn = RUP_GRD(2ms);

    /* RCVRUNBLANK & RCVRBLANK mechanism is normally taken care 
       of by the RF pulse itself. However RCVN entry point does
       not have any RF pulse so we wrap data acqusition window. */

    /* Unblank receiver before Data Acquisition window */
    RCVRUNBLANK(rcvrbl,unblank_delay_rcvn,);
  
    /* Data Acquisiton with 2K filter */
    ACQUIREDATA(rcvn_fid, tsamp_delay_rcvn, , ,DABNORM);

    postemp  = RUP_GRD(tsamp_delay_rcvn + rcvn_tdaq);
    postemp2 = RUP_GRD(unblank_delay_rcvn + postemp);

    ATTENUATOR(rcvn_attenkey, postemp);

    /* Blank receiver after Data Acquisition is completed */
    RCVRBLANK(rcvrbl2,postemp2,);
  
    SEQLENGTH(seqrcvn, rcvn_tr, seqrcvn);

    return SUCCESS;
}

/*
 *  CFHpulsegen
 *  
 *  Type: Private Function
 *  
 *  Description:
 *  
 */
STATUS
CFHpulsegen( INT posstart )
{
    INT postemp;
    INT cfh_te;   /* Run at operator TE unless below min cfh te.
                     Then run at min cfh te */
    INT cfh_te2;   /* vmx 05/02/95 YO */
    INT cfh_acq_window_pos;   /* vmx 05/02/94 YO */
    INT tsamp_delay_cfh;
    INT start_time = 0;
    INT start_time_CS = 0;
    INT start_time_MT = 0;
    INT pos_rf2;
    INT newcfh_tr;
    INT min_ssp;

    /* variables for presscfh */
    INT pos_rf3;

    tsamp_delay_cfh = RUP_GRD(1ms);

    start_time = posstart;
    start_time_MT = posstart;
    start_time_CS = posstart;

    /* Check for CFH Inversion inclusion */
    if (PSD_ON == PSir)
    {
        short *temp_wave_space; /* temporary waveform space for rf scaling */
        short *wave_space; /* temporary waveform space for rf scaling */

        /* adiabatic pulse */
        SLICESELZEXT2(rf0cfh, posstart, pw_rf0cfh, 
               ((presscfh_ctrl == PRESSCFH_NONE) ? opslthick : presscfh_ir_slthick ), 
                      flip_rf0cfh, cyc_rf0cfh, 0,
                      1, NULL, res_rf0cfh, shNvrg5b.rho, RF0_CFH_SLOT,
                      TYPNDEF, cfhloggrd);

        EXTWAVE2(THETA, omegarf0cfh, posstart, pw_omegarf0cfh,-1,
                res_rf0cfh, shNvrg5b.pha, 0.0, RF0_CFH_SLOT);

        TRAPEZOID(YGRAD, gyrf0kcfh, pend(&gzrf0cfh,"gzrf0cfh",0) +
                  pw_gyrf0kcfha, 0, TYPNDEF, cfhloggrd);
    
        /* Setting up WAIT duration */
        /* 1ms is added at the end so that ssp sequencer has no overlap 
           as 'frq', and RF amp unblank pulses are played. */
        dur_invse = RUP_GRD(cfh_ti - pw_rf0cfh/2 - pw_gyrf0kcfha - pw_gyrf0kcfh 
                            - pw_gyrf0kcfhd  - pw_gzrf1cfha - pw_rf1cfh/2 - 1ms);
#ifdef PSD_CFH_CHEMSAT
        if ((cs_sat == PSD_ON) && PScs_sat)
        {
            /* GEHmr03577 : Subtract specir_delay in case of SPECIAL. */
#ifdef PSD_CFH_CHEMSAT_SPECIAL
            dur_invse -= RUP_GRD(cs_sattime - specir_delay);
#else
            dur_invse -= RUP_GRD(cs_sattime);
#endif
        }
#endif
%ifdef PSD_CFH_MT
        if ((opmt == PSD_ON) && PSmt)
        {
            dur_invse -= RUP_GRD(mt_time);
        }
%endif /* PSD_CFH_MT */
        dur_invse = RUP_GRD(dur_invse);

        WAIT(ZGRAD,zticfh,pend(&gyrf0kcfhd,"gyrf0kcfhd",0),GRAD_UPDATE_TIME);
        WAIT(RHO,rticfh,pend(&gyrf0kcfhd,"gyrf0kcfhd",0),GRAD_UPDATE_TIME);
        WAIT(XGRAD,xticfh,pend(&gyrf0kcfhd,"gyrf0kcfhd",0),GRAD_UPDATE_TIME);
        WAIT(YGRAD,yticfh,pend(&gyrf0kcfhd,"gyrf0kcfhd",0),GRAD_UPDATE_TIME);
        WAIT(SSP,sticfh,pend(&gyrf0kcfhd,"gyrf0kcfhd",0),GRAD_UPDATE_TIME);
    
        setperiod(dur_invse,&zticfh,0);
        setperiod(dur_invse,&rticfh,0);
        setperiod(dur_invse,&xticfh,0);
        setperiod(dur_invse,&yticfh,0);
        setperiod(dur_invse,&sticfh,0);
    
        /* Change start time for the 90 180 sequence; pw_gzrf1cfha added to 
           offset start_time calc in rf1cfh call */
        start_time    = pmid(&gzrf0cfh,"gzrf0cfh",0) + cfh_ti - pw_rf1cfh/2;
        start_time_MT = pend(&gyrf0kcfhd,"gyrf0kcfhd",0) + dur_invse + 300us;
        start_time_CS = pend(&gyrf0kcfhd,"gyrf0kcfhd",0) + dur_invse + 300us;

        amp_gyrf0kcfh = ia_gyrf0kcfh;
    }
    /* End of IR sequence check */
  
  
%ifdef PSD_CFH_MT
    if ((opmt == PSD_ON) && PSmt)
    {
        MTPG(start_time_MT, &mtcfh_index);
        mtcfh_index -= 1;
        start_time_CS += mt_time;

        if ((PSir != PSD_ON) && (oppseq != PSD_IR || ssfse_ir_on == PSD_OFF))
           start_time += mt_time;
    }
%endif /* PSD_CFH_MT */

#ifdef PSD_CFH_CHEMSAT
    if ((cs_sat == PSD_ON) && PScs_sat)
    {
        ChemSatPG(start_time_CS, &cscfh_satindex);
        cscfh_satindex -= 1;
    }
    if( (PSir != PSD_ON) && ( ((oppseq!=PSD_IR)
                               || (ssfse_ir_on == PSD_OFF)) && PScs_sat ) )
    {
        /* MRIge30640 - already caught in SLICESELZ2 call below! */
        /* GEHmr03577 : Subtract specir_delay in case of SPECIAL. */
#ifdef PSD_CFH_CHEMSAT_SPECIAL
        start_time += (cs_sattime - specir_delay);
#else
        start_time += cs_sattime;
#endif
    }
#endif
  
    /*  90 slice sel pulse  */
    SLICESELZ2( rf1cfh, RF1_CFH_SLOT, RUP_GRD(start_time+pw_gzrf1cfha),
                pw_rf1cfh, 
                ((presscfh_ctrl == PRESSCFH_NONE) ? thk_rf1cfh : presscfh_slthick ), 
                flip_rf1cfh, cyc_rf1cfh, TYPNDEF, cfhloggrd );

    min_ssp = RUP_GRD(-rfupa + rfupd + RFUNBLANK_LENGTH + RFFREQ_LENGTH);

    if(PSfield_strength <= B0_5000)
    {

        cfh_te = (0.5 * pw_rf1cfh + PSoff90 + pw_rf2cfh
                  + IMax(2, min_ssp,
                         (pw_gzrf1cfhd + pw_gzrf2lcfha
                          + pw_gzrf2lcfh + pw_gzrf2lcfhd))
                  + pw_gzrf2rcfha + pw_gzrf2rcfh + pw_gzrf2rcfhd
                  + (cfh_tdaq  * cfh_ec_position));

        cfh_te2 = ((IMax(2, min_ssp,
                         (pw_gzrf2lcfha + pw_gzrf2lcfh + pw_gzrf2lcfhd
                          + pw_gzrf2rcfha + pw_gzrf2rcfh + pw_gzrf2rcfhd))
                    + pw_rf2cfh) / 2
                   + (cfh_tdaq * cfh_ec_position));
        cfh_te2 *= 2;

        cfh_te = IMax(3, cfh_te, cfh_te2, min180te);

        pos_rf2 = RDN_GRD((int)(pmid(&gzrf1cfh,"gzrf1cfh", 0)
                                + 0.5 * cfh_te - 0.5 * pw_rf2cfh));
    }
    else
    {
        if( presscfh_ctrl == PRESSCFH_NONE )
        {
            cfh_te =  0.5 * (pw_rf1cfh + pw_rf2cfh) + PSoff90
                           + IMax(2, min_ssp,
                                  (pw_gzrf1cfhd + pw_gzrf2lcfha
                                   + pw_gzrf2lcfh + pw_gzrf2lcfhd));
            if (opspf == 0)
            {
                cfh_te = 2*IMax(3, cfh_te, pw_gxrf2cfhd, pw_gzrf2rcfha + pw_gzrf2rcfh + pw_gzrf2rcfhd);
            }
            else
            {
                cfh_te = 2*IMax(3, cfh_te, pw_gyrf2cfhd, pw_gzrf2rcfha + pw_gzrf2rcfh + pw_gzrf2rcfhd);
            }

            /* MRIge90312 - use 50ms TE for cfh */
            cfh_te = IMax(3, cfh_te, min180te, eff_cfh_te);

            pos_rf2 = RDN_GRD((int)(pmid(&gzrf1cfh,"gzrf1cfh", 0) - PSoff90 
                                    + 0.5 * cfh_te - 0.5 * pw_rf2cfh));
        }
        else
        {
            int temp_time = 0;

            cfh_te = IMax(2, presscfh_minte, eff_cfh_te);
            if( cfh_steam_flag != PSD_ON )
            {
                int echo1te = 0;

                echo1te =  0.5 * (pw_rf1cfh + pw_rf2cfh) + PSoff90 +
                    IMax(2, min_ssp, (pw_gzrf1cfhd + pw_gzrf2lcfha + pw_gzrf2lcfh + pw_gzrf2lcfhd))
                    + presscfh_wait_rf12;

                echo1te = 2*IMax(3, echo1te, pw_gxrf2cfhd, pw_gzrf2rcfha + pw_gzrf2rcfh + pw_gzrf2rcfhd);

                echo1te = IMax(2, echo1te, min180te);
                pos_rf2 = RDN_GRD((int)(pmid(&gzrf1cfh,"gzrf1cfh", 0) - PSoff90 
                                        + 0.5 * echo1te - 0.5 * pw_rf2cfh));

                temp_time = IMax(2, pw_gyrf3cfha, pw_gzrf3lcfha + pw_gzrf3lcfh + pw_gzrf3lcfhd);
                cfh_te = IMax(2, 2*temp_time + echo1te + pw_rf3cfh, cfh_te);

                pos_rf3 = RDN_GRD((int)(pos_rf2 + 0.5*pw_gxrf2cfh + 0.5*cfh_te - 0.5*pw_rf3cfh));
            }
            else
            {
                int mix_time = 0;

                temp_time = IMax(2, 2*min_ssp,  pw_gzrf2rcfha + pw_gzrf2rcfh + pw_gzrf2rcfhd 
                                 + steam_pg_gap + pw_gzrf3lcfha + pw_gzrf3lcfh + pw_gzrf3lcfhd);
                mix_time = RUP_GRD(0.5 * (pw_gxrf2cfh + pw_gyrf3cfh) + temp_time);

                pos_rf2 = RUP_GRD( pmid(&gzrf1cfh,"gzrf1cfh",0) + 0.5*cfh_te - 0.5*pw_gxrf2cfh );
                pos_rf3 = RUP_GRD( pos_rf2 + 0.5*pw_gxrf2cfh + mix_time - 0.5*pw_gyrf3cfh );
            }
        }
    }

    if (rfpulseInfo[RF2_CFH_SLOT].change==PSD_ON)
    {
        res_rf2cfh = rfpulseInfo[RF2_CFH_SLOT].newres;
    }

    SINC(RHO, rf2cfh, pos_rf2, pw_rf2cfh, a_rf2cfh,,,,, cfhloggrd);
    if( presscfh_ctrl != PRESSCFH_NONE ){ /* for presscfh_ctrl */
          SINC(RHO, rf3cfh, pos_rf3, pw_rf3cfh, a_rf3cfh,,,,, cfhloggrd);  
    } 

    if(PSdebugstate)	/* vmx 05/02/95 YO */
    {
	printf("CFH : TE = %d\n", cfh_te);
	printf("CFH : Mid Position of rf2cfh = %d\n", (int)(pos_rf2+pw_rf2cfh/2));
    }

    /* FOV selective gradients */
    if (opspf == 0 || presscfh_ctrl != PRESSCFH_NONE)
    {
        TRAPEZOID(XGRAD, gxrf2cfh, pbegall(&rf2cfh,0) - psd_rf_wait, 
                  0, TYPNDEF, cfhloggrd);
    }
    else
    {
        TRAPEZOID(YGRAD, gyrf2cfh, pbegall(&rf2cfh,0) - psd_rf_wait, 
                  0, TYPNDEF, cfhloggrd);
    }
    
    /* Z crushers */
    TRAPEZOID(ZGRAD, gzrf2lcfh, pbegall(&rf2cfh,0)  
              - (pw_gzrf2lcfh + pw_gzrf2lcfhd) - psd_rf_wait, 0, TYPNDEF, cfhloggrd);
    TRAPEZOID(ZGRAD, gzrf2rcfh, pendall(&rf2cfh,0) + pw_gzrf2rcfha 
              - psd_rf_wait, 0, TYPNDEF, cfhloggrd);

   if( presscfh_ctrl != PRESSCFH_NONE ){
        TRAPEZOID(YGRAD, gyrf3cfh, pbegall(&rf3cfh,0) - psd_rf_wait, 
                      0, TYPNDEF, cfhloggrd);

        /* Z crushers */
        TRAPEZOID(ZGRAD, gzrf3lcfh, pbegall(&rf3cfh,0)  
                  - (pw_gzrf3lcfh + pw_gzrf3lcfhd) - psd_rf_wait, 0, TYPNDEF, cfhloggrd);
        TRAPEZOID(ZGRAD, gzrf3rcfh, pendall(&rf3cfh,0) + pw_gzrf3rcfha 
                  - psd_rf_wait, 0, TYPNDEF, cfhloggrd);

        /* steam_flag */
        if( cfh_steam_flag == PSD_ON ){
            INT pos_g1cfh = 0;
            pos_g1cfh = RUP_GRD( pbeg(&gzrf2lcfh, "gzrf2lcfh", 0) - pw_gzrf2lcfha 
                                 - pw_gy1cfh - pw_gy1cfhd );
            TRAPEZOID(YGRAD, gy1cfh, pos_g1cfh, 0, TYPNDEF, cfhloggrd);
            pos_g1cfh = RUP_GRD( pend(&gzrf3rcfh, "gzrf3rcfh", 0) + pw_gzrf3rcfhd + pw_gx1cfha );
            TRAPEZOID(XGRAD, gx1cfh, pos_g1cfh, 0, TYPNDEF, cfhloggrd);
        }
    }

    /* Data Acquisiton with .5K/.25k filter */
    if(PSfield_strength <= B0_5000)	/* vmx 05/02/94 */
    {
	cfh_acq_window_pos = RUP_GRD( (int)(pmid(&gzrf1cfh,"gzrf1cfh", 0)
                                            + cfh_te - (cfh_tdaq * cfh_ec_position)));
    }
    else
    {
        if( presscfh_ctrl == PRESSCFH_NONE ) {
	    cfh_acq_window_pos = RUP_GRD(pendall(&gzrf2rcfh,0) + tsamp_delay_cfh);
        } else {
            if( cfh_steam_flag != PSD_ON ){
	        cfh_acq_window_pos = RUP_GRD(pendall(&gzrf3rcfh,0) + tsamp_delay_cfh);
            }else{
                cfh_acq_window_pos = RUP_GRD(pendall(&gx1cfh,0) + tsamp_delay_cfh);
            }
        }
    }

    if(PSdebugstate)	/* vmx 05/02/95 YO */
    {
	printf("CFH : Start of data window = %d\n", cfh_acq_window_pos);
    }

    ACQUIREDATA( cfh_fid, cfh_acq_window_pos, , , DABNORM );
    /* vmx 05/02/05 YO */
    /* Assert the ESSP flag on the rf1cfh  and rf2cfh pulse */
    attenflagon(&rf1cfh, 0);  
    attenflagon(&rf2cfh, 0);
    if(presscfh_ctrl != PRESSCFH_NONE)
        attenflagon(&rf3cfh, 0); /* for presscfh */

    postemp = RUP_GRD(cfh_acq_window_pos + cfh_tdaq + pw_gykcfha);
    ATTENUATOR(cfh_attenkey, postemp);
    TRAPEZOID(YGRAD, gykcfh, postemp, 0, TYPNDEF, cfhloggrd);

    /*  If the TE is so long that the readout and killer are pushed out beyond
        the default cfh_tr, cfh_tr must be increased.  Setting cfh_tr to the
        end of the killer + 10ms (time_ssi should never be more than 10ms) 
        should do the trick.   */

    newcfh_tr = RUP_GRD( (((pendall(&gykcfh,0)+10ms)>cfh_tr) ? (pendall(&gykcfh,0)+10ms) : cfh_tr) );

    if(PSdebugstate)	/* vmx 05/02/95 YO */
    {
	printf("CFH : TR = %d\n", newcfh_tr);
    }

    SEQLENGTH(seqcfh, newcfh_tr, seqcfh);

    return SUCCESS;
}


/*
 *  CoilSwitchPG
 *  
 *  Type: Private Function
 *  
 *  Description: Creates a ssp sequence which can set RF HUB index on
 *  RFHUBSEL. The sequence length needs to change depending upon
 *  setrcvportimm flag. If setrcvportimm needs to be called we need to
 *  provide additional time before starting to acquire as there is time
 *  delay in setting HW. So we add a wait pulse whos pulsewidth will be
 *  decided based on setrcvportimm flag. We also need a 'delay' sequence
 *  as explained in MRIhc14300.
 *  
 */
STATUS
CoilSwitchPG( void )
{
    INT PosContRFHubSel;

    /* SSP Packet for setting the hub index corresponding to the desired
     * coil configuration */
    short dcontrfhubsel[4] = {
        SSPDS,
        SSPOC | RFHUBSEL,
        SSPD,
        SSPDS
    };

    /* SSP Packet for changing receiver input */
    short dcontrfsel[4] = {     
        SSPDS,
        SSPOC | RRFSEL,
        SSPD | RFAUX,
        SSPDS
    };

    PosContRFHubSel = 15us + delay_rfhubsel;

    /* SSP sequence for changing RF Hub index for coil switch */
    SSPPACKET(contrfhubsel, PosContRFHubSel, pw_contrfhubsel, dcontrfhubsel, 0);

    /* SSP Sequence for changing receiver input */
    SSPPACKET(contrfsel, pendallssp(&contrfhubsel, 0), pw_contrfsel, dcontrfsel, 0);

    /* Insert a wait pulse to allow us to change the actual TR when
       a setrcvportimm() call is necessary */
    WAIT(SSP, csw_wait, pendallssp(&contrfsel, 0), SSP_UPDATE_TIME);

    csw_tr = 15us + delay_rfhubsel + pw_contrfhubsel + pw_contrfsel
        + SSP_UPDATE_TIME + csw_time_ssi;

    if( csw_tr < 1ms ) {
        /* Switch time needs to be long enough for RF Hub to switch the coils.
           This is much less than 1ms. */
        csw_tr = RUP_GRD(1ms);
    }
    SEQLENGTH(seqcsw, RUP_GRD(csw_tr - csw_time_ssi), seqcsw);

    /* MRIhc14300: Short wait pulse before setrcvportimm to avoid race
       condition with SCP */
    SEQLENGTH(seqcswWaitBefore, RUP_GRD(csw_wait_before), seqcswWaitBefore);

    return SUCCESS;
}

/*
 *  PSpulsegen
 *  
 *  Type: Public Function
 *  
 *  Description:
 *  
 */
STATUS
PSpulsegen( void )
{
    INT posstart;

    /* MRIge82455 */
    psc_vol_index = 0;

    posstart = RUP_GRD(IMax(2, pw_gzrf1mps1a, pw_gzrf1cfla) + 1ms);

    /***********************************************************************
     * MPS1/APS1 SECTION
     ***********************************************************************/

    PS1pulsegen( posstart );

    /***********************************************************************
     * CFL SECTION
     ***********************************************************************/

    CFLpulsegen( posstart );

    /***********************************************************************
     * RCVN SECTION
     ***********************************************************************/

    /* GEHmr03545 */
    RCVNpulsegen( posstart );

    /***********************************************************************
     * CFH SECTION
     ***********************************************************************/

    CFHpulsegen( posstart );

    /***********************************************************************
     * CoilSwitch SECTION
     ***********************************************************************/

    CoilSwitchPG( );

    return SUCCESS;
}   /* end PSpulsegen() */


/*
 *  FTGpulsegen
 *  
 *  Type: Public Function
 *  
 *  Description:
 *  
 */
STATUS
FTGpulsegen( void )
{
    INT ftgl_tr = 2s;
    INT PosGz1FTG;          /* Slice dephaser location   */
    INT PosReadoutWindow;   /* Readout window location   */
    INT PosReadoutWindow2;  /* Readout window location   */
    INT postemp;
    INT prescan_start;

    ftgl_tr = ftgtr;

    /* MRIge30645 */
    /* Need to change posstart to eliminate dwnld failures with .8 G/cm. */
    prescan_start = RUP_GRD(pw_gzrf1ftga + 1ms);

    /* Z-axis slice selective with x readout gradient for 1-d projection
       after theta2 pulse; positioning of signals after theta3 pulse */

    /* Theta1 selective pulse */
    SLICESELZ2(rf1ftg, RF1_FTG_SLOT, prescan_start, pw_rf1ftg, FTGopslthickz1, flip_rf1ftg,
               cyc_rf1ftg, TYPNDEF, ps1loggrd);

    /* Theta1 rephaser - split into two parts, 2nd part after rf2tg */
    /*                   is defined.                                */
    PosGz1FTG = pendall(&gzrf1ftg, 0) + pw_gz1ftga;

    TRAPEZOID( ZGRAD, gz1ftg, PosGz1FTG, 0, TYPNDEF, ps1loggrd );


    /* Theta2 selective pulse */
    postemp = (int) (pend(&rf1ftg,"gzrf1ftg",0)-pw_rf1ftg/2+FTGtau1-pw_rf2ftg/2);

    SLICESELZ2(rf2ftg, RF2_FTG_SLOT, RUP_GRD(postemp), pw_rf2ftg,
               FTGopslthickz2, flip_rf2ftg, cyc_rf2ftg, TYPNDEF, ps1loggrd);

    /* Theta2 rephaser - split into two parts: gz2tg and gz2btg (below) */
    PosGz1FTG = pendall(&gzrf2ftg, 0)+pw_gz2ftga;

    TRAPEZOID( ZGRAD, gz2ftg, PosGz1FTG, 0, TYPNDEF, ps1loggrd );


    /* Theta3 selective pulse */
    postemp = (int) (pend(&rf1ftg,"gzrf1ftg",0)-pw_rf1ftg/2+FTGtau2-pw_rf3ftg/2);

    SLICESELZ2(rf3ftg, RF3_FTG_SLOT, RUP_GRD(postemp), pw_rf3ftg,
               FTGopslthickz3, flip_rf3ftg, cyc_rf2ftg, TYPNDEF, ps1loggrd);

    /* Theta3 rephaser */
    PosGz1FTG = pendall(&gzrf3ftg, 0)+pw_gz3ftga;

    TRAPEZOID( ZGRAD, gz3ftg, PosGz1FTG, 0, TYPNDEF, ps1loggrd );

    /*----------------------------------------------------------*/
    /* Readout windows and dephasers                            */ 
    /*----------------------------------------------------------*/

    postemp = (int) (pbeg(&rf2ftg,"gzrf2ftga",0)-pw_gx1ftg-pw_gx1ftgd);

    TRAPEZOID( read_axis, gx1ftg, RUP_GRD(postemp), 0, TYPNDEF, ps1loggrd );

    postemp = (int) (pendall(&gzrf2ftg,0) + pw_gx1bftga);

    TRAPEZOID( read_axis, gx1bftg, RUP_GRD(postemp), 0, TYPNDEF, ps1loggrd);

    PosReadoutWindow=RUP_GRD((int)(pend(&rf1ftg,"gzrf1ftg",0) - pw_rf1ftg/2 + 2*FTGtau1 - pw_gxw1ftg/2));

    /* HD--Error Check For Gradient Overlapp. If gradients gx1bftg
     * and gxw1ftg overlap then shift the start of gxw1ftg after end of
     * gx1bftg 
     */
    if ( pendall(&gx1bftg,0) >= (PosReadoutWindow - pw_gxw1ftga)){
        PosReadoutWindow = pendall(&gx1bftg,0) + pw_gxw1ftga ;
    }

    TRAPEZOID( read_axis, gxw1ftg, PosReadoutWindow, 0, TYPNDEF, ps1loggrd );

    PosReadoutWindow=RUP_GRD((int)(pend(&gxw1ftg,"gxw1ftgd",0))+pw_postgxw1ftga);

    TRAPEZOID( read_axis, postgxw1ftg, PosReadoutWindow, 0, TYPNDEF, ps1loggrd );

    PosReadoutWindow =  RUP_GRD((int)(pend(&gxw1ftg, "gxw1ftga", 0)));

    ACQUIREDATA( echo1ftg, PosReadoutWindow+psd_grd_wait, , , DABNORM);

    /* Second part of theta2 rephaser */
    PosGz1FTG = pbegall(&rf3ftg, 0)-(pw_gz2bftg + pw_gz2bftgd + pw_gzrf3ftga);

    TRAPEZOID( ZGRAD, gz2bftg, PosGz1FTG, 0, TYPNDEF, ps1loggrd);

    /* Another refocusing pulse to insure S1 forms tau1 ms after
       center of rf3.  This is the time at which we want to
       the S1 signal to refocus:  */
    PosReadoutWindow =  RUP_GRD((int)(pendall(&rf3ftg, 0) + pw_gx2ftga));

    TRAPEZOID( read_axis, gx2ftg, PosReadoutWindow, 0, TYPNDEF, ps1loggrd);

    /* Second readout window */
    PosReadoutWindow2 = RUP_GRD((int)(pmidall(&rf3ftg, 0) + FTGtau1 - pw_gxw2ftgleft));

    /* HD--Error Check For Gradient Overlapp. If gradients gx2ftg and
     * gxw2ftg overlap then shift the start of gxw2ftg after end of
     * gx2ftg 
     */
    if ( pendall(&gx2ftg,0) >= (PosReadoutWindow2 - pw_gxw2ftga)){
        PosReadoutWindow2 = pendall(&gx2ftg,0) + pw_gxw2ftga ;
    }

    TRAPEZOID( read_axis, gxw2ftg, PosReadoutWindow2, 0, TYPNDEF, ps1loggrd );

    if (FTGtestpulse == 1)
    {
        PosReadoutWindow =  RUP_GRD((int)(pbegall(&rf3ftg, 0)+pw_gx2test + pw_gx2testd));
        TRAPEZOID( read_axis, gx2test, PosReadoutWindow, 0, TYPNDEF, ps1loggrd);
    }    

    PosReadoutWindow =  RUP_GRD((int)(pend(&gxw2ftg, "gxw2ftga", 0)));

    ACQUIREDATA(echo2ftg, PosReadoutWindow+psd_grd_wait, , , DABNORM);

    ATTENUATOR(ftg_attenkey, RUP_GRD(pbegall(&gxw2ftg,0) + 1ms + pw_gxw2ftg));

    SEQLENGTH(seqftg, ftgl_tr, seqftg);

    return SUCCESS;
}   /* end FTGpulsegen() */


/*
 *  XTGpulsegen
 *  
 *  Type: Public Function
 *  
 *  Description:
 *  
 */
STATUS
XTGpulsegen( void )
{
    INT xtgl_tr = 2s;
    INT PosGz1XTG;          /* Slice dephaser location   */
    INT PosReadoutWindow;   /* Readout window location   */
    INT postemp;
    INT prescan_start;

    xtgl_tr = xtgtr;

    /* MRIge30645 */
    /* Need to change posstart to eliminate dwnld failures with .8 G/cm. */
    prescan_start = RUP_GRD(pw_gzrf1xtga + 1ms);

    /* Theta1 selective pulse */
    SLICESELZ2(rf1xtg, RF1_XTG_SLOT, prescan_start, pw_rf1xtg, XTGopslthick,
               flip_rf1xtg, cyc_rf1xtg, TYPNDEF, ps1loggrd);

    /* Theta1 rephaser - split into two parts, 2nd part after rf2tg */
    /*                   is defined.                                */
    postemp = RUP_GRD(pend(&rf1xtg,"rf1xtg",0)+pw_gykxtgla);
    TRAPEZOID( killer_axis, gykxtgl, postemp, 0, TYPNDEF, ps1loggrd );
    
    postemp = RUP_GRD(pendall(&gykxtgl,0));
    OFFRESPULSE_STX(rf3xtg, postemp, pw_bsrf, flip_rf3xtg);

    PosGz1XTG = RUP_GRD(pend(&rf1xtg,"gzrf1xtg",0)-pw_rf1xtg/2+XTGtau1-
                        pw_rf2xtg/2-pw_gzrf2xtga-pw_gz1xtgd-pw_gz1xtg);
    TRAPEZOID( ZGRAD, gz1xtg, PosGz1XTG, 0, TYPNDEF, ps1loggrd );

    /* Theta2 selective pulse */
    postemp = (int) (pend(&rf1xtg,"gzrf1xtg",0)-pw_rf1xtg/2+XTGtau1-
                     pw_rf2xtg/2);

    SLICESELZ2(rf2xtg, RF2_XTG_SLOT, RUP_GRD(postemp), pw_rf2xtg,
               XTGopslthick, flip_rf2xtg, cyc_rf2xtg, TYPNDEF, ps1loggrd);

    /* Theta2 rephaser - split into two parts: gz2tg and gz2btg (below) */
    PosGz1XTG = pendall(&gzrf2xtg, 0)+pw_gz2xtga;

    TRAPEZOID( ZGRAD, gz2xtg, PosGz1XTG, 0, TYPNDEF, ps1loggrd );

    postemp = RUP_GRD(pendall(&gz2xtg, 0));
    OFFRESPULSE_STX(rf4xtg, postemp, pw_bsrf, flip_rf4xtg);

    postemp = RUP_GRD(pendall(&rf4xtg, 0)+pw_gykxtgra);
    TRAPEZOID(killer_axis, gykxtgr, postemp, 0, TYPNDEF, ps1loggrd);
    
    postemp = RUP_GRD(pendall(&rf2xtg,0)-pw_rf2xtg/2+XTGtau1-pw_gxw1xtg/2
                      -pw_gxw1xtga-pw_gx1bxtgd-pw_gx1bxtg);
    TRAPEZOID( read_axis, gx1bxtg, RUP_GRD(postemp), 0, TYPNDEF, ps1loggrd);

    PosReadoutWindow=RUP_GRD((int)(pend(&rf1xtg,"gzrf1xtg",0) - pw_rf1xtg/2 + 
                                   2*XTGtau1 - pw_gxw1xtg/2));

    /* HD--Error Check For Gradient Overlapp. If gradients gx1bftg
     * and gxw1ftg overlap then shift the start of gxw1ftg after end of
     * gx1bftg 
     */
    if ( pendall(&gx1bxtg,0) >= (PosReadoutWindow - pw_gxw1xtga)){
        PosReadoutWindow = pendall(&gx1bxtg,0) + pw_gxw1xtga ;
    }

    TRAPEZOID( read_axis, gxw1xtg, PosReadoutWindow, 0, TYPNDEF, ps1loggrd );

    PosReadoutWindow =  RUP_GRD((int)(pend(&gxw1xtg, "gxw1xtga", 0)));

    ACQUIREDATA( echo1xtg, PosReadoutWindow+psd_grd_wait, , , DABNORM);

    /* position for ATTENUATOR */
    postemp = RUP_GRD(pbegall(&gxw1xtg,0) + 1ms + pw_gxw1xtg);

    ATTENUATOR(xtg_attenkey, postemp);

    SEQLENGTH(seqxtg, xtgl_tr, seqxtg);

    return SUCCESS;
}   /* end XTGpulsegen() */


/*
 *  ASpulsegen
 *  
 *  Type: Public Function
 *  
 *  Description:
 *  
 */
STATUS
ASpulsegen( void )
{
    INT PosTemp;

    /***************************************
     * Z board
     ***************************************/
    /* Z gradient slice select */
    PosTemp = RUP_GRD(td0as + tleadas - rfupa + pw_gzrf1asa);
    SLICESELZ2( rf1as, RF1_AUTOSHIM, PosTemp, pw_rf1as, asslthick, flip_rf1as,
                cyc_rf1as, TYPNDEF, asloggrd ) ;

    /* Assert the ESSP flag on the rf1as pulse */
    attenflagon(&rf1as, 0);

    /* Z gradient rephaser */
    PosTemp = pendall(&gzrf1as, 0) + pw_gz1asa;
    TRAPEZOID( ZGRAD, gz1as, pendall(&gzrf1as, 0)+pw_gz1asa, 0, TYPNDEF,
               asloggrd );

    /***************************************
     * X board
     ***************************************/
    PosTemp = RUP_GRD(pmid(&gzrf1as,"gzrf1as",0)-off90as+te_as-pw_gxwas/2);
    TRAPEZOID( XGRAD, gxwas, PosTemp, 0, TYPNDEF, asloggrd );

    PosTemp = RUP_GRD(PosTemp+psd_grd_wait);
    ACQUIREDATA(echo1as, PosTemp, , , DABNORM);
    attenflagon(&echo1as,0);

    /* frequency dephaser */
    PosTemp = RUP_GRD(pbegall(&gxwas,0) - (pw_gx1as+pw_gx1asd));
    TRAPEZOID( XGRAD, gx1as, PosTemp, 0, TYPNDEF, asloggrd );

    /*****************************************
     * Attenuator lock
     *****************************************/
    PosTemp = RUP_GRD( pend(&gxwas, "gxwas",0) );
    ATTENUATOR(attenuator_keyas, PosTemp);

    /*****************************************
     * Y board
     *****************************************/
    /* HSI - changed SINUSOID to TRAPEZOID2 */
    /* encode */
    PosTemp = RUP_GRD(pend(&gz1asd,"gz1asd",0));
    TRAPEZOID2(YGRAD, gy1as, PosTemp, TRAP_ALL_SLOPED, , , endview_scaleas, asloggrd);

    /* rewind */
    PosTemp= RUP_GRD(pend(&gxwas,"gxwas",0));
    TRAPEZOID2(YGRAD, gy1ras, PosTemp, TRAP_ALL_SLOPED, , , endview_scaleas, asloggrd);


    /*******************
     * X and Z Killers
     *******************/
    PosTemp= RUP_GRD(pend(&gxwasd,"gxwasd",0) + pw_gxkasa);
    TRAPEZOID(XGRAD, gxkas, PosTemp, 0, TYPNDEF, asloggrd);

    PosTemp= RUP_GRD(pend(&gxwasd,"gxwasd",0) + pw_gzkasa);
    TRAPEZOID(ZGRAD, gzkas, PosTemp, 0, TYPNDEF, asloggrd);

    /**************
     * dixon shifts
     **************/
    PosTemp = RUP_GRD(td0as + tleadas - rfupa);
    CONST(XGRAD, xdixon, PosTemp, GRAD_UPDATE_TIME, 0.0, asloggrd);
    CONST(YGRAD, ydixon, PosTemp, GRAD_UPDATE_TIME, 0.0, asloggrd);

    PosTemp = pend(&gz1asd,"gz1asd",0);
    CONST(ZGRAD, zdixon, PosTemp, GRAD_UPDATE_TIME, 0.0, asloggrd);

    /* just pad the ssp somewhere beyond the rf unblank */
    rfdisable_add = YES;
    PosTemp = RUP_RF(pend(&rf1as,"rf1as",0) + rfupd + 12);
    CONST(SSP, sdixon, PosTemp, GRAD_UPDATE_TIME, 0.0, asloggrd);

    PosTemp = RUP_RF(pbeg(&gzkas,"gzkas",0));
    CONST(SSP, sdixon2, PosTemp, GRAD_UPDATE_TIME+dix_timeas, 0.0, asloggrd);

    rfdisable_add = NO;
    SEQLENGTH(seqaushim, RUP_GRD((int)(tr_as - time_ssias)), seqaushim);
    attenflagon(&seqaushim, 0);

    /***********************************************************
     * Pass Packet sequence
     ***********************************************************/
    PASSPACK(pass_aushim, RUP_GRD(TR_PASS3D-1ms));
    SEQLENGTH(seqpassas, RUP_GRD(TR_PASS3D), seqpassas);

    return SUCCESS;
}   /* end ASpulsegen() */


@rsp PSeplist
/*********************************************************************
 *                    PRESCAN.E RSP SECTION                          *
 *                           PSeplist                                *
 *                                                                   *
 * Additional list of entry points for Prescan.                      *
 *********************************************************************/
         "cfl",
         "cfh",
         "mps1",
         "aps1", 
         "autoshim",
         "fasttg",
         "rcvn",
         "expresstg",
         0	/* 0 is needed for the parser */


@rsp PScore
/*********************************************************************
 *                    PRESCAN.E RSP SECTION                          *
 *                            PScore                                 *
 *                                                                   *
 * Write here the functional code for the real time processing (Tgt  *
 * side). You may declare standard C variables, but of limited types *
 * short, int, long, float, double, and 1D arrays of those types.    *
 *********************************************************************/

long PSrsptrigger[MAX_PSC_VQUANT]={0};   /* prescan trigger */ /* vmx 10/13/94 YI */
long finalpscrot[MAX_PSC_VQUANT][9]={{0}};

%ifndef VOXTG
/*
 *  mps1
 *  
 *  Type: Public Function
 *  
 *  Description:
 *  
 */
STATUS
mps1( void )
{
    rspent = L_MPS1;
    strcpy(psdexitarg.text_arg, "MPS1");

/* begin aps1_mod changes (GE) */
    PSinit(PSrot_mod);
    PSmps1(2);
    rspexit();

    return SUCCESS;
}   /* end mps1() */


/*
 *  aps1
 *  
 *  Type: Public Function
 *  
 *  Description:
 *  
 */
STATUS
aps1( void )
{
    rspent=L_APS1;
    strcpy(psdexitarg.text_arg, "APS1");

/* begin aps1_mod changes (GE) */
    PSinit(PSrot_mod);
    PSmps1(1);
    rspexit();

    return SUCCESS;
}   /* end aps1() */
%endif


/*
 *  cfl
 *  
 *  Type: Public Function
 *  
 *  Description:
 *  
 */
STATUS
cfl( void )
{
    rspent=L_CFL;
    strcpy(psdexitarg.text_arg, "CFL");

    PSinit(PSrot);
    PScfl();
    rspexit();

    return SUCCESS;
}   /* end cfl() */


/*
 *  rcvn
 *  
 *  Type: Public Function
 *  
 *  Description:
 *  
 */
STATUS
rcvn( void )
{
    rspent=L_RCVN;
    strcpy(psdexitarg.text_arg, "RCVN");

    PSinit(PSrot);
    PSrcvn();
    rspexit();

    return SUCCESS;
}   /* end rcvn() */


%ifndef SPECTROSCOPY
/*
 *  cfh
 *  
 *  Type: Public Function
 *  
 *  Description:
 *  
 */
STATUS
cfh( void )
{
    rspent=L_CFH;
    strcpy(psdexitarg.text_arg, "CFH");

    PSinit(rsp_PSrot);
    PScfh();
    rspexit();

    return SUCCESS;
}   /* end cfh() */
%endif

/*
 *  fasttg
 *
 *  Type: Public Function
 *  
 *  Description:
 *  
 */
STATUS
fasttg( void )
{
    rspent=L_FTG;
    PSinit(PSrot_mod);
    PSfasttg();
    rspexit();

    return SUCCESS;
}   /* end fasttg() */

/*
 *  expresstg
 *
 *  Type: Public Function
 *  
 *  Description:
 *  
 */
STATUS
expresstg( void )
{
    rspent=L_XTG;
    PSinit(PSrot_mod);
    PSexpresstg();
    rspexit();

    return SUCCESS;    
}   /* end expresstg() */

/*
 *  autoshim
 *
 *  Type: Public Function
 *  
 *  Description:
 *  
 */
STATUS
autoshim( void )
{
    rspent=L_AUTOSHIM;
    strcpy(psdexitarg.text_arg, "Autoshim");

    ASautoshim();
    rspexit();

    return SUCCESS;
}   /* end autoshim() */


/*
 *  PSmps1
 *  
 *  Type: Public Function
 *  
 *  Description:
 *  
 */
STATUS
PSmps1( INT mps1nex )
{
    INT acq_type;
    SHORT temp_short;
    INT slice_freq;
    FLOAT receive_freq;
    float rsp_PStloc=0;
    float rsp_PSrloc=0;
    int old_psc_vol_index = -1; /* Last prescan volume for SWIFT switching */

    printdbg("Starting PSmps1",PSdebugstate);
    boffset(off_seqmps1);

    scopeon(&seqmps1);	/* Activate scope for core */
    syncon(&seqmps1);	/* Activate sync for core */

    rsp_PStloc = PStloc_mod;
    rsp_PSrloc = PSrloc_mod;
    /* begin aps1_mod changes (GE) */

    slice_freq = ( GAM * a_gzrf1mps1 * rsp_PStloc/
                   (10* TARDIS_FREQ_RES) );
    /* factor 10 is because rsptloc is in mm */
    setfrequency(slice_freq, &rf1mps1,0);

    slice_freq = ( GAM * a_gzrf2mps1 * rsp_PStloc
                   / (10* TARDIS_FREQ_RES) );
    /* factor 10 is because rsptloc is in mm */
    setfrequency(slice_freq, &rf2mps1,0);

    receive_freq = 2 * 16000 * rsp_PSrloc / mpsfov;

    /* end aps1_mod changes (GE) */ 

    setfrequency((int)((PSfreq_offset[rspent] + receive_freq)/TARDIS_FREQ_RES),
		 &echo1mps1, 0);

    if(PSdebugstate)
    {
        printf("\nAPS1/MPS1 Slthick = %f FOV = %f Xmit location = %f, Rcv location  = %f\n", 
                thickPS_mod,mpsfov, PStloc_mod, PSrloc_mod);
        printf("\nAPS1/MPS1 Xmit Freq = %i, Rcv Freq  = %f\n", slice_freq, receive_freq);
        printf("%s\n","APS1/MPS1 Rotation Matrix:");
        printf("\t %6ld %6ld %6ld \n", PSrot_mod[0][0], PSrot_mod[0][1], PSrot_mod[0][2]);
        printf("\t %6ld %6ld %6ld \n", PSrot_mod[0][3], PSrot_mod[0][4], PSrot_mod[0][5]);
        printf("\t %6ld %6ld %6ld \n", PSrot_mod[0][6], PSrot_mod[0][7], PSrot_mod[0][8]);
        fflush(stdout);
    }

    if( PSdebugstate )
    {
        printf("PSmps1: psc_vol_index = %d\n",psc_vol_index );
    }


    
    for (view = -1; view <= 30000; view++)
    {
        for (excitation = 1; excitation <= mps1nex; excitation++)
    	{
            if (view > 0)

      	    {
                if (excitation == mps1nex)
      	      	{
                    attenlockoff(&attenuator_keymps1);
      	      	}
                else
      	      	{
                    attenlockon(&attenuator_keymps1);
      	      	}
      	    }

            if ((view > 0) && (excitation >= 1))
      	    {
                acq_type = (int)DABON;
      	    }
            else
      	    {
                acq_type = (int)DABOFF;
      	    }

            loaddab(&echo1mps1, (int)1,(int)1,(int)1,(int)1,
                    (TYPDAB_PACKETS)acq_type, PSD_LOAD_DAB_ALL); 

            /* For SWIFT(PSmps1): We are trying not to use the CV opswift.
               Instead use psc_vol_index.
               */
            if( (psc_vol_index > 0) && (opvquant > 1))
            {
                if( old_psc_vol_index != psc_vol_index )
                {
                    if( noswitch_slab_psc == PSD_OFF )
                    {
                        rsp_PStloc = rsp_info[psc_vol_index-1].rsptloc;
                        rsp_PSrloc = rsp_info[psc_vol_index-1].rsprloc;
                    }
                    else
                    {
                        rsp_PStloc = rsp_info[PStest_slab-1].rsptloc;
                        rsp_PSrloc = rsp_info[PStest_slab-1].rsprloc;
                    }

                    if( noswitch_coil_psc == PSD_OFF )
                    {
                        if (FAILURE == CoilSwitchSetCoil(coilInfo_tgt[psc_vol_index-1], 0))
                        {
                            return FAILURE;
                        }
                        boffset(off_seqmps1);
                    }
                        
                    {
                        char tempstr[200]={0};
                        sprintf(tempstr,"PSmps1: psc_vol_index=%d,rsp_PStloc=%f\n",psc_vol_index,rsp_PStloc);
                        printdbg(tempstr,PSdebugstate);
                    }
                    old_psc_vol_index = psc_vol_index;
                }
            }

            slice_freq = ( GAM * a_gzrf1mps1 * rsp_PStloc/
                    (10* TARDIS_FREQ_RES) );
            /* factor 10 is because rsptloc is in mm */
            setfrequency(slice_freq, &rf1mps1,0);

            slice_freq = ( GAM * a_gzrf2mps1 * rsp_PStloc
                    / (10* TARDIS_FREQ_RES) );
            /* factor 10 is because rsptloc is in mm */
            setfrequency(slice_freq, &rf2mps1,0);

            receive_freq = 2 * 16000 * rsp_PSrloc / mpsfov;

            /* end aps1_mod changes (GE) */

            setfrequency((int)((PSfreq_offset[rspent] + receive_freq)/TARDIS_FREQ_RES),
                    &echo1mps1, 0);

            /* Makesure psc_vol_index is in range */
            if( (psc_vol_index > 0 && psc_vol_index <= opvquant) && (opvquant > 1)) 
            {
                if(noswitch_slab_psc == PSD_OFF)
                    startseq((short)(psc_vol_index-1), (SHORT)MAY_PAUSE);
                else
                    startseq((short)(PStest_slab-1),(SHORT)MAY_PAUSE);
            }
            else
                startseq((short)0, (SHORT)MAY_PAUSE);

            syncoff(&seqmps1);

            /* Chopper logic */
            getiamp(&temp_short, &rf1mps1,0);
            setiamp((-temp_short),&rf1mps1,0);
	}
    }

    printdbg("Returning from PSmps1",PSdebugstate);

    return SUCCESS;
}   /* end PSmps1() */


/*
 *  PScfl
 *  
 *  Type: 
 *  
 *  Description:
 *  
 */
STATUS
PScfl( void )
{
    INT acq_type; /* enable or disable data acquisiton */
    INT slice_freq; /* transmit frequency for the prescan slice */
    SHORT temp_short; /* temp variable */
    float rsp_PStloc;
    int old_psc_vol_index = -1; /* Last prescan volume for SWIFT switching */

    printdbg("Entering PScfl", PSdebugstate);
    boffset(off_seqcfl); 
    scopeon(&seqcfl);
    syncon(&seqcfl); 
    attenlockoff(&cfl_attenkey);
  
    rsp_PStloc = PStloc;
    slice_freq = GAM * a_gzrf1cfl * PStloc/
        (10* TARDIS_FREQ_RES); /* factor 10 is because rsptloc is in mm */
  
    setfrequency(slice_freq, &rf1cfl,0);
    setfrequency((int)(PSfreq_offset[rspent]/TARDIS_FREQ_RES),
                 &cfl_fid, 0);

    if( PSdebugstate )
    {
        printf("PScfl: psc_vol_index = %d\n",psc_vol_index );
    }

    for (view = -cfl_dda+1 ; view <= 30000; view ++)
    {
        for (excitation = 1; excitation <= cfl_nex; excitation ++)
        {
            if (view > 0) 
            {
                acq_type = (int)DABON;

                if(excitation == cfl_nex)
                    attenlockoff(&cfl_attenkey);
                else
                    attenlockon(&cfl_attenkey);
            } 
            else
            {
                acq_type = (int)DABOFF;
            }

            /* DAB packet is not used in prescan.  loaddab
               is used just to turn on or off the receiver */

            loaddab(&cfl_fid, (int)1,(int)1,(int)1,(int)1,
                    (TYPDAB_PACKETS)acq_type, PSD_LOAD_DAB_ALL);

            /* For SWIFT(PScfl): We are trying not to use the CV opswift.
               Instead use psc_vol_index.
             */
            if( (psc_vol_index > 0) && (opvquant > 1) )
            {
                if( old_psc_vol_index != psc_vol_index )
                {
                    if( noswitch_slab_psc == PSD_OFF )
                    {
                        rsp_PStloc = rsp_info[psc_vol_index-1].rsptloc;
                    }
                    else
                    {
                        rsp_PStloc = rsp_info[PStest_slab-1].rsptloc;
                    }

                    if( noswitch_coil_psc == PSD_OFF )
                    {
                        if (FAILURE == CoilSwitchSetCoil(coilInfo_tgt[psc_vol_index-1], 0))
                        {
                            return FAILURE;
                        }
                        boffset(off_seqcfl); 
                    }

                    {
                        char tempstr[200]={0};
                        sprintf(tempstr,"PScfl: psc_vol_index=%d,rsp_PStloc=%f\n",psc_vol_index,rsp_PStloc);
                        printdbg(tempstr,PSdebugstate);
                    }
                    old_psc_vol_index = psc_vol_index;
                }
            }

            slice_freq = GAM * a_gzrf1cfl * rsp_PStloc/
                (10* TARDIS_FREQ_RES); /* factor 10 is because rsptloc is in mm */

            setfrequency(slice_freq, &rf1cfl,0);

            /* Makesure psc_vol_index is in range */
            if( (psc_vol_index > 0 && psc_vol_index <= opvquant) && (opvquant > 1)) 
            {
                if(noswitch_slab_psc == PSD_OFF)
                    startseq((short)(psc_vol_index-1), (SHORT)MAY_PAUSE);
                else
                    startseq((short)(PStest_slab-1),(SHORT)MAY_PAUSE);
            }
            else
                startseq((short)0, (SHORT)MAY_PAUSE);

            /*printdbg("S", PSdebugstate);*/
            syncoff(&seqcfl);

            if(PSdebugstate && view == 0)
            {
                printf("\n CFL:  Xmit Location = %f Receive Location = %f\n ", PStloc, 0.0 );
                printf ("%s\n","CFL : Rotation Matrix");
                printf("\t %6ld %6ld %6ld \n", PSrot[0][0], PSrot[0][1], PSrot[0][2]);
                printf("\t %6ld %6ld %6ld \n", PSrot[0][3], PSrot[0][4], PSrot[0][5]);
                printf("\t %6ld %6ld %6ld \n", PSrot[0][6], PSrot[0][7], PSrot[0][8]);
                fflush(stdout);
            }

            getiamp(&temp_short, &rf1cfl, 0);
            setiamp(-temp_short, &rf1cfl, 0);

            if (view < 1)
            {
                break; /* Skip excitation loop for disdaqs */
            }

        }  /* end excitation */
    }  /* end view */

    return SUCCESS;
}   /* end PScfl() */


/*
 *  PSrcvn
 *  
 *  Type: 
 *  
 *  Description:
 *  
 */
STATUS
PSrcvn( void )
{
    INT acq_type; /* enable or disable data acquisiton */
    int old_psc_vol_index = -1; /* Last prescan volume for SWIFT switching */

    printdbg("Entering PSrcvn", PSdebugstate);

    /* GEHmr03545 */
    if ( rcvn_flag != PSD_OFF )
    {
        boffset( off_pre_rcvn );
        startseq((SHORT)0, (SHORT)MAY_PAUSE);
    }

    boffset(off_seqrcvn); 
    scopeon(&seqrcvn);
    syncon(&seqrcvn); 
    attenlockoff(&rcvn_attenkey);
  
    setfrequency((int)(PSfreq_offset[rspent]/TARDIS_FREQ_RES),
                 &rcvn_fid, 0);

    if( PSdebugstate )
    {
        printf("PSrcvn: psc_vol_index = %d\n",psc_vol_index );
    }

    for (view = -rcvn_dda+1 ; view <= rcvn_loops ; view ++)
    {
        for (excitation = 1; excitation <= rcvn_nex; excitation ++)
        {
            if (view > 0) 
            {
                acq_type = (int)DABON;

                if(excitation == rcvn_nex)
                    attenlockoff(&rcvn_attenkey);
                else
                    attenlockon(&rcvn_attenkey);
            } 
            else
            {
                acq_type = (int)DABOFF;
            }

            /* DAB packet is not used in prescan.  loaddab
               is used just to turn on or off the receiver */

            loaddab(&rcvn_fid, (int)1,(int)1,(int)1,(int)1,
                    (TYPDAB_PACKETS)acq_type, PSD_LOAD_DAB_ALL);

            /* For SWIFT(PSrcvn): We are trying not to use the CV opswift.
               Instead use psc_vol_index.
             */
            if( (psc_vol_index > 0) && (opvquant > 1) )
            {
                if( old_psc_vol_index != psc_vol_index )
                {
                    /* No need to switch slab*/
                    if( noswitch_coil_psc == PSD_OFF )
                    {
                        if (FAILURE == CoilSwitchSetCoil(coilInfo_tgt[psc_vol_index-1], 0))
                        {
                            return FAILURE;
                        }
                        boffset(off_seqrcvn); 
                    }

                    {
                        char tempstr[200]={0};
                        sprintf(tempstr,"PSrcvn: psc_vol_index=%d\n",psc_vol_index);
                        printdbg(tempstr,PSdebugstate);
                    }
                    old_psc_vol_index = psc_vol_index;
                }
            }
            startseq((short)0,(short)MAY_PAUSE);

            syncoff(&seqrcvn);

            if (view < 1)
            {
                break; /* Skip excitation loop for disdaqs */
            }

        }  /* end excitation */
    }  /* end view */

    return SUCCESS;
}   /* end PSrcvn() */


%ifndef SPECTROSCOPY
/*
 *  PScfh
 *  
 *  Type: Public Function
 *  
 *  Description:
 *  
 */
STATUS
PScfh( void )
{
    INT acq_type;      /* enable or disable data acquisiton */
    INT slice_freq;    /* transmit frequency for the prescan slice */
    INT slice2_freq;   /* transmit frequency for the prescan slice */
    INT slice3_freq;   /* presscfh: transmit frequency for the prescan slice */
    INT ir_slice_freq; /* IR sequence transmit freq. for the prescan slice */
    SHORT temp_short;  /* temp variable */
    long init_cfh_deadtime; /* initial deadtime of the seqcfh sequence */
    long new_cfh_deadtime;  /* updated cfh deadtime */

    float rsp_PStloc;    /* RSP transmit location */
    float rsp_PSrloc;    /* RSP receive location */
    float rsp_PSphasoff; /* RSP phase off location */
    int old_psc_vol_index = -1; /* Last prescan volume for SWIFT switching */
    
    showfp = 0;
    slice3_freq = 0; /* initialize */

    printdbg("Entering PScfh", PSdebugstate);
    boffset(off_seqcfh); 
    scopeon(&seqcfh);
    syncon(&seqcfh); 
    attenlockoff(&cfh_attenkey);

    if (psc_vol_index == 0) {

        slice_freq = cfh_rf1freq;
        setfrequency(slice_freq, &rf1cfh,0);

        slice2_freq = cfh_rf2freq;
        setfrequency((int)slice2_freq, &rf2cfh,0);
        setfrequency((int)(PSfreq_offset[rspent]/TARDIS_FREQ_RES), &cfh_fid, 0);
    }

#ifdef PSD_CFH_CHEMSAT

    if ((cs_sat == PSD_ON) && PScs_sat)
    {
        if (PSir != PSD_ON)
        {
            setiamp(ia_rfcssatcfh, &rfcssat, cscfh_satindex);
        }
        else
        {
            setiamp(0, &rfcssat, cscfh_satindex);
        }
    }
    cstun = 0;

#endif
%ifdef PSD_CFH_MT
    if ((opmt == PSD_ON) && PSmt)
    {
        setiamp(ia_rfmtcfh, &rfmt, mtcfh_index);
        mttun = 0;
    }
%endif /* PSD_CFH_MT */


    if (PSD_ON == PSir)
    {
        /* Inversion */
        /* factor 10 is because rsptloc is in mm */
        ir_slice_freq = GAM * a_gzrf0cfh * PStloc/(10 * 
                                                   TARDIS_FREQ_RES); 
        setfrequency(ir_slice_freq, &rf0cfh,0);
    }

    /* Setting tislice to PSslice_num+1 so user sees slices start at 1*/
    tislice=PSslice_num+1;
    tislice_start=PSslice_num+1;

    /* Setting titime to opti; results in every entry into CFH having 
       titime=opti*/
    titime = cfh_ti/1000;
    noir = 0;

    /* Finding initial deadtime (titime = opti) of seqcfh */
    getperiod( &init_cfh_deadtime, &seqcfh, 0 );

    if( PSdebugstate )
    {
        printf("PScfh: psc_vol_index = %d\n",psc_vol_index );
    }


    for (view = -cfh_dda+1 ; view <= 30000; view ++)
    {
        /* Modified for MULTI VOLUME Prescan - AP */
        if( presscfh_ctrl != PRESSCFH_NONE )
        {
            int cfh_slab_index = 0;

            if (oppscvquant == 1){
                cfh_slab_index = 0;
            }else{
                cfh_slab_index = noswitch_slab_psc ? PStest_slab : psc_vol_index;
                cfh_slab_index = IMax(2, 0, cfh_slab_index-1);
            }
            if( presscfh_ctrl != PRESSCFH_SHIMVOL ){
                rsp_PStloc = presscfh_info[cfh_slab_index].oppsctloc;
                rsp_PSrloc = presscfh_info[cfh_slab_index].oppscrloc;
                rsp_PSphasoff = presscfh_info[cfh_slab_index].oppscphasoff;
            }else{
                rsp_PStloc = rsp_psc_info[cfh_slab_index].rsppsctloc;
                rsp_PSrloc = rsp_psc_info[cfh_slab_index].rsppscrloc;
                rsp_PSphasoff = rsp_psc_info[cfh_slab_index].rsppscphasoff;
            }

            slice_freq = GAM * a_gzrf1cfh * rsp_PStloc / (10 * TARDIS_FREQ_RES);
            slice2_freq = GAM * rsp_PSrloc * a_gxrf2cfh / (10 * TARDIS_FREQ_RES);
            slice3_freq = GAM * rsp_PSphasoff * a_gyrf3cfh / (10 * TARDIS_FREQ_RES);

            if( PSD_ON == PSir ) {
                ir_slice_freq = ( GAM * a_gzrf0cfh * rsp_PStloc
                                  / (10 * TARDIS_FREQ_RES) );
                setfrequency(ir_slice_freq, &rf0cfh,0);
            }
        }  else if ( psc_vol_index > 0 ) {
            /* Only for Vibrant and Swift with PRESSCFH_NONE*/

            int cfh_slab_index=0;
            int cfh_switch_coil=0;
            /* First check for shim, then slab */
            if( oppscvquant >= 1 )
            {
                /*1 shim vol: No coil switch. Disregard psc_vol_index*/
                /*2 shim vol: No coil switch.*/
                cfh_switch_coil = 0;

                if (oppscvquant == 1) 
                    cfh_slab_index = 0;
                else 
                    cfh_slab_index = noswitch_slab_psc ? (PStest_slab-1):(psc_vol_index-1);

                rsp_PStloc = rsp_psc_info[cfh_slab_index].rsppsctloc;
                rsp_PSrloc = rsp_psc_info[cfh_slab_index].rsppscrloc;
                rsp_PSphasoff = rsp_psc_info[cfh_slab_index].rsppscphasoff;
            } else { /*No shim vol */ 
                cfh_slab_index = psc_vol_index - 1;
                if ( opswift == PSD_ON ) {
                    cfh_switch_coil = 1;
                    cfh_slab_index = noswitch_slab_psc ? (PStest_slab-1):(psc_vol_index-1);
                }
                rsp_PStloc = rsp_info[cfh_slab_index].rsptloc;
                rsp_PSrloc = rsp_info[cfh_slab_index].rsprloc;
                rsp_PSphasoff = rsp_info[cfh_slab_index].rspphasoff;
            }

            if( cfh_switch_coil )
            {
                if( old_psc_vol_index != psc_vol_index )
                {
                    if( noswitch_coil_psc == PSD_OFF )
                    {
                        if (FAILURE == CoilSwitchSetCoil(coilInfo_tgt[psc_vol_index-1], 0))
                        {
                            return FAILURE;
                        }
                        boffset(off_seqcfh); 
                    }

                    {
                        char tempstr[200]={0};
                        sprintf(tempstr,"PScfh: psc_vol_index=%d, cfh_slab_index=%d\n",psc_vol_index,cfh_slab_index);
                        printdbg(tempstr,PSdebugstate);
                    }
                    old_psc_vol_index = psc_vol_index;
                }
            }
            slice_freq = ( GAM * a_gzrf1cfh * rsp_PStloc
                           / (10 * TARDIS_FREQ_RES) );
            if((opcoax != 0) && cfh_newmode) {
                slice2_freq = ( GAM * (opspf ? rsp_PSphasoff * a_gyrf2cfh
                                       : rsp_PSrloc * a_gxrf2cfh )
                                / (10 * TARDIS_FREQ_RES) );
            } else {
                slice2_freq = 0;
            }

            if( PSD_ON == PSir ) {
                ir_slice_freq = ( GAM * a_gzrf0cfh * rsp_PStloc
                                  / (10 * TARDIS_FREQ_RES) );
                setfrequency(ir_slice_freq, &rf0cfh,0);
            }
        } else {

            /* This is for the 3D scan volume - psc_vol_index = 0 */
            rsp_PStloc = PStloc;
            rsp_PSrloc = PSrloc;
            rsp_PSphasoff = PSphasoff;

            slice_freq = cfh_rf1freq;
            slice2_freq = cfh_rf2freq;
            if( PSD_ON == PSir ) {
                ir_slice_freq = GAM * a_gzrf0cfh * PStloc/(10 * 
                                                           TARDIS_FREQ_RES); 
                setfrequency(ir_slice_freq, &rf0cfh,0);
            }
        }  /* end of psc_vol_index */

        setfrequency(slice_freq, &rf1cfh,0);
        setfrequency((int)slice2_freq, &rf2cfh,0);
        if( presscfh_ctrl != PRESSCFH_NONE )
        {
            setfrequency((int)slice3_freq, &rf3cfh,0); /* for presscfh_ctrl */
        }

        setfrequency((int)(PSfreq_offset[rspent]/TARDIS_FREQ_RES), &cfh_fid, 0);

        if(PSdebugstate)
        {
            INT pscrot_idx;

            pscrot_idx = psc_vol_index;
            if (psc_vol_index > 0)
            {
                pscrot_idx = psc_vol_index-1;
            }
            printf("\n %d rsp_PStloc (tx, rec, phase) -> %0f %0f %0f\n",
                   psc_vol_index, rsp_PStloc, rsp_PSrloc, rsp_PSphasoff);
            printf ("CFH : Rotation Matrix for Prescan Volume Index %d\n", psc_vol_index);
            printf("\t %6ld %6ld %6ld \n", rsp_PSrot[pscrot_idx][0], rsp_PSrot[pscrot_idx][1], rsp_PSrot[pscrot_idx][2]);
            printf("\t %6ld %6ld %6ld \n", rsp_PSrot[pscrot_idx][3], rsp_PSrot[pscrot_idx][4], rsp_PSrot[pscrot_idx][5]);
            printf("\t %6ld %6ld %6ld \n", rsp_PSrot[pscrot_idx][6], rsp_PSrot[pscrot_idx][7], rsp_PSrot[pscrot_idx][8]);

            fflush(stdout);
        }


        /* INVERSION RECOVERY CODE */
        if (PSD_ON == PSir)
        {
            /* Check for IR pulse being turned off */
            if (noir==1)
            {
                setiamp(0,&rf0cfh,0);
            }
            else
            {
                setiamp(ia_rf0cfh,&rf0cfh,0);
            }
            /* End check for IR pulse being on/off */

            if (tislice != (PSslice_num+1))
            {
                /* user has changed slice, so change frequencies */
                /* Check for proper range */
                if ((tislice < 1) || (tislice > opslquant)) 
                {
                    tislice=tislice_start;
                    PSslice_num=tislice-1;
                } 
                else
                {
                    PSslice_num=tislice-1;
                    /* changing the slice */

                    /* Calculation of new slice loc */
                    new_slice_loc = PStloc - ((opslspace + opslthick) *
                                              (tislice_start - tislice));
                    /* Inversion */
                    /* factor 10 is because rsptloc is in mm */
                    ir_slice_freq = GAM * a_gzrf0cfh * new_slice_loc/(10 * 
                                                                      TARDIS_FREQ_RES); 

                    setfrequency(ir_slice_freq, &rf0cfh,0);

                    /* Spin Echo */
                    /* factor 10 is because rsptloc is in mm */
                    slice_freq = GAM * a_gzrf1cfh * new_slice_loc/(10 * 
                                                                   TARDIS_FREQ_RES); 

                    setfrequency(slice_freq, &rf1cfh,0);
                    /* end of changing the slice */
                }
            }

            /* Changing the inversion time realtime */
            if (titime<50)
            {
                titime=50;
            }

            /* Need this check to reset titime to opti. Download error were
               happening because of invalid titimes */
            if(titime > 300)
            {
                titime = cfh_ti/1000;
            }

            titime_us = (titime*1000);
            new_dur=dur_invse + (titime_us - cfh_ti);
            new_dur = RUP_GRD(new_dur);
            setperiod(new_dur,&zticfh,0);
            setperiod(new_dur,&xticfh,0);
            setperiod(new_dur,&yticfh,0);
            setperiod(new_dur,&rticfh,0);
            setperiod(new_dur,&sticfh,0);

            /* Change the deadtime of seqcfh to preserve the cfh_tr value */
            /* If the seqcfh deadtime < 5 then set to the minimum of 4us */
            new_cfh_deadtime = init_cfh_deadtime + (cfh_ti - titime_us);
            if (new_cfh_deadtime < 5)
            {
                new_cfh_deadtime = 4;
            }

            setperiod( new_cfh_deadtime, &seqcfh, 0 );
            /* End deadtime change */

            /* End inversion time change */

            /* END INVERSION CORE CODE */
        } 
#ifdef PSD_CFH_CHEMSAT
        if( (cs_sat == PSD_ON) && PScs_sat )
        {
            CsSatMod((int)(cscfh_satindex+1));
        }
#endif
%ifdef PSD_CFH_MT
        if( (opmt == PSD_ON) && PSmt )
        {
            MTMod((int)(mtcfh_index+1));
        }
%endif /* PSD_CFH_MT */
        /* Check should play cs/mt or stir for manual cfh */
        if( PSir )
        {
            StIRMod();
        }

        for( excitation = 1; excitation <= cfh_nex; excitation ++ )
        {
            if( view > 0 ) 
            {
                acq_type = (int)DABON;

                if( excitation == cfh_nex )
                {
                    attenlockoff(&cfh_attenkey);
                }
                else
                {
                    attenlockon(&cfh_attenkey);
                }
            } 
            else
            {
                acq_type = (int)DABOFF;
            }

            /* DAB packet is not used in prescan.  loaddab
               is used just to turn on or off the receiver */
            loaddab(&cfh_fid, (int)1,(int)1,(int)1,(int)1,
                    (TYPDAB_PACKETS)acq_type, PSD_LOAD_DAB_ALL);

            /* Makesure psc_vol_index is in range */
            if( (psc_vol_index > 0 && psc_vol_index <= opvquant) && (opvquant > 1)) 
            {
                if(noswitch_slab_psc == PSD_OFF)
                    startseq((short)(psc_vol_index-1), (SHORT)MAY_PAUSE);
                else
                    startseq((short)(PStest_slab-1),(SHORT)MAY_PAUSE);
            }
            else
                startseq((short)0, (SHORT)MAY_PAUSE);

            /*printdbg("S", PSdebugstate);*/
            syncoff(&seqcfh); 

            getiamp(&temp_short, &rf1cfh, 0);
            setiamp(-temp_short, &rf1cfh, 0);

            if (view < 1)
            {
                break; /* Skip excitation loop for disdaqs */
            }
        }  /* end excitation */
    }  /* end view */

    return SUCCESS;
}   /* end PScfh() */


/* disable STIR cfh if cs/mttun is on */
void
StIRMod(void)
{
    int do_mttun = 0;
    int do_cstun = 0;

%ifdef PSD_CFH_MT
    do_mttun = mttun;
%endif /* PSD_CFH_MT */

#ifdef PSD_CFH_CHEMSAT
    if(PSD_ON == PScs_sat)
    {
        do_cstun = cstun;
    }
    else
    {
        do_cstun = PSD_OFF;
    }
#endif /* PSD_CFH_CHEMSAT */

    if( do_cstun || do_mttun || showfp ) 
    {
        /* disable stir cfh pulse */
        rfoff(&rf0cfh, 0);
        setiampt(0, &gyrf0kcfh, 0);
    }
    else {
        rfon(&rf0cfh, 0);
        setiampt(amp_gyrf0kcfh, &gyrf0kcfh, 0);
    }
    return;
}
%endif /* !SPECTROSCOPY */

/*  begin aps1_mod changes (GE) */
/*
 *  PSinit
 *  
 *  Type: Public Function
 *  
 *  Description:
 *  
 */
STATUS
PSinit( long (*PSrotmat)[9] )
{   
    /* Range check for psc_vol_index */
    if ( psc_vol_index < 0 )
    {
        if(PSdebugstate)
            printf("WRONG psc_vol_index =%d\n",psc_vol_index);
        psc_vol_index = 0;
    }

    setrfconfig((short) 5);	/* only activate rho1 */
    setssitime(100);	/* set ssi counter to 400 us. */
    rspqueueinit(200);	/* initialize to 200 entries */

    if( presscfh_cgate && rspent == L_CFH ){
        PSrsptrigger[0] = TRIG_ECG;
    } else{
        PSrsptrigger[0] = PStrigger;
    }
    setrotatearray((short)1, *PSrotmat);
    settriggerarray((short)1, PSrsptrigger);
    
    /* Always use scan rot matrix for SWIFT*/
    if(opswift == PSD_ON)
    {
        int i=0;
        setrotatearray((short)opvquant,rsprot[0]);
        for(i=0;i<opvquant;i++)
            PSrsptrigger[i] = PStrigger;

        settriggerarray((short)opvquant,PSrsptrigger);
    }

    return SUCCESS;
}   /* end PSinit() */

/* end aps1_mod changes (GE) */

/*
 *  PSfasttg
 *  
 *  Type: Public Function
 *  
 *  Description:
 *  
 */
STATUS
PSfasttg( void )
{
    printdbg("Greetings from FastTG", debugstate);
    rspent = L_FTG;
    rspdda = ftg_dda;
    rspbas = 0;
    rspvus = 30000;
    rspgy1 = 0;
    rspnex = 1;
    rspesl = pre_slice;
    rspasl = 0;
    rspslq = 1;
    rspsct = -1;

    strcpy(psdexitarg.text_arg, "FastTG");
    
    FastTGCore( PStloc_mod,
                (int)rspdda,
                (int)rspvus,
                (int)rspnex,
                (int)debugstate);
 
    printdbg("Normal End of FastTG", debugstate);
    rspexit();

    return SUCCESS;
}   /* end PSfasttg() */


/*
 *  PSexpresstg
 *  
 *  Type: Public Function
 *  
 *  Description:
 *  
 */
STATUS
PSexpresstg( void )
{
    printdbg("Greetings from eXpress TG", debugstate);
    rspent = L_XTG;
    rspdda = xtg_dda;
    rspbas = 0;
    rspvus = 30000;
    rspgy1 = 0;
    rspnex = 1;
    rspesl = pre_slice;
    rspasl = 0;
    rspslq = 1;
    rspsct = -1;

    strcpy(psdexitarg.text_arg, "eXpressTG");
    
    eXpressTGCore( PStloc_mod,
                   (int)rspdda,
                   (int)rspvus,
                   (int)rspnex,
                   (int)debugstate);
 
    printdbg("Normal End of eXpressTG", debugstate);
    rspexit();

    return SUCCESS;
}   /* end PSexpresstg() */


/*
 *  FastTGCore
 *  
 *  Type: Public Function
 *  
 *  Description:
 *    slice_loc: location of prescan slice  in mm
 *    slice_num: slice number to be excited
 *    ftg_disdaqs: number of disdaq pairs
 *    ftg_views: # of max views to be executed in ftg
 *    ftg_nex: # of excitations in ftg
 *    ftg_chop: No chop if = 2
 *    ftg_debug: debug state
 */
STATUS
FastTGCore( DOUBLE slice_loc,
            INT ftg_disdaqs,
            INT ftg_views,
            INT ftg_nex,
            INT ftg_debug )
{
    INT slice_freq1;    /* transmit frequency for the prescan slices */
    INT slice_freq2;
    INT slice_freq3;
    INT rcv_freq;
    float rsp_PStloc=0;
    float rsp_PSrloc=0;
    int old_psc_vol_index = -1; /* Last prescan volume for SWIFT switching */

    if (FTGacq1 == 1)
    {
        ftg_acq1 = (int)DABON;
    }
    else
    {
        ftg_acq1 = (int)DABOFF;
    }
      
    if (FTGacq2 == 1)
    {
        ftg_acq2 = (int)DABON;
    }
    else
    {
        ftg_acq2 = (int)DABOFF;
    }
      
    printdbg("Entering FastTGCORE", (SHORT)ftg_debug);
    boffset(off_seqftg);
    scopeon(&seqftg);
    syncon(&seqftg);
 
    /* Initialize */
    rsp_PStloc = slice_loc;
    rsp_PSrloc = PSrloc_mod;

    attenlockoff(&ftg_attenkey);
    slice_freq1 = GAM * a_gzrf1ftg * rsp_PStloc/
        (10* TARDIS_FREQ_RES); /* factor 10 is because rsptloc is in mm */
    slice_freq2 = GAM * a_gzrf2ftg * rsp_PStloc/
        (10* TARDIS_FREQ_RES); /* factor 10 is because rsptloc is in mm */
    slice_freq3 = GAM * a_gzrf3ftg * rsp_PStloc/
        (10* TARDIS_FREQ_RES); /* factor 10 is because rsptloc is in mm */
    rcv_freq = 2*1000*FTGecho1bw*rsp_PSrloc / FTGfov;

    setfrequency(slice_freq1, &rf1ftg,0);
    setfrequency(slice_freq2, &rf2ftg,0);
    setfrequency((INT)((rcv_freq + PSfreq_offset[rspent])/TARDIS_FREQ_RES),
                 &echo1ftg, 0);
    setfrequency(slice_freq3, &rf3ftg,0);
    setfrequency((INT)((rcv_freq + PSfreq_offset[rspent])/TARDIS_FREQ_RES),
                 &echo2ftg, 0);
 
    if(PSdebugstate)
    {
        printf("\nFTG Slthick = %f FOV = %f Xmit location = %f, Rcv location  = %f\n", 
                FTGslthk,FTGfov, slice_loc, PSrloc_mod);
        printf("\nFTG Xmit Freq = %f, Rcv Freq  = %f\n", (float)slice_freq1, (float)rcv_freq);
        printf("%s\n","FTG Rotation Matrix:");
        printf("\t %6ld %6ld %6ld \n", PSrot_mod[0][0], PSrot_mod[0][1], PSrot_mod[0][2]);
        printf("\t %6ld %6ld %6ld \n", PSrot_mod[0][3], PSrot_mod[0][4], PSrot_mod[0][5]);
        printf("\t %6ld %6ld %6ld \n", PSrot_mod[0][6], PSrot_mod[0][7], PSrot_mod[0][8]);
        fflush(stdout);
    }

    if( PSdebugstate )
    {
        printf("FTG: psc_vol_index = %d\n",psc_vol_index );
    }

    for (view = (-ftg_disdaqs + 1) ; view <= ftg_views; view ++)
    {
        for (excitation = 1; excitation <= ftg_nex; excitation ++)
        {
            if (view > 0)
            {
                if (excitation == 1)
                    attenlockon(&ftg_attenkey);
                if (excitation == ftg_nex)
                    attenlockoff(&ftg_attenkey);
                ftg_acq2 = (int)DABON;
            }
            else
            {
                ftg_acq2 = (int)DABOFF;
            }
 
            /* DAB packet is not used in prescan.  loaddab
               is used just to turn on or off the receiver */
            loaddab(&echo1ftg, (INT)1,(INT)1,(INT)1,(INT)1,
                    (TYPDAB_PACKETS)ftg_acq1, PSD_LOAD_DAB_ALL);
            loaddab(&echo2ftg, (INT)1,(INT)1,(INT)1,(INT)1,
                    (TYPDAB_PACKETS)ftg_acq2, PSD_LOAD_DAB_ALL);

            /* For SWIFT(fasttg): We are trying not to use the CV opswift.
               Instead use psc_vol_index.
               */
            if( (psc_vol_index > 0) && (opvquant > 1) )
            {
                if( old_psc_vol_index != psc_vol_index )
                {
                    /* No need to switch coil*/
                    if( noswitch_slab_psc == PSD_OFF )
                    {
                        rsp_PStloc = rsp_info[psc_vol_index-1].rsptloc;
                        rsp_PSrloc = rsp_info[psc_vol_index-1].rsprloc;
                    }
                    else
                    {
                        rsp_PStloc = rsp_info[PStest_slab-1].rsptloc;
                        rsp_PSrloc = rsp_info[PStest_slab-1].rsprloc;
                    }

                    {
                        char tempstr[200]={0};
                        sprintf(tempstr,"fasttgcore: psc_vol_index=%d,rsp_PStloc=%f\n",psc_vol_index,rsp_PStloc);
                        printdbg(tempstr,PSdebugstate);
                    }
                    old_psc_vol_index = psc_vol_index;
                }
            }

            slice_freq1 = GAM * a_gzrf1ftg * rsp_PStloc/
                (10* TARDIS_FREQ_RES); /* factor 10 is because rsptloc is in mm */
            slice_freq2 = GAM * a_gzrf2ftg * rsp_PStloc/
                (10* TARDIS_FREQ_RES); /* factor 10 is because rsptloc is in mm */
            slice_freq3 = GAM * a_gzrf3ftg * rsp_PStloc/
                (10* TARDIS_FREQ_RES); /* factor 10 is because rsptloc is in mm */
            rcv_freq = 2*1000*FTGecho1bw*rsp_PSrloc / FTGfov;

            setfrequency(slice_freq1, &rf1ftg,0);
            setfrequency(slice_freq2, &rf2ftg,0);
            setfrequency((INT)((rcv_freq + PSfreq_offset[rspent])/TARDIS_FREQ_RES),
                    &echo1ftg, 0);
            setfrequency(slice_freq3, &rf3ftg,0);
            setfrequency((INT)((rcv_freq + PSfreq_offset[rspent])/TARDIS_FREQ_RES),
                    &echo2ftg, 0);


            /* Makesure psc_vol_index is in range */
            if( (psc_vol_index > 0 && psc_vol_index <= opvquant) && (opvquant > 1)) 
            {
                if(noswitch_slab_psc == PSD_OFF)
                    startseq((short)(psc_vol_index-1), (SHORT)MAY_PAUSE);
                else
                    startseq((short)(PStest_slab-1),(SHORT)MAY_PAUSE);
            }
            else
                startseq((short)0, (SHORT)MAY_PAUSE);

            printdbg("S", (SHORT)ftg_debug);
            syncoff(&seqftg);
 
 
            if (view < 1)
            {
                break; /* Skip excitation loop for disdaqs */
            }
 
        }  /* end excitation */
 
    }  /* end view */

    return SUCCESS;
}   /* end FastTGCore() */


/*
 *  eXpressTGCore
 *  
 *  Type: Public Function
 *  
 *  Description:
 *    slice_loc: location of prescan slice  in mm
 *    slice_num: slice number to be excited
 *    xtg_disdaqs: number of disdaq pairs
 *    xtg_views: # of max views to be executed in xtg
 *    xtg_nex: # of excitations in xtg
 *    xtg_chop: No chop if = 2
 *    xtg_debug: debug state
 */
STATUS
eXpressTGCore( DOUBLE slice_loc,
            INT xtg_disdaqs,
            INT xtg_views,
            INT xtg_nex,
            INT xtg_debug )
{
    INT slice_freq1;    /* transmit frequency for the prescan slices */
    INT slice_freq2;
    INT rcv_freq;
    SHORT temp_short;
    INT rf3xtg_freq, rf4xtg_freq;
    float rsp_PStloc=0;
    float rsp_PSrloc=0;
    int old_psc_vol_index = -1; /* Last prescan volume for SWIFT switching */

    if (XTGacq1 == PSD_ON)
    {
        xtg_acq1 = (int)DABON;
    }
    else
    {
        xtg_acq1 = (int)DABOFF;
    }
            
    printdbg("Entering eXpressTGCORE", (SHORT)xtg_debug);
    boffset(off_seqxtg);
    scopeon(&seqxtg);
    syncon(&seqxtg);
 
    /* Initialize */
    rsp_PStloc = slice_loc;
    rsp_PSrloc = PSrloc_mod;

    attenlockoff(&xtg_attenkey);
    slice_freq1 = GAM * a_gzrf1xtg * rsp_PStloc/
        (10* TARDIS_FREQ_RES); /* factor 10 is because rsptloc is in mm */
    slice_freq2 = GAM * a_gzrf2xtg * rsp_PStloc/
        (10* TARDIS_FREQ_RES); /* factor 10 is because rsptloc is in mm */
    rcv_freq = 2*1000*XTGecho1bw*rsp_PSrloc / XTGfov;

    setfrequency(slice_freq1, &rf1xtg,0);
    setfrequency(slice_freq2, &rf2xtg,0);
    setfrequency((INT)((rcv_freq + PSfreq_offset[rspent])/TARDIS_FREQ_RES),
                 &echo1xtg, 0);
 
    if(PSdebugstate)
    {
        printf("\nXTG Slthick = %f FOV = %f Xmit location = %f, Rcv location  = %f\n", 
                XTGopslthick, XTGfov, slice_loc, PSrloc_mod);
        printf("\nXTG Xmit Freq = %f, Rcv Freq  = %f\n", (float)slice_freq1, (float)rcv_freq);
        printf("%s\n","XTG Rotation Matrix:");
        printf("\t %6ld %6ld %6ld \n", PSrot_mod[0][0], PSrot_mod[0][1], PSrot_mod[0][2]);
        printf("\t %6ld %6ld %6ld \n", PSrot_mod[0][3], PSrot_mod[0][4], PSrot_mod[0][5]);
        printf("\t %6ld %6ld %6ld \n", PSrot_mod[0][6], PSrot_mod[0][7], PSrot_mod[0][8]);
        fflush(stdout);
        printf("XTG: psc_vol_index = %d\n",psc_vol_index );
    }

    rf4xtg_freq = (int)(xtg_offres_freq/TARDIS_FREQ_RES);
    rf3xtg_freq = rf4xtg_freq;

    for (view = (-xtg_disdaqs + 1) ; view <= xtg_views; view ++)
    {
        if( PSdebugstate )
        {
            printf("XTG: view = %d\n",view);
        }

        if(view == 1)
        {
            getiamp(&temp_short, &rf4xtg, 0);
            if(temp_short < 0.0)
            {
                setiamp(-1.0*temp_short, &rf4xtg, 0);
                getiamp(&temp_short, &rf3xtg, 0);
                setiamp(-1.0*temp_short, &rf3xtg, 0);
            }
            getiamp(&temp_short, &rf1xtg, 0);
            if(temp_short < 0.0)
            {
                setiamp(-1.0*temp_short, &rf1xtg, 0);
            }
        }

        for (excitation = 1; excitation <= xtg_nex; excitation ++)
        {
            if (view > 0)
            {
                if (excitation == 1)
                    attenlockon(&xtg_attenkey);
                if (excitation == xtg_nex)
                    attenlockoff(&xtg_attenkey);
                if(XTGacq1 == PSD_ON)
                {
                    xtg_acq1 = (int)DABON;
                }
            }
            else
            {
                xtg_acq1 = (int)DABOFF;
            }
 
            /* DAB packet is not used in prescan.  loaddab
               is used just to turn on or off the receiver */
            loaddab(&echo1xtg, (INT)1,(INT)1,(INT)1,(INT)1,
                    (TYPDAB_PACKETS)xtg_acq1, PSD_LOAD_DAB_ALL);

            /* For SWIFT(expresstg): We are trying not to use the CV opswift.
               Instead use psc_vol_index.
               */
            if( (psc_vol_index > 0) && (opvquant > 1) )
            {
                if( old_psc_vol_index != psc_vol_index )
                {
                    /* No need to switch coil*/
                    if( noswitch_slab_psc == PSD_OFF )
                    {
                        rsp_PStloc = rsp_info[psc_vol_index-1].rsptloc;
                        rsp_PSrloc = rsp_info[psc_vol_index-1].rsprloc;
                    }
                    else
                    {
                        rsp_PStloc = rsp_info[PStest_slab-1].rsptloc;
                        rsp_PSrloc = rsp_info[PStest_slab-1].rsprloc;
                    }

                    {
                        char tempstr[200]={0};
                        sprintf(tempstr,"eXpresstgcore: psc_vol_index=%d,rsp_PStloc=%f\n",psc_vol_index,rsp_PStloc);
                        printdbg(tempstr,PSdebugstate);
                    }
                    old_psc_vol_index = psc_vol_index;
                }
            }

            slice_freq1 = GAM * a_gzrf1xtg * rsp_PStloc/
                (10* TARDIS_FREQ_RES); /* factor 10 is because rsptloc is in mm */
            slice_freq2 = GAM * a_gzrf2xtg * rsp_PStloc/
                (10* TARDIS_FREQ_RES); /* factor 10 is because rsptloc is in mm */
            rcv_freq = 2*1000*XTGecho1bw*rsp_PSrloc / XTGfov;
            
            setfrequency(slice_freq1, &rf1xtg,0);
            setfrequency(slice_freq2, &rf2xtg,0);
            setfrequency((INT)((rcv_freq + PSfreq_offset[rspent])/TARDIS_FREQ_RES),
                         &echo1xtg, 0);

            /* phase and freq cycling  */
            if((view > 0) && (view%2 == 1)) /* freq cycling */
            {
                setfrequency(rf4xtg_freq, &rf4xtg, 0);
                setfrequency((-1.0*rf3xtg_freq), &rf3xtg, 0);
            }
            else
            {
                setfrequency((-1.0*rf4xtg_freq), &rf4xtg, 0);
                setfrequency(rf3xtg_freq, &rf3xtg, 0);
            }
            
            if((view > 1) && (view%2 == 1)) /* phase cycling */
            {
                getiamp(&temp_short, &rf4xtg, 0);
                setiamp(-1.0*temp_short, &rf4xtg, 0);
                getiamp(&temp_short, &rf3xtg, 0);
                setiamp(-1.0*temp_short, &rf3xtg, 0);
            }

            /* Makesure psc_vol_index is in range */
            if( (psc_vol_index > 0 && psc_vol_index <= opvquant) && (opvquant > 1)) 
            {
                if(noswitch_slab_psc == PSD_OFF)
                    startseq((short)(psc_vol_index-1), (SHORT)MAY_PAUSE);
                else
                    startseq((short)(PStest_slab-1),(SHORT)MAY_PAUSE);
            }
            else
                startseq((short)0, (SHORT)MAY_PAUSE);

            printdbg("S", (SHORT)xtg_debug);
            syncoff(&seqxtg);
 
            if (view < 1)
            {
                break; /* Skip excitation loop for disdaqs */
            }
 
            /* Chopper logic */
            getiamp(&temp_short, &rf1xtg,0);
            setiamp((-temp_short),&rf1xtg,0);

        }  /* end excitation */
 
    }  /* end view */

    return SUCCESS;
}   /* end eXpressTGCore() */


/*
 *  ASautoshim
 *  
 *  Type: Public Function
 *  
 *  Description:
 *  
 */
STATUS
ASautoshim( void )
{
    SHORT temp_short;
    INT dix_shift;

#ifdef OFFSET_AS
    INT phase_step[3];
    INT phase_sign[3];
    INT yres_phase;
#endif /* OFFSET_AS */

    long asrottemp[3][9]; /* vmx 3/6/95 YI */
    INT *rf1_freq;
    INT *receive_freq1;
    SHORT viewtable[513];
    long trigger_temp[3]; /* vmx 10/13/94 YI */
    INT acquire_echo1;
    INT dab_view,dab_op;
    SHORT disdaqs;
    FLOAT tempGAM;

    printdbg("Greetings from autoshim", debugstate);
    boffset(off_seqaushim);

    disdaqs = as_dda;
    tempGAM = GAM;
    GAM = GAMMA_PROTON; /*only shim on proton*/

    /**************************************************************
      so here's how we loop through this entry point:

        pass = rspslq;               each slice is a pass
        for all slices
           for all views              
             (first echo)           note: these aren't really 2echos,
             reset dixon shift to 0       we just pretend they are.
             do the disdaqs               there's really a separate
             do the baselines             excitation for each 'echo'.
             collect the data 'echo'
             (second echo)
             do dixon shift
             do the disdaqs
             do the baselines
             collect the data 'echo'
           next view
           decrement pass
           send pass packet
        next slice
        send end of scan packet
     ******************************************************************/

    setrfconfig((short) 5);

    /* Allocate memory for various arrays.
     * An extra 2 locations are saved in case the user wants to do
     * some tricks. */
    rf1_freq = (int *)AllocNode((as_slquant + 2)*sizeof(int));
    receive_freq1 = (int *)AllocNode((as_slquant + 2)*sizeof(int));

    rf1_freq[0] = astloc1*GAM * a_gzrf1as /(10 * TARDIS_FREQ_RES);
    rf1_freq[1] = astloc2*GAM * a_gzrf1as /(10 * TARDIS_FREQ_RES);
    rf1_freq[2] = astloc3*GAM * a_gzrf1as /(10 * TARDIS_FREQ_RES);

    asrottemp[0][0] = hostToRspRotMat( asrot0 );  asrottemp[0][1] = hostToRspRotMat( asrot1 ); asrottemp[0][2] = hostToRspRotMat( asrot2 );  
    asrottemp[0][3] = hostToRspRotMat( asrot3 );  asrottemp[0][4] = hostToRspRotMat( asrot4 ); asrottemp[0][5] = hostToRspRotMat( asrot5 );
    asrottemp[0][6] = hostToRspRotMat( asrot6 );  asrottemp[0][7] = hostToRspRotMat( asrot7 ); asrottemp[0][8] = hostToRspRotMat( asrot8 );

    asrottemp[1][0] = hostToRspRotMat( asrot9 );  asrottemp[1][1] = hostToRspRotMat( asrot10 ); asrottemp[1][2] = hostToRspRotMat( asrot11 ); 
    asrottemp[1][3] = hostToRspRotMat( asrot12 ); asrottemp[1][4] = hostToRspRotMat( asrot13 ); asrottemp[1][5] = hostToRspRotMat( asrot14 );
    asrottemp[1][6] = hostToRspRotMat( asrot15 ); asrottemp[1][7] = hostToRspRotMat( asrot16 ); asrottemp[1][8] = hostToRspRotMat( asrot17 );

    asrottemp[2][0] = hostToRspRotMat( asrot18 ); asrottemp[2][1] = hostToRspRotMat( asrot19 ); asrottemp[2][2] = hostToRspRotMat( asrot20 ); 
    asrottemp[2][3] = hostToRspRotMat( asrot21 ); asrottemp[2][4] = hostToRspRotMat( asrot22 ); asrottemp[2][5] = hostToRspRotMat( asrot23 );
    asrottemp[2][6] = hostToRspRotMat( asrot24 ); asrottemp[2][7] = hostToRspRotMat( asrot25 ); asrottemp[2][8] = hostToRspRotMat( asrot26 );

    scalerotmats(asrottemp, &asloggrd, &phygrd, 3, asobl_debug);

    if(PSdebugstate)
    {
        printf("\n%d astlocs (1,2,3) -> %0f %0f %0f\n",
               as_slquant, astloc1, astloc2, astloc3);
        printf("\nrf1_freq (1,2,3) -> %d %d %d\n",
               rf1_freq[0], rf1_freq[1], rf1_freq[2]);
    }

    /* AutoShim Changes to Center Image always - HH- Sept 21, 2004 */
    /* Use asrot0 , asrot11 and asrot20 to decide on read and phase directions on the 3 planes */
    /* Slice 1 = Axial    - Z slice (astloc1) - X/Y read (astloc2/astloc3) - Y/X Phase (astloc3/astloc2) */
    /* Slice 2 = Sagittal - X slice (astloc2) - Z/Y read (astloc1/astloc3) - Y/Z Phase (astloc3/astloc1) */
    /* Slice 3 = Coronal  - Y slice (astloc3) - Z/X read (astloc1/astloc2) - X/Z Phase (astloc2/astloc1) */

/* Old Code */
    receive_freq1[0] = (PSfreq_offset[rspent] +
                        ((2 * echo1bwas*1000 / (asfov))
                         * asrloc1)) / TARDIS_FREQ_RES;
    receive_freq1[1] = (PSfreq_offset[rspent] +
                        ((2 * echo1bwas*1000 / (asfov))
                         * asrloc2)) / TARDIS_FREQ_RES;
    receive_freq1[2] = (PSfreq_offset[rspent] +
                        ((2 * echo1bwas*1000 / (asfov))
                         * asrloc3)) / TARDIS_FREQ_RES;

#ifdef OFFSET_AS
    /* New code Begin */

    receive_freq1[0] = (PSfreq_offset[rspent]
                        + ((2 * echo1bwas*1000 / (asfov)) * 
                           (asrottemp[0][0] != 0) ? astloc2 : astloc3)) / TARDIS_FREQ_RES;
    receive_freq1[1] = (PSfreq_offset[rspent]
                        + ((2 * echo1bwas*1000 / (asfov)) * 
                           (asrottemp[1][2] != 0) ? astloc1 : astloc3)) / TARDIS_FREQ_RES;
    receive_freq1[2] = (PSfreq_offset[rspent]
                        + ((2 * echo1bwas*1000 / (asfov)) * 
                           (asrottemp[2][2] != 0) ? astloc1 : astloc2)) / TARDIS_FREQ_RES;

    phase_step[0] = .5 + fabs(FS_2PI*(asrottemp[0][0] != 0)? astloc3 : astloc2 / asfov);
    phase_step[1] = .5 + fabs(FS_2PI*(asrottemp[1][2] != 0)? astloc3 : astloc1 / asfov);
    phase_step[2] = .5 + fabs(FS_2PI*(asrottemp[2][2] != 0)? astloc2 : astloc1 / asfov);
 
    if ( ((asrottemp[0][0] != 0) ? astloc3 : astloc2) >= 0.0 )
    {
        phase_sign[0] = 1.0;
    }
    else
    {
        phase_sign[0] = -1.0;
    }

    if ( ((asrot[1][2] != 0) ? astloc3 : astloc1) >= 0.0 )
    {
        phase_sign[1] = 1.0;
    }
    else
    {
        phase_sign[1] = -1.0;
    }

    if ( ((asrottemp[2][2] != 0) ? astloc2 : astloc1) >= 0.0 )
    {
        phase_sign[2] = 1.0;
    }
    else
    {
        phase_sign[2] = -1.0;
    }

    /* New code End */
#endif /* OFFSET_AS */

    trigger_temp[0] = TRIG_INTERN;
    trigger_temp[1] = TRIG_INTERN;
    trigger_temp[2] = TRIG_INTERN;

    setupphasetable(viewtable, TYPNORM,(int)asyres);


    /* Set ssi time.  This is time from eos to start of sequence   */
    /* interrupt in internal triggering.  The minimum time is 50us */
    /* plus 2us*(number of waveform and instruction words modified */
    /* in the update queue).                                       */
    setssitime((LONG)time_ssias/HW_GRAD_UPDATE_TIME);

    settriggerarray(as_slquant, trigger_temp);

    if(PSdebugstate)
    {
        printf("\n AUTOSHIM:  FOV = %f astloc1 = %f, astloc2 = %f  astloc3  = %f\n", asfov, astloc1,astloc2,astloc3); 
        printf("\n AUTOSHIM:  asrloc1 = %f, asrloc2 = %f  asrloc3  = %f\n", asrloc1,asrloc2,asrloc3); 
        printf("%s\n", "AUTOSHIM:  Slice 1 Rotation Matrix");
        printf("\t %6ld  %6ld  %6ld\n", asrottemp[0][0],asrottemp[0][1],asrottemp[0][2]);
        printf("\t %6ld  %6ld  %6ld\n", asrottemp[0][3],asrottemp[0][4],asrottemp[0][5]);
        printf("\t %6ld  %6ld  %6ld\n", asrottemp[0][6],asrottemp[0][7],asrottemp[0][8]);
        printf("%s\n", "AUTOSHIM:  Slice 2 Rotation Matrix");
        printf("\t %6ld  %6ld  %6ld\n", asrottemp[1][0],asrottemp[1][1],asrottemp[1][2]);
        printf("\t %6ld  %6ld  %6ld\n", asrottemp[1][3],asrottemp[1][4],asrottemp[1][5]);
        printf("\t %6ld  %6ld  %6ld\n", asrottemp[1][6],asrottemp[1][7],asrottemp[1][8]);
        printf("%s\n", "AUTOSHIM:  Slice 3 Rotation Matrix");
        printf("\t %6ld  %6ld  %6ld\n", asrottemp[2][0],asrottemp[2][1],asrottemp[2][2]);
        printf("\t %6ld  %6ld  %6ld\n", asrottemp[2][3],asrottemp[2][4],asrottemp[2][5]);
        printf("\t %6ld  %6ld  %6ld\n", asrottemp[2][6],asrottemp[2][7],asrottemp[2][8]);

#ifdef OFFSET_AS
        printf("\n AUTOSHIM SLICE1:   slice loc= %f read loc= %f phase loc= %f \n", 
               astloc1, (asrottemp[0][0] != 0)   astloc2 : astloc3,
               (asrottemp[0][0] != 0) ? astloc3 : astloc2); 
        printf("\n AUTOSHIM SLICE2:   slice loc= %f read loc= %f phase loc= %f \n", 
               astloc1, (asrottemp[1][2] != 0) ? astloc1 : astloc3,
               (asrottemp[1][2] != 0) ? astloc3 : astloc1); 
        printf("\n AUTOSHIM SLICE3:   slice loc= %f read loc= %f phase loc= %f \n", 
               astloc1, (asrottemp[2][2] != 0) ? astloc1 : astloc2,
               (asrottemp[2][2] != 0) ? astloc2 : astloc1); 
#endif /* OFFSET_AS */
        fflush(stdout);
    }

    setrotatearray(as_slquant,asrottemp[0]);

    /* Reset scope and attenuator lock*/
    attenlockon(&attenuator_keyas);
    scopeon(&seqaushim);

    /* Reset all the synchronizations  - no need to use one in pass */
    syncon(&seqaushim);
    syncoff(&seqpassas);


    for (slice = 0; slice < as_slquant; slice++)
    {
        dab_op = DABSTORE;

        setiamp((short)(ia_rf1as*flip_pctas), &rf1as, 0);

        for (view=(-(disdaqs+asbaseline)+1); view<= asyres; view++)
        {
            if (view<=0)
            {
                /* for the baselines in each slice, set the phase encode
                   amplitude to the first view */
                /* HSI change to setiampt from setiamp */
                setiampt(-viewtable[1], &gy1as, 0);
                setiampt(viewtable[1], &gy1ras, 0);
#ifdef OFFSET_AS
                /* New code Begin */

                setiphase(0,&echo1as,0);

                /* New code End */
#endif /* OFFSET_AS */
            }
            else
            {
                if (view > 0) 
                {
                    /* HSI change to setiampt from setiamp */
                    setiampt(-viewtable[view], &gy1as, 0);
                    setiampt(viewtable[view], &gy1ras, 0);
#ifdef OFFSET_AS

                    /* New code Begin */

                    yres_phase = -phase_sign[slice]*(((view-1)* phase_step[slice] + 3L*FS_PI)%FS_2PI-FS_PI);
                    setiphase(yres_phase, &echo1as, 0);

                    /* New code End */
#endif /* OFFSET_AS */
                }
            }
            for (excitation=1; excitation <= 1; excitation++)
            {
                /* Scope Trigger */
                if ((rspsct == slice) || (rspsct == -1))
                {
                    scopeon(&seqaushim);
                }
                else
                {
                    scopeoff(&seqaushim);
                }
                /*MRIge33520 - VB*/
                if ((view > -asbaseline)&&(excitation > 0))
                {
                    acquire_echo1 = (int)DABON;
                }
                else
                {
                    acquire_echo1 = (int)DABOFF;
                }

                /* Load Transmit and Receive frequencies */
                setfrequency(rf1_freq[slice], &rf1as, 0);
                setfrequency(receive_freq1[slice], &echo1as, 0);

                if (view > 0)
                {
                    dab_view = view;
                    if (excitation == 1)
                    {
                        dab_op = 0;
                    }
                    else
                    {
                        dab_op = 3 - 2*(excitation % 2);
                    }
                }
                else
                {
                    dab_view = 0;
                }

                /* set up the dixon shift for every other seq */
                for(dixon=0;dixon<=1;++dixon)
                {
                    dix_shift = (dixon % 2)*RUP_GRD((dix_timeas))+GRAD_UPDATE_TIME;
                    setperiod(dix_shift,&xdixon,0);
                    setperiod(dix_shift,&ydixon,0);
                    setperiod(dix_shift,&zdixon,0);
                    setperiod(dix_shift,&sdixon,0);
                    setperiod((int)(dix_shift - dixon*dix_timeas),&sdixon2,0);
                    /* All DAB Info is Set. Pretend even views are 2nd echo. */
                    /* Load up dab packet!                                   */
                    loaddab(&echo1as,(short)0,dixon,dab_op,dab_view,(TYPDAB_PACKETS)acquire_echo1, PSD_LOAD_DAB_ALL);

                    /*printdbg("S", debugstate);*/
                    startseq((short)slice, (short)MAY_PAUSE);

                    syncoff(&seqaushim);
                }

                if (view == (-asbaseline+1))
                    dab_op = 1; /* add baseviews */

                getiamp(&temp_short, &rf1as, 0);
                setiamp(-temp_short, &rf1as, 0);

            } /* excitation */

        }  /* view */


        boffset(off_seqpassas);
        if (slice == (as_slquant-1)  ) /* Last pass */
        {
            /* Set DAB pass packet to end of scan */
            setwamp(SSPD + DABPASS + DABSCAN, &pass_aushim, 2);
            printdbg("End of Scan and Pass", debugstate);
        }
        else
        {
            /* Set DAB pass packet to end of pass */
            setwamp(SSPD + DABPASS, &pass_aushim, 2);
            printdbg("End of Pass", debugstate);
        }


        startseq((short)0, (short)MAY_PAUSE);
        boffset(off_seqaushim);

    }

    /* Reset the rotation matrix */
    setrotatearray((short)opslquant,rsprot[0]);

    GAM = tempGAM;

    printdbg("Normal End of autoshim", debugstate);

    return SUCCESS;
} /* End of ASautoshim */

/* CoilSwitchSetCoil
 *
 *  Description: MRIhc15304
 *    This sets the RF HUB index for the coil by changing the data on an ssp
 *    pulse and/or with the sethubindeximm function.  Also calls
 *    setrcvportimm, if requested.
 *
 *  Parameters:
 *  (I: for input parameters, O: for output parameters)
 *  
 *  (O) STATUS return variable - Did function complete successfully.
 *  (I) const COIL_INFO - coil info structure of coil to switch to
 *  (I) const INT setRcvPortFlag - Flag indicating that setrcvportimm needs
 *        to be executed.  This needs to be set whenever switching to or
 *        from the BODY coil.
 *  
 *  Globals:
 *  (I) txCoilInfo
 *
 */

STATUS
CoilSwitchSetCoil( const COIL_INFO coil,
                   const INT setRcvPortFlag)
{
    SHORT device = 0;

    if( setRcvPortFlag || (COIL_SWITCH_RSP_SETHUBINDEXIMM & cfcoilswitchmethod) )
    {
        int wait_rspimm = 0;

        /* MRIhc14300: When switching coils, play a delay sequence to
           wait for scan prep to complete and scanning to start before the
           first setrcvportimm() & sethubindeximm().  These functions may not
           be called until scanning starts and the first startseq will not
           return until scanning starts. After scanning has started this will
           merely add an additional short delay to the switch time */

        boffset(off_seqcswWaitBefore);
        startseq((short)0, (SHORT)MAY_PAUSE);

        /* Need additional delay for setrcvpowerimm & sethubindeximm
         * to take effect.  The delay time must be set to guarantee
         * the completion of both setrcvportimm and sethubindeximm */

        if( COIL_SWITCH_RSP_SETHUBINDEXIMM & cfcoilswitchmethod )
        {
            wait_rspimm = csw_wait_sethubindeximm;
        }

        if( setRcvPortFlag )
        {
            wait_rspimm = IMax(2, wait_rspimm, csw_wait_setrcvportimm);
        }

        setperiod(wait_rspimm, &csw_wait, 0);
    }
    else
    {
        /* No additional delay needed when not calling setrcvportimm */
        setperiod(SSP_UPDATE_TIME, &csw_wait, 0);
    }

    /* Setup hub index switching packet */
    device = 0;
    if( COIL_SWITCH_SSP_HUB_INDEX & cfcoilswitchmethod )
    {
        device = RDC;
        /* Set hub index on SSP packet */
        setwamp( (SHORT)(SSPD | (HUBIND + coil.hubIndex)), &contrfhubsel, (LONG)2 );
    }
    setwamp( (SHORT)(SSPDS | device), &contrfhubsel, (LONG)0 );

    /* Set receiver port & receiver input */
    if( setRcvPortFlag || (COIL_SWITCH_SSP_RECEIVER_INPUT & cfcoilswitchmethod) )
    {
        SHORT coil_port = PSD_RP_BODY;
        SHORT rcv_input = RFBODYI;

        /* Calculate coil port & receiver input */
        switch (coil.rxCoilType) 
        {
        case RX_COIL_BODY:
        default:
            /* Transmit & receive with body coil */
            coil_port = PSD_RP_BODY;
            rcv_input = RFBODYI;
            break;
        case RX_COIL_LOCAL:
            {
                /* Assume there is only one transmit coil.  If two transmit
                   coils, the primary will be used */
                n32 txCoilType = TX_INDEX_NONE; 
                if (TX_INDEX_NONE != coil.txIndexPri)
                {
                    txCoilType = txCoilInfo_tgt[coil.txIndexPri].txCoilType;
                }
                else if (TX_INDEX_NONE != coil.txIndexSec)
                {
                    txCoilType = txCoilInfo_tgt[coil.txIndexSec].txCoilType;
                }
                else
                {
                    printf("CoilSwitchSetCoil: No transmit coil defined!\n");
                    return FAILURE;
                }

                if (TX_COIL_LOCAL == txCoilType)
                {
                    /* Local transmit coil */
                    coil_port = PSD_RP_HEAD;
                    rcv_input = RHEADI;
                }
                else
                {
                    /* Surface coil */
                    coil_port = PSD_RP_SURFACE;
                    rcv_input = RFAUX;
                }
            }
            break;
        }

        /* Setup receiver input switching packet */
        device = 0;
        if( COIL_SWITCH_SSP_RECEIVER_INPUT & cfcoilswitchmethod )
        {
            device = RDC;
            setwamp( (SHORT)(SSPD | rcv_input), &contrfsel, (LONG)2 );
        }
        setwamp( (SHORT)(SSPDS | device), &contrfsel, (LONG)0 );

        /* Set receiver port using RSP function call when switching to/from 
         * body coil */
        if (setRcvPortFlag) 
        {
#ifdef PSD_HW
            setrcvportimm( (SHORT)coil_port );
#endif /* PSD_HW */
        }
    }

    /* Select coil using RSP function on MGD Rx chain */
    if(COIL_SWITCH_RSP_SETHUBINDEXIMM & cfcoilswitchmethod)
    {
        sethubindeximm( coil.hubIndex );
    }
   
    boffset( off_seqcsw );
    startseq( (short)0, (SHORT)MAY_PAUSE );
    return SUCCESS;
}

@rspvar PSrspvar
/*********************************************************************
 *                    PRESCAN.E RSPVAR SECTION                       *
 *                             PSrspvar                              *
 *                                                                   *
 * Declare here the real time variables that can be viewed and modi- *
 * fied while the Tgt PSD process is running. Only limited standard  *
 * C types are provided: short, int, long, float, double, and 1D     *
 * arrays of those types.                                            *
 *                                                                   *
 * NOTE: Do not declare all real-time variables here because of the  *
 *       overhead required for viewing and modifying them.           *
 *********************************************************************/

float dipir_ratio;
int dixon,slice,view,excitation;
int dur_invse; /* Initial duration for the WAIT pulse before the 90 pulse */
               /* based on OPTI */
int new_dur;   /* New duration time calulated with titime and cfh_ti delta */
int titime; /* input TI time on milliseconds; initialized to OPTI */
int titime_us;  /* input TI time in microseconds */
int tislice; /* holds new slice value that is changed during CFH */
int tislice_start; /* initial slice number passed into CFH */
float new_slice_loc; /* holds new slice location in cm */
int noir; /* Flag for IR pulse to be turned off */
int cscfh_satindex = 1; /* default chemsat occurence in cfh */
int mtcfh_index = 1; /* default mt occurence in cfh */
int amp_gyrf0kcfh;
int psc_vol_index; /* index of the prescan volume */
int showfp;

/****  FastTG RSPvar    ****/
int ftg_acq1, ftg_acq2; /* flags for data acquisiton, windows 1 and 2
                           1=on, 0=off */

/****  eXpressTG RSPvar    ****/
int xtg_acq1;           /* flag for data acquisiton windows 1=on, 0=off */

/****  AutoShim RSPvar  ****/
float asrot0,asrot1,asrot2,asrot3,asrot4,asrot5,asrot6,asrot7,asrot8,
    asrot9,asrot10,asrot11,asrot12,asrot13,asrot14,asrot15,asrot16,
    asrot17,asrot18,asrot19,asrot20,asrot21,asrot22,asrot23,asrot24,
    asrot25, asrot26;
float asrloc1,asrloc2,asrloc3;
float astloc1,astloc2,astloc3;
short as_slquant;
int PSdebugrotmat=0;

/* SWIFT debug */
int swift_debug = 0;
/************************ END OF PRESCAN.E ****************************/
