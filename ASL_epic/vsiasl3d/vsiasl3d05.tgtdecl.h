/*
 *  vsiasl3d05.tgtdecl.h
 *
 *  Do not edit this file. It is automatically generated by EPIC.
 *
 *  Date : Feb 28 2018
 *  Time : 11:06:49
 */

#ifndef h_vsiasl3d05_tgtdecl_h
#define h_vsiasl3d05_tgtdecl_h

/* id tags - one for MROR-WS, one for pulse sequence, 0 based */
int cont_wsid, cont_psdid;
/* time stamp to recon for MROR header */
int cont_time;
/* variable communicating to recon that real time is on */
int cont_synch;
/* stops the cont imaging after acquisition of the last pass */
int cont_stop, cont_stop_usr;
/* x, y, and z offsets from original prescription - mm units, 1/10 mm
   accuracy, range TBD */
float cont_xoffset, cont_yoffset, cont_zoffset;
/* alpha, beta, and gamma offsets from original prescription - degree
   units, 1/10 degree accuracy, range 0-360 degrees */
float cont_alpha, cont_beta, cont_gamma;

/* variable for checking scan plane offset (cont_x/y/zofset )change */
int cont_sp_changed;

/* Fix for spr MRIge46008 making sure the following prescan vars*/
/* are available for all psds */
/* The following are set by prescan upon successful completion: */
int pscR1;  /* R1 receive attenuation (dB) */
int pscR2;  /* R2 receive attenuation (barrel shift) */
int pscCF;  /* CF center frequency (Hz) */
int pscTG;  /* TG transmit gain (dB) */

int maxTGAtOffset = MAX_SYS_TG;  /* TG limit for off-center */

#ifdef PSD_HW
/* MRIge81678 - FGRET Timestamp */
int perf_tdel_id = 0;

/* 
 * If promo_rescantime_id==2, the rescan time communication between psd and recon 
 * will be initialized, the value of rescan time will then be sent to host for display
 */
int promo_rescantime_id = 0;
#endif /* PSD_HW */

/*@Start***********************************************************/
/* GEMSBG Include File
 * Copyright (C) 2000 The General Electric Company
 *
 *      Include File Name:  dabrecord   
 *      Developer: Jeff Hopkins
 *
 * $Source: dabrecord.h $
 * $Revision: 1.0 $  $Date: 6/19/00 17:04:15 $
 *
 */

/*@Synopsis 
	Contains arrays and variables for recording dab information
*/     

/*@Description
 * ******************************************
   dabrecord.h
   Author: Jeff Hopkins
   Date:   12/05/89
   
   Description: Contains constants, arrays, and counters
   for keeping track of the contents of DAB packets as
   they are passed in. This is the epic.h version for
   the ipgexport (rsp?) section. This needs to match
   one-to-one with dabrecord_pgen.h which gives the
   external links for the loaddab functions.

   Revision Date    Author     Desc.
   ________________________________________
   07/19/2000       JAH        Created view_record variables.

**********************************************/

/*JAH: for recording view order table*/
int view_record[DABRECLEN];
int dabview_counter = 0;
int fse_dabview_flag = 0;
int record_views = 0;



int iv, ifr, isl, kzcount;
int tmp;
int i, k;
int bangn;
int trig, dtype;

short trigonpkt[3] = {0, SSPOC+DREG, SSPD+DSSPD4};
short trigoffpkt[3] = {0, SSPOC+DREG, SSPD};
short trigonwd, trigoffwd;

/* variables needed for prescan */
short chopamp;
int view, slice, dabop, excitation, seqCount, ec, rf1Phase, seqCount;
int rspent, rspdda, rspbas, rspvus, rspgy1, rspasl;
int rspesl, rspchp, rspnex, rspslq, rspsct;
int dabmask;
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
int dixon,as_slice,as_view,as_excitation;
int dur_invse; /* Initial duration for the WAIT pulse before the 90 pulse */
               /* based on OPTI */
int new_dur;   /* New duration time calulated with psctitime and cfh_ti delta */
int psctitime; /* input TI time on milliseconds; initialized to OPTI */
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
float asrot0,asrot1,asrot2,asrot3,asrot4,asrot5,asrot6,asrot7,asrot8;
float asrot9,asrot10,asrot11,asrot12,asrot13,asrot14,asrot15,asrot16,asrot17;
float asrot18,asrot19,asrot20,asrot21,asrot22,asrot23,asrot24,asrot25,asrot26;
float astloc1,asrloc1,asploc1;
float astloc2,asrloc2,asploc2;
float astloc3,asrloc3,asploc3;
float asdim1,asdim2,asdim3;
short as_slquant = 3;
short as_index = 1;
int PSdebugrotmat=0;

/* SWIFT debug */
int swift_debug = 0;
/************************ END OF PRESCAN.E ****************************/
extern PSD_EXIT_ARG psdexitarg;


  WF_PULSE gzrf1a = INITPULSE;
  WF_PULSE gzrf1  = INITPULSE;
  WF_PULSE gzrf1d = INITPULSE;
  WF_PULSE rf1 = INITPULSE;

  WF_PULSE gzrf1ra = INITPULSE;
  WF_PULSE gzrf1r = INITPULSE;
  WF_PULSE gzrf1rd = INITPULSE;

  SEQUENCE_ENTRIES  off_tipdown_core;
  WF_PULSE tipdown_core;
#if defined(HOST_TGT)
  int idx_tipdown_core;   /* sequence entry index */
#endif

  WF_PULSE gz180crush1a = INITPULSE;
  WF_PULSE gz180crush1 = INITPULSE;
  WF_PULSE gz180crush1d = INITPULSE;

  WF_PULSE gzrf2a = INITPULSE;
  WF_PULSE gzrf2  = INITPULSE;
  WF_PULSE gzrf2d = INITPULSE;
  WF_PULSE rf2 = INITPULSE;

  WF_PULSE gz180crush2a = INITPULSE;
  WF_PULSE gz180crush2 = INITPULSE;
  WF_PULSE gz180crush2d = INITPULSE;

  WF_PULSE gzphase1a = INITPULSE;
  WF_PULSE gzphase1 = INITPULSE;
  WF_PULSE gzphase1d = INITPULSE;

  WF_PULSE mapx = INITPULSE;

  WF_PULSE mapy = INITPULSE;

  WF_PULSE mapz = INITPULSE;

  WF_PULSE mapt = INITPULSE;


	    WF_PULSE gx = INITPULSE;


	    WF_PULSE gy = INITPULSE;

  WF_PULSE maps1 = INITPULSE;

  WF_PULSE echo1 = INITPULSE;

  WF_PULSE mapx2 = INITPULSE;

  WF_PULSE mapy2 = INITPULSE;

  WF_PULSE mapz2 = INITPULSE;

  WF_PULSE mapt2 = INITPULSE;


	    WF_PULSE gx2 = INITPULSE;


	    WF_PULSE gy2 = INITPULSE;

  WF_PULSE maps2 = INITPULSE;

  WF_PULSE echo2 = INITPULSE;

  WF_PULSE maps3 = INITPULSE;

  WF_PULSE gzphase2a = INITPULSE;
  WF_PULSE gzphase2 = INITPULSE;
  WF_PULSE gzphase2d = INITPULSE;

  WF_PULSE gxpre = INITPULSE;

  WF_PULSE gypre = INITPULSE;

  SEQUENCE_ENTRIES  off_readout_core;
  WF_PULSE readout_core;
#if defined(HOST_TGT)
  int idx_readout_core;   /* sequence entry index */
#endif

  WF_PULSE gzrfdummya = INITPULSE;
  WF_PULSE gzrfdummy  = INITPULSE;
  WF_PULSE gzrfdummyd = INITPULSE;
  WF_PULSE rfdummy = INITPULSE;

  SEQUENCE_ENTRIES  off_dummrycore;
  WF_PULSE dummycore;
#if defined(HOST_TGT)
  int idx_dummycore;   /* sequence entry index */
#endif

  WF_PULSE endpass = INITPULSE;

  SEQUENCE_ENTRIES  off_pass;
  WF_PULSE pass;
#if defined(HOST_TGT)
  int idx_pass;   /* sequence entry index */
#endif

  
  EXTERN_FILENAME fileloc_BS0rf =  "sech_7360.rho";
  WF_PULSE BS0rf = INITPULSE;

  
  EXTERN_FILENAME fileloc_BS0rf_theta =  "sech_7360.theta";
  WF_PULSE BS0rf_theta = INITPULSE;

  WF_PULSE gzBS0rfspoilera = INITPULSE;
  WF_PULSE gzBS0rfspoiler = INITPULSE;
  WF_PULSE gzBS0rfspoilerd = INITPULSE;

  SEQUENCE_ENTRIES  off_preBScore;
  WF_PULSE preBScore;
#if defined(HOST_TGT)
  int idx_preBScore;   /* sequence entry index */
#endif


	    WF_PULSE vsitag1 = INITPULSE;


	    WF_PULSE vsitag1_theta = INITPULSE;


	    WF_PULSE gztag1 = INITPULSE;

  SEQUENCE_ENTRIES  off_astcore;
  WF_PULSE astcore;
#if defined(HOST_TGT)
  int idx_astcore;   /* sequence entry index */
#endif


	    WF_PULSE vsictl1 = INITPULSE;


	    WF_PULSE vsictl1_theta = INITPULSE;


	    WF_PULSE gzctl1 = INITPULSE;

  SEQUENCE_ENTRIES  off_controlcore;
  WF_PULSE controlcore;
#if defined(HOST_TGT)
  int idx_controlcore;   /* sequence entry index */
#endif

  SEQUENCE_ENTRIES  off_nothingcore;
  WF_PULSE nothingcore;
#if defined(HOST_TGT)
  int idx_nothingcore;   /* sequence entry index */
#endif

  SEQUENCE_ENTRIES  off_vsi_gapcore;
  WF_PULSE vsi_gapcore;
#if defined(HOST_TGT)
  int idx_vsi_gapcore;   /* sequence entry index */
#endif

  
  EXTERN_FILENAME fileloc_BS1rf =  "sech_7360.rho";
  WF_PULSE BS1rf = INITPULSE;

  
  EXTERN_FILENAME fileloc_BS1rf_theta =  "sech_7360.theta";
  WF_PULSE BS1rf_theta = INITPULSE;

  
  EXTERN_FILENAME fileloc_BS2rf =  "sech_7360.rho";
  WF_PULSE BS2rf = INITPULSE;

  
  EXTERN_FILENAME fileloc_BS2rf_theta =  "sech_7360.theta";
  WF_PULSE BS2rf_theta = INITPULSE;


	    WF_PULSE ASrf_mag = INITPULSE;


	    WF_PULSE ASrf_theta = INITPULSE;


	    WF_PULSE ASrf_grad = INITPULSE;

    WF_PULSE rf0 = INITPULSE;

  WF_PULSE gz0a = INITPULSE;
  WF_PULSE gz0 = INITPULSE;
  WF_PULSE gz0d = INITPULSE;

  SEQUENCE_ENTRIES  off_tdelaycore;
  WF_PULSE tdelaycore;
#if defined(HOST_TGT)
  int idx_tdelaycore;   /* sequence entry index */
#endif

  WF_PULSE trigon = INITPULSE;

  WF_PULSE trigoff = INITPULSE;

  WF_PULSE TRdelay = INITPULSE;

  SEQUENCE_ENTRIES  off_tadjustcore;
  WF_PULSE tadjustcore;
#if defined(HOST_TGT)
  int idx_tadjustcore;   /* sequence entry index */
#endif

  WF_PULSE waitStart = INITPULSE;

  SEQUENCE_ENTRIES  off_waitpass;
  WF_PULSE waitpass;
#if defined(HOST_TGT)
  int idx_waitpass;   /* sequence entry index */
#endif

  WF_PULSE waitEnd = INITPULSE;

  SEQUENCE_ENTRIES  off_waitend;
  WF_PULSE waitend;
#if defined(HOST_TGT)
  int idx_waitend;   /* sequence entry index */
#endif

    WF_PULSE rf1mps1 = INITPULSE;

  WF_PULSE gyrf1mps1a = INITPULSE;
  WF_PULSE gyrf1mps1 = INITPULSE;
  WF_PULSE gyrf1mps1d = INITPULSE;

  WF_PULSE gy1mps1a = INITPULSE;
  WF_PULSE gy1mps1 = INITPULSE;
  WF_PULSE gy1mps1d = INITPULSE;

  WF_PULSE gzrf1mps1a = INITPULSE;
  WF_PULSE gzrf1mps1 = INITPULSE;
  WF_PULSE gzrf1mps1d = INITPULSE;

  WF_PULSE gz1mps1a = INITPULSE;
  WF_PULSE gz1mps1 = INITPULSE;
  WF_PULSE gz1mps1d = INITPULSE;

  WF_PULSE gx1mps1a = INITPULSE;
  WF_PULSE gx1mps1 = INITPULSE;
  WF_PULSE gx1mps1d = INITPULSE;

  WF_PULSE gzrf2mps1a = INITPULSE;
  WF_PULSE gzrf2mps1  = INITPULSE;
  WF_PULSE gzrf2mps1d = INITPULSE;
  WF_PULSE rf2mps1 = INITPULSE;

  WF_PULSE gzrf2lmps1a = INITPULSE;
  WF_PULSE gzrf2lmps1 = INITPULSE;
  WF_PULSE gzrf2lmps1d = INITPULSE;

  WF_PULSE gzrf2rmps1a = INITPULSE;
  WF_PULSE gzrf2rmps1 = INITPULSE;
  WF_PULSE gzrf2rmps1d = INITPULSE;

  WF_PULSE gxwmps1a = INITPULSE;
  WF_PULSE gxwmps1 = INITPULSE;
  WF_PULSE gxwmps1d = INITPULSE;

  WF_PULSE echo1mps1 = INITPULSE;

  WF_PULSE attenuator_keymps1 = INITPULSE;

  SEQUENCE_ENTRIES  off_seqmps1;
  WF_PULSE seqmps1;
#if defined(HOST_TGT)
  int idx_seqmps1;   /* sequence entry index */
#endif

  WF_PULSE gzrf1cfla = INITPULSE;
  WF_PULSE gzrf1cfl  = INITPULSE;
  WF_PULSE gzrf1cfld = INITPULSE;
  WF_PULSE rf1cfl = INITPULSE;

  WF_PULSE gz1cfla = INITPULSE;
  WF_PULSE gz1cfl = INITPULSE;
  WF_PULSE gz1cfld = INITPULSE;

  WF_PULSE cfl_fid = INITPULSE;

  WF_PULSE cfl_attenkey = INITPULSE;

  WF_PULSE gykcfla = INITPULSE;
  WF_PULSE gykcfl = INITPULSE;
  WF_PULSE gykcfld = INITPULSE;

  SEQUENCE_ENTRIES  off_seqcfl;
  WF_PULSE seqcfl;
#if defined(HOST_TGT)
  int idx_seqcfl;   /* sequence entry index */
#endif

  WF_PULSE gxkrcvna = INITPULSE;
  WF_PULSE gxkrcvn = INITPULSE;
  WF_PULSE gxkrcvnd = INITPULSE;

  WF_PULSE gykrcvna = INITPULSE;
  WF_PULSE gykrcvn = INITPULSE;
  WF_PULSE gykrcvnd = INITPULSE;

  WF_PULSE gzkrcvna = INITPULSE;
  WF_PULSE gzkrcvn = INITPULSE;
  WF_PULSE gzkrcvnd = INITPULSE;

  WF_PULSE rcvn_wait = INITPULSE;

  SEQUENCE_ENTRIES  off_pre_rcvn;
  WF_PULSE pre_rcvn;
#if defined(HOST_TGT)
  int idx_pre_rcvn;   /* sequence entry index */
#endif

  WF_PULSE rcvrbl = INITPULSE;
  short rcvrbl_pack[4];

  WF_PULSE rcvn_fid = INITPULSE;

  WF_PULSE rcvn_attenkey = INITPULSE;

  WF_PULSE rcvrbl2 = INITPULSE;
  short rcvrbl2_pack[4];

  SEQUENCE_ENTRIES  off_seqrcvn;
  WF_PULSE seqrcvn;
#if defined(HOST_TGT)
  int idx_seqrcvn;   /* sequence entry index */
#endif

  EXTERN_FILENAME grad_zrf0cfh= "NULL";
  EXTERN_FILENAME rf_rf0cfh = "shNvrg5b.rho";

  WF_PULSE gzrf0cfha = INITPULSE;
  WF_PULSE gzrf0cfh  = INITPULSE;
  WF_PULSE gzrf0cfhd = INITPULSE;
  WF_PULSE rf0cfh = INITPULSE;


  EXTERN_FILENAME fileloc_omegarf0cfh =  "shNvrg5b.pha";
  WF_PULSE omegarf0cfh = INITPULSE;

  WF_PULSE gyrf0kcfha = INITPULSE;
  WF_PULSE gyrf0kcfh = INITPULSE;
  WF_PULSE gyrf0kcfhd = INITPULSE;

  WF_PULSE zticfh = INITPULSE;

  WF_PULSE rticfh = INITPULSE;

  WF_PULSE xticfh = INITPULSE;

  WF_PULSE yticfh = INITPULSE;

  WF_PULSE sticfh = INITPULSE;

  WF_PULSE gzrf1cfha = INITPULSE;
  WF_PULSE gzrf1cfh  = INITPULSE;
  WF_PULSE gzrf1cfhd = INITPULSE;
  WF_PULSE rf1cfh = INITPULSE;

    WF_PULSE rf2cfh = INITPULSE;

    WF_PULSE rf3cfh = INITPULSE;

    WF_PULSE rf4cfh = INITPULSE;

  WF_PULSE gxrf2cfha = INITPULSE;
  WF_PULSE gxrf2cfh = INITPULSE;
  WF_PULSE gxrf2cfhd = INITPULSE;

  WF_PULSE gyrf2cfha = INITPULSE;
  WF_PULSE gyrf2cfh = INITPULSE;
  WF_PULSE gyrf2cfhd = INITPULSE;

  WF_PULSE gzrf2lcfha = INITPULSE;
  WF_PULSE gzrf2lcfh = INITPULSE;
  WF_PULSE gzrf2lcfhd = INITPULSE;

  WF_PULSE gzrf2rcfha = INITPULSE;
  WF_PULSE gzrf2rcfh = INITPULSE;
  WF_PULSE gzrf2rcfhd = INITPULSE;

  WF_PULSE gyrf3cfha = INITPULSE;
  WF_PULSE gyrf3cfh = INITPULSE;
  WF_PULSE gyrf3cfhd = INITPULSE;

  WF_PULSE gzrf3lcfha = INITPULSE;
  WF_PULSE gzrf3lcfh = INITPULSE;
  WF_PULSE gzrf3lcfhd = INITPULSE;

  WF_PULSE gzrf3rcfha = INITPULSE;
  WF_PULSE gzrf3rcfh = INITPULSE;
  WF_PULSE gzrf3rcfhd = INITPULSE;

  WF_PULSE gy1cfha = INITPULSE;
  WF_PULSE gy1cfh = INITPULSE;
  WF_PULSE gy1cfhd = INITPULSE;

  WF_PULSE gx1cfha = INITPULSE;
  WF_PULSE gx1cfh = INITPULSE;
  WF_PULSE gx1cfhd = INITPULSE;

  WF_PULSE gzrf4cfha = INITPULSE;
  WF_PULSE gzrf4cfh = INITPULSE;
  WF_PULSE gzrf4cfhd = INITPULSE;

  WF_PULSE isi_slice1 = INITPULSE;

  WF_PULSE rot_slice1 = INITPULSE;

  WF_PULSE isi_slice2 = INITPULSE;

  WF_PULSE rot_slice2 = INITPULSE;

  WF_PULSE gzrf4lcfha = INITPULSE;
  WF_PULSE gzrf4lcfh = INITPULSE;
  WF_PULSE gzrf4lcfhd = INITPULSE;

  WF_PULSE gzrf4rcfha = INITPULSE;
  WF_PULSE gzrf4rcfh = INITPULSE;
  WF_PULSE gzrf4rcfhd = INITPULSE;

  WF_PULSE cfh_fid = INITPULSE;

  WF_PULSE cfh_attenkey = INITPULSE;

  WF_PULSE gykcfha = INITPULSE;
  WF_PULSE gykcfh = INITPULSE;
  WF_PULSE gykcfhd = INITPULSE;

  SEQUENCE_ENTRIES  off_seqcfh;
  WF_PULSE seqcfh;
#if defined(HOST_TGT)
  int idx_seqcfh;   /* sequence entry index */
#endif

  WF_PULSE contrfhubsel = INITPULSE;

  WF_PULSE contrfsel = INITPULSE;

  WF_PULSE csw_wait = INITPULSE;

  SEQUENCE_ENTRIES  off_seqcsw;
  WF_PULSE seqcsw;
#if defined(HOST_TGT)
  int idx_seqcsw;   /* sequence entry index */
#endif

  SEQUENCE_ENTRIES  off_seqcswWaitBefore;
  WF_PULSE seqcswWaitBefore;
#if defined(HOST_TGT)
  int idx_seqcswWaitBefore;   /* sequence entry index */
#endif

  WF_PULSE gzrf1ftga = INITPULSE;
  WF_PULSE gzrf1ftg  = INITPULSE;
  WF_PULSE gzrf1ftgd = INITPULSE;
  WF_PULSE rf1ftg = INITPULSE;

  WF_PULSE gz1ftga = INITPULSE;
  WF_PULSE gz1ftg = INITPULSE;
  WF_PULSE gz1ftgd = INITPULSE;

  WF_PULSE gzrf2ftga = INITPULSE;
  WF_PULSE gzrf2ftg  = INITPULSE;
  WF_PULSE gzrf2ftgd = INITPULSE;
  WF_PULSE rf2ftg = INITPULSE;

  WF_PULSE gz2ftga = INITPULSE;
  WF_PULSE gz2ftg = INITPULSE;
  WF_PULSE gz2ftgd = INITPULSE;

  WF_PULSE gzrf3ftga = INITPULSE;
  WF_PULSE gzrf3ftg  = INITPULSE;
  WF_PULSE gzrf3ftgd = INITPULSE;
  WF_PULSE rf3ftg = INITPULSE;

  WF_PULSE gz3ftga = INITPULSE;
  WF_PULSE gz3ftg = INITPULSE;
  WF_PULSE gz3ftgd = INITPULSE;

  WF_PULSE gx1ftga = INITPULSE;
  WF_PULSE gx1ftg = INITPULSE;
  WF_PULSE gx1ftgd = INITPULSE;

  WF_PULSE gx1bftga = INITPULSE;
  WF_PULSE gx1bftg = INITPULSE;
  WF_PULSE gx1bftgd = INITPULSE;

  WF_PULSE gxw1ftga = INITPULSE;
  WF_PULSE gxw1ftg = INITPULSE;
  WF_PULSE gxw1ftgd = INITPULSE;

  WF_PULSE postgxw1ftga = INITPULSE;
  WF_PULSE postgxw1ftg = INITPULSE;
  WF_PULSE postgxw1ftgd = INITPULSE;

  WF_PULSE echo1ftg = INITPULSE;

  WF_PULSE gz2bftga = INITPULSE;
  WF_PULSE gz2bftg = INITPULSE;
  WF_PULSE gz2bftgd = INITPULSE;

  WF_PULSE gx2ftga = INITPULSE;
  WF_PULSE gx2ftg = INITPULSE;
  WF_PULSE gx2ftgd = INITPULSE;

  WF_PULSE gxw2ftga = INITPULSE;
  WF_PULSE gxw2ftg = INITPULSE;
  WF_PULSE gxw2ftgd = INITPULSE;

  WF_PULSE gx2testa = INITPULSE;
  WF_PULSE gx2test = INITPULSE;
  WF_PULSE gx2testd = INITPULSE;

  WF_PULSE echo2ftg = INITPULSE;

  WF_PULSE ftg_attenkey = INITPULSE;

  SEQUENCE_ENTRIES  off_seqftg;
  WF_PULSE seqftg;
#if defined(HOST_TGT)
  int idx_seqftg;   /* sequence entry index */
#endif

    WF_PULSE rf1xtg = INITPULSE;

  WF_PULSE gyrf1xtga = INITPULSE;
  WF_PULSE gyrf1xtg = INITPULSE;
  WF_PULSE gyrf1xtgd = INITPULSE;

  WF_PULSE gzrf1xtga = INITPULSE;
  WF_PULSE gzrf1xtg = INITPULSE;
  WF_PULSE gzrf1xtgd = INITPULSE;

  WF_PULSE gykxtgla = INITPULSE;
  WF_PULSE gykxtgl = INITPULSE;
  WF_PULSE gykxtgld = INITPULSE;

        WF_PULSE rf3xtg    = INITPULSE;

        /* for RF shimming */
        WF_PULSE phs_rf3xtg      = INITPULSE;

  WF_PULSE gz1xtga = INITPULSE;
  WF_PULSE gz1xtg = INITPULSE;
  WF_PULSE gz1xtgd = INITPULSE;

  WF_PULSE gzrf2xtga = INITPULSE;
  WF_PULSE gzrf2xtg  = INITPULSE;
  WF_PULSE gzrf2xtgd = INITPULSE;
  WF_PULSE rf2xtg = INITPULSE;

  WF_PULSE gz2xtga = INITPULSE;
  WF_PULSE gz2xtg = INITPULSE;
  WF_PULSE gz2xtgd = INITPULSE;

        WF_PULSE rf4xtg    = INITPULSE;

        /* for RF shimming */
        WF_PULSE phs_rf4xtg      = INITPULSE;

  WF_PULSE gykxtgra = INITPULSE;
  WF_PULSE gykxtgr = INITPULSE;
  WF_PULSE gykxtgrd = INITPULSE;

  WF_PULSE gx1bxtga = INITPULSE;
  WF_PULSE gx1bxtg = INITPULSE;
  WF_PULSE gx1bxtgd = INITPULSE;

  WF_PULSE gxw1xtga = INITPULSE;
  WF_PULSE gxw1xtg = INITPULSE;
  WF_PULSE gxw1xtgd = INITPULSE;

  WF_PULSE echo1xtg = INITPULSE;

  WF_PULSE xtg_attenkey = INITPULSE;

  SEQUENCE_ENTRIES  off_seqxtg;
  WF_PULSE seqxtg;
#if defined(HOST_TGT)
  int idx_seqxtg;   /* sequence entry index */
#endif

  WF_PULSE gzrf1asa = INITPULSE;
  WF_PULSE gzrf1as  = INITPULSE;
  WF_PULSE gzrf1asd = INITPULSE;
  WF_PULSE rf1as = INITPULSE;

  WF_PULSE gz1asa = INITPULSE;
  WF_PULSE gz1as = INITPULSE;
  WF_PULSE gz1asd = INITPULSE;

  WF_PULSE gxwasa = INITPULSE;
  WF_PULSE gxwas = INITPULSE;
  WF_PULSE gxwasd = INITPULSE;

  WF_PULSE echo1as = INITPULSE;

  WF_PULSE gx1asa = INITPULSE;
  WF_PULSE gx1as = INITPULSE;
  WF_PULSE gx1asd = INITPULSE;

  WF_PULSE attenuator_keyas = INITPULSE;

  WF_PULSE gy1asa = INITPULSE;
  WF_PULSE gy1as = INITPULSE;
  WF_PULSE gy1asd = INITPULSE;

  WF_PULSE gy1rasa = INITPULSE;
  WF_PULSE gy1ras = INITPULSE;
  WF_PULSE gy1rasd = INITPULSE;

  WF_PULSE gxkasa = INITPULSE;
  WF_PULSE gxkas = INITPULSE;
  WF_PULSE gxkasd = INITPULSE;

  WF_PULSE gzkasa = INITPULSE;
  WF_PULSE gzkas = INITPULSE;
  WF_PULSE gzkasd = INITPULSE;

  WF_PULSE xdixon = INITPULSE;

  WF_PULSE ydixon = INITPULSE;

  WF_PULSE zdixon = INITPULSE;

  WF_PULSE sdixon = INITPULSE;

  WF_PULSE sdixon2 = INITPULSE;

  SEQUENCE_ENTRIES  off_seqaushim;
  WF_PULSE seqaushim;
#if defined(HOST_TGT)
  int idx_seqaushim;   /* sequence entry index */
#endif

  WF_PULSE pass_aushim = INITPULSE;

  SEQUENCE_ENTRIES  off_seqpassas;
  WF_PULSE seqpassas;
#if defined(HOST_TGT)
  int idx_seqpassas;   /* sequence entry index */
#endif

  WF_PULSE dDDIQ = INITPULSE;

  SEQUENCE_ENTRIES  off_seqIQControl;
  WF_PULSE seqIQControl;
#if defined(HOST_TGT)
  int idx_seqIQControl;   /* sequence entry index */
#endif

  WF_PULSE rf1rs = INITPULSE;

  WF_PULSE gzrf1rsa = INITPULSE;
  WF_PULSE gzrf1rs = INITPULSE;
  WF_PULSE gzrf1rsd = INITPULSE;

  WF_PULSE gxkbsrsa = INITPULSE;
  WF_PULSE gxkbsrs = INITPULSE;
  WF_PULSE gxkbsrsd = INITPULSE;

  WF_PULSE gz1rsa = INITPULSE;
  WF_PULSE gz1rs = INITPULSE;
  WF_PULSE gz1rsd = INITPULSE;

  WF_PULSE rfbrs = INITPULSE;

  WF_PULSE thetarfbrs = INITPULSE;

  WF_PULSE gzkbsrsa = INITPULSE;
  WF_PULSE gzkbsrs = INITPULSE;
  WF_PULSE gzkbsrsd = INITPULSE;

  WF_PULSE gxwrsa = INITPULSE;
  WF_PULSE gxwrs = INITPULSE;
  WF_PULSE gxwrsd = INITPULSE;

  WF_PULSE echo1rs = INITPULSE;

  WF_PULSE gx2rsa = INITPULSE;
  WF_PULSE gx2rs = INITPULSE;
  WF_PULSE gx2rsd = INITPULSE;

  WF_PULSE gy2rsa = INITPULSE;
  WF_PULSE gy2rs = INITPULSE;
  WF_PULSE gy2rsd = INITPULSE;

  WF_PULSE gxw2rsa = INITPULSE;
  WF_PULSE gxw2rs = INITPULSE;
  WF_PULSE gxw2rsd = INITPULSE;

  WF_PULSE gx1rsa = INITPULSE;
  WF_PULSE gx1rs = INITPULSE;
  WF_PULSE gx1rsd = INITPULSE;

  WF_PULSE gy1rrsa = INITPULSE;
  WF_PULSE gy1rrs = INITPULSE;
  WF_PULSE gy1rrsd = INITPULSE;

  WF_PULSE gy1rsa = INITPULSE;
  WF_PULSE gy1rs = INITPULSE;
  WF_PULSE gy1rsd = INITPULSE;

  WF_PULSE gzkrsa = INITPULSE;
  WF_PULSE gzkrs = INITPULSE;
  WF_PULSE gzkrsd = INITPULSE;

  WF_PULSE gxkrsa = INITPULSE;
  WF_PULSE gxkrs = INITPULSE;
  WF_PULSE gxkrsd = INITPULSE;

  WF_PULSE attenuator_keyrs = INITPULSE;

  SEQUENCE_ENTRIES  off_seqrs;
  WF_PULSE seqrs;
#if defined(HOST_TGT)
  int idx_seqrs;   /* sequence entry index */
#endif

  WF_PULSE pass_rs = INITPULSE;

  SEQUENCE_ENTRIES  off_seqpassrs;
  WF_PULSE seqpassrs;
#if defined(HOST_TGT)
  int idx_seqpassrs;   /* sequence entry index */
#endif

  WF_PULSE rf1dtg = INITPULSE;

  WF_PULSE gzrf1dtga = INITPULSE;
  WF_PULSE gzrf1dtg = INITPULSE;
  WF_PULSE gzrf1dtgd = INITPULSE;

  WF_PULSE gxkbsdtga = INITPULSE;
  WF_PULSE gxkbsdtg = INITPULSE;
  WF_PULSE gxkbsdtgd = INITPULSE;

  WF_PULSE gz1dtga = INITPULSE;
  WF_PULSE gz1dtg = INITPULSE;
  WF_PULSE gz1dtgd = INITPULSE;

  WF_PULSE rfbdtg = INITPULSE;

  WF_PULSE thetarfbdtg = INITPULSE;

  WF_PULSE gzkbsdtga = INITPULSE;
  WF_PULSE gzkbsdtg = INITPULSE;
  WF_PULSE gzkbsdtgd = INITPULSE;

  WF_PULSE gxwdtga = INITPULSE;
  WF_PULSE gxwdtg = INITPULSE;
  WF_PULSE gxwdtgd = INITPULSE;

  WF_PULSE echo1dtg = INITPULSE;

  WF_PULSE gx2dtga = INITPULSE;
  WF_PULSE gx2dtg = INITPULSE;
  WF_PULSE gx2dtgd = INITPULSE;

  WF_PULSE gy2dtga = INITPULSE;
  WF_PULSE gy2dtg = INITPULSE;
  WF_PULSE gy2dtgd = INITPULSE;

  WF_PULSE gxw2dtga = INITPULSE;
  WF_PULSE gxw2dtg = INITPULSE;
  WF_PULSE gxw2dtgd = INITPULSE;

  WF_PULSE gx1dtga = INITPULSE;
  WF_PULSE gx1dtg = INITPULSE;
  WF_PULSE gx1dtgd = INITPULSE;

  WF_PULSE gy1rdtga = INITPULSE;
  WF_PULSE gy1rdtg = INITPULSE;
  WF_PULSE gy1rdtgd = INITPULSE;

  WF_PULSE gy1dtga = INITPULSE;
  WF_PULSE gy1dtg = INITPULSE;
  WF_PULSE gy1dtgd = INITPULSE;

  WF_PULSE gzkdtga = INITPULSE;
  WF_PULSE gzkdtg = INITPULSE;
  WF_PULSE gzkdtgd = INITPULSE;

  WF_PULSE gxkdtga = INITPULSE;
  WF_PULSE gxkdtg = INITPULSE;
  WF_PULSE gxkdtgd = INITPULSE;

  WF_PULSE attenuator_keydtg = INITPULSE;

  SEQUENCE_ENTRIES  off_seqdtg;
  WF_PULSE seqdtg;
#if defined(HOST_TGT)
  int idx_seqdtg;   /* sequence entry index */
#endif

  WF_PULSE pass_dtg = INITPULSE;

  SEQUENCE_ENTRIES  off_seqpassdtg;
  WF_PULSE seqpassdtg;
#if defined(HOST_TGT)
  int idx_seqpassdtg;   /* sequence entry index */
#endif

#endif /* h_vsiasl3d05_tgtdecl_h */

