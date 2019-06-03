/*@Start***********************************************************/
/* spiral sequence.
 *   s2.e 11/9/92   cleaned up @rsp of spiral.e, in order to add movies.
 *   s3.e 11/17/92  added map making.
 *   s4.e 11/27/92  added auto-prescan (from gaw's bl9.e).
 *   s5.e 11/27/92  flow-comped spin-echo sequence.
 *   s6.e 11/30/92  put flow-comp lobe on right side of 180.
 *                  fixed numerous things.  s1-s5.e are now out of date.
 *                  turned on "Gradient Echo" plasma option.
 *                  added fastrecon switch and passthrough filter
 *                  (with > 1 ms increase in "gradient-echo TE").
 *                  added option for longer, minimum-phase pulse.
 *         1/11/93  Changed the plasma around quite a bit.
 *                  Added small-FOV gradients for surface- and head-coil scans.
 *                  Fixed nextra so that physics setup time is either
 *                  3 heartbeats (cardiac gating) or 1.5 seconds.
 *         1/19/93  Revamped cardiac gating to fix bugs and add features.
 *         2/17/93  Added empirical rf1 amplitude fudge factor of 1.4.
 *   s7.e  2/23/93  Added support for Kevin King's tardis routine.
 *                  Got rid of "correct" setting of xres and yres,
 *                  because they need to be 256.  Turned off fids.
 *                  Still calculates gres, but doesn't currently use it.
 *                  Added slice stepping and slice offset.
 *                  Wired up tdel button.
 *   s8.e  3/18/93  Changed to fast, 127-tap filter for short TR's.
 *                  Removed fastrecon cv.  Added fastfilt cv.
 *   sf1.e 4/18/93  Beginnings of a 5.X fluoro/activation sequence.
 *                  Cut down to one gradient set (20-cm, 20-interleaves).
 *                  Kissed off the maps.  Allowed consecutive images
 *                  using nframes.  Removed ncviews stuff.  Changed
 *                  to sf1.*.in parameter file.
 *   C. Meyer, 1992, 1993
 *
 *   sf3.e adding spatial sat pulse and multislice D.M.S
 *   2/22/94 adding multiple gradient sets, DMS
 *   4/20/94 adding fastfilt for shorter tr
 *
 *	sg4.e	3/95	ghg	it does the te and tr right
 *	sg5.e	4/7/95	ghg	relax the 256 image limit by using echoes
 *	     	4/15/95	ghg	add rf spoiling
 *	sg6.e	5/1/95	ghg	add b0 map as first echo shift
 *	sg7.e	5/14/95	ghg	add autoshim
 *	sg8.e	6/9/95	ghg	rf pulse is a sinc2
 *		7/3/95		can do spin echo, sinc1.
 *	sg9.e	7/4/95	ghg	uses sg* gradients.
 *		7/12/95		selects between 20 and 32 cm fov waves
 *		8/9/95		fix create image bug, no longer allocates
 *				can start sequence by aux trig
 *    sprl501	9/29/95	ghg	First shot at 5.5 version
 *    sprl502	10/27/95	fix prescan to chop
 *    sprl503	11/17/95	fix gram duty cycle calc for cardiac
 *    sprl504	11/18/95	fix triggering for gating and aux
 *    sprl505	12/10/95	use filter_gen
 *    sprl507	2/2/96		use fast receiver.  All data is in
 *				single slice
 *    sprl508	2/8/96		back to usual multislice data set
 *    sprl509	3/20/96		add ext trig opuser
 *		4/16/96		add data store opuser
 *    sprl510	7/17/96		move scalerotmats to inner loop
 *    sprl511	11/24/96	new data format: 1 echo, all data in rhnframes
 *   				add support for grev != 5
 *    sprl512	12/18/96	add grev user cv
 *    sprl513	5/23/97		support for grad waves with higher dutycycle
 *		7/16/97		remove artifical limit on nframes
 *    sprl514	2/2/98		support for grev 7
 *		3/25/98		don't do B0 mapping if spin echo
 *    sprl801	5/9/98		first LX version
 *    sprl803	8/23/98		gradient files contain a little hdr at end
 *				to give res_gx, tsp, etc.
 *		9/4/98		for CV1 version, can use setfilter.c with
 *				last argument changed
 *    sprl804	10/1/98		add fast_rec init code for prescan and
*				scan with boffset
*    sprl805	10/26/98	does concat acquisition if rhfrsize >8K
*    sprl806	10/30/98	puts second acq in next view, rather
*				than echo
*    sprl808	1/18/99		put spoilers on the readout axes to
*				fix the hardware hysteresis prob
*    sprl809	1/23/99		real time spiral gen
*    sprl810	2/27/99		add filter delay stuff
*    sprl811	3/3/99		save grad rotator for maxwell correction
*    sprl812	3/31/99		try to get clock to start with aux_trig
*                              also fix aps2 for doconcat
*		4/13/99		ihtr = optr to get the raw header right
*    sprl813	5/8/99		allow slice clustering, remove
*				non-interleave option
*		10/8/99		fix gslew problem
*		10/23/99	rationalize params vs gradients
*		1/20/00		install capability for TTL output pulse
*				on CAP board, J10
*    sprl816	3/10/00		variable bandwidth
*				also remove recon size cv
*    sprl830	5/27/00		make 8.3 version
*				add nextra opuser
*    sprl831	7/27/00		eliminate concat acq because rhfrsize works
*    		8/5/00		put in RT hooks (scancore)
*    		8/14/00		make slow receiver prescan work
*    sprl832	10/14/00	allow new ext trig mode 2
*    sprlio832	10/22/00	spiral out, in, combiner
*    sprlio834 2/10/01         install bmapnav: don't move dacq for first frame
*		12/30/01	don't use bmapnav if gtype != 0
*    sprlio835 5/28/02         put in ss rf pulse
*    sprlio836 7/5/02         	fix prephaser sinusoid amplitude
*		7/15/02		trap over-range in phase for ss pulse
*    sprlio837 7/16/02         make an error message for trap in host
*		8/6/02		add check for rawsize vs cftpssize
*    sprlio838 8/7/02          play out ss pulse on theta board by
*				integrating the frequency offsets
*		8/12/02		remove the nex button
*    sprlio839 9/4/02          fix bug in aps2 for gtype=2
*    sprlio840 9/15/02         add code to switch ACGD transition mode
*              9/23/02         move acgd stuff to prescan
*              11/13/02        acgd_tr = 50ms + 0.036 ms/image
*    sprlio841 11/24/02        use fat sat instead of ss pulse
*		3/21/03		add code to discriminate daqdeloff
*				based on gradient type
*              3/25/03         piimages = 0 to not allocate images
*    sprlio843 11/25/03        make bmapnav work for gtype = 2-
*				move only the spiral-out for first frame
*    sprlio1043 11/25/03       first 10-11x version
*	  	 1/27/04        add opuser to control physio acquisition
*    sprlio1045 9/26/04        offset receiver FOV: search thetrec
*    sprlio1146 1/29/05        add GP3 Gradient mode control
*		 2/23/05	fix a little bug in daqdeloff switch
*		 3/19/05	make optr, ihtr max 24s
*
*    sprlio1246 4/29/05        12.0 version- remove GP3 include and physio
*				cv- included in epic.h
*		 5/5/05		remove numrec from rhrawsize
*		 10/23/05	allow smaller bandwidths
*    sprlio1247 11/14/05       make pw_rf1 9ms if opslthick < 1 mm
*                              This allows thinner slices (600 um).
*               12/18/05       fixed phys_record_flag to work (prescan inline
		*                              was setting it to 0).
*               1/8/06      	Added realtime opuser cv
*    sprlio1248 5/25/06        put in rcvr delay = tsp/2 us.
*    sprlio1249 8/23/06        tighten up tmin timing calc: tlead, psd_rf_wait
*               1/4/08         zone_cntl = 0
*    sprlio1250 1/12/08        check rawsize to exceed 2^31.
*    sprlio2050 2/21/08        make 20x version- don't check rawsize, remove
*				zone control
*    sprlio2051 3/29/08        redo the spoilers to avoid PNS.  Not used if nl=1.
*		 4/5/08		fix bmapnav = 0 (add maps3)
*    sprlio2053 5/11/08        don't do phys record except in scan() entry
*		 7/1/08		add file trig mode- start scan on StartScan exists (opuser2)
*		 8/13/08	add physiochansel CV to select only PPGTrig and RESPData
*               9/6/08         in realtime mode, put wait at end for recon to finish
*               12/19/08       drop slew rate to 150 (sigh)
*    sprlio2054 12/21/08       make obl_method = 1 to shorten ramps
*    sprlio2055 2/23/09        read optional daqdel.off file
*    sprlio2056 3/22/09        set prescan filter to same as scan
*		 3/27/09	add little delay for cluster to deadtime
*    sprlio2057 4/13/09        make prescan filter with oprbw bw, but length psfrsize
*    sprlio2058 7/21/09        20M4- remove phys_record_shannelsel as c
*               8/7/09         make clock run for ext trig = 2
*    sprlio2059 11/5/09        limit time for disdaqs to 6 sec in prescan
*    sprlio2060 11/11/09       change cyc_rf1 to 4  (and a_rf1)
*    sprlio2062 4/4/10        	now fix the quantization error in optr/opslquant
*    sprlio2063 9/3/10        	recompile with Prescan.e that eliminates
*				the noise cal crusher (replaces with wait)
*				add opuser12 to allow short rf1 pulse
*    sprlio2263 11/8/10        for 22x add }; after PSeplist inline
*    sprlio2264 12/1/10        add slice loc phase correction
*    sprlio2265 8/1/11         add variable density code from catie
*		 10/6/11	fix clustered acquisition- add 452us to tmin
*    sprlio2266 12/10/11       add sequential slice order mode
*
*
* (c) Board of Trustees, Leland Stanford Junior University 1993-2011.
*      Gary H. Glover
*
*    pcasl3d06	6/28/12		Using Stanford spiral as the base for importing PCASL sequence
*				from the Signa into the MR750 platform
*				real time support
*				3D acquitision
*    pcasl3d07	12/15/12	Adding two Background suppression pulses during the post inversion delay
* 				and one before labeling
*    pcasl3d08 6/11/13		turning BS off during early time points to get M0 maps and better
*				selection of R1 and R2 gains 
*
*    pcasl3d09			this one included the SPINS (or nautilus) readout
*
*    pcasl3d10   1/29/16	back to stack of spirals (derived from pcasl3d08).  implementing 
*         			GRAPPA in the kz direction only.
*
*   vsiasl3d03.e 9/27/16	major change:  replacing PCASL train with VSAI pulses
*				these are read externally from text files.  
*				essentially merging pcasl3d10 and vsiasl2d05.
*/
/*@End*************************************************************/

@inline epic.h
@inline intwave2.h

@global
#include <math.h>
#ifndef HW_IO
#include <stdio.h>
#include <stdlib.h>
#else /* HW_IO */
#include <stdioLib.h>
#endif /* !HW_IO */
#include <string.h>

#include "stddef_ep.h"
#include "epicconf.h"
#include "pulsegen.h"
#include "em_psd_ermes.in"
#include "epic_error.h"
#include "support_func.h"
#include "filter.h"
#include "epicfuns.h"
#include "epic_loadcvs.h"
#include "InitAdvisories.h"
#include "psdiopt.h"
#include "epic_iopt_util.h"

#include "grad_rf_sprlio.globals.h"

#define GRESMAX       21000     /* number of points grad max */
#define RES_GRAMP       100     /* number of points grad rampdown  */
#define TSP             2.0us   /* slow rcvr sampling 500kHz */
#define MAX_BW		250
#define IRINT(f) ((((f)-(int)(f)) < 0.5) ? ((int) (f)) : (1+(int)(f)))
#define MAX(a,b) ( (a > b) ? a : b )
#define MIN(a,b) ( (a < b) ? a : b )
#define ABS(a) ( (a < 0) ? -a : a )

	/*-------------------------------------------------------------*/
	/* LHG define constants for spin labeling version */

#define GAMMA_H1 26754          /* in (rad/s)/Gauss */
#define LHG_DEBUG 1
	/*-------------------------------------------------------------*/


@inline Prescan.e PSglobal
	int debugstate = 1;

@cv
@inline loadrheader.e rheadercv
@inline Prescan.e PScvs

int nl = 1 with {1,128,,,"number of interleaves",}; /* number of interleaves */
int scluster = 0 with {0,1,,,"1=cluster slices together in time, 0=not",};
int nextra = 2 with {0,,,,"number of disdaqs",};
int nframes = 1 with {1,,,,"number of time frames",};
int total_views;	/*  total shots to gather per slice */
int gating = TRIG_INTERN with {,,,,"1=line,3=ecg,5=aux,7=intern,9=nline",};
int psdseqtime;     /* sequence repetition time */
int psdtottime;     /* actual TR = psdseqtime*opslquant */
/* int timessi=400us with {0,,400us,INVIS,"time from eos to ssi in intern trig",};*/
int timessi=120us with {0,,400us,INVIS,"time from eos to ssi in intern trig",};
int trerror;		/* error in TR total */
int nerror;		/* num slices  to take up the TR error */
int gtype = 0 with {0,99,1,VIS, "trajectory: (0) spiral out, (1) spiral in, (2) both",};
/* float gslew = 150.0 with {0.0,3000.0,1,VIS, "readout gradient max slew, mT/m/ms",};  */
float gslew = 150.0 with {0.0,3000.0,1,VIS, "readout gradient max slew, mT/m/ms",};
float gamp = 2.3 with {0.0,50.0,1,VIS, "readout gradient max amp",};
float gfov = 24.0 with {0.0,100.0,1,VIS, "readout gradient design FOV, cm",};
int gres = 0 with {0,32768,1,VIS, "gradient waveform resolution",};
float gmax = 2.4 with {0.0,10.0,1,VIS, "bw-limited readout gradient amp",};
float rtimescale = 2.1 with {0.1,10.0,1,VIS, "grad risetime scale relative to cfsrmode",};
int nramp = RES_GRAMP with {0,4096,1,VIS, "spiral gradient ramp res",};
float agxpre = 0.;      /* calculated x dephasing amplitude for spiral in-out */
float agypre = 0.;      /* calculated y dephasing amplitude for spiral in-out */
float prefudge = 1.004 with {0.0,2.0,1,VIS, "prephaser fudge fac ",};
int pwgpre = 1ms;       /* calculated dephasing width for spiral in-out */
int tadmax;     /* maximum readout duration  */
int tlead = 160;      /* initial lead time for rf pulse prep*/
float daqdeloff = 0.0us; /* fudge factor, tuned at SRI  */
int daqdel = 128us; /*  gradient delay vs readout */
int thetdeloff = 0us;
int espace = 140us with {0,,1,VIS, "space between echoes",};
int readpos;        /* position of readout gradients, based on oppseq */
int daqpos;        /* position of readout window, based on oppseq */
int minte;
int seqtr = 0 with {0,,1,VIS, "total time to play seq",};
int endtime = 500ms with {0,,,,"time at end of seq in rt mode",};

int vdflag = 1 with {0,1,,VIS, "variable-density flag",};
float alpha = 3.6 with {1.01,200,,VIS, "variable-density parameter alpha",};
float kmaxfrac = 0.5 with {0.05,0.95,,VIS, "fraction of kmax to switch from constant to variable density spiral",};

int slord;

float satBW = 440.0 with {0.0,10000.0,1,VIS, "sat bw, Hz",};
float satoff = -520.0 with {-10000.0,1000.0,1,VIS, "sat center freq, Hz",};
int pwrf0;
float arf0;
int 	fuzz, fatsattime;

int pwrf1 = 6400us;
int cycrf1 = 4;

int domap = 1 with {0,2,1,VIS, "1=do B0 map & recon, 2 = only do B0 map",};
int mapdel = 2ms;   /* variable delay for field map */
int bmapnav = 1 with {0,1,1,VIS, "1=do nav cor for bmap",};

int off_fov = 1 with {0,1,1,VIS, "1 for rcvr offset fov, 0 for fftshift",};

float tsp = 4.0us with {1.0,12.75,,, "A/D sampling interval",};
float bandwidth = 125.0 with {2.0,250.0,1,VIS, "CERD rec low pass freq, kHz",};
float decimation = 1.0;

int filter_aps2 = 1;
int psfrsize = 512 with {0,,,INVIS, "num data acq during APS2",};
int queue_size = 1024 with {0,,,INVIS, "Rsp queue size",};

/* int spgr_flag = 1 with {0,256,,VIS,"Spoiling: 0=no, 1=GEMS, 2=rand, >2 =quad slope",}; */
int seed = 16807 with {0,16807,,VIS,"Spoiled Grass seed value",};
int nbang = 0 with {0,,,,"total rf shots",};
float zref = 1.0 with {,,,VIS, "refocussing pulse amp factor",};

int trigloc = 1ms with {0,,0,VIS, "location of ext trig pulse after RF",};
int triglen = 4ms with {0,,0,VIS, "duration of ext trig pulse",};
int trigfreq = 1 with {1,128,0,VIS, "num views between output trigs",};
int maketrig = 1 with {0,1,0,VIS, "trigger output port (1) or not (0)",};

int ndisdaq = 4;                 /* For Prescan: # of disdaqs ps2*/
int pos_start = 0 with {0,,,INVIS, "Start time for sequence. ",};
int obl_debug = 0; /* for obloptimize, I guess */
int obl_method = 1 with {0,1,0,INVIS,
	"On(=1) to optimize the targets based on actual rotation matrices",};

/***************************************************************************/

/* LHG 9.27.16  ***********************
   CVs for the spin tagging part .....
 ***************************************/
/* int opslquant = 16 with {8,128,,INVIS, "# pixels along z",};*/
int 	t_tag;              /*duration of tagging pulse*/
int 	t_delay;            /*delay after tagging before Aquisition*/
int 	B1tag_offset;       /*frequency offset of tagging pulse*/
int 	astseqtime;         /*duration of a single cycle of ast part of the sequence */
int 	astseqtime2;        /*duration of the total ast part of the sequence */
int 	t_adjust;           /*delay before ast pulse to make up the rest of the TR */
int 	delta1 = 0us;       /* time between flow spoiling gradients */
int 	isLabel=0;           /* variable used to toggle the RF */
int 	mycrush=0;          /* 1 = flow crushing ON   0 = OFF */
int 	crushtype=1;        /* 1 = crushers straddle 180    else   bipolar pair after 180 */
int 	mycontrol=0;        /* 0 = toggle RF   1= flip gradients */
int 	xres=64;            /* used in place of opxres to  fix a simulation error */  

int	vsi_controlmode = 5 with {1,,,VIS, "Control Mode: (1) flip the magnitude  (2) no RF  (3) negative vel target (4) flip grads (5) no VelSel grads, (6) monopolar grads", } ;  

int	vsi_rfdelay = 50; /* allow a little time between RF and gradient */
float 	vsi_flip ;  /* flip angle of each Hanning pulse */
int  	myramptime = 200; /* 120; *** lhg 9/21/12   ramp times for the pcasl pulses */
double 	vsi_Gmax = 2.0; /* default gradient max aplitude in VSI pulse train */
double 	vsi_Gcontrol = 2.0 ; /* default gradient max aplitude in control pulse train */
double	vel_target ; /* the decceleration rate the we wish to use for labeleing */
int	vsi_Ncycles = 1;  /* number of VSI pulse bursts */
int 	M0frames = 2;  /* number of images withoug labeling (at the begining) */
int	isOdd = 0;

/* stuff for peak B1 calculations */
double	vsi_RFmax = 117 ; /* mGauss. Peak B1 of the pulse train 
				for a 180 deg, 0.4 ms  hard pulse 
				B1max should be 292.50 mG*/ 
double 	rf1_RFmax = 147;   /* the excitation pulse is a hanning windowed 7 lobe sinc
				for 5 ms: should be B1max = 375.16 mG for a 180 deg flip*/ 
/* this pulse is r13d.rho :   still need to calculate this number, so this is likely wrong !!  */

double	dummy_RFmax = 200;    /* mGauss.  Peak B1 of reference pulse rfdummy in rf_grad_sprilo.h
				this is the largest pulse, and therefore the 1.0*/

/*LHG 10.2015 */
int	vsi_train_len = 250;  /* number of points in the pulse.  
				Note that the external files are named according to the number of points */
int 	Npoints;	   /* the number of points that were read from the RF pulse file.  Should match vsi_train_len */
double 	vsi_velocity = 0;   /* velocity of the spins to target for inversion */
int	vsi_timegap = 50000;   /* gap between VSI pulses if we are doing mulitple pulsees */

/* arterial suppression by BIR pulses */
double 	BIR_Gmax = 1.0; /* default gradient max aplitude in VSI pulse train */
int	BIR_len = 7000;  /* number of points in the iArterial suppression BIR8 pulse. */ 

int  	multiFlipFlag = 0;  /* LHG 10.9.14 - Flag for VSI flip angle optimization  */  
int  	multiphsFlag = 0;  /* LHG 10.3.12 - Flag for phase correction optimization mode */  
int  	nfr_phs = 2; /* Number of frameas collected at each phase correction increment */

/* CV's for 3D imaging JFN 09.29.10     */
int 	spgr_flag = 3 with {0,256,0,VIS,"Spoiling: 0=no, 1=GEMS, 2=rand, 3=RF-spoiling (117), >3 =quad slope",};
int 	rfamp_flag = 0; /* with {0,99,0,VIS,"Flip angle schedule: 0=constant, 99=Gai, >0 = power law",}; */
int 	rf_spoil_seed = 117 with {0,180,,VIS,"RF-spoiling phase increment (degrees)",};

/* LHG 6.29.12 - these are the variables for the kz phase encoding gradient plus some other stuff(?)*/	
float 	target;
int 	rtime, ftime;
float 	zfov;
float 	ampGzstep;
float 	areaGzstep;
int 	iampGzstep;
int 	pwGzstep = 500 ; /*plateau of the kz encoding gradient (us) */
int 	dopresat = 0 with {0,1,,INVIS, "pre-saturate 3D readout?",};
int	doZgrappa = 0;

float  	rf1maxTime = (456.0 / 500.0);  /* the max of the pulse happens at 456.  the whole pulse has 500 points) */
double  area_gzrf1;
double  area_gzrf1r;
float   rf1_bw;
float   slab_fraction=0.85;  /* excite only a fraction of the desired z FOV*/

/* LHG 12/14/12 : Variables for the backgroun Suppression pulses */
int	BStime = 0;    /* total time needed for background suppresson block */
int	BS1_time = 1100000; /* delay between label and first BS inversion pulse , default to 1100 ms after tag*/
int 	BS2_time = 300000;  /* delay between first and second BS inversion pulses, default to 300 ms after the first pulse */
float	rfscalesech;	/* ratio of areas: autoprescan sinc to SECH pulse */
int	doBS = 1;	/* background suppression pulses */
int	doArtSup = 0;	/*arterial suppression pulses */
int	t_preBS;	/* this is the duration of the core containing the first BS pulse (beofre tagging */
int	AStime = 50000; /* time between arterial suppression pulse and fat sat pulse */

/* some variables to adjust timing - leftovers from signa */
float 	slwid180 = 1.2 with {0.0,4.0,,,"180 slice width as fraction of 90",}; /*this is so that the 180 stuff will work */
int	tpre=0;
int	tdel=0;
int	mytpre = 0;
int 	tcs;  /* duration of the chem sat pulses */
int	cyc_rf1 = 4; /* this seems like it was used by the stanford code, but was missing from declarations */
int 	echoshift = 0us with {,,,,"180 advance for spin echo T2* wting",};

/* int textra_tadjust = 24000;   old number - LHG 9.20.07 */
int 	textra_tadjust = 100;
int 	textra_astcore = 200;
int 	textra_controlcore = 200;
int 	textra_nothingcore = 200;
int 	textra_delaycore = 300;
int 	textra_vsi_gapcore = 100;

int textra = 0;  /* adjustment to psdeqtime... don't really need it any more */

/* cv's for chem sat pulse */
int 	fat_chemshift = -440;

/*LHG 2017.01.06:  creating variables for VSI pulse: */
float 	flip_vsitag1 = 180;
int 	wg_vsitag1 =  TYPRHO1 with {0, WF_MAX_PROCESSORS*2-1,
                                           TYPRHO1, VIS, , };


/****************************************/

@ipgexport
@inline Prescan.e PSipgexport
long savrot[TRIG_ROT_MAX][9];   /* copy of rotation matrices */
RF_PULSE_INFO rfpulseInfo[RF_FREE];
int Gx[GRESMAX];
int Gy[GRESMAX];	 /*  for spiral waves  */

/* LHG 9.27.16  VSI pulses data arrays */
int vsi_pulse_mag[50000];
int vsi_pulse_phs[50000];
int vsi_pulse_ctl_phs[50000];
int vsi_pulse_grad[50000];
int vsi_pulse_ctl_grad[50000];
/*LHG 10.4.16 BIR pulses for arterial suppression */
int BIR_mag[50000];
int BIR_phs[50000];
int BIR_grad[50000];


@host
#include "stdio.h"
#include "sar_pm.h"
#include "grad_rf_sprlio_test.h"
#include "fudgetargets.c"

static char supfailfmt[] = "Support routine %s exploded!";
FILTER_INFO echo1_filt;
FILTER_INFO aps2_filt;

int genspiral(float D, int N, float Tmax, float dts);
int genspiralcdvd(float D, int N, float Tmax, float dts, float alpha, float kmaxfrac);
int gram_duty(void);

/*LHG 9.27.16 */
int read_vsi_pulse(int *vsi_pulse_mag, int *vsi_pulse_phs, int *vsi_pulse_grad, int vsi_train_len);
int calc_vsi_phs_from_velocity (int* vsi_pulse_mag, int* vsi_pulse_phs, int* vsi_pulse_grad,	float vel_target, int vsi_train_len, double vsi_Gmax);



@inline Prescan.e PShostVars

abstract("vel. sel. inversion ASL  with 3D, Spiral acquisition ");
psdname("vsiasl3d03");

int cvinit()
{
	fprintf(stderr,"\nStarting cvint");
	EpicConf();
	inittargets(&loggrd, &phygrd);
	fudgetargets(&loggrd, &phygrd, rtimescale);
	if (obloptimize(&loggrd, &phygrd, scan_info, exist(opslquant), exist(opplane),
				exist(opcoax), obl_method, obl_debug, &opnewgeo, cfsrmode)==FAILURE)
		return FAILURE;

@inline Prescan.e PScvinit
#include "cvinit.in"	/* Runs the code generated by macros in preproc.*/

	cvmax(optr,24s);
	cvmax(ihtr,24s);

	/* restrict sequence type to gradient echo */
	cvdef(oppseq,2);
	cvmin(oppseq,2);
	cvmax(oppseq,2);

	/* don't make the user even worry about selecting the number of echoes
	   or NEX  */
	piechnub = 0;
	pinexnub = 0;

	/* turn on variable bandwidth button */
	pircbnub = 6;
	pircb2nub = 0;
	pircbval2 = 15.875;
	pircbval3 = 31.75;
	pircbval3 = 62.5;
	pircbval3 = 125.0;
	pircbval4 = 200.0;
	cvmin(oprbw, 15.875);
	cvmax(oprbw, 250);
	cvdef(oprbw, 200);
	oprbw = 200;

	return SUCCESS;

}

int cveval()
{
	int tmptr, entry;

	if (_psd_rf_wait.fixedflag == 0)  { /* sets psd_grd_wait and psd_rf_wait */
		if (setsysparms() == FAILURE)  {
			epic_error(use_ermes,"Support routine setsysparams failed",
					EM_PSD_SUPPORT_FAILURE,1, STRING_ARG,"setsysparms");
			return FAILURE;
		}
	}

	if (obloptimize(&loggrd, &phygrd, scan_info, exist(opslquant), exist(opplane),
				exist(opcoax), obl_method, obl_debug, &opnewgeo, cfsrmode)==FAILURE)
		return FAILURE;

	if (existcv(opcgate) && (opcgate == PSD_ON))

	{
		tmptr = RUP_GRD((int)((float)(exist(ophrep))
					*(60.0/exist(ophrate))* 1e6));
		pitrnub = 0;
	}
	else
	{
		tmptr = optr;
		pitrnub = 6;
		if (oppseq == 2) /* gradient echo */
		{
			pitrval2 = 80ms;
			pitrval3 = 250ms;
			pitrval4 = 500ms;
			pitrval5 = 1s;
			pitrval6 = 2s;
		}
		else
		{
			pitrval2 = 300ms;
			pitrval3 = 500ms;
			pitrval4 = 1s;
			pitrval5 = 1.5s;
			pitrval6 = 2s;
		}
	}
	pisctim1 = nframes*nl*tmptr;
	pisctim2 = pinexval2*pisctim1;
	pisctim3 = pinexval3*pisctim1;
	pisctim4 = pinexval4*pisctim1;
	pisctim5 = pinexval5*pisctim1;
	pisctim6 = pinexval6*pisctim1;

	opflip = 90.0;
	cvdef(opflip, 90.0);
	if (oppseq == 2) /* gradient echo */
	{
		pifanub = 6;
		pifaval2 = 10;
		pifaval3 = 30;
		pifaval4 = 40;
		pifaval5 = 60;
		pifaval6 = 90;
	}
	else /* spin echo */
	{
		pifanub = 0; /* turn off flip angle buttons */
	}

	/* label TE buttons */
	pite1nub = 63; /* apparently a bit mask */
	pite1val2 = PSD_MINIMUMTE;
	pite1val3 = 20ms;
	pite1val4 = 30ms;
	pite1val5 = 40ms;
	pite1val6 = 50ms;

	/* label FOV buttons */
	cvdef(opfov, 200);
	opfov = 200;
	pifovval2 = 200;
	pifovval3 = 220;
	pifovval4 = 240;
	pifovval5 = 360;
	pifovval6 = 480;

	/* turn on user cv page */
	piuset = use0+use1+use2+use3+use4+use5+use6+use10+use11+use12+use13+use14+use15+use16;
	pititle = 1;
	cvdesc(pititle, "Spiral User CV page");

	piuset += use20+use21+use22+use23+use24;
	/*
	   pititle = 2;
	   cvdesc(pititle, "PCASL User CV page");
	 */

	cvdesc(opuser0, "# of interleaves");
	cvdef(opuser0, 1.0);
	opuser0 = 1;
	cvmin(opuser0, 1.0);
	cvmax(opuser0, 128.0);
	nl = opuser0;

	cvdesc(opuser1, "number of temporal frames");
	cvdef(opuser1, 1);
	opuser1 = 4;
	cvmin(opuser1, 1);
	cvmax(opuser1, 16384);
	nframes = opuser1;

	cvdesc(opuser2, "Ext trig: (0)none, start (1)scan, (2)frame, (3)by file");
	cvdef(opuser2, 0);
	opuser2 = 0;
	cvmin(opuser2, 0);
	cvmax(opuser2, 3);
	gating = (opuser2==1 || opuser2==2)? TRIG_AUX : TRIG_INTERN;        /* aux trig */

	cvdesc(opuser3, "Recon script_number or none(0)");
	cvdef(opuser3, 0);
	opuser3 = 0;
	cvmin(opuser3, 0);
	cvmax(opuser3, 999);

	cvdesc(opuser4, "Cluster slice acquisition (1) or not(0)");
	cvdef(opuser4, 0);
	opuser4 = 0;
	cvmin(opuser4, 0);
	cvmax(opuser4, 1);
	scluster = opuser4;

	cvdesc(opuser5, "Number of extra shots before data acq");
	cvdef(opuser5, 2);
	opuser5 = 2;
	cvmin(opuser5, 0);
	cvmax(opuser5, 100);
	nextra = opuser5;

	cvdesc(opuser6, "gtype: (0) spiral out, (1) spiral in, (2) both");
	cvdef(opuser6, 2);
	opuser6 = 0;
	cvmin(opuser6, 0);
	cvmax(opuser6, 2);
	gtype = opuser6;
	bmapnav = (gtype == 1)? 0:1;		/* for now */

	cvdesc(opuser10, "Record physio data (1) or not(0)");
	cvdef(opuser10, 0);
	opuser10 = 0;
	cvmin(opuser10, 0);
	cvmax(opuser10, 1);

	/*	
		cvdesc(opuser11, "Realtime mode (1) or not(0)");
		cvdef(opuser11, 0);
		opuser11 = 0;
		cvmin(opuser11, 0);
		cvmax(opuser11, 1);
	 */

	cvdesc(opuser11,"Collect Fiel Map scan? (0=No, 1=Yes");
	cvmin(opuser11,0);
	cvmax(opuser11,1);
	cvdef(opuser11,0);
	opuser11 = 0;

	domap = opuser11;


	cvdesc(opuser12, "short rf pulse (1) or not(0)");
	cvdef(opuser12, 0);
	opuser12 = 0;
	cvmin(opuser12, 0);
	cvmax(opuser12, 1);
	if(opuser12 == 1)  {
		cycrf1 = 2;
		pwrf1 = 3200;
	}  else  {
		cycrf1 = 4;
		pwrf1 = 6400;
	}

	cvdesc(opuser13, "vd (1) or cd (0)");
	cvdef(opuser13, 0);
	opuser13 = 0;
	cvmin(opuser13, 0);
	cvmax(opuser13, 1);
	vdflag = opuser13;

	cvdesc(opuser14, "alpha");
	cvdef(opuser14, 3.6);
	opuser14 = 3.6;
	cvmin(opuser14, 1.01);
	cvmax(opuser14, 10);
	alpha = opuser14;

	cvdesc(opuser15, "kmaxfrac");
	cvdef(opuser15, 0.5);
	opuser15 = 0.5;
	cvmin(opuser15, 0.05);
	cvmax(opuser15, 0.95);
	kmaxfrac = opuser15;
	
	cvdesc(opuser16, "Background Suppression? (0=no, 1=yes)");
	cvdef(opuser16, 0);
	opuser16 = 0;
	cvmin(opuser16, 0.);
	cvmax(opuser16, 1.0);
	doBS = opuser16;

	
	cvdesc(opuser17, "Target Velocity of the blood to be labeled. (cm/s)");
	cvdef(opuser17, 0);
	opuser17 = 0;
	cvmin(opuser17, -300);
	cvmax(opuser17, 300);
	vel_target = opuser17;

	
	cvdesc(opuser19,"Gradient amplitude for VEL Selective Inversion");
	cvmin(opuser19,0);
	cvmax(opuser19,4);
	cvdef(opuser19,3);
	opuser19 = 2.0;
	vsi_Gmax = opuser19;
	vsi_Gcontrol = vsi_Gmax;
	
	
	cvdesc(opuser20,"Arterial Suppression Pulse (0=NO  1=YES)");
	cvmin(opuser20,0);
	cvmax(opuser20, 1);
	cvdef(opuser20,0);
	opuser20 = 0;
	doArtSup = opuser20;
	
	
	cvdesc(opuser21,"Delay between tagging period and imaging time (ms)");
	cvmin(opuser21,10);
	cvmax(opuser21, 5000);
	cvdef(opuser21,1400);
	opuser21 = 1400;
	t_delay = RUP_GRD(opuser21 * 1000);

	cvdesc(opuser22,"Delay between Arterial Suppression and imaging (ms)");
	cvmin(opuser22,10);
	cvmax(opuser22, 200);
	cvdef(opuser22,100);
	opuser22 = 100;
	AStime = RUP_GRD(opuser22*1000);
	

	cvdesc(opuser23,"Choice of labeling pulse");
	cvmin(opuser23,0);
	cvmax(opuser23, 999999);
	cvdef(opuser23,15999);
	opuser23 = 12360;
	vsi_train_len = opuser23;
	
	
	cvdesc(opuser24,"Do GRAPPA-Z (0=NO  1=YES)");
	cvmin(opuser24,0);
	cvmax(opuser24, 1);
	cvdef(opuser24,0);
	opuser24 = 0;
	doZgrappa = opuser24;

	
	cvdesc(opuser26,"VSI control image type (1) flip the magnitude  (2) no RF  (3) flip phase (4) flip grads (5) no VelSel grads");
	cvmin(opuser26,0);
	cvmax(opuser26, 10);
	cvdef(opuser26,6);
	opuser26 = 6;
	vsi_controlmode = opuser26;	


	slord = TYPNORMORDER;

	piamnub = 7;
	pixresnub = 3;
	cvmin(opxres, 16);
	cvdef(opxres, 64);
	opxres = 64;
	pixresval2 = 64;
	pixresval3 = 96;
	pixresval4 = 128;
	piyresnub = 0;
	

	/* Duration of whole labeling period : */
	/* LHG 11.7.14:   now includes a gap between five VSI bursts */
	t_tag = vsi_Ncycles*(vsi_timegap + astseqtime );
	t_tag = astseqtime ;
	astseqtime = pw_vsitag1 + mytpre;
	

	/* ************************************************
	   RF Scaling
	   Scale Pulses to the peak B1 in whole seq.
	 ********************************************** */
	/* 
	   maxB1Seq = 0.0;
	   for (entry=0; entry < MAX_ENTRY_POINTS; entry++) {
	   if (peakB1(&maxB1[entry], entry, RF_FREE, rfpulse) == FAILURE) {
	   epic_error(use_ermes,"peakB1 failed",EM_PSD_SUPPORT_FAILURE,1,STRING_ARG,"peakB1");
	   return FAILURE;
	   }
	   if (maxB1[entry] > maxB1Seq)
	   maxB1Seq = maxB1[entry];
	   }
	 */
	if( findMaxB1Seq(&maxB1Seq, maxB1, MAX_ENTRY_POINTS, rfpulse, RF_FREE) == FAILURE ) {
		epic_error(use_ermes,supfailfmt,EM_PSD_SUPPORT_FAILURE,EE_ARGS(1),STRING_ARG, "findMaxB1Seq");
		return FAILURE;
	}

	/*---------------------------------------------------------*/
	/* LHG 11/22/2010 set up the slab select trapezoid stuff here for 3D*/
	pw_rf1 = 5000;
	res_rf1 = 500;
	/* calculate the bandwidth of the min phase pulse: LHG 11/22/2010 */
	/* the rf3dmin pulse has 18 lobes */
	rf1_bw = 1000000.0 * 18.0 / (float)(pw_rf1);
	pw_gzrf1 = pw_rf1;
	pw_gzrf1a = myramptime;
	pw_gzrf1d = myramptime;
	pw_gzrf1r = 200;   /* LHG 11.30.12 - changing the duration of the refocuser */
	pw_gzrf1ra = myramptime;
	pw_gzrf1rd = myramptime;
	zfov = (double)opslquant*opslthick/10;

	if(doZgrappa) slab_fraction = 0.75;

	a_gzrf1 =  rf1_bw / (slab_fraction*zfov*GAMMA_H1/(2*3.14159)) ;
	area_gzrf1 = (double)a_gzrf1 * (pw_gzrf1 + pw_gzrf1d);
	area_gzrf1r = -(double)area_gzrf1 * (1-rf1maxTime);
	a_gzrf1r = (double)area_gzrf1r/(pw_gzrf1r + pw_gzrf1rd);

	/* ----------------------------------------------------------------------------------------*/
	/* Calculate the  z phase-encode gradient params */
	ampGzstep = 1 / ((pwGzstep + myramptime) *1e-6 * GAMMA_H1/(2*PI) * zfov);  /* z phase-ecncode readout gradient for this zfov G/cm*/

	if (doZgrappa) 
		ampGzstep *= 2;
	
	areaGzstep = ampGzstep * (pwGzstep + myramptime) ;        /* max area of z encoding gradient in G/cm*usec */

	/* get the specs of the z gradient */
	gettarget(&target, ZGRAD, &loggrd);
	getramptime(&rtime, &ftime, ZGRAD, &loggrd);

	/* previous code from JFN : 
	 *  Use ampwpgrad macro to design a trapezoid given a desired area and the hardware specs 
	 if (amppwgrad(	gzphasearea, target, 
	 0.0, 0.0, 
	 rtime, 
	 MIN_PLATEAU_TIME, 
	 &a_gzphase, &pw_gzphasea,  
	 &pw_gzphase,  &pw_gzphased) 
	 == FAILURE)
	 return FAILURE; 
	 */
	/* calculate the amplitude of the max. Gz encode gradient in DAC units */
	iampGzstep = 2 * ((ampGzstep / target) * MAX_PG_IAMP / 2.0);

	/* Now stick those values into the actual trapezoid: */
	a_gzphase1 = -ampGzstep;
	a_gzphase2 = ampGzstep;

	pw_gzphase1= pwGzstep;
	pw_gzphase1a= myramptime;
	pw_gzphase1d= myramptime;
	pw_gzphase2= pwGzstep;
	pw_gzphase2a= myramptime;
	pw_gzphase2d= myramptime;

	fprintf(stderr, "predownload(): ampGzstep = %.3f, iampGzstep = %d \n", ampGzstep, iampGzstep);
	/* ----------------------------------------------------------------------------------------*/
	pw_gz0=4000;
	pw_gz0a=400;
	pw_gz0d=400;
	a_gz0 = 2.0;

	/* LHG 12/20/12
	   The crushers for the Background suppression 90 are the same as for the fat sat pulses */	
	pw_gzBS0=4000;
	pw_gzBS0a=400;
	pw_gzBS0d=400;
	a_gzBS0 = 2.0;
	a_gzBS0 = 0;  /* LHG 1/9/13:  simulations indicate that this doesnt' help much. */

	/* Calculations for Background suppression pulses: */
	/* 
	   % LHG 12/14/12        Determining the amplitude of the SECH pulses 
	   % First find the right sech pulse file!!

	   s1 =loadrho('myhsec1t.rho');
	   s1th = loadrho('myhsec1t.theta'); 
	   s1=s1/max(s1);
	   s1 = s1.*exp(-i*s1th);
	   area_sech = sum(abs(s1))* 8.0 /length(s1)

	   % Area of the actual prescan pulse:  a Hanning windowed sinc
	   % when a_rf2mps=1, that means you get a 180

	   t=linspace(-1600,1600,3200)';
	   s1=sinc(t/800).*hanning(3200);
	   plot(s1)
	   area_sinc=sum(abs(s1(:))) * 3.2 / length(s1);

	   rfscalesech = area_sinc/area_sech
	   ====> 0.2153
	 */

	rfscalesech = 0.2153;
	rfscalesech = 1.0;  /* Empirically on sphere phantom, I didn't get an inversion unutl 0.65) */


	pw_BS1rf = 8000;
	a_BS1rf = rfscalesech * (float)doBS ; /* the sech is used to generate a 180 */ 
	res_BS1rf = 2000; /* the number of points in the file */
	res_BS1rf = 500; /* the number of points in the file */

	pw_BS2rf = pw_BS1rf;;
	a_BS2rf = a_BS1rf;
	res_BS2rf = res_BS1rf;

	pw_BS3rf = pw_BS1rf;;
	a_BS3rf = a_BS1rf;
	res_BS3rf = res_BS1rf;

	/*remove the pre label inversion pulse in the BackSUpp train */
	a_BS3rf = 0;
	
	/*---------------------------*/
	/* LHG 6.29.12 - replace their tmin calculation and calculation of psdseqtime:
	 old: 
	 tmin = RUP_GRD(tlead +pw_rf0 + pw_gz0 + 2*pw_gz0a + pw_gzrf1a +
	 pw_rf1/2 + opte + pw_gx + daqdel + mapdel + pw_gzspoil +
	 2*pw_gzspoila);
	 */
	if (doBS==0){
		BS1_time = 10000;
		BS2_time = 10000;
	}

	BStime = BS1_time + BS2_time + pw_BS1rf + pw_BS2rf + pw_BS0rf + pw_gzBS0 + 2*pw_gzBS0a; 

	fatsattime = pw_rf0 + 2*pw_gz0d +pw_gz0 + 2*timessi;
	fuzz = 10000;

	/* if (doBS)  t_delay = t_delay - BStime;  LHG 10/4/16 reomivng this */

	astseqtime2 = t_tag + t_delay;  /*duration of tagging + postabeling delay */
	
	tmin = RUP_GRD(tlead +tdel+
		tpre + 0.5*pw_gzrf1 + pw_gzrf1a +
		opte + pw_gx + mapdel +
		pw_gzspoil + 2*pw_gzspoila +
		timessi);

	psdseqtime =  RUP_GRD(tmin) ;
	
	seqtr = tmin * opslquant;
	
	fprintf(stderr, "\npsdseqtime, seqtr, tmin, optr=  %d   %d  %d  %d", psdseqtime,  seqtr,tmin,optr);
	
	t_adjust = RUP_GRD(optr - psdseqtime*opslquant - t_tag - t_delay - tpre - pw_BS3rf);                 /* LHG 11.29.10 */
	fprintf(stderr, "\nt_adjust = %d , astseqtime= %d \n", t_adjust, astseqtime);
	/*---------------------------*/
	

	fprintf(stderr,"\nEnded cveval");



@inline Prescan.e PScveval

	return SUCCESS;
}

int cvcheck()
{
	return SUCCESS;
}

int predownload()
{
	int pdi, pdj;
	int i;
	float max_rbw;
	FILE *fpin;
	/* int  BIR_len = 7000;*/

	fprintf(stderr,"\nPredownload stuff:");

	/* Load the BIR-8 pulses for Arterial Suppression */
	Npoints = read_vsi_pulse(BIR_mag, BIR_phs, BIR_grad, BIR_len);
	
	/* LHG 9.27.16 : Timing calculation from the VSI sequence... */
	/* Npoints and vsi_train_len should come out the same ... good sanity check */
	Npoints = read_vsi_pulse(vsi_pulse_mag, vsi_pulse_phs, vsi_pulse_grad, vsi_train_len);
	fprintf(stderr,"\nNpoints: %d ... n. segments = %d ", Npoints, (Npoints*4 - 4000)/2000);

	/* note that in some cases, we may want a different phase in the control pulses */
	for (i=0; i<Npoints ; i++){
		vsi_pulse_ctl_phs[i] = vsi_pulse_phs[i];
	}
	/* in some cases, we want the abs of the gradient for the control pulses (redundant code ) */
	for (i=0; i<Npoints ; i++){
		vsi_pulse_ctl_grad[i] = vsi_pulse_grad[i];
		if (vsi_controlmode==6)
			vsi_pulse_ctl_grad[i] = abs(vsi_pulse_grad[i]);
	}


	/*But ... in the label case, we update phase wave form for velocity selectivity */
	calc_vsi_phs_from_velocity (vsi_pulse_mag, vsi_pulse_phs,  vsi_pulse_grad, vel_target, vsi_train_len, vsi_Gmax);  
	
	astseqtime = pw_vsitag1 + mytpre ;


	opnecho = 1;
	cvdef(opnecho, 1);
	cvmin(opnecho, 1);
	cvmax(opnecho, 1);

	/* image header variables set for correct annotation */
	ihflip = opflip;
	ihnex = opnex;
	ihtr = optr;

	/*  make the end dead time = TR */
	endtime = optr;

	/*  figure out some params depending on system/grads  */
	/* gamp = loggrd.xfs;  NEED to reduce to avoid PNS  */
	/*  gslew = cfsrmode;   can't use it all  */

	if(fpin=fopen("daqdeloff.sprl", "r")) {     /* look for ext. file */
		fscanf(fpin, "%f", &daqdeloff);
		fclose(fpin);
	}

	/*  try to make tr check */

	/*---------------------------*/
	/* LHG 6.29.12 - replace their tmin calculation and calculation of psdseqtime:
old: 
tmin = RUP_GRD(tlead +pw_rf0 + pw_gz0 + 2*pw_gz0a + pw_gzrf1a +
pw_rf1/2 + opte + pw_gx + daqdel + mapdel + pw_gzspoil +
2*pw_gzspoila);
	 */
	if (doBS==0){
		BS1_time = 10000;
		BS2_time = 10000;
	}


	BStime = BS1_time + BS2_time + pw_BS1rf + pw_BS2rf + pw_BS0rf + pw_gzBS0 + 2*pw_gzBS0a; 

	fatsattime = pw_rf0 + 2*pw_gz0d +pw_gz0 + 2*timessi;
	fuzz = 10000;

	/*if (doBS)  t_delay = t_delay - BStime;  LHG 10/4/16 : removing this: */

	astseqtime2 = t_tag + t_delay;  /*duration of tagging + postabeling delay */
	
	tmin = RUP_GRD(tlead +tdel+
		tpre + 0.5*pw_gzrf1 + pw_gzrf1a +
		opte + pw_gx + mapdel +
		pw_gzspoil + 2*pw_gzspoila +
		timessi);

	psdseqtime =  RUP_GRD(tmin) ;

	seqtr = tmin * opslquant;

	fprintf(stderr, "\npsdseqtime, seqtr, tmin, optr=  %d   %d  %d  %d", psdseqtime,  seqtr,tmin,optr);
	
	t_adjust = RUP_GRD(optr - psdseqtime*opslquant - t_tag - t_delay - tpre - pw_BS3rf);                 /* LHG 11.29.10 */
	fprintf(stderr, "\nt_adjust = %d , astseqtime= %d \n", t_adjust, astseqtime);
	/*---------------------------*/

	if(gtype==1)  tmin -= pw_gx;
	if(gtype==2)  tmin += espace;
	if(nl>1) tmin += pw_gxspoil + 2*pw_gxspoila; 
	/* tmin += pw_gxspoil + 2*pw_gxspoila; */
	seqtr = tmin*opslquant;
	/* lhg 9.21.12:  include time for labeling in calculation of psdseqtime*/
	/* psdseqtime = (scluster==0)? RUP_GRD(optr/opslquant) : tmin + 452us;	/* add some time */

	/* LHG 7.19.17 */ 
	psdseqtime =  RUP_GRD(tmin) ;

	psdtottime = psdseqtime*opslquant;
	trerror = psdtottime - optr;
	nerror = fabs((double)trerror)/GRAD_UPDATE_TIME;
	printf("\nseqtr, tmin, psdseqtime, psdtottime, trerror=  %d  %d  %d  %d  %d\n", seqtr,tmin,psdseqtime,psdtottime,trerror);
	if(seqtr > optr)  {
		epic_error(use_ermes,"Oops! optr must be > %.1f ms.\n", EM_PSD_SUPPORT_FAILURE,1,FLOAT_ARG,seqtr/1000.0);
		return FAILURE;
	}

	/* adjust grads, tsp, to reflect bandwidth  */

	gfov = opfov*.1;
	calcvalidrbw((double)oprbw,&bandwidth,&max_rbw,&decimation,OVERWRITE_OPRBW,0);
	gmax = 2e3*bandwidth/(GAM*gfov);    /* max grad allowed for BW */
	gamp = MIN(cfxfs, gmax);
	tsp = TSP*decimation;

	/* readout gradients - calculate 'em here  */
	fprintf(stderr, "\nCalling genspiral  ... ");

	tadmax = tsp*GRESMAX;		/* max readout */
	if(vdflag==0) {
		if(!(genspiral(gfov, opxres, tadmax*1.e-6, 4.e-6))) {
			epic_error(use_ermes,"Oops! genspiral failed\n",0,0);
			return FAILURE;
		}
	}
	else{
		if(!(genspiralcdvd(gfov, opxres, tadmax*1.e-6, 4.e-6, alpha, kmaxfrac))) {
			epic_error(use_ermes,"Oops! genspiral CDVD failed\n",0,0);
			return FAILURE;
		}
	}
	if(gram_duty() == FAILURE) return FAILURE;

#include "predownload.in"

	/* set up RF pulse  */
	/* a_rf1 = cyc_rf1*(6400us/pw_rf1)*opflip/360; */
	/* 
	   LHG 11/29/12        Determining the amplitude of the pcasl pulse 
	   if the prescan calibration was done with rf1 (rf3d.rho -the time shifted 9 cycle sinc)
	   then a 90 degree flip  wass obtained when the  amplitude of rf1 == 0.9.

	   s1 =loadrho('rf3d.rho'); 
	   s1=s1/max(s1);
	   area_rf3d = sum(abs(s1))*5.0/length(s1)

	   % Area of the actual prescan pulse:  a Hanning windowed sinc
	   % when a_rf2mps=1, that means you get a 180

	   t=linspace(-1600,1600,3200)';
	   s1=sinc(t/800).*hanning(3200);
	   plot(s1)
	   area_sinc=sum(abs(s1(:))) * 3.2 / length(s1);

	   rfscale180 = area_sinc/area_rf3d
	   ====>  1.3125

EMPIRICALLY: we know that for this pulse, if pw_rf1=5ms, then a_rf1=0.9 gives you a 90 
so don't change the duration of pw_rf1 without adjusting this calculation ...*/

	/*
	a_rf1 = opflip *1.3125/180;
	ia_rf1 = a_rf1 * max_pg_iamp;
	*/
	a_rf1 = opflip/180 * rf1_RFmax / dummy_RFmax;
	a_vsitag1 = vsi_RFmax / dummy_RFmax;
	a_vsictl1 = vsi_RFmax / dummy_RFmax;
	ia_rf1 = a_rf1 * max_pg_iamp;
	ia_vsitag1 = a_vsitag1 * max_pg_iamp;
	ia_vsictl1 = a_vsictl1 * max_pg_iamp;


	/* set up sat pulse */
	satBW = (cffield==15000)? 220 : 440;
	satoff = (cffield==15000)? -260 : -520;  /* offset  down a little */
	arf0 = 0.5*satBW/1250.0;
	pwrf0 = RUP_GRD((int)(4s*cyc_rf0/satBW));
	fprintf(stderr,"\nsatBW, satoff, arf0, pwrf0 = %f  %f  %f  %d\n", satBW, satoff,
			arf0, pwrf0);

	pw_rf0 = pwrf0;
	a_rf0 = arf0;		

	/* LGH 12/20/12:  
	   The Backgroup suppression 90deg pulse is the same as the fat sat pulse, but without the freq. offset 
	   (it's also a sinc with same duration and amplitude)
	 */
	pw_BS0rf = pw_rf0;
	a_BS0rf = a_rf0 * (float)doBS;
	a_BS0rf = 0; /* LHG 1/1/13: siimulations indicate that this pulse is not useful */ 		

	/*  spoilers to return the readout grads to a constant state  */

	if(nl>1)  {		/* used to be only if interleaving (nl>1)*/
		pw_gxspoila = 2*cfrmp2xfs;
		pw_gxspoild = 2*cfrmp2xfs;
		pw_gxspoil = 500us;
		a_gxspoil = 0.5*loggrd.yfs;
		pw_gyspoila = 2*cfrmp2yfs;
		pw_gyspoild = 2*cfrmp2yfs;
		pw_gyspoil = 500us;
		a_gyspoil = 0.5*loggrd.yfs;
	}

	nbang = (nextra+nl*nframes)*opslquant;
	ndisdaq = MIN(nextra, 6e6/optr);	/* for prescan only  */

	/* set up clock */
	if ((exist(opcgate) == PSD_ON) && existcv(opcgate))
	{
		pidmode = PSD_CLOCK_CARDIAC;
		piviews = nextra+nl*nframes;
		piclckcnt = ophrep;
		pitscan = (float)(optr)*nbang/opslquant;
		pitslice = psdseqtime;
	}
	else
	{
		pidmode = PSD_CLOCK_NORM;
		pitslice = psdseqtime;
		pitscan = (float)(optr)*nbang/opslquant;
	}

	/* initialize slice sel spoiler gradient. */

	/*  LHG 12/14/12 : replace these for the values we had in the old Signa 
	    pw_gzspoil = 800us;
	    a_gzspoil = .5*loggrd.zfs;
	 */
	pw_gzspoil = 1000;
	a_gzspoil = 0.7*loggrd.zfs;
	pw_gzspoila = cfrmp2zfs;
	pw_gzspoild = cfrmp2zfs;


	/*-------------------------------------------------*/
	minte = pw_rf1/2+pw_gzrf1d + pw_gzrf1r+2*pw_gzrf1ra;
	minte += 2*pw_gzphase1a + pw_gzphase1;

	/* LHG 6.29.12 - these are the flow spoilers: */
	if(mycrush==1){
		fprintf(stderr, "\nSetting Flow Crusher parameters ...");
		pw_gzsp3 = 3000;
		pw_gzsp3a = 300us;
		pw_gzsp3d = pw_gzsp3a ;
		pw_gzsp4 = pw_gzsp3;
		pw_gzsp4a = pw_gzsp3a;
		pw_gzsp4d = pw_gzsp3a;
		a_gzsp3 = 4.0;
		a_gzsp4 = -a_gzsp3;
/*
		pw_gxsp3 = pw_gzsp3;
		pw_gxsp3a = pw_gzsp3a;
		pw_gxsp3d = pw_gzsp3a;
		pw_gxsp4 = pw_gzsp3;
		pw_gxsp4a = pw_gzsp3a;
		pw_gxsp4d = pw_gzsp3a;
		a_gxsp3 = 0.0;
		a_gxsp4 = -a_gxsp3;

		pw_gysp3 = pw_gzsp3;
		pw_gysp3a = pw_gzsp3a;
		pw_gysp3d = pw_gzsp3a;
		pw_gysp4 = pw_gzsp3;
		pw_gysp4a = pw_gzsp3a;
		pw_gysp4d = pw_gzsp3a;
		a_gysp3 = 0.0;
		a_gysp4 = -a_gysp3;

*/
		delta1 = 50us;

		/* LHG 2/1/11  adjust minte to account for the flow crushers */
		minte += 2*pw_gzsp3 + 4*pw_gzsp3a + delta1 + 50 ;
		fprintf(stderr,"\nmin. TE after flow crushers: %d \n", minte);
	}

	/*-------------------------------------------------*/

	if(gtype>0) minte += pw_gx + pwgpre;
	cvmin(opte, minte);
	cvdef(opte, minte);
	if ((exist(opautote) == PSD_MINTE)||(exist(opautote) == PSD_MINTEFULL))
		opte = minte;
	ihte1 = opte;

@inline loadrheader.e rheaderinit

	rhimsize = 64;
	while(rhimsize < opxres)  rhimsize *= 2;
	if(off_fov==1)  {
		rhxoff = 0;
		rhyoff = 0;
	}
	else  {
		rhxoff = rhimsize*scan_info[0].oprloc/opfov;
		rhyoff = rhimsize*scan_info[0].opphasoff/opfov;
	}
	rhuser0 = gfov; 	/* grad design fov in cm  */
	rhuser1 = nframes; 	/* number of temporal frames */
	rhuser3 = opxres;		/*  equiv resolution  */
	rhuser4 = nl;   	/* number of interleaves  */
	rhuser5 = gtype;   	/* gradient wavewform rev  */
	rhuser6 = gamp; 	/* grad design ampl  */
	rhuser7 = gslew; 	/* grad slew rate  */
	rhuser9 = bmapnav;  /* if set, keeps dacq fixed in time for ifr=0 */
	rhuser10 = 0;	/* ngap between concat acq's */
	rhuser11 = 0;
	rhuser12 = 0;  	/* receiver center freq, kHz */
	rhuser13 = tsp; 	/* a/d sampling time (us) */
	rhuser14 = (nframes>1 && oppseq==2)? domap : 0;  /* do it if GRE */
	rhuser15 = mapdel;  /* echo difference (us) */
	rhuser21 = alpha;	/* for vd trajectory  */
	rhuser22 = kmaxfrac;
	rhuser23 = vdflag;
	rhuser30 = opslthick;
	rhuser31 = opslspace;

	rhbline = 0;
	rhnecho = (gtype==2)? 2:1;
	rhnslices = opslquant;
	rhrcctrl = 1;    		/* lx wants to create images.  */
	rhexecctrl = 2;  		/* just save the raw data */
	/* rhdacqctrl = 0;		LHG: 11/29/13 - take this out in DV 23 */
	rhrecon = opuser3;		/* invoke son of recon */

	/*  LHG 9/8/14 - disable this: 	
	    if(opuser11==1.0)           // turn on/off RDS for realtime 
	    rhtype1 |= (0x00080000 | 0x00100000);
	    else
	    rhtype1 &=  0x0007FFFF;
	 */


	acqs = 1;
	slquant1 = rhnslices;

	/* for straight sequential order of slices. */
	if (!orderslice(slord, rhnslices, slquant1, gating))
		epic_error(use_ermes,"orderslice call failed",0,0);

	/* initialize copy of original rotation matrices */
	for (pdi = 0; pdi < rhnslices; pdi++)
		for (pdj = 0; pdj < 9; pdj++)
			savrot[pdi][pdj] = rsprot[pdi][pdj];
	scalerotmats(rsprot, &loggrd, &phygrd, rhnslices, 0);

	/* save stuff for maxwell correction */
	rhmaxcoef1a = rsprot[0][0]/(float)cfxfull;	/* save x rotator */
	rhmaxcoef1b = rsprot[0][1]/(float)cfxfull;
	rhmaxcoef2a = rsprot[0][3]/(float)cfyfull; /* y  */
	rhmaxcoef2b = rsprot[0][4]/(float)cfyfull;
	rhmaxcoef3a = rsprot[0][6]/(float)cfzfull; /* z  */
	rhmaxcoef3b = rsprot[0][7]/(float)cfzfull;

	rhdab0s = cfrecvst;
	rhdab0e = cfrecvend;

	if(entrytabinit(entry_point_table, (int)ENTRY_POINT_MAX)
			== FAILURE) {
		epic_error(use_ermes,"Can't initialize entry point table.",0,0);
		return FAILURE;
	}

	/* set up receiver */

	initfilter();

	cvmax(rhfrsize, 32768);		/* for now  */
	rhfrsize = (res_gx-RES_GRAMP)*4us/tsp;      /* num points sampled */
	rhfrsize = 4*(rhfrsize/4);          /* wants to be divisible by 4 */
	total_views=2*((nl*nframes+1)/2);  /* has to be an even number */
	cvmax(rhnframes, total_views);
	rhnframes = total_views;

	if (calcfilter( &echo1_filt,bandwidth,rhfrsize,OVERWRITE_OPRBW ) == FAILURE) {
		epic_error(use_ermes,"%s failed",EM_PSD_SUPPORT_FAILURE,1,STRING_ARG,"calcfilter");
		return FAILURE;
	}
	setfilter( &echo1_filt, SCAN );
	filter_echo1 = echo1_filt.fslot;

	if (calcfilter( &aps2_filt,bandwidth,psfrsize,OVERWRITE_OPRBW ) == FAILURE) {
		epic_error(use_ermes,"%s failed",EM_PSD_SUPPORT_FAILURE,1,STRING_ARG,"calcfilter");
		return FAILURE;
	}
	setfilter( &aps2_filt, PRESCAN );
	filter_aps2 = aps2_filt.fslot;

	rhrawsize = 2*rhptsize*rhfrsize*(rhnframes+1)*rhnslices*rhnecho;
	numrecv = rhdab0e-rhdab0s+1;
	if((float)rhrawsize*numrecv > cftpssize-2e7)  {      /* reserve 20 MB */
		epic_error(use_ermes,"Oops! Tooo much memory requested.\n",0,0);
		return FAILURE;
	}

	daqdel = psd_grd_wait + daqdeloff + 0.5;

	strcpy(entry_point_table[L_SCAN].epname, "scan");
	entry_point_table[L_SCAN].epfastrec = 0;
	entry_point_table[L_SCAN].epstartrec = rhdab0s;
	entry_point_table[L_SCAN].ependrec = rhdab0e;
	entry_point_table[L_SCAN].epfilter = (unsigned char)echo1_filt.fslot;
	entry_point_table[L_SCAN].epprexres = rhfrsize;
	entry_point_table[L_SCAN].epxmtadd = txCoilInfo[getTxIndex(coilInfo[0])].coilAtten;
	entry_point_table[L_APS2] =
		entry_point_table[L_MPS2] =
		entry_point_table[L_SCAN];      /* copy scan into APS2 & MPS2 */
	strcpy(entry_point_table[L_APS2].epname,"aps2");
	strcpy(entry_point_table[L_MPS2].epname,"mps2");
	entry_point_table[L_APS2].epfilter = (unsigned char)aps2_filt.fslot;
	entry_point_table[L_MPS2].epfilter = (unsigned char)aps2_filt.fslot;
	entry_point_table[L_APS2].epprexres = psfrsize;
	entry_point_table[L_MPS2].epprexres = psfrsize;

	trigfreq = opslquant; /* default for TTL trig is every TR */

	/*  Some prescan stuff  */

	pislquant = opslquant;

@inline Prescan.e PSfilter
@inline Prescan.e PSpredownload

		phys_record_flag = opuser10;        /* put here so it isn't overwritten in predownload */
	phys_record_channelsel = 14;        /* PG wave,trig & RESP */

	return SUCCESS;
} /* End-Of-Predownload */

@inline Prescan.e PShost
#include "genspiral5.e"
#include "genspiralcdvd.e"
#include "gram_duty.e"
/*LHG 9.17.16 */
#include "read_vsi_pulse.e"

@rsp

int pre = 0; /* prescan flag */
short thamp;

CHAR *entry_name_list[ENTRY_POINT_MAX] = { "scan", "mps2", "aps2",
@inline Prescan.e PSeplist
};

WF_PULSE thetrec = INITPULSE;   /* fov offset */
WF_PULSE *thetrecintl;
WF_HW_WAVEFORM_PTR *thetrecintlp;
WF_PULSE thetrec2 = INITPULSE;   /* fov offset */
WF_PULSE *thetrecintl2;
WF_HW_WAVEFORM_PTR *thetrecintlp2;
int *fxmit, *frec;
int *rfphastab;
int *slcphastab;
int *slordtab;
short ibuf[GRESMAX]; 	/*  for intwave.h  */

/* LHG 9/26/12/ Tables for kz encoding (grad amplitudes and phases) */ 
int *rfphastab,*rfamptab;

@rspvar


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
@inline Prescan.e PSrspvar	/* For Prescan */
extern PSD_EXIT_ARG psdexitarg;

@pg

#include <epic_loadcvs.h>
long deadtime_core;

/*----------------------------------
  LHG 6.29.12
  --------------------------------*/
long 	deadtime_astcore;
long 	deadtime_controlcore;
long 	deadtime_nothingcore;
long	deadtime_tadjustcore;
long	deadtime_preBScore;
long	deadtime_tdelaycore;
long	deadtime_vsi_gapcore;

/* LHG 6.29.12:  declaration of pulse waveforms for the kz-phase encode pulses 
   WF_PULSE gzphase1a = INITPULSE; WF_PULSE gzphase1 = INITPULSE; WF_PULSE gzphase1d = INITPULSE;
   WF_PULSE gzphase2a = INITPULSE; WF_PULSE gzphase2 = INITPULSE; WF_PULSE gzphase2d = INITPULSE;
/* -----------------------------------------------------*/

STATUS pulsegen(void)
{
	int waitloc, waitloc2;
	char tstr[40];
	short *ts;
	int j, jj;
	float rdx, rdy;
	float cphi, sphi, x;

	sspinit(psd_board_type);

	/* LHG 12/14/12 Chem sat pulse (rf0) has been moved to the delay core */ 

	/*  LHG 6.29.12  - replacing 90 degree pulse with min phase pulse and new slice sel grads by JFN */
	fprintf(stderr, "\npulsegen: generating RF1 pulse from rf3d.rho (min phase) ... ");

	EXTWAVE(RHO, 
			rf1,
			RUP_GRD(tlead + pw_gzrf1a + psd_rf_wait),
			pw_rf1,
			a_rf1,
			res_rf1,
			rf3d.rho,
			,loggrd); 


	/* this may be an alternative ...
	   SLICESELZEXT2(	rf1,
	   RUP_GRD(pend(&gz0d, "gz0d",0) + pw_gzrf1a + psd_rf_wait),
	   pw_rf1,
	   opslthick*opslquant,
	   opflip,0,1,
	   NULL, res_rf1,
	   rf3d.rho,
	   RF1_SLOT,
	   TYPNDEF,
	   ,loggrd);
	 */

	fprintf(stderr, " ... done.");

	fprintf(stderr, "\npulsegen: generating slice select gradient  ... ");

	TRAPEZOID(ZGRAD,
			gzrf1,
			RUP_GRD(tlead + pw_gzrf1a)  ,
			0,
			0,loggrd);

	fprintf(stderr, " at time= %d ... done.", RUP_GRD(tlead + pw_gzrf1a));


	/* -----------------------------------------------------*/
	/*  LHG 6.29.12:   original  code for first excitation pulse:

	    SLICESELZ(rf1, pend(&gz0d,"gz0d",0)+pw_gzrf1a, (opslthick>=1.0)? pwrf1:9000, opslthick, opflip,cycrf1,, loggrd);
	    TRAPEZOID(ZGRAD, gz2, RUP_GRD(pend(&gzrf1d, "gzrf1d", 0)+ pw_gz2a), -.5*zref*a_gzrf1*(pw_rf1+pw_gzrf1d),1, loggrd);
	 */
	/* -----------------------------------------------------*/

	/*waitloc = pend(&gzrf1d,"gzrf1d",0);   /* the end slice select gradient */
	readpos = RUP_GRD(pmid(&rf1,"rf1",0) + opte); /* echo time in GRE option */

	daqpos = readpos;

	if(gtype>0)  {
		readpos -= gres*4us;
		daqpos -= (gres - nramp)*4us;
		SINUSOID(XGRAD, gxpre, readpos-pwgpre, pwgpre, prefudge*agxpre, 0, 0., 0.5, 0., loggrd);
		SINUSOID(YGRAD, gypre, readpos-pwgpre, pwgpre, prefudge*agypre, 0, 0., 0.5, 0., loggrd);
	}

	/* if oppseq != 1 ... gradient-echo case requires a refocuser for the slice select Z grad.
	   Note that in sprlio, it's called gz2.  Here it's gzrf1r */
	TRAPEZOID(ZGRAD,
			gzrf1r,
			pend(&gzrf1d, "gzrf1d",0) + pw_gzrf1ra,
			0,
			0,loggrd);

	fprintf(stderr, "\npulsegen: generating z phase-encode gradient... ");
	TRAPEZOID(ZGRAD,
			gzphase1,
			RUP_GRD(pend(&gzrf1rd, "gzrf1rd",0) + pw_gzphase1a),
			0, 0,loggrd);
	fprintf(stderr, "... done. " );


	/* waitloc = pendall(&gzrf1r,0); */
	waitloc = RUP_GRD(pendall(&gzphase1,0)); 

	if(gtype<2)  {
		fprintf(stderr, "\npulsegen: generating map waits ... ");
		WAIT(XGRAD, mapx, waitloc, 4us);
		WAIT(YGRAD, mapy, waitloc, 4us);
		WAIT(ZGRAD, mapz, waitloc, 4us);
		WAIT(THETA, mapt, waitloc, 4us);
		fprintf(stderr, "... done.  ");
	}

	if(mycrush==1)  /* Now put in the flow crushers for the GRE case */
	{

		fprintf(stderr, "\npulsegen: generating flow crusher gradients... ");
		/*
		TRAPEZOID(XGRAD,
				gxsp3,
				RUP_GRD(pend(&gzphase1d,"gzrf1rd",0)+pw_gxsp3a + 50 ),
				0,0,loggrd);
		TRAPEZOID(XGRAD,
				gxsp4,
				RUP_GRD(pend(&gxsp3d,"gxsp3d",0)+pw_gxsp4a + delta1),
				0,0,loggrd);
		fprintf(stderr, "\npulsegen: ... done with X ");
		TRAPEZOID(YGRAD,
				gysp3,
				RUP_GRD(pend(&gzphase1d,"gzrf1rd",0)+pw_gxsp3a + 50 ),
				0,0,loggrd);
		TRAPEZOID(YGRAD,
				gysp4,
				RUP_GRD(pend(&gysp3d,"gysp3d",0)+pw_gysp4a + delta1),
				0,0,loggrd);

		fprintf(stderr, "\npulsegen: ... done with Y ");
		*/
		TRAPEZOID(ZGRAD,
				gzsp3,
				RUP_GRD(pend(&gzphase1d,"gzphase1d",0) + pw_gzsp3a + 50) ,
				0,0,loggrd);
		TRAPEZOID(ZGRAD,
				gzsp4,
				RUP_GRD(pend(&gzsp3d,"gzsp3d",0) + pw_gzsp4a + delta1),
				0,0,loggrd);

		waitloc = pendall(&gzsp4,0);
		fprintf(stderr, "\npulsegen: ... done with Z ");
	}	

	/*-------------------------------------------------------*/


	fprintf(stderr, "\npulsegen: reading in the spiral gradient waveforms ... ");
	INTWAVE(XGRAD, 
			gx, 
			readpos, 
			gamp, gres, gres*4, Gx, 
			1, 
			loggrd);
	INTWAVE(YGRAD, 
			gy, 
			readpos, 
			gamp, gres, gres*4, Gy, 
			1, 
			loggrd);

	WAIT(SSP, maps1, waitloc, 4us);
	/* LHG 10/16/12 
	   ACQUIREDATA(echo1, daqpos+daqdel,,, ); */
	ACQUIREDATA(echo1, readpos+daqdel,,, ); 

	pulsename(&thetrec, "thetrec");
	createreserve(&thetrec, THETA, gres);
	createinstr(&thetrec, RUP_RF(readpos+daqdel+thetdeloff), pw_gx, max_pg_iamp);

	waitloc2 = pend(&gx,"gx",0);

	/* spiral out as echo 2 */

	if(gtype==2)  {
		WAIT(XGRAD, mapx2, waitloc2, 4us);
		WAIT(YGRAD, mapy2, waitloc2, 4us);
		WAIT(ZGRAD, mapz2, waitloc2, 4us);
		WAIT(THETA, mapt2, pend(&thetrec,"thetrec",0), 4us);
		readpos = RUP_GRD(pend(&gx,"gx",0) + espace);
		INTWAVE(XGRAD, gx2, readpos, gamp, gres, gres*4, Gx, -1, loggrd);
		INTWAVE(YGRAD, gy2, readpos, gamp, gres, gres*4, Gy, -1, loggrd);
		WAIT(SSP, maps2, waitloc2, 4us);
		ACQUIREDATA(echo2, readpos+daqdel,,, );

		pulsename(&thetrec2, "thetrec2");
		createreserve(&thetrec2, THETA, gres);
		createinstr(&thetrec2, RUP_RF(readpos+daqdel+thetdeloff), pw_gx, max_pg_iamp);

		waitloc2 = pend(&gx2,"gx2",0);
	}

	WAIT(SSP, maps3, waitloc2+daqdel, mapdel+4us);

	/* Fov offset  */
	ts = AllocNode(gres*sizeof(short));
	thetrecintl = (WF_PULSE *) AllocNode(nl*sizeof(WF_PULSE));
	thetrecintlp = (WF_HW_WAVEFORM_PTR *) AllocNode(nl*sizeof(WF_HW_WAVEFORM_PTR));
	if(gtype==2)  {
		thetrecintl2 = (WF_PULSE *) AllocNode(nl*sizeof(WF_PULSE));
		thetrecintlp2 = (WF_HW_WAVEFORM_PTR *) AllocNode(nl*sizeof(WF_HW_WAVEFORM_PTR));
	}
	rdx = rsp_info[0].rsprloc;
	rdy = rsp_info[0].rspphasoff;
	for (i = 0; i < nl; i++) {
		cphi = cos(2.0*3.14159265*i / (double) nl);
		sphi = sin(2.0*3.14159265*i / (double) nl);
		x = off_fov*2.*FS_PI*gamp*GAM*GRAD_UPDATE_TIME*1e-6/(10.0*max_pg_wamp);
		sprintf(tstr, "thetrecint_%d", i);
		pulsename(&thetrecintl[i], tstr);
		createreserve(&thetrecintl[i], THETA, gres);
		/* The demodulation waveform necessary for this interleaf */
		if(gtype==0)  {
			ts[0] = 0;
			for (j = 1; j < gres; j++)  {
				ts[j] = (short) (ts[j-1] + x*((cphi*Gx[j] - sphi*Gy[j])*rdx +
							(cphi*Gy[j] + sphi*Gx[j])*rdy)) & ~WEOS_BIT;
			}
		}  else  {
			ts[gres-1] = 0;
			for (j = gres-2; j >=0; j-=1)  {
				ts[j] = (short) (ts[j+1] - x*((cphi*Gx[j] - sphi*Gy[j])*rdx +
							(cphi*Gy[j] + sphi*Gx[j])*rdy)) & ~WEOS_BIT;
			}
		}
		ts[gres - 1] |= WEOS_BIT;
		movewaveimm(ts, &thetrecintl[i], (int) 0, gres, TOHARDWARE);
		thetrecintlp[i] = thetrecintl[i].wave_addr;
		if(gtype==2)  {
			sprintf(tstr, "thetrecint2_%d", i);
			pulsename(&thetrecintl2[i], tstr);
			createreserve(&thetrecintl2[i], THETA, gres);
			ts[0] = 0;
			for (j = 1; j < gres; j++)  {
				jj = gres -j - 1;
				ts[j] = (short) (ts[j-1] - x*((cphi*Gx[jj] - sphi*Gy[jj])*rdx +
							(cphi*Gy[jj] + sphi*Gx[jj])*rdy)) & ~WEOS_BIT;
			}
			ts[gres - 1] |= WEOS_BIT;
			movewaveimm(ts, &thetrecintl2[i], (int) 0, gres, TOHARDWARE);
			thetrecintlp2[i] = thetrecintl2[i].wave_addr;
		}
	}
	FreeNode(ts);
	setwave(thetrecintlp[0], &thetrec, 0);      /* initial load */
	if(gtype==2)  setwave(thetrecintlp2[0], &thetrec2, 0);

	/*------------------------------------------------------*/
	/* z phase-encode rewinder JFN 09.29.10 
	   pulsename(&gzphase2a, "gzphase2a");
	   pulsename(&gzphase2, "gzphase2");
	   pulsename(&gzphase2d, "gzphase2d");
	   trapezoid(ZGRAD, "gzphase2", 
	   &gzphase2, &gzphase2a, &gzphase2d,
	   pw_gzphase, pw_gzphasea, pw_gzphased, 
	   -ia_gzphase,
	   0, 0, 0, 0,
	   RUP_GRD(pend(&gx,"gx",0)),
	   TRAP_ALL, &loggrd);
	   ------------------------------------------------------*/
	fprintf(stderr, "\npulsegen: generating z phase-encode rewinder... ");
	TRAPEZOID(ZGRAD,
			gzphase2,
			RUP_GRD(pend(&gx, "gx",0) + pw_gzphase2a),
			0,
			0,loggrd);
	fprintf(stderr, "...done ");

	/*------------------------------------------------------*/

	/* Make room for the z rewinder:
	   TRAPEZOID(ZGRAD, gzspoil, RUP_GRD(waitloc2+pw_gzspoila),0,TYPNDEF, loggrd);
	 */
	fprintf(stderr, "\npulsegen: generating Spoilers... ");
	TRAPEZOID(ZGRAD, 
			gzspoil, 
			RUP_GRD(pendall(&gzphase2, 0)+pw_gzspoila),
			0,
			TYPNDEF, loggrd);

	if(nl>1)  {  /* used to be only if interleaving (nl>1)*/
		TRAPEZOID(XGRAD, 
				gxspoil, 
				RUP_GRD(pend(&gzspoild,"gzspoild",0)+pw_gxspoila),
				0,
				TYPNDEF, loggrd);
		TRAPEZOID(YGRAD, 
				gyspoil, 
				RUP_GRD(pend(&gzspoild,"gzspoild",0)+pw_gxspoila),
				0,
				TYPNDEF, loggrd);
	}

	/* fprintf(stderr, "\npw_gzphase=%d , pw_gzphasea = %d", pw_gzphase1, pw_gzphase1a);  */
	fprintf(stderr, "\npw_gx =%d , pw_gy = %d", pw_gx, pw_gy);
	fprintf(stderr, "\npsdseqtime, seqtr, tmin, optr=  %d   %d  %d  %d\n", psdseqtime,  seqtr,tmin,optr);

	SEQLENGTH(seqcore,psdseqtime-timessi+textra ,seqcore);
	getperiod(&deadtime_core, &seqcore, 0);

	/* dummy core to calibrate RF pulses */
	fprintf(stderr, "\npulsegen: generating RFDUMMY pulse with SLICESELZ ... ");
	SLICESELZ(rfdummy, 
		RUP_GRD( psd_rf_wait + tlead),
		pwrf1, 
		5.0, 
		180,
		cycrf1,
		, loggrd);
	fprintf(stderr, " ... done.");
	SEQLENGTH(dummycore, 50ms, dummrycore);


	/*  pass packet sequence (pass).  */
	PASSPACK(endpass, 49ms);
	SEQLENGTH(pass, 50ms, pass);

	/* LHG 1/14/13:  add an extra inversion pulse to be played before the labeling train: */
	fprintf(stderr, "\npulsegen: generating BS3  ... ");
	EXTWAVE(RHO,
			BS3rf,
			RUP_RF(psd_rf_wait),
			pw_BS3rf, a_BS3rf, res_BS3rf,
			/*/usr/g/bin/myhsec1t.rho,,loggrd);*/
		sech_7360.rho,,loggrd);

	EXTWAVE(THETA,
			BS3rf_theta,
			RUP_RF(psd_rf_wait),
			pw_BS3rf, 1.0, res_BS3rf,
			/* /usr/g/bin/myhsec1t.theta,,loggrd);*/
		sech_7360.theta,,loggrd);

	fprintf(stderr, "\n ... done , pw_BS3rf = %d", pw_BS3rf);
	t_preBS = pw_BS3rf + 600;

	SEQLENGTH(preBScore, RUP_RF(t_preBS - timessi),preBScore);
	fprintf(stderr, " \n...  ");
	getperiod(&deadtime_preBScore, &preBScore, 0);
	fprintf(stderr, " \n... finisehd the BS core that goes in front of the labeling train ");

	/* LHG 9.27.16 Here is the VSAI tagging pulse */
	/* VSAI labeling core */
        /* these should be in DAC units */
	INTWAVE(RHO,
			vsitag1,
			RUP_GRD(psd_rf_wait),
			vsi_RFmax / dummy_RFmax,
			vsi_train_len,
			vsi_train_len*4,
			vsi_pulse_mag,  
			1,
			loggrd);

	
        INTWAVE(THETA,
                        vsitag1_theta,
                        RUP_GRD(psd_rf_wait) ,
                        1.0,
                        vsi_train_len,
			vsi_train_len*4,
                        vsi_pulse_phs, 
                        1,
                        loggrd);



        fprintf(stderr, "\n ... vsitag1 done ." );


        INTWAVE(ZGRAD,
                        gztag1,
                        0 ,
                        vsi_Gmax,
                        vsi_train_len,
			vsi_train_len*4,
                        vsi_pulse_grad,
                        1,
                        loggrd);

        fprintf(stderr, "\n ... gztag1 done ." );

        fprintf(stderr, "\nastseqtime=  %d , textra_astcore= %d  "
                , astseqtime, textra_astcore);

        fprintf(stderr, "\n ... astcore (label) done ." );

	SEQLENGTH(astcore, astseqtime + textra_astcore  ,astcore);
	getperiod(&deadtime_astcore, &astcore, 0);
	

	/* VSAI control pulse core */

        INTWAVE(RHO,
                        vsictl1,
                        RUP_GRD(psd_rf_wait)  ,
			vsi_RFmax / dummy_RFmax,
                        vsi_train_len,
			vsi_train_len*4,
                        vsi_pulse_mag,
			1,
                        loggrd);


        INTWAVE(THETA,
                        vsictl1_theta,
                        RUP_GRD(psd_rf_wait),
                        1.0,
                        vsi_train_len,
			vsi_train_len*4,
                        vsi_pulse_ctl_phs, 
                        1,
                        loggrd);


	fprintf(stderr, "\n ... vsictl1 done ." );

        INTWAVE(ZGRAD,
                        gzctl1,
                        0 ,
                        vsi_Gcontrol,
                        vsi_train_len,
			vsi_train_len*4,
                        vsi_pulse_ctl_grad, /* these should be in DAC units */
                        1,
                        loggrd);


        fprintf(stderr, "\n ... gzctl1 done ." );
        fprintf(stderr, "\n ... controlcore (control pulses) done ." );

	SEQLENGTH(controlcore, astseqtime + textra_astcore  ,controlcore);
	getperiod(&deadtime_controlcore, &controlcore, 0);



	/* tagging delay between tagging pulse and image acquisition */
	/* LHG- 12/14/12 -  add the background suppression pulses (sech) */

	fatsattime = pw_rf0 + 2*pw_gz0d +pw_gz0 + 2*timessi;
	fuzz = 10000;
	/* 
	   fprintf(stderr, "\npulsegen: generating PID wait with BS pulses ... ");
	   WAIT(SSP, 
	   astdelay1, 
	   RUP_GRD(24us + fuzz), 
	   RUP_GRD(t_delay - BS1_time - BS2_time - pw_BS1rf/2 - 24us - 2*fuzz ));
	 */

	/* LHG 9.28.16 */
        fprintf(stderr, "\n ... nothingcore (Do nothing) done ." );
	SEQLENGTH(nothingcore, astseqtime + textra_astcore  ,nothingcore);
	getperiod(&deadtime_nothingcore, &nothingcore, 0);

	/* LHG 6.9.17 */
        fprintf(stderr, "\n ... vsi_gap core (Do nothing between VSI pulses) done ." );
	SEQLENGTH(vsi_gapcore, vsi_timegap + textra_vsi_gapcore  ,vsi_gapcore);
	getperiod(&deadtime_vsi_gapcore, &vsi_gapcore, 0);




	fprintf(stderr, "\npulsegen: generating BS0 (sat) pulse  ... ");
	SINC2(RHO, 
			BS0rf,  
			RUP_GRD( psd_rf_wait ), 
			pw_BS0rf, 
			a_BS0rf,
			,0.5,,,loggrd);
	fprintf(stderr," start: %d	end: %d", 
			RUP_GRD( psd_rf_wait ), 
			RUP_GRD( psd_rf_wait + pw_BS0rf )); 

	fprintf(stderr, "\npulsegen: generating BS crushers after the 90   ... ");
	TRAPEZOID(ZGRAD, 
			gzBS0, 
			RUP_GRD(pw_BS0rf + pw_gzBS0a), 
			0, 
			0, loggrd);
	fprintf(stderr," start: %d	end: %d", 
			RUP_GRD(pw_BS0rf + pw_gzBS0a), 
			RUP_GRD(pw_BS0rf + pw_gzBS0a + pw_gzBS0a)); 

	fprintf(stderr, "\npulsegen: generating BS1  ... ");
	EXTWAVE(RHO,
			BS1rf,
			RUP_RF(BS1_time - pw_BS1rf/2 + psd_rf_wait),
			pw_BS1rf, a_BS1rf, res_BS1rf,
			/*/usr/g/bin/myhsec1t.rho,,loggrd);*/
		sech_7360.rho,,loggrd);

	EXTWAVE(THETA,
			BS1rf_theta,
			RUP_RF(BS1_time - pw_BS1rf/2 + psd_rf_wait),
			pw_BS1rf, 1.0, res_BS1rf,
			/* /usr/g/bin/myhsec1t.theta,,loggrd);*/
		sech_7360.theta,,loggrd);

	fprintf(stderr," start: %d	end: %d", 
			RUP_RF(BS1_time - pw_BS1rf/2 + psd_rf_wait),
			RUP_RF(BS1_time - pw_BS1rf/2 + pw_BS1rf  + psd_rf_wait ));
	/*
	   fprintf(stderr, "\npulsegen: generating wait between BS pulses  ... ");
	   WAIT(SSP, 
	   astdelay2, 
	   RUP_GRD(t_delay - BS1_time - BS2_time +  pw_BS1rf/2 + 2*fuzz ), 
	   RUP_GRD(BS1_time -  pw_BS1rf - fuzz ));
	 */
	fprintf(stderr, "\npulsegen: generating BS2  ... ");
	EXTWAVE(RHO,
			BS2rf,
			RUP_RF(BS1_time +  BS2_time - pw_BS2rf/2 + psd_rf_wait),
			pw_BS2rf, a_BS2rf, res_BS2rf,
			/* /usr/g/bin/myhsec1t.rho,,loggrd); */
		sech_7360.rho,,loggrd);


	EXTWAVE(THETA,
			BS2rf_theta,
			RUP_RF(BS1_time +  BS2_time - pw_BS2rf/2 + psd_rf_wait),
			pw_BS2rf, 1.0, res_BS2rf,
			/* /usr/g/bin/myhsec1t.theta,,loggrd);  */
		sech_7360.theta,,loggrd);

	fprintf(stderr," start: %d	end: %d", 
			RUP_RF(BS1_time +  BS2_time - pw_BS2rf/2 + psd_rf_wait),
			RUP_RF(BS1_time +  BS2_time - pw_BS2rf/2 + psd_rf_wait + pw_BS2rf));
	/*
	   fprintf(stderr, "\npulsegen: generating wait between BS2 and fatsat pulses  ... ");
	   WAIT(SSP, 
	   astdelay3, 
	   RUP_GRD(t_delay - BS2_time + pw_BS2rf/2 + fuzz), 
	   RUP_GRD(BS2_time - pw_BS2rf/2 - fatsattime - 2*fuzz));
	 */

	fprintf(stderr, "\npulsegen: generating ASrf_* Arterial Suppression BIR pulses  ... ");
	/* Arterial saturation pulses (preloaded BIR-8 pulses from a file) */
	INTWAVE(RHO,
			ASrf_mag,
			RUP_GRD(t_delay - fatsattime - BIR_len*4 - AStime + psd_rf_wait ), 
			(float)doArtSup ,  
			BIR_len,
			BIR_len*4,
			BIR_mag ,  
			1,
			loggrd);

	
        INTWAVE(THETA,
                        ASrf_theta,
			RUP_GRD(t_delay - fatsattime - BIR_len*4 - AStime +psd_rf_wait ), 
			(float)doArtSup,   
			BIR_len,
			BIR_len*4,
			BIR_phs ,  
                        1,
                        loggrd);

	
        INTWAVE(ZGRAD,
                        ASrf_grad,
			RUP_GRD(t_delay - fatsattime - BIR_len*4 - AStime  ), 
			(float)doArtSup * BIR_Gmax ,  
			BIR_len,
			BIR_len*4,
			BIR_grad ,  
                        1,
                        loggrd);



        fprintf(stderr, "\n ... ASrf done ." );



	/* LHG - 12/14/12 -  fat sat pulse from the scan core 
	 * regular interleaved sequence (core).
	 */

	fprintf(stderr, "\npulsegen: generating Fat Sat pulse  ... ");
	SINC2(RHO, 
			rf0,  
			RUP_GRD(t_delay - fatsattime + psd_rf_wait ), 
			pw_rf0, 
			a_rf0,
			,0.5,,,loggrd);
	fprintf(stderr," start: %d	end: %d", 
			RUP_GRD(t_delay - fatsattime + psd_rf_wait ), 
			RUP_GRD(t_delay - fatsattime + psd_rf_wait + pw_rf0 )); 

	fprintf(stderr, "\npulsegen: generating Fat Sat gradients  ... ");
	TRAPEZOID(ZGRAD, 
			gz0, 
			RUP_GRD(t_delay - fatsattime + pwrf0 + pw_gz0a), 
			0, 
			0, loggrd);
	fprintf(stderr," start: %d	end: %d", 
			RUP_GRD(t_delay - fatsattime + pwrf0 + pw_gz0a), 
			RUP_GRD(t_delay - fatsattime + pwrf0 + 2*pw_gz0a + pw_gz0));

	fprintf(stderr, "\n ...  done ." );


	SEQLENGTH(tdelaycore,t_delay - timessi, tdelaycore);
	getperiod(&deadtime_tdelaycore, &tdelaycore, 0);
	fprintf(stderr, "\n ...  done ." );


	/* delay after image acquisition to fill out the TR*/
	t_adjust =  RUP_RF(optr - t_preBS - t_tag - t_delay - psdseqtime *opslquant );
	fprintf(stderr, "\npulsegen: Creating TRdelay adjust core ... ");
	fprintf(stderr,"\n optr: %d, t_preBS: %d, t_tag: %d,  t_delay: %d,  psdseqtime: %d, opslquant: %d, t_adjust:  %d", 
			optr, t_preBS, t_tag, t_delay, psdseqtime,opslquant,t_adjust );

	/*  make a little trig on DABOUT6 J12  */
	if(maketrig)  {
		trigonwd = SSPD + DABOUT6;
		trigoffwd = SSPD;
		trigonpkt[0] = SSPDS + EDC;
		trigoffpkt[0] = SSPDS + EDC;
		trigonpkt[2] = trigonwd;
		trigoffpkt[2] = trigoffwd;
		SSPPACKET(trigon, trigloc+waitloc, 3us, trigonpkt, 0);
		SSPPACKET(trigoff, trigloc+waitloc+triglen, 3us, trigoffpkt, 0);
	}


	WAIT(SSP, 
		TRdelay,
		RUP_GRD(24us),  
		/*RUP_RF(optr - t_preBS - t_tag - t_delay - psdseqtime *opslquant ) );*/
		5ms);
	/* Note that there is a 3ms fudge factor.  Come back and find where the error is!! */
	fprintf(stderr, "\ntimessi: %d ", timessi );

	SEQLENGTH(tadjustcore, RUP_RF(t_adjust - timessi ) ,tadjustcore);
	fprintf(stderr, "\n ...  " );
	getperiod(&deadtime_tadjustcore, &tadjustcore, 0);
	fprintf(stderr, "\n ...  done ." );

	/*-------------------------------------------------*/

	/*  wait_for_scan_to_start sequence */
	WAIT(SSP, waitStart, 24us, 1ms);
	SEQLENGTH(waitpass, 10ms, waitpass);

	/*  wait for scanend for real time */
	WAIT(SSP, waitEnd, 24us, 2us);
	SEQLENGTH(waitend, endtime, waitend);

@inline Prescan.e PSpulsegen	/*  prescan sequences  */

		buildinstr();              	/* load the sequencer memory */
	return SUCCESS;

} /* end of pulsegen */

@inline Prescan.e PSipg
/* end of @pg */

@rsp

STATUS scancore(void);
int doleaf(int leafn, int framen, int slicen, int* trig, int* bangn, int dabop, int dtype);
/*---------------------------------------------------------
  LHG 6.29.12 : My function declarations
  ----------------------------------------------------*/
void doast(int* trig, int isLabel);
void dodelay(int* trig); 
void doadjust(int* trig);
/*---------------------------------------------------------*/



@inline Prescan.e PScore

/*manual prescan */
int mps2() {
	pre = 2;
	scancore();
	rspexit();
	return SUCCESS;
}
/*auto prescan */
int aps2() {
	pre = 1;
	scancore();
	rspexit();
	return SUCCESS;
}
/*Actual scan*/
int scan()
{
	pre = 0;
	scancore();
	rspexit();
	return SUCCESS;
}

void get_spgr_phase(int ntab, int *rfphase);


STATUS scancore()
{
	int counter;
#define RESERVED_IPG_MEMORY (0xbfff00)
	int *comm_buffer;		/* to talk to real time recon */
	FILE *fpin;                 /* for StartScan file */
	float phi;
	int i, j, n;

	int 	*vsi_iphase;
	int 	*vsi_iphase_control;
	double 	*rftag_weights;
	double 	vsi_rfphase=0;
	int	multphs_ctr=0;

	printf("\nEntering scancore   pre = %d\n", pre);

	/* ---------------------------------------------------------------*/
	/* LHG 12/24/08:  initialize the interactive Real Time vars: */
	/* ---------------------------------------------------------------*/
	int myrfamp ; /* LHG 6.29.12 used for the variable amplitude 3D excitation pulses */
/*
	setperiod(astseqtime, &astcore, 0);
	setperiod(astseqtime, &controlcore, 0);
	setperiod(astseqtime, &nothingcore, 0);
*/
	printf("\ncalculated  a_rf1= %f , opflip= %f: , ia_rf1=  %d", a_rf1, opflip, ia_rf1);
        
        printf("\n TARDIS_FREQ_RES: %d", TARDIS_FREQ_RES);
        printf("\n AST VARIABLES: ");
        printf("\n Gradient max amplitude (a_gztag1):   %f", a_gztag1);
        printf("\n Gradient Duration  (pw_gztag1):  %d", pw_gztag1);
        printf("\n-------------------\n ");
        fflush(0);
        /* ---------------------------------------------------------------*/
        

	/* set the pointer comm_buffer to the 64 ints at the top of ipg memory */
	comm_buffer = (int *) RESERVED_IPG_MEMORY;
	/* set the first entry to the address, as a handshake to recon client */
	comm_buffer[0] = (int) (&(comm_buffer[0]));

	if(maketrig==1)  {
		setwamp(trigoffwd, &trigon, 2);  /* no trigs yet */
		counter = 0;
	}
	setrfconfig(ENBL_RHO1 + ENBL_THETA);
	setssitime((LONG)timessi/GRAD_UPDATE_TIME);
	rspqueueinit(queue_size);
	scopeon(&seqcore);
	syncon(&seqcore);
	syncoff(&pass);
	setrotatearray((SHORT)opslquant, rsprot[0]);
	settriggerarray((SHORT)opslquant, rsptrigger);

	dabmask = PSD_LOAD_DAB_ALL;
	if(pre)  {
		setrfltrs(filter_aps2, &echo1);
		if(rhnecho==2) setrfltrs(filter_aps2, &echo2);
	}
	else  {
		setrfltrs(filter_echo1, &echo1);
		if(rhnecho==2) setrfltrs(filter_echo1, &echo2);
	}

	/* fix the clock for aux trig */

	if (!pre && gating==TRIG_AUX)  {
		setscantimemanual();
		setscantimestop();
		setscantimeimm(pidmode, pitscan, piviews,
				pitslice, opslicecnt);
	}

	/* Allocate memory for RF pulse param tables */ 
	rfphastab = 	(int *) AllocNode(nbang*sizeof(int));
	fxmit = 	(int *) AllocNode(opslquant*sizeof(int));
	frec  = 	(int *) AllocNode(opslquant*sizeof(int));
	slcphastab  = 	(int *) AllocNode(opslquant*sizeof(int));
	slordtab  = 	(int *) AllocNode((opslquant+1)*sizeof(int));
 	rfamptab = 	(int *) AllocNode(opslquant*sizeof(int));	

	fprintf(stderr,"\nSetting up slice sel xmit Phase table (for RF spoiling) spgr_flag= %d ...", spgr_flag);
	get_spgr_phase(nbang, rfphastab);

	fprintf(stderr,"\nSetting up the slice sel xmit frequency table ...");
	setupslices(fxmit, rsp_info, opslquant, a_gzrf1, 1.0, opfov, TYPTRANSMIT);
	/* for (i=0; i<opslquant; i++) fprintf(stderr,"\nXmit freq:  fxmit[%d] = %d ", i, fxmit[i]);	*/


	/*LHG 9/26/12 - setting up the amplitudes of the slab select pulses */
	fprintf(stderr,"\nSetting up the slice sel xmit Amplitude table ...");
	get_rfamp(opslquant, rfamptab);

	for (isl = 0; isl < opslquant; isl++)
		rsp_info[isl].rsprloc = 0;
	setupslices(frec, rsp_info, opslquant, 0.0, 1.0, 2.0, TYPREC);
	/*for (i=0; i<opslquant; i++) fprintf(stderr,"\nreceive freq:  frec[%d] = %d ", i, frec[i]);	  */

	/* make table for loaddab slice num */
	fprintf(stderr,"\nSetting up the loaddab table ...");

	for (i=0; i<opslquant; i++)
		slordtab[i] = i;

	/*skip this part:  in a 3D sequence the slice order is not determined here*/
	if(0){
		switch((int)opuser16)  {
			case 0:			/* sequential */
				for (i=0; i<opslquant; i++)
					slordtab[i] = i;
				break;
			case 1:			/* interleaved */
				j = 0;
				for (i=0; i<opslquant; i+=2)  {
					slordtab[j++] = i;
					slordtab[j+(opslquant-1)/2] = i+1;
				}
				break;
		}
	}
	/*------------------------------------------------------------*/
	/* phase shift due to sliceselect freq shift  */
	for (isl = 0; isl < opslquant; isl++)  {
		/* LHG 10/7/2012   ... I don't know why there was a factor of 4 there in the first place!
		   phi = -4.0*cyc_rf1*rsp_info[isl].rsptloc*(1.0 + psd_rf_wait/(float)pw_rf1)/opslthick; 
		 */
		phi = -1.0*cyc_rf1*rsp_info[isl].rsptloc*(1.0 + psd_rf_wait/(float)pw_rf1)/opslthick;
		slcphastab[isl] = FS_PI*phi;
		/*fprintf(stderr,"\nisl slc loc phaseshifts = %d %f %f\n", isl, rsp_info[isl].rsptloc, phi); */
	}

	/* Frequency for Fat Sat pulse */
	setfrequency(satoff/TARDIS_FREQ_RES, &rf0, 0);

	if(rhuser14>0 && !pre)  {	/* first shot is delayed for map */
		fprintf(stderr,"\nSetting up the map delays ...");
		if(gtype<2)  {
			setperiod(mapdel+4us, &mapx, 0);
			setperiod(mapdel+4us, &mapy, 0);
			setperiod(mapdel+4us, &mapz, 0);
			setperiod(mapdel+4us, &mapt, 0);
			setperiod(mapdel+4us, &maps1, 0);  /*LHG 10/16/12  */
		}  else  {
			setperiod(mapdel+4us, &mapx2, 0);
			setperiod(mapdel+4us, &mapy2, 0);
			setperiod(mapdel+4us, &mapz2, 0);
			setperiod(mapdel+4us, &mapt2, 0);
		}
		if(bmapnav==0)  {         /* a/d window follows readout grads */
			if(gtype==2)
				setperiod(mapdel+4us, &maps2, 0);
			else
				setperiod(mapdel+4us, &maps1, 0);
			setperiod(4us, &maps3, 0);
		}
		fprintf(stderr,"Done");
	}

	bangn = 0;		/* init spoiler counter */
	dabop = DABSTORE;


	/* -------------------------------------------------------*/
	/*LHG 6.29.12 : setting the currents for the flow crushers*/
	if (mycrush) {
	/*
		setiampt( max_pg_iamp*(a_gxsp3/loggrd.xfs), &gxsp3, 0);
		setiampt( max_pg_iamp*(a_gxsp4/loggrd.xfs), &gxsp4, 0);

		setiampt( max_pg_iamp*(a_gysp3/loggrd.yfs), &gysp3, 0);
		setiampt( max_pg_iamp*(a_gysp4/loggrd.yfs), &gysp4, 0);

	*/
		setiampt( max_pg_iamp*(a_gzsp3/loggrd.zfs), &gzsp3, 0);
		setiampt( max_pg_iamp*(a_gzsp4/loggrd.zfs), &gzsp4, 0);
	}

	/* --------------------------------------------------------*/
	/* prescan loop. */
	/* --------------------------------------------------------*/
	ifr = 0;

	fprintf(stderr,"\nDone with setup stuff, getting into loops ...");

	if(pre)  {
		fprintf(stderr,"\n...PRE-SCAN loop ...");
		if(gtype==0 || gtype==2)
			setperiod(4us, &maps1, 0);
		else
			setperiod(pw_gx-psfrsize*tsp, &maps1, 0);  /* offset window to k=0 */
		for (iv = 0; iv < ndisdaq; iv++) {	/*  disdaqs  */
			trig = (gating==TRIG_AUX)? TRIG_INTERN : gating;

			/* rt_updates(); */
			doadjust(&trig);
			doast(&trig, isLabel);
			dodelay(&trig);

			/* center-out z phase-encode order */
			for (isl = 0; isl < opslquant/2; isl++) {
				doleaf(0,0, opslquant/2 -1-isl, &trig, &bangn, dabop, DABOFF);
				doleaf(0,0, opslquant/2 +isl,   &trig, &bangn, dabop, DABOFF);
			}
			getiamp(&chopamp, &rf1, 0);
			setiamp(-chopamp, &rf1, 0);
		}
		for (iv = 0; iv < 10000; iv++) {	/* get data  */
			trig = (gating==TRIG_AUX)? TRIG_INTERN : gating;

			/* rt_updates();*/
			doadjust(&trig);
			doast(&trig, isLabel);
			dodelay(&trig);

			if (dopresat) {
				printf("\nDoing a dummy bang");
				doleaf(0, 0, opslquant/2-1, &trig, &bangn, dabop, DABOFF);
				bangn-=1;
			}

			/* center-out z phase-encode order */
			for (isl = 0; isl < opslquant/2; isl++) {
				dtype = (pre==1 || isl==0)? DABON:DABOFF;
				doleaf(0,0, opslquant/2-isl-1, &trig, &bangn, dabop, dtype);
				doleaf(0,0, opslquant/2+isl,   &trig, &bangn, dabop, dtype);
			}
			getiamp(&chopamp, &rf1, 0);
			setiamp(-chopamp, &rf1, 0);
		}
		rspexit();
	}
	/* ---------------------------------------------------------*/
	/* disdaq loop. */
	/* ---------------------------------------------------------*/
	fprintf(stderr,"\n....DISDAQ loop ...");

	if(opuser2==3)  {           /* wait for StartScan to exist */
		boffset(off_waitpass);
		settrigger(TRIG_INTERN, 0);
		while((fpin=fopen("/usr/g/mrraw/StartScan", "r"))==NULL)
			startseq(0, MAY_PAUSE);     /* spin in a 10 ms wait loop  */
	}

	trig = gating;
	counter = 0;
	for (iv = 0; iv < nextra; iv++)  {

		/* if(maketrig==1 && (counter++)%trigfreq==0)*/
		if(maketrig==1 )
			setwamp(trigonwd, &trigon, 2);
		else if(maketrig==1)
			setwamp(trigoffwd, &trigon, 2);

		/* rt_updates();*/
		doadjust(&trig);
		doast(&trig, isLabel);
		dodelay(&trig);
		if (dopresat) {
			printf("\nDoing a dummy bang");
			doleaf(0, 0, opslquant/2-1, &trig, &bangn, dabop, DABOFF);
			bangn-=1;
		}
		/* center-out z phase-encode order */
		for (isl = 0; isl < opslquant/2; isl++) {
			doleaf(0,0, opslquant/2-isl-1, &trig, &bangn, dabop, DABOFF);
			doleaf(0,0, opslquant/2+isl,   &trig, &bangn, dabop, DABOFF);
		}
		/* ---------------------------------------------------------*/
		if(opuser2 == 2) trig = TRIG_AUX;  	/* wait for trig */

		if (iv==0)
			setwamp(trigoffwd, &trigon, 2); /* end TTL trigger out code by SJP ?*/

	}

	/* ---------------------------------------------------------*/
	/* scan loop */
	/* ---------------------------------------------------------*/
	fprintf(stderr,"\n...SCAN loop ...");
	counter = 0;
	isOdd = 0;

	/* LHG : 9.27.16 : we don't want the field map when we're adjusting the phase of pcasl pulses*/ 
	if (multiphsFlag) 
	{
		fprintf(stderr,"\nMulti-phase calibration  Mode = %d", domap);
		domap=0;
	}
	fprintf(stderr,"\ndomap = %d", domap);
	fprintf(stderr,"\nM0frames = %d", M0frames);


	for (ifr = 0; ifr < nframes; ifr++)  {

		/* LHG: 1/29/16: implementing GRAPPA along the z direction */
		isOdd = !isOdd ;  /* the first frame (ie, when ifr == 0) is not odd */

		if(rhuser14>0 && ifr == 1)  {	/*  reset  */
			if(gtype<2)  {
				setperiod(4us, &mapx, 0);
				setperiod(4us, &mapy, 0);
				setperiod(4us, &mapz, 0);
				setperiod(4us, &mapt, 0);
				setperiod(4us, &maps1, 0);
			}  else  {
				setperiod(4us, &mapx2, 0);
				setperiod(4us, &mapy2, 0);
				setperiod(4us, &mapz2, 0);
				setperiod(4us, &mapt2, 0);
				setperiod(4us, &maps2, 0);
			}
			if(bmapnav==0)  {
				setperiod(mapdel+4us, &maps3, 0);
			}
		}
		/* alternate between control and tag after collecting field map (if field map is used)*/
		if( ifr>1 || domap==0)
			isLabel = !(isLabel);

		for (iv = 0; iv < nl; iv++)  {
			/* rt_updates(); */
			/* if(maketrig==1 && (counter++)%trigfreq==0)*/
			if(maketrig==1 )
				setwamp(trigonwd, &trigon, 2);
			else if(maketrig==1)
				setwamp(trigoffwd, &trigon, 2);
		
			doadjust(&trig);
			doast(&trig, isLabel);
			dodelay(&trig);

			if (dopresat) {
				printf("\nDoing a dummy bang");
				doleaf(0, 0, opslquant/2-1, &trig, &bangn, dabop, DABOFF);
				bangn-=1;
			}

			/* center-out z phase-encode order */
			kzcount=0;

			if (ifr < M0frames ) 
				fprintf(stderr,"collecting M0 images and the GRAPPA-Z cal data at the odd frames: ifr=%d isOdd=%d", ifr, isOdd); 

			for (isl = 0; isl < opslquant/2; isl++) {
				/* eg - if opslquant is 18 : 		*/
				/* "slices" 9 to 17 */
				doleaf(iv, ifr, opslquant/2+isl,   &trig, &bangn, dabop, DABON);
				/* "slices" 8 to 0 			*/
				doleaf(iv, ifr, opslquant/2-isl-1, &trig, &bangn, dabop, DABON);
			}

			if(opuser2 == 2) trig = TRIG_AUX;  	/* wait for trig */

		}
		/*LHG 10/3/12 phase correction calibration 
		if (multiphsFlag && ifr>M0frames-1) {
			multphs_ctr++;
			if (multphs_ctr == nfr_phs){
				myphase_increment += 0.4;
				fprintf(stderr,"\nincrementing the  phase correction: %f", myphase_increment);
				get_pcasl_phases(pcasl_iphase, myphase_increment, nreps);
				multphs_ctr=0;
			}
		}
		*/

		comm_buffer[1] = ifr;           /* talk to RT grecons */
	}      /* end frames */

	/*  LHG 9/8/14 - disable this: 	
	    if(opuser11 == 1.0)  {      // wait around for grecon to finish /
	    boffset(off_waitend);
	    settrigger(TRIG_INTERN, 0);
	    startseq(0, MAY_PAUSE);   // fat lady sings 
	    }
	 */



	/* tell 'em it's over */
	boffset(off_pass);
	setwamp(SSPD+DABPASS+DABSCAN, &endpass, 2);
	settrigger(TRIG_INTERN, 0);
	startseq(0, MAY_PAUSE);	/* fat lady sings */

	return SUCCESS;
}

/*-----------------------------------------------------
  LHG 5/11/01

  function doast()

  This function sets up the Arterial Spin Labeling frequencies, amplitudes, ...
  -----------------------------------------------------------------------*/
void doast(int* trig, int isLabel)
{

	int slicen, td, n, m;
	int rotmatx[1][9];
	float phi, cphi, sphii, DABphase;
	slicen=0;

	fprintf(stderr,"\ndoast ...");

	/* This is where we do the first background suppression pulse: */	
	fprintf(stderr,"\ndoing the first BS pulse before labeling ...");
	boffset(off_preBScore);
	startseq(0, MAY_PAUSE);

	/*-----------------------------------------------------
	  LHG 5/24/01
	  always use same rotation for the tag
	  ( un-do the slice rotation )
	  --------------------------------------------------------*/
	/* This Assures that gztag will always occur on the z axis only
	   even for an oblique Rx */

	rotmatx[0][0] = 0;
	rotmatx[0][1] = 32767;
	rotmatx[0][2] = 0;
	rotmatx[0][3] = -32767;
	rotmatx[0][4] = 0;
	rotmatx[0][5] = 0;
	rotmatx[0][6] = 0;
	rotmatx[0][7] = 0;
	rotmatx[0][8] = 32767;

	scalerotmats(rotmatx, &loggrd, &phygrd, 1, 0);
	/* setrotate(rotmatx,slicen); */
	setrotate(rotmatx,0);
	/*-----------------------------------------------------------*/

	/* Set up gradient amplitudes  for VSAI train */
	setiamp((int)(max_pg_iamp*(vsi_Gmax/loggrd.zfs)), &gztag1, 0);
	setiamp(ia_vsitag1, &vsitag1,0);

	switch (vsi_controlmode){
		case 1:
			/* case 1: alternate the sign of the flip angle */
			setiamp( -ia_vsitag1 , &vsictl1, 0);
			break;

		case 2:
			/* case 2: trurn off the RF pulses */
			setiamp( 0 , &vsictl1, 0);
			break;
		case 3:
			/* case 3: Phase of control pulses targets the negative velocities */
			break;
		case 4:
			/* case 4: gradients are reversed (also target negative velocities */
			setiamp( (int)(max_pg_iamp*(-vsi_Gmax/loggrd.zfs)), &gzctl1, 0);
			break;	
		case 5:
			/* case 5: control pulses are NOT velselective. everything gets inverted */
			setiamp( 0 , &gzctl1, 0);
			break;
		case 6:
			/* case 6: the control gradient is the absolute value of the gradient (as in paper by Qin)*/
			for (i=0; i<Npoints ; i++){
				vsi_pulse_ctl_grad[i] = abs(vsi_pulse_grad[i]);
			}
			break;
	}


	/* LHG: 11.7.14:   we will do N_cycles cycles of pulses */
	for(m=0; m<vsi_Ncycles; m++)
	{
		
		fprintf(stderr," calling vsi_gapcore  ...");
		boffset(off_vsi_gapcore);
		startseq(0, MAY_PAUSE);
		
		/* LHG 7/11/13:  while collecting the M0 images, we want the label icompletely off */
		if (ifr < M0frames)
		{
			fprintf(stderr, " calling nothingcore  ...");
			boffset(off_nothingcore);
			startseq(0, MAY_PAUSE);
			fprintf(stderr,"done");
		}
		else
		{
			if (isLabel)
			{
				fprintf(stderr, " astcore ...");
				boffset(off_astcore);
				startseq(0, MAY_PAUSE);
				fprintf(stderr,"done");

			}
			else
			{

				fprintf(stderr, "Control core ...");
				boffset(off_controlcore);
				startseq(0, MAY_PAUSE);
				fprintf(stderr,"done");

			}
		}
	}
	/* reset offset and go for it */
	if(*trig == TRIG_AUX)
	{                /*  make clock start  */

		setscantimeimm(pidmode,
				astseqtime,
				piviews,
				pitslice,
				opslicecnt);
		setscantimestart();
	}

	*trig = TRIG_INTERN;
}
/*--------------------------------------------------*/


void dodelay(   int* trig)
{
	/*--------------------------------------------------------*/
	/* LH 12/14/08 - this is where the t_delay gets updated on the fly  */
	/*--------------------------------------------------------*/
	fprintf(stderr,"\ndodelay ...");
	/* setperiod(rt_t_delay  ,
	   &astdelay1,0);
	 */
	/*setperiod(rt_t_delay + tpre + tdel + textra_delaycore,  
	  &tdelaycore,0); 
	/*--------------------------------------------------------*/

	/* turn the Back. Supp.  off - LHG 6/11/13 */
	if (ifr<M0frames){
		setiamp(0, &BS1rf, 0);	
		setiamp(0, &BS2rf, 0);	
		setiamp(0, &BS3rf, 0);	
	}else
	{
		/* turning them back on - use DAC units */
		setiamp(ia_BS1rf, &BS1rf, 0);	
		setiamp(ia_BS2rf, &BS2rf, 0);	
		setiamp(ia_BS3rf, &BS3rf, 0);	

	}

	fprintf(stderr," calling tdelaycore");	
	boffset(off_tdelaycore);
	startseq(0, MAY_PAUSE);
	

	/* reset offset and go for it */
	if(*trig == TRIG_AUX)
	{                /*  make clock start  */

		setscantimeimm(pidmode,
				2600, /* I don't think this ever is used here... */
				piviews,
				pitslice,
				opslicecnt);
		setscantimestart();
	}

	*trig = TRIG_INTERN;
}


void doadjust(  int* trig)
{
	fprintf(stderr,"\ndoadjust ...");
	/*--------------------------------------------------------*/ 
	/* LH 12/14/08 - this is where the t_delay gets updated on the fly  */
	/*--------------------------------------------------------*/
	/* setperiod( rt_t_adjust,
	   &TRdelay,0); */
	/*setperiod( t_adjust + tlead + tdel + tpre + textra_tadjust,  
	  &tadjustcore,0);
	/*--------------------------------------------------------*/


	boffset(off_tadjustcore);
	startseq(0, MAY_PAUSE);

	/* reset offset and go for it */
	if(*trig == TRIG_AUX)
	{                /*  make clock start  */

		setscantimeimm(pidmode,
				2600, /* I don't think this ever is used here... */
				piviews,
				pitslice,
				opslicecnt);
		setscantimestart();
	}

	*trig = TRIG_INTERN;
}


/*
 * function doleaf()
 *
 * this function does most of the things connected with a single
 * spiral excitation.  in keeping with standard epic practices,
 * it secretly uses a variety of global variables, including
 * opslquant, savrot, nl, nframes, fxmit, frec, echo1, omrf1, core
 * and gating.
 *
 * arguments
 *     leafn  -- interleaf number.       0 <= leafn  < nl.
 *     framen -- temporal frame number.  0 <= framen < nframes.
 *     slicen -- slice number.           0 <= slicen < opslquant.
 * 		 LHG: in the 3D version, this is now the kz phase encode position.
 *     bangn  -- rf bang number.         0 <= bangn  < nbang.
 *     dabop  -- dabop for loaddab.
 *     dtype  -- type of data acquisition (DABON or DABOFF).
 *
 */

int doleaf(int leafn, int framen, int slicen, int* trig, int* bangn, int dabop, int dtype)
{
	int k, viewn;
	int echon;
	long rotmatx[1][9];
	float phi, cphi, sphi;
	int phase;
	int myxmitfreq, myrecfreq;
	float zpeamp;
	int i_zpeamp;

	/*fprintf(stderr,"\ndoleaf..."); */
	if (1)
	{
		/* Slab selection for 3D AQ, so it's all different from sprlio */

		/* set the frequency of xmitter and receiver.  The offset is right in between the two middle slicesn*/
		myxmitfreq = (int)((fxmit[opslquant/2] + fxmit[opslquant/2-1])/2);
		setfrequency(myxmitfreq, &rf1, 0); 

		myrecfreq = (int)((frec[opslquant/2] + frec[opslquant/2-1])/2);
		setfrequency(myrecfreq , &echo1, 0);  
		/* setfrequency(0 , &echo1, 0); */ 


		/* set xmit  amplitude according to precalcualted table by JFN 2011 
		   already in DAC units*/

		/*setiamp(rfamptab[slicen], &rf1, 0);	LHG fix for ampl. order 7/11/13 */
		/*kzcount = (int)(fmod((float)*bangn,(float)opslquant));*/
		setiamp(rfamptab[kzcount], &rf1, 0);
		/*fprintf(stderr,"\nkzcount ... rf numer: %d , bangn ..%d , isl ...%d  ", kzcount, *bangn,isl );*/
		kzcount++;
		if (kzcount >= opslquant)
			kzcount=0;


		/* set the phase of the xmitter and recvr according to precalculated table */
		setiphase(rfphastab[*bangn], &echo1, 0);
		/* set up receiver phase to account for shift in the z-direction */
		phase = rfphastab[*bangn] + slcphastab[slicen];
		phase = phase % FS_2PI;
		if(phase>FS_PI) phase -= FS_2PI;
		if(phase<-FS_PI) phase += FS_2PI;   /*  +/- pi  */
		setiphase(phase, &rf1, 0);
		setiphase(phase, &rf0, 0);

	}
	else
	{ /* old sprlio code: */

		setfrequency(fxmit[slicen], &rf1, 0);

		phase = rfphastab[*bangn] + slcphastab[slicen];
		phase = phase % FS_2PI;
		if(phase>FS_PI) phase -= FS_2PI;
		if(phase<-FS_PI) phase += FS_2PI;   /*  +/- pi  */
		setiphase(phase, &rf1, 0);

		setfrequency(frec[slicen], &echo1, 0);
		setiphase(rfphastab[*bangn], &echo1, 0);
		if(rhnecho==2)   {
			setfrequency(frec[slicen], &echo2, 0);
			setiphase(rfphastab[*bangn], &echo2, 0);
		}
	} /* end old sprlio code: */

	/* set Gz amplitude for kz encodings */
	/* phase encode in 3D:  the middle slice should be the center of k-space (kz = 0) */
	if (framen < M0frames &&  isOdd && doZgrappa) 
		/* in GRAPPA, we collect data in the kz gaps during the calbration scans */ 
		zpeamp =  ((float)slicen - (float)opslquant/2.0 + 0.5 ) * (float)iampGzstep ; 
	else	
		zpeamp =  ((float)slicen - (float)opslquant/2.0 ) * (float)iampGzstep ; 

	i_zpeamp = (short)zpeamp;
	i_zpeamp = (short)(2*i_zpeamp/2.0);  /* force amplitudes to be even. necessary? */
	setiampt(i_zpeamp, &gzphase1, 0);
	setiampt( -i_zpeamp, &gzphase2, 0);

	/* I think this is for the double echo case, which we don't use */
	if(rhnecho==2)   {
		setfrequency(myrecfreq, &echo2, 0);
		setiphase(rfphastab[*bangn], &echo2, 0);
	}

	/* Incrementing the bang counter    */
	(*bangn)++;
	*bangn = *bangn%nbang;      /*  for mps2  */

	/* set up timing. */
	tmp = deadtime_core;
	if (scluster==1)  {
		if (slicen == opslquant-1)
			tmp = optr - opslquant*psdseqtime;
	}
	else
		if(slicen < nerror)		/* spread the error ... */
			tmp = deadtime_core - trerror/nerror;	/* ... while keeping the sign */
	setperiod(RUP_GRD(tmp), &seqcore, 0);

	settrigger(*trig, slicen);


	/* set up dab */
	viewn = framen*nl + leafn + 1;
	echon = 0;
	if(gtype==2 && pre)  {	/* no first echo acq */
		loaddab(&echo1, slordtab[slicen], echon, dabop, viewn, DABOFF, dabmask);
	}
	else  {
		loaddab(&echo1, slordtab[slicen], echon, dabop, viewn, dtype, dabmask);
	}
	if(rhnecho==2)  {
		loaddab(&echo2, slordtab[slicen], echon+1, dabop, viewn, dtype, dabmask);
	}

	/* LHG 6.29.12:  this code was needed for RT in the Signa platform:
	   if (doRT)
	   loadhsdab(&echo1,
	   (LONG)slicen,
	   (LONG)echon,
	   (LONG)dabop,
	   (LONG)0,
	   (LONG)1,
	   (LONG)1,
	   (LONG)1,
	   (LONG)viewn,
	   (TYPDAB_PACKETS)DABON,
	   (LONG)hsdabmask);

	   else
	   loaddab(&echo1, slicen, echon, dabop, viewn, dtype, dabmask);
	 */


	/* select interleaf */


	phi = 2.0*3.14159265*(leafn)/(float)nl;
	cphi = cos(phi); sphi = sin(phi);
	for (k = 0; k < 9; k += 3)
	{
		rotmatx[0][k] = IRINT(cphi*savrot[slicen][k]+sphi*savrot[slicen][k+1]);
		rotmatx[0][k+1] = IRINT(-sphi*savrot[slicen][k]+cphi*savrot[slicen][k+1]);
		rotmatx[0][k+2] = savrot[slicen][k+2];
	}
	scalerotmats(rotmatx, &loggrd, &phygrd, 1, 0);
	setrotate(rotmatx[0],slicen);

	setwave(thetrecintlp[leafn], &thetrec, 0);  /* fov offset for receiver */
	if(gtype==2) setwave(thetrecintlp2[leafn], &thetrec2, 0);

	/*  spoiler gradients on readout axis don't rotate; they are here to
	    fix an evil hysteresis problem in the grads  */
	if(nl>1){	
		ia_gxspoil = (cphi + sphi)*(a_gxspoil/loggrd.xfs)*max_pg_iamp;
		ia_gyspoil = (cphi - sphi)*(a_gyspoil/loggrd.yfs)*max_pg_iamp;
		setiampt(ia_gxspoil, &gxspoil, 0);
		setiampt(ia_gyspoil, &gyspoil, 0);
	}

	/* reset offset and go for it */
	boffset(off_seqcore);
	startseq(slicen, MAY_PAUSE);
	if(*trig == TRIG_AUX && opuser2 >= 1)  {	/*  make clock start  */
		setscantimestart();
	}
	*trig = TRIG_INTERN;
	/* fprintf(stderr,"..doleaf done.");  */
	return SUCCESS;

}   /* end of function doleaf() */

void get_spgr_phase(ntab, rfphase)
	int ntab;
	int rfphase[];
{

	int IA = 141;         /*  24 bit overflow  */
	int IC = 28411;
	int IM = 134456;
	int i, jran;
	float x;
	float frfphase,ftmp;

	switch (spgr_flag)  {
		case 0:           /* no spoiling  */
			for (i=0;i<ntab;i++)
				rfphase[i] = 0;
			break;
		case 1:           /*  GEMS algorithm  */
			rfphase[0] = 0;
			for (i=1;i<ntab;i++)
				rfphase[i] = ((int)((float)rfphase[i-1] + (float)i*seed + 3L*FS_PI) %
						FS_2PI)-FS_PI;
			break;
		case 2:       /* random number generator (Num. Recipes, p.209  */
			jran = seed;
			for (i=0;i<ntab;i++)  {
				jran = (jran*IA + IC) % IM;
				rfphase[i] = (int)((float)(FS_2PI + 1)*(float)jran/(float)IM) - FS_PI;
			}
			break;
		case 3:           /* LHG add RF spoiling with linear phase increment from JFN 3/11/11 */
			rfphase[0] = 0;
			frfphase = 0;
			for (i=1;i<ntab;i++) {
				frfphase += (float) rf_spoil_seed * M_PI / 180.0 * (float) i;
				frfphase = atan2 (sin(frfphase), cos(frfphase));      /* wrap phase to (-pi,pi) range */
				rfphase[i] = (int)(frfphase/M_PI * (float)FS_PI);
			}
		default:
			break;
	}  /* switch  */

	/* 
	   fprintf(stderr, "\nPhase Table: ");
	   for (i=0;i<ntab;i++)  
	   fprintf(stderr,"\nget_spgr_phase: rfphase[%d] = %d / ", i, rfphase[i]);
	 */
	return;

}	/*  end get_spgr_phase  */


void get_rfamp(ntab, rfamp)
	int ntab;
	int rfamp[];
{
	int i;
	float a0,a,sina,E1;

	E1 = exp(-(float)psdseqtime/(float)1300ms);
	fprintf(stderr,"\nget_rfamp: E1 = %f", E1);

	switch (rfamp_flag)  {
		case 0:           /* constant flip angle */
			for (i=0;i<ntab;i++)
				rfamp[i] = ia_rf1;
			break;
		case 99:          /* Gai schedule */
			a0 = opflip/180.0*M_PI;
			a = a0;
			rfamp[0] = (int) (a0/(M_PI/2.0)*max_pg_iamp);
			for (i=1;i<ntab;i++) {
				sina = sin(a0)*tan(a)/(E1*sin(a0)+(1-E1)*tan(a));
				if (sina < 1)
					a = asin(sina);
				else
					a = (float) (M_PI/2.0);
				rfamp[i] = (int) (a/(M_PI/2.0)*max_pg_iamp);
			}
			break;
		default:           /* some power */
			for (i=0;i<ntab;i++) {
				rfamp[i] = opflip + (int) (pow(i,rfamp_flag)*(90-opflip)/pow(ntab-1,rfamp_flag));
				/* fprintf(stderr,"\nget_rfamp: rfamp[%d] = %d / ", i, rfamp[i]); */
				rfamp[i] = (int) ((float)rfamp[i]/90.0*max_pg_iamp);
				/*fprintf(stderr,"%d", rfamp[i]);*/
			}
			break;
	}  /* switch  */

}       /*  end get_rfamp */



@pg
/********************************************
 * dummylinks
 *
 * This routine just pulls in routines from
 * the archive files by making a dummy call.
 ********************************************/
void dummylinks()
{
	epic_loadcvs("thefile"); /* for downloading CVs */
}


