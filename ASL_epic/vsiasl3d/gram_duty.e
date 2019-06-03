/* This routine checks gram heating limitations, using 
*  the slewing duty cycle limit.
*
*       rev 00  9/24/02         GH Glover
*	rev 1	2/22/08		for dvmr- doesn't need to do anything
*/

int gram_duty(void)
{
  float gradduty;		/* readout fraction of TR */
  float allowed_duty;
  int newtr;

  if(cfgradamp==8920)      	/* DVMR 5 G/cm SR200 */
    allowed_duty = 1.0;		/* pedal to the metal */
  else
    allowed_duty = 0.65;	/* Maybe ACGD's or SGD's */

  gradduty = opslquant*(float)pw_gx/optr;
  if(gtype == 2) gradduty *= 2.0;	/* in-out waveform  */
  printf("grad slew_frac, allowed = %g  %g\n", gradduty, allowed_duty);

  if(gradduty > allowed_duty)  {
      newtr = optr*allowed_duty/gradduty;
      epic_error(use_ermes,"Oops! Grams (slew) want optr > %.1f ms.\n",
	  EM_PSD_SUPPORT_FAILURE,1,FLOAT_ARG,newtr/1000.0);
      return FAILURE;
  }
    
  return SUCCESS;

}
