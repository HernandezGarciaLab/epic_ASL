/*   multi- shot spiral design
*    uses Duyn's approximate slewrate limited design
*    augmented with archimedian gmax limit
*    inputs (args)
*        D = FOV, cm
*        N = matrix size
*	 Tmax = longest acquisition allowed, s
*	 dts = output sample spacing, s
*    inputs (CVs)
*        nl = number of interleaves
*        gamp = design grad max, G/cm
*        gslew = design slew rate, mT/m/ms
*	 nramp = number of rampdown points
*	 gtype = trajectory type (0 = spiral out, 1 = in-out)
*    outputs
*        Gx, Gy
*        gres
*        pwgpre
*        ag[x,y]pre
*    time is in sec
*
*		rev 0 12/26/98	original
*		rev 1 4/15/99	little better calc of ts
*		rev 3 1/10/00	can do spiral in-out (gtype = 1)
*		rev 4 4/10/00	can do spiral in, too (gtype = 2)
*		      5/27/00	fix for 8.3 (replace serror)
*		rev 5 10/22/00	1 = sprl in, 2 = single shot in-out
*/

int genspiral(float D, int N, float Tmax, float dts)
{
  float gamma = 2*M_PI*GAM;

  float dt;
  float S0;
  int i, j, n;
  float Ts, T, ts, t, theta, dthdt, thetas, x, y;
  float a1, a2, beta, Gmax, bx, by, c, s;
  float gx[2*GRESMAX], gy[2*GRESMAX];
  float gabs, gmax;
  float q;
  int nint;
  float areax, areay, area;
#define LIMIT(x, xmin, xmax)   ( (x<xmax)? ((x>xmin)? x:xmin):xmax )


  q = 5;		/* nice number  */
  S0 = gslew*100;
  dt = dts*.5;

  printf("entering genspiral  N, D, nl = %d  %6.2f cm  %d\n", N, D, nl);

  nint = nl;
  if(gtype > 2)  {
    epic_error(use_ermes,"wrong trajectory type",0,0);
    return FAILURE;
  }

/*  slew-rate limited approximation  */

  Ts = .666667/nint*sqrt(pow((double)(M_PI*N), (double)3.0)/(gamma*D*S0));
  if(Ts > Tmax)  {
    epic_error(use_ermes,"slew limited readout too long\n",0,0);
    return FAILURE;
  }
  a2 = N*M_PI/(nint*pow(Ts, .666667));
  a1 = 1.5*S0/a2;
  beta = S0*gamma*D/nint;
  Gmax = a1*pow(Ts, .333333);
  gmax = 0;
  i = 0;
  for (t = 0; t<=Ts; t+=dt)  {
    x = pow(t, 1.333333);
    theta = .5*beta*t*t/(q + .5*beta/a2*x);
    y = q+.5*beta/a2*x;
    dthdt = beta*t*(q+.166667*beta/a2*x)/(y*y);
    c = cos(theta);
    s = sin(theta);
    gx[i] = nint/(D*gamma)*dthdt*(c - theta*s);
    gy[i] = nint/(D*gamma)*dthdt*(s + theta*c);
    gabs = hypot(gx[i], gy[i]);
    if(gabs>=gamp)  {
      if(gmax==0)  gmax = hypot(gamp/theta, gamp);
      if(gabs>gmax) break;
    }
    ts = t;
    thetas = theta;
    i++;
  }

/*  gmax limited approximation  */

  if(Gmax > gamp)  {
    T = ((M_PI*N/nint)*(M_PI*N/nint) - thetas*thetas)/(2*gamma*gamp*D/nint) + ts;
    if(T > Tmax)  {
      epic_error(use_ermes,"gmax limited readout too long\n",0,0);
      return FAILURE;
    }
    for (t=ts+dt; t<=T; t+=dt)  {
      theta = sqrt(thetas*thetas + 2*gamma*gamp*D*(t-ts)/nint);
      c = cos(theta);
      s = sin(theta);
      gx[i] = gamp*(c/theta - s);
      gy[i++] = gamp*(s/theta + c);
    }
  }

/*  decimate by 2 to get back to 4us sampling */

  n = 0;
  areax = areay = 0;
  for (j=0; j<i; j+=2)  {
    gx[n] = LIMIT(gx[j], -gamp, gamp);
    gy[n] = LIMIT(gy[j], -gamp, gamp);
    areax += gx[n];
    areay += gy[n++];
  }

/*  add ramp down to zero  */

  bx = gx[n-1];
  by = gy[n-1];
  for (j=0; j<nramp; j++)  {
    c = 1 - j/(float)nramp;
    gx[n] = bx*c;
    gy[n] = by*c;
    areax += gx[n];
    areay += gy[n++];
  }

/*  final scaling and fill output buffer  */
 
  switch(gtype)  {
    case 0:  	/* conventional spiral out */
      for (j=0; j<n; j++)  {
        Gx[j] = 2*((int)(.5*gx[j]*MAX_PG_WAMP/gamp));
        Gy[j] = 2*((int)(.5*gy[j]*MAX_PG_WAMP/gamp));
      }
      gres = n;
      break;
    case 1:  	/* spiral in */
    case 2:  	/* spiral in-out */
      for (j=0; j<n; j++)  {
        Gx[n-j-1] = -2*((int)(.5*gx[j]*MAX_PG_WAMP/gamp));
        Gy[n-j-1] = -2*((int)(.5*gy[j]*MAX_PG_WAMP/gamp));
      }
      gres = n;
      areax *= 4;  	/* *dt  */
      areay *= 4;
      area = MAX(fabs(areax), fabs(areay));
      pwgpre = RUP_GRD((int)(area*.5*M_PI/gamp));
      if(pwgpre<1000) pwgpre = 1000;	/* 1 ms min */
      agxpre = .5*M_PI*areax/pwgpre;
      agypre = .5*M_PI*areay/pwgpre;
      break;
  }

  printf("Tslew, Ttot (ms), tot points =  %f  %f  %d\n", ts*1000, T*1000, gres);

  return SUCCESS;

}

