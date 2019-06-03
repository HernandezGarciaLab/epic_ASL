/* variable-density spiral
 * (that is constant-density up to a specified fraction of kmax)
 *
 * inputs (args)
 *       D = FOV, cm
 *       N = matrix size
 *	 Tmax = longest acquisition allowed, s
 *	 dts = output sample spacing, s
 *       alpha = vd parameter (> 1)
 *       kmaxfrac = fraction of kmax at which to switch from constant to vd 
 *                   (0 < kmaxfrac < 1) 
 *    inputs (CVs)
 *       gamp = design grad max, G/cm
 *       gslew = design slew rate, mT/m/ms
 *	 nramp = number of rampdown points
 *	 gtype = trajectory type (0 = spiral out, 1 = in, 2 = in/out)
 *    outputs (these vars defined in higher-level code)
 *       Gx, Gy
 *       gres
 *       pwgpre
 *       ag[x,y]pre
 *    time is in sec
 *
 *           rev 0 5/25/09 original
 *           rev 1 9/26/09 allow VD to hit amplitude-limited range
 *           rev 2 4/5/10 allow CD to hit amplitude-limited range
 *
 * author:catie chang, catie@stanford.edu
 * 
 * the constant-density spiral code is borrowed from genspiral5.e
 */


int genspiralcdvd(float D, int N, float Tmax, float dts, float alpha, float kmaxfrac)
{

  float gamma = 2*M_PI*GAM;
  float dt;
  float S0;
  int i, j, n;
  float kmax, r1;
  float areax, areay, area;
  float gx[8*GRESMAX], gy[8*GRESMAX];
  float kx[8*GRESMAX], ky[8*GRESMAX];
  int nn;
  FILE *fdebug;

 
  /* announce 
  printf("entering genspiralcdvd, D=%f, N=%d, alpha=%f, kmaxfrac=%f\n",D,N,alpha,kmaxfrac);
*/
  
  /* isection */
  kmax = M_PI*N/D;
  r1 = kmaxfrac*kmax;
  dt = dts*.25;  /* upsample by 4 for design */


  /* constant-density spiral */

  float Ts, T, ts, t, theta, dthdt, thetas, x, y;
  float a1, a2, beta, Gmax, bx, by, c, s;
  float gabs, gmax;
  float q;
  int nint;
  float lambda1;
 
 
  q = 5;		/* nice number  */
  S0 = gslew*100;
  nint = 1;  /* for now, allow one interleave only */

  /*  slew-rate limited approximation  */
  Ts = .666667/nint*sqrt(pow((double)(M_PI*N), (double)3.0)/(gamma*D*S0));

  lambda1 = 1/D;
  a2 = N*M_PI/(nint*pow(Ts, .666667));
  a1 = 1.5*S0/a2;
  beta = S0*gamma*D/nint;
  Gmax = a1*pow(Ts, .333333);
  gmax = 0;
  i = 0;
  nn = 0; /* will count total points in upsample g and k  */
  for (t = 0; t<=Ts; t+=dt)  {
    x = pow(t, 1.333333);
    theta = .5*beta*t*t/(q + .5*beta/a2*x);
    if(theta>(r1/lambda1))  /* then start the VD spiral */
      break;
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
    nn++; /* keep track */
  }
  printf("cd: finished slew-lim at amplitude gabs=%f, length nn=%d\n",gabs,nn);
  
  /*  gmax limited approximation  */
  if(Gmax > gamp)  {
    T = ((M_PI*N/nint)*(M_PI*N/nint) - thetas*thetas)/(2*gamma*gamp*D/nint) + ts;
    for (t=ts+dt; t<=T; t+=dt)  {
      theta = sqrt(thetas*thetas + 2*gamma*gamp*D*(t-ts)/nint);
      if(theta>(r1/lambda1)) { /* then start the VD spiral */
	break;
      }
      c = cos(theta);
      s = sin(theta);
      gx[i] = gamp*(c/theta - s);
      gy[i++] = gamp*(s/theta + c);
      nn++;  /* keep track */
    }
  } 
  printf("cd: finished amp-lim at amplitude gabs=%f, length nn=%d\n",gabs,nn);


  /* compute k-space (integrate) */
  x = y = 0;
  float kmaxabs = 0;
  for(j=0; j<nn; j++) {
    x += gamma*dt*gx[j];
    kx[j] = x;
    y += gamma*dt*gy[j];
    ky[j] = y;
    kmaxabs = hypot(kx[j],ky[j]);
  }
  printf("cd: max k-space amplitude reached = %f\n",kmaxabs);


  /* variable-density */

  int nturn;
  float lambda2, omega, beta2, gconst;
  float tau_d, tau_o, DeltaTau, DeltaR, ro, to;
  float tau, delta;
  float kx2_next, ky2_next, gx2_next, gy2_next;

  nturn = N/2;
  lambda2 = 2*M_PI*N/(2*D);
  omega = 2*M_PI*nturn;

  tau_d = pow(lambda1*omega/(lambda2*alpha),1/(alpha-1));
  tau_o = (r1/(omega*lambda1));
  DeltaTau = tau_o - tau_d;
  DeltaR = r1 - lambda2*pow(tau_d,alpha);

  /* slew-rate limited approximation */
  /*generate tau(t) iteratively from the differential equation*/
  ro = DeltaR/lambda2;
  to = DeltaTau;
  beta2 = S0*gamma/(lambda2*pow(omega,2));
  j=0; 
  tau = kmaxfrac; 

  while(tau<=1) {   
    delta = pow(ro + pow(tau-to,alpha),-0.5) * sqrt(beta2) * dt;
    tau += delta;

    /* compute functions on current tau */
    kx2_next = lambda2*(ro + pow(tau-to,alpha))*cos(omega*tau);
    ky2_next = lambda2*(ro + pow(tau-to,alpha))*sin(omega*tau);
    gconst = (lambda2/gamma) * sqrt(beta2/(ro+pow(tau-to,alpha)));
    gx2_next = gconst * (alpha*pow(tau-to,alpha-1)*cos(omega*tau) - (ro+pow(tau-to,alpha))*omega*sin(omega*tau));
    gy2_next = gconst * (alpha*pow(tau-to,alpha-1)*sin(omega*tau) + (ro+pow(tau-to,alpha))*omega*cos(omega*tau));
    
    /* break if done with trajectory (kmax) */
    if (hypot(kx2_next,ky2_next)>kmax) {
      printf("vd, slew-lim: max k-space amplitude reached = %f\n",kmax);
      break;
    }
    
    /* break if we exceed g0 */
    gabs = hypot(gx2_next, gy2_next);
    if(gabs > gamp) {
      printf("vd, slew-lim: max gradient amplitude reached = %f\n",gabs);
      break;
    }
    
    /* assign */
    kx[nn]=kx2_next;
    ky[nn]=ky2_next;
    gx[nn]=gx2_next;
    gy[nn]=gy2_next;
    nn++;    
  }

  /* amplitude-limited approximation */
  while(tau<=1) {    
    delta = pow(pow(alpha,2)*pow(tau-to,2*(alpha-1)) + pow(ro+pow(tau-to,alpha),2)*pow(omega,2),-0.5)*gamma*gamp/lambda2*dt;
    tau += delta;

    /* compute functions on current tau */
    kx2_next = lambda2*(ro + pow(tau-to,alpha))*cos(omega*tau);
    ky2_next = lambda2*(ro + pow(tau-to,alpha))*sin(omega*tau);
    gconst = gamp*pow(pow(alpha,2)*pow(tau-to,2*(alpha-1)) + pow(ro+pow(tau-to,alpha),2)*pow(omega,2),-0.5);
    gx2_next = gconst * (alpha*pow(tau-to,alpha-1)*cos(omega*tau) - (ro+pow(tau-to,alpha))*omega*sin(omega*tau));
    gy2_next = gconst * (alpha*pow(tau-to,alpha-1)*sin(omega*tau) + (ro+pow(tau-to,alpha))*omega*cos(omega*tau));

   
    /* break if done with trajectory (kmax) */
    if (hypot(kx2_next,ky2_next)>kmax) {
      printf("vd, amp-lim: max k-space amplitude reached = %f\n",kmax);
      break;
    }

    /* assign */
    kx[nn]=kx2_next;
    ky[nn]=ky2_next;
    gx[nn]=gx2_next;
    gy[nn]=gy2_next;
    nn++;
  }

/*
  fdebug = fopen("grads_psd.txt","w");
  for(i=0; i<nn; i++) {
    fprintf(fdebug,"%f %f\n",gx[i],gy[i]);
  }
  fclose(fdebug);
*/

  /* downsample (by 4) to dts */
  n = 0; /* maintain count (length) */
  for (j=0; j<nn; j+=4)  {
    gx[n] = gx[j];
    gy[n] = gy[j];
    kx[n] = kx[j];
    ky[n] = ky[j];
    n++;
  }

  /* cut off at appropriate place (kmax) */
  for(i=0;i<n;i++) {
    kmaxabs = hypot(kx[i],ky[i]);
    if(kmaxabs>kmax)
      break;
  }
  n = i;  /* update length */
 
  /* calculate area */
  areax = areay = 0;
  for (j=0; j<n; j++)  {
    areax += gx[j];
    areay += gy[j];
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
      printf("gtype=%d, gres=%d\n",gtype,gres);
      break;
    case 1:  	/* spiral in */
    case 2:  	/* spiral in-out */
      for (j=0; j<n; j++)  {
        Gx[n-j-1] = -2*((int)(.5*gx[j]*MAX_PG_WAMP/gamp));
        Gy[n-j-1] = -2*((int)(.5*gy[j]*MAX_PG_WAMP/gamp));
      }
      gres = n;
      printf("gtype=%d, gres=%d\n",gtype,gres);
      areax *= 4;  	/* dt  */
      areay *= 4;
      area = MAX(fabs(areax), fabs(areay));
      pwgpre = RUP_GRD((int)(area*.5*M_PI/gamp));
      if(pwgpre<1000) pwgpre = 1000;	/* 1 ms min */
      agxpre = .5*M_PI*areax/pwgpre;
      agypre = .5*M_PI*areay/pwgpre;
      break;
  }

  printf("Tslew, Ttot (ms), tot points =  %f  %f  %d\n", ts*1000, n*dts*1000, gres);
  return SUCCESS;

}
