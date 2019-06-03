
/*
sph2cart 
converts arrays of spherical coordinates into cartesian coordinates
*/
int sph2cart(
	float*	pfx, 
	float*	pfy, 
	float*	pfz, 
	float*	pfr, 
	float*	pftheta, 
	float*	pfphi, 
	int	Npts)
{
int c;
float r, cth,sth, cphi, sphi;
fprintf(stderr, "\ntranslating spherical to cartesian ...");
for (c=0;c<Npts;c++)
{
	sth = sin(pftheta[c]);
	cth = cos(pftheta[c]);
	sphi = sin(pfphi[c]);
	cphi = cos(pfphi[c]);
	r = pfr[c];
 
	pfx[c] = r*sphi*cth; 
	pfy[c] = r*sth*sphi; 
	pfz[c] = r*cphi; 

/*		
	fprintf(stderr, "\nCART: %f \t %f \t %f \t SPHERICAL %f \t %f \t %f \n", 
		pfx[c], pfy[c], pfz[c], pfr[c], pftheta[c], pfphi[c]);
*/	
}
return 1;

} 
	
/* 
genspins3d
Generate 3 spiral waveforms for a 3D yarn trajectory 
	float *gx,		pointer to gradient waveform data  
	float *gy, 
	float *gz,
	int   iGsize,		length of gradient waveform data buffers 
	int   iNshots_theta,	Number of interleaves along the azimuthal angle
	int   iNshots_phi,	Number of interleaves along the elevation angle
	float fFOV,		image space Field of View
	float fXres,		image space resolution
	float fDeltaT		waveform sample time (granularity)
*/
int genspins3d(
	float *gx, 
	float *gy, 
	float *gz,
	int   iGsize, 
	float   Nshots_theta,
	float   Nshots_phi,
	float fFOV,
	float fXres,
	float fDeltaT
	 )
{

float gamma = 267.5 * 1e2 ; /* rad/s/Gauss  */
float Kmax, deltaK;
float Ntheta, Nphi;
float r[iGsize+1];
float t[iGsize+1];
float theta[iGsize+1];
float phi[iGsize+1];
float phi2;
float Nturns;  /*how many times turns of the spiral are needed in plane to get adequate deltaK*/
float rampSlopex, rampSlopey, rampSlopez;
float x, y, z;
float areax, areay, areaz;
float kx[iGsize];
float ky[iGsize];
float kz[iGsize];
int i,j;
float dr, dtheta, dphi, dphi2;
float  maxslewx, maxslewy, maxslewz;
FILE *fid;
FILE *fid2;
FILE *fid3;
int  NRAMPDN = 100;
int  NRAMPUP = 20;

Kmax = 1.0/(2.0*fXres);
deltaK = 1.0/(fFOV);
dr = Kmax/iGsize;
Nturns = Kmax/deltaK;


/*
dtheta = atan(deltaK/Kmax);
dphi = atan(deltaK/Kmax);
*/
/* alternative calculation using the spacing between the turns */
dtheta = (Nturns*2.0*M_PI)/((float)iGsize*Nshots_theta);
dphi = (Nturns*2.0*M_PI)/((float)iGsize*Nshots_phi);

fprintf(stderr, "\nNshots_theta = %f", Nshots_theta);
fprintf(stderr, "\nNshots_phi = %f", Nshots_phi);
fprintf(stderr, "\n Generating yarn threads ...\n FOV = %f,  xres=  %f,   deltaK=%f, Kmax=%f, Nturns= %f, dtheta= %f , dphi=%f \n", 
	fFOV, fXres, deltaK, Kmax, Nturns, dtheta, dphi);

/* make the spherical coordinates  */
for (i=0; i<iGsize; i++)
{
	t[i] = i*fDeltaT;
	r[i] = i*dr;
	theta[i] = i*dtheta; 
	phi[i] = i*dphi ; 
	/* phi = 0;*/  /*keep it 2D for now */

}
fprintf(stderr,"\nTranslating from Spherical to Cartesian k-space traj ...");
sph2cart(kx, ky, kz, r,theta, phi, iGsize);
fprintf(stderr,"...done");

/* The Gradients are the derivative of the kspace trajectory */
areax=0.0;
areay=0.0;
areaz=0.0;

fprintf(stderr,"\nCalculating grads from k-space traj ...");
gx[0] = 0.0; gy[0] = 0.0; gz[0] = 0.0;
for (i=1; i<iGsize-NRAMPDN; i++)
{
	gx[i] = 2.0*M_PI*(kx[i] - kx[i-1])/fDeltaT/gamma;	
	gy[i] = 2.0*M_PI*(ky[i] - ky[i-1])/fDeltaT/gamma;	
	gz[i] = 2.0*M_PI*(kz[i] - kz[i-1])/fDeltaT/gamma;	
	
	areax+=gx[i];
	areay+=gy[i];
	areaz+=gz[i];

}

/* The first NRAMPUP points are used to ramp the up gradients gentlyso they don't blow up */
rampSlopex = gx[NRAMPUP-1]/NRAMPUP;
rampSlopey = gy[NRAMPUP-1]/NRAMPUP;
rampSlopez = gz[NRAMPUP-1]/NRAMPUP;

fprintf(stderr, "\nFixing the Ramp up: %f , %f , %f ", rampSlopex, rampSlopey, rampSlopez);

for (i=1; i<NRAMPUP; i++)
{
        gx[i] = gx[i-1] + rampSlopex;
        gy[i] = gy[i-1] + rampSlopey;
        gz[i] = gz[i-1] + rampSlopez;
}

/* The last NRAMPDN points are used to ramp the gradients back down to zero */
rampSlopex = gx[iGsize-NRAMPDN-1]/NRAMPDN;
rampSlopey = gy[iGsize-NRAMPDN-1]/NRAMPDN;
rampSlopez = gz[iGsize-NRAMPDN-1]/NRAMPDN;

fprintf(stderr, "\nFixing the Ramp Down: %f , %f , %f ", rampSlopex, rampSlopey, rampSlopez);

for (i=iGsize-NRAMPDN; i<iGsize; i++)
{
	gx[i] = gx[i-1] - rampSlopex; 	
	gy[i] = gy[i-1] - rampSlopey; 	
	gz[i] = gz[i-1] - rampSlopez; 	

	areax+=gx[i];
	areay+=gy[i];
	areaz+=gz[i];

}
gx[iGsize]=0.0; gy[iGsize] = 0.0; gx[iGsize]=0.0;

/* dump to text file for later recon */
fid = fopen("grad.txt","w");
fid2 = fopen("ktraj.txt","w");
fid3 = fopen("ktraj_sph.txt","w");
for (i=0; i<iGsize; i++){
	fprintf(fid , "%f\t%f\t%f\n", gx[i], gy[i], gz[i]);
	fprintf(fid2, "%f\t%f\t%f\n", kx[i], ky[i], kz[i]);
	fprintf(fid3, "%f\t%f\t%f\n", r[i], theta[i], phi[i]);
}
fclose (fid);
fclose (fid2);
fclose (fid3);


/*  Check the slew rates */
maxslewx = 0;
maxslewy = 0;
maxslewz = 0;
for (i=0; i<iGsize-2; i++)
{
	if (fabs(gx[i+1] - gx[i]) > maxslewx) 
		maxslewx = gx[i+1] - gx[i];
	if (fabs(gy[i+1] - gy[i]) > maxslewy) 
		maxslewy = gy[i+1] - gy[i];
	if (fabs(gz[i+1] - gz[i]) > maxslewz)
		maxslewz = gz[i+1] - gz[i];

}

fprintf(stderr, "\n\nGradient moments  are: %f %f %f", 
	areax, areay, areaz);
fprintf(stderr, "\n\nMaximum Slew RATES  are: %f %f %f", 
	maxslewx, maxslewy, maxslewz);


fprintf(stderr,"\n ..... end genspiral3d \n");

return 1;

}
/* 

rotGradWaves 
rotate the gradient waveforms by a matrix 

*/
int rotGradWaves(
	float* 	gx, 
	float* 	gy, 
	float* 	gz, 
	float* 	gx2, 
	float* 	gy2, 
	float* 	gz2, 
	int   		glen,
	float*		mat)
{
int c;
for(c=0; c<glen; c++)
{
	gx2[c] = mat[0]*gx[c] + mat[1]*gy[c] + mat[2]*gz[c];
	gy2[c] = mat[3]*gx[c] + mat[4]*gy[c] + mat[5]*gz[c];
	gz2[c] = mat[6]*gx[c] + mat[7]*gy[c] + mat[8]*gz[c];

}

fprintf(stderr,"RotMat:  ");
for (c=0; c<9 ; c++) fprintf(stderr,"%f \t", mat[c]);
fprintf(stderr,"\n");


return 1;
}



/* 
euler2mat:   generates a rotation matrix from Euler angles
takes from Varian at  psg/GradientBase.cpp.... see ::calc_obl_matrix
*/
int euler2mat(float ang1,float ang2,float ang3,float *tm)
{
    float D_R;
    float sinang1,cosang1,sinang2,cosang2,sinang3,cosang3;
    float m11,m12,m13,m21,m22,m23,m31,m32,m33;
    float im11,im12,im13,im21,im22,im23,im31,im32,im33;
    float tol = 1.0e-14;

    /* Convert the input to the basic mag_log matrix */
    D_R = M_PI / 180;
    D_R = 1.0; /*Varian works in degrees, this program is is radians, though */

    cosang1 = cos(D_R*ang1);
    sinang1 = sin(D_R*ang1);

    cosang2 = cos(D_R*ang2);
    sinang2 = sin(D_R*ang2);

    cosang3 = cos(D_R*ang3);
    sinang3 = sin(D_R*ang3);

    m11 = (sinang2*cosang1 - cosang2*cosang3*sinang1);
    m12 = (-1.0*sinang2*sinang1 - cosang2*cosang3*cosang1);
    m13 = (sinang3*cosang2);

    m21 = (-1.0*cosang2*cosang1 - sinang2*cosang3*sinang1);
    m22 = (cosang2*sinang1 - sinang2*cosang3*cosang1);
    m23 = (sinang3*sinang2);

    m31 = (sinang1*sinang3);
    m32 = (cosang1*sinang3);
    m33 = (cosang3);

    if (fabs(m11) < tol) m11 = 0;
    if (fabs(m12) < tol) m12 = 0;
    if (fabs(m13) < tol) m13 = 0;
    if (fabs(m21) < tol) m21 = 0;
    if (fabs(m22) < tol) m22 = 0;
    if (fabs(m23) < tol) m23 = 0;
    if (fabs(m31) < tol) m31 = 0;
    if (fabs(m32) < tol) m32 = 0;
    if (fabs(m33) < tol) m33 = 0;

    /* Generate the transform matrix for mag_log ******************/

    /*HEAD SUPINE*/
        im11 = m11;       im12 = m12;       im13 = m13;
        im21 = m21;       im22 = m22;       im23 = m23;
        im31 = m31;       im32 = m32;       im33 = m33;

    /*Transpose intermediate matrix and return***********
    *tm11 = im11;     *tm21 = im12;     *tm31 = im13;
    *tm12 = im21;     *tm22 = im22;     *tm32 = im23;
    *tm13 = im31;     *tm23 = im32;     *tm33 = im33;
*/
    /*Transpose intermediate matrix and return***********/
    tm[0] = im11;     tm[1] = im12;     tm[2] = im13;
    tm[3] = im21;     tm[4] = im22;     tm[5] = im23;
    tm[6] = im31;     tm[7] = im32;     tm[8] = im33;


   return 1;
}



