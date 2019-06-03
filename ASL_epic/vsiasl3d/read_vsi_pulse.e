/* Reads in text files with the magnitude and phase of the VSI pulses 
inputs are the array pointers where the data will be stored:
	int* vsi_pulse_mag, 
        int* vsi_pulse_phs,
        int* vsi_pulse_grad,

	vsi_train_len:  this number is the number of points in the file.  
	The program assumes that the pulses are named like
	       sprintf(strVSIpulsename, "/usr/g/bin/luis/myRFpulses/myVSI_%d.rho.txt" , vsi_train_len);

output is the number of points that we read from the file
*/

int read_vsi_pulse( 
	int* vsi_pulse_mag, 
	int* vsi_pulse_phs,
	int* vsi_pulse_grad,
	int vsi_train_len)
{

        char    strVSIpulsename[256];
        FILE*   fpin;
        int     Npoints, i;
        int     DACMAX = 32766;
        int   	tmp;

        fprintf(stderr, "\n Reading VSI  pulses ... ");

        /*  Magnitude of RF */

        sprintf(strVSIpulsename, "/usr/g/bin/luis/myRFpulses/myVSI_%d.rho.txt" , vsi_train_len);
        fpin=fopen(strVSIpulsename, "r");
        if(fpin!=0) {     /* look for ext. file */
                i=0;
                while ( fscanf(fpin, "%d\n",  &tmp) != EOF)
                {
			tmp = 2*(tmp/2); 
                        vsi_pulse_mag[i++] = tmp ;
                        /* fprintf(stderr, "  mag[%d] = %d  ##  " ,  i-1  , vsi_pulse_mag[i-1]); */

                }
                fclose(fpin);

                Npoints = i;
                fprintf(stderr, "\n Read %s ... %d points " , strVSIpulsename, Npoints);
        }
        else{
                fprintf(stderr, "\n* didn't find %s ... ! ", strVSIpulsename);
        }


        /*  Phase of RF */

        sprintf(strVSIpulsename, "/usr/g/bin/luis/myRFpulses/myVSI_%d.theta.txt" , vsi_train_len);
        fpin=fopen(strVSIpulsename, "r");
        if(fpin!=0) {     /* look for ext. file */
                i=0;
                while ( fscanf(fpin, "%d\n",  &tmp) != EOF)
                {
			tmp = 2*(tmp/2); 
                        vsi_pulse_phs[i++] = (int)tmp;
                }
                fclose(fpin);

                Npoints = i;
                fprintf(stderr, "\n Read %s ... %d points " , strVSIpulsename, Npoints);
        }
        else{
                fprintf(stderr, "\n* didn't find %s ... ! ", strVSIpulsename);
        }

        /*  Gradients */

        sprintf(strVSIpulsename, "/usr/g/bin/luis/myRFpulses/myVSI_%d.grad.txt" , vsi_train_len);
        fpin=fopen(strVSIpulsename, "r");
        if(fpin!=0) {     /* look for ext. file */
                i=0;
                while ( fscanf(fpin, "%d\n",  &tmp) != EOF)
                {
			tmp = 2*(tmp/2); 
                        vsi_pulse_grad[i++] = (int)tmp;
                }
                fclose(fpin);

                Npoints = i;
                fprintf(stderr, "\n Read %s ... %d points " , strVSIpulsename, Npoints);
        }
        else{
                fprintf(stderr, "\n* didn't find %s ... ! ", strVSIpulsename);
        }

        /*  Phase of RF */
        return Npoints;

}

/* This is where we calcualte the phase of the RF pulses in order to get velocity selectivity */
int calc_vsi_phs_from_velocity (
		int* vsi_pulse_mag,
		int* vsi_pulse_phs,
		int* vsi_pulse_grad,
		float vel_target,
		int vsi_train_len,
		double	vsi_Gmax)
{
	/* GAMMA_H1 26754 in (rad/s)/Gauss */
	double phase_val;
	double grad_val;
	double pos=0.0;
	double dt = 4e-6;
	double delta_phs = 0.0;
	int 	i;
	double 	n = 0.0;
	int	DACMAX = 32766;
	int	tmp;

	for (i=1; i<vsi_train_len; i++)
	{
		/* from DAC units to radians */
		phase_val = M_PI * (double)(vsi_pulse_phs[i]) /  (double)FS_PI  ; 		
		/* from DAC units to G/cm */
		grad_val = vsi_Gmax * (double)(vsi_pulse_grad[i]) / (double)DACMAX ;

		/* increment the phase */
		phase_val -=  delta_phs;
		/* calc the phase gained by moving spins during THIS interval */
		pos += vel_target*dt; 
		delta_phs += GAMMA_H1 * grad_val * pos * dt ;

		/* from radians to DAC ... unwrap first  , then make them even numbers only. */	
		phase_val = atan2( sin(phase_val), cos(phase_val));
		tmp = (int)(phase_val / M_PI * FS_PI);
		vsi_pulse_phs[i] = 2*(tmp/2);

		/*
		if ((vsi_pulse_grad[i] ==0) && (vsi_pulse_grad[i-1] == 0) )
			pos = 0.0;
		*/
	
                fprintf(stderr, " \nphs[%d] = %d , pos = %f , grad_val=%f , phase_val=%f , delta_phs = %f; " , 
				 i  , vsi_pulse_phs[i], pos, grad_val, phase_val, delta_phs); 
	}

	for (i=0; i<vsi_train_len; i++)
	{
		if (vsi_pulse_mag[i] == 0) 
			vsi_pulse_phs[i] = 0;

	}

}	


