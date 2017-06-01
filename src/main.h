using namespace std;

float compute_abs_dist(float x, float y, float z){
	return sqrt(x*x + y*y + z*z);
}

float snapshot_to_redshift(int snapshot){
	float z_initial = 200;
	float z_final = 0;
	float n_steps = 500;
	double p_initial = 1/(1+z_initial);
	double p_final = 1/(1+z_final);
	double Delta = (p_final - p_initial)/n_steps;
	double z_snap = -1.0 + 1.0/( p_initial + (snapshot+1)*Delta );
	return z_snap;
}

float comput_comov_dist(float redshift, float Omega_M) {
	// compute distance, given redshift

	float result; 

    if (redshift==0.0){
    	result = 0.0;
    }

    else{ 
	    int dummy;
	    double x[1000], y[1000];
	    float c_light = 2.9979e5; //km/s
   
	    for(dummy=0;dummy<1000;dummy++){
	      	x[dummy] = 0.0 + redshift/999 * dummy;
	      	y[dummy] = 1.0/( sqrt( Omega_M*pow(1+x[dummy],3) + (1.0-Omega_M) ) );
	    }

    	gsl_interp_accel *acc = gsl_interp_accel_alloc ();
 		gsl_spline *spline = gsl_spline_alloc (gsl_interp_cspline, 1000);
   		gsl_spline_init (spline, x, y, 1000);
   		result = c_light/100.0 * gsl_spline_eval_integ (spline, x[0], x[dummy-1], acc); // comoving distance in Mpc/h
   		gsl_spline_free (spline);
   		gsl_interp_accel_free (acc);
    }
    
    return result;
}

string int_to_string(int integer){
	stringstream ss;
	ss << integer;
	return ss.str();
}