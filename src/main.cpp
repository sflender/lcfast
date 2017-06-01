#include <iostream>
#include <mpi.h>
#include <vector>

#include "GenericIODefinitions.hpp"
#include "GenericIOReader.h"
#include "GenericIOMPIReader.h"
#include "GenericIOPosixReader.h"
#include "dtk/all.hpp"
#include <gsl/gsl_spline.h>
#include <gsl/gsl_integration.h>
#include "fitsio.h"
#include "main.h"

using namespace std;

int main(int argc, char** argv){

	dtk::AutoTimer t;

	MPI_Init(&argc,&argv);
	MPI_Comm comm = MPI_COMM_WORLD;
	int nranks, my_rank;
	MPI_Comm_size ( comm, &nranks );	
 	MPI_Comm_rank ( comm, &my_rank );
	if(my_rank==0) cout<< "I am rank "<<my_rank<<" of " <<nranks<<endl;

	if (argc!=2){
		if(my_rank==0) cout<<"usage: lcfast <param_file>"<<endl; return -1;
	}
	string input_file = argv[1];

	// get input parameters
 	dtk::Param prm(input_file);
 	float Omega_M = prm.get<float>("Omega_M");
 	float Lbox = prm.get<float>("Lbox");
 	int N_snap = prm.get_length("list_snapshots");
 	int* snapshot = new int[N_snap];
 	prm.get_array<int>("list_snapshots", snapshot);
	string fileroot = prm.get<string>("fileroot");
	string file_divider = prm.get<string>("file_divider");
	string file_end = prm.get<string>("file_end");
	string output_file = prm.get<string>("output_file");    

    // each rank should estimate the total lightcone count in order to allocate an array to store the data:
    gio::GenericIOMPIReader reader;
    reader.SetCommunicator(MPI_COMM_SELF);
    reader.SetFileName(fileroot + int_to_string(snapshot[N_snap-1]) + file_divider + int_to_string(snapshot[N_snap-1]) + file_end);
    reader.OpenAndReadHeader();
    long int N_halos_final = reader.GetNumberOfElements(); // number of halos in the final snapshot
    reader.Close();
    // estimated lightcone count (in one octant!)
	long int N_lc_est = (1.0/8.0) * (4.0/3.0) * M_PI * pow( comput_comov_dist(snapshot_to_redshift(snapshot[0]), Omega_M) , 3) * N_halos_final / pow(Lbox,3);

	float *lc_x = new float[N_lc_est];
	float *lc_y = new float[N_lc_est];
	float *lc_z = new float[N_lc_est];
	float *lc_mass = new float[N_lc_est];
	float *lc_redshift = new float[N_lc_est];

 	int lc_count_tot = 0;
 	long int N_halos;

 	// now construct the lightcone
 	for(int i=1; i<N_snap; i++){
 		
		string filename = fileroot + int_to_string(snapshot[i]) + file_divider + int_to_string(snapshot[i]) + file_end;

		if(my_rank==0) cout<<"...now reading file "<<filename<<" from snapshot "<<snapshot[i]<<" at redshift "<<snapshot_to_redshift(snapshot[i])<<endl;
 		float dmin = comput_comov_dist(snapshot_to_redshift(snapshot[i]), Omega_M);
 		float dmax = comput_comov_dist(snapshot_to_redshift(snapshot[i-1]), Omega_M);
		if(my_rank==0) cout<<"   ...dmin,dmax = "<<dmin<<", "<<dmax<<" Mpc/h\n";

 		float *x, *y, *z, *mass;
 		dtk::read_gio_quick(filename, "sod_halo_min_pot_x", x, N_halos);
  		dtk::read_gio_quick(filename, "sod_halo_min_pot_y", y, N_halos);
   		dtk::read_gio_quick(filename, "sod_halo_min_pot_z", z, N_halos);
      	dtk::read_gio_quick(filename, "sod_halo_mass", mass , N_halos);

		int maxrep = ceil(dmax/Lbox); 
		if(my_rank==0) cout<<"   ...maxrep = "<<maxrep<<endl;

 		int lc_count = 0;

	 	for (int h=0; h<N_halos; h++){
	 		for (int i=0; i<maxrep; i++){
	 			for (int j=0; j<maxrep; j++){
	 				for (int k=0; k<maxrep; k++){
	 					float d_obs_halo = compute_abs_dist(x[h]+i*Lbox, y[h]+j*Lbox, z[h]+k*Lbox);
	 					if (dmin<d_obs_halo and d_obs_halo<dmax){
	 						lc_x[lc_count_tot] = x[h]+i*Lbox;
	 						lc_y[lc_count_tot] = y[h]+j*Lbox;
	 						lc_z[lc_count_tot] = z[h]+k*Lbox;
	 						lc_mass[lc_count_tot] = mass[h];
	 						lc_count++;
	 						lc_count_tot++;
	 						if(lc_count_tot>N_lc_est){
	 							"ERROR: lc_count_tot>N_lc_est \n"; return -1;
	 						}
	 					}
	 				}
	 			}
	 		}
	 	}

	 	delete [] x,y,z,mass;

 		if(my_rank==0) cout<<"   ...lightcone count in snapshot "<<snapshot[i]<<": "<<lc_count<<"\n";
	}

 	if(my_rank==0) cout<<"...total count : "<<lc_count_tot<<", estimated: "<< N_lc_est<<endl;



 	//-------------//
 	if(my_rank==0) cout<<"...computing redshifts\n";

 	int interp_size = 1000;
    double *d_array = new double[interp_size+1]; // array of d values for interpolation
    double *zeta = new double[interp_size+1]; // array of zeta-values (zeta = z/d) for interpolation
    float zhigh = 2.0 * snapshot_to_redshift(snapshot[0]);
    float zlow = snapshot_to_redshift(snapshot[N_snap-1]); // interpolation range

    float Delta_z = (zhigh-zlow)/interp_size;

    for (int i=0; i<interp_size+1; i++){
        float z = zlow + i*Delta_z;
        d_array[i] = comput_comov_dist(z, Omega_M);
        zeta[i]= z/d_array[i];
        if (d_array[i]==0){zeta[i]=0;}
    }

    {
        gsl_interp_accel *acc = gsl_interp_accel_alloc ();
        gsl_spline *spline = gsl_spline_alloc (gsl_interp_cspline, interp_size);
        gsl_spline_init (spline, d_array, zeta, interp_size);

        for (int h=0; h<lc_count_tot; h++){
        	float comov_dist = compute_abs_dist(lc_x[h],lc_y[h],lc_z[h]);
        	lc_redshift[h] = comov_dist * gsl_spline_eval(spline, comov_dist, acc);
        }
        

        gsl_spline_free (spline);
        gsl_interp_accel_free (acc);
    }

	//----------//

 	if(my_rank==0) cout<<"...writing output file "<<output_file<<endl;

  	gio::GenericIO gio(MPI_COMM_WORLD, output_file);
  	gio.setNumElems(lc_count_tot);
  	gio.addVariable("x",lc_x,gio::GenericIO::VarIsPhysCoordX);
  	gio.addVariable("y",lc_y,gio::GenericIO::VarIsPhysCoordY);
  	gio.addVariable("z",lc_z,gio::GenericIO::VarIsPhysCoordZ);
  	gio.addVariable("M200red",lc_mass);
  	gio.addVariable("redshift",lc_redshift);
	gio.write();

  	if(my_rank==0) cout<<"time: "<<t<<endl;
	if(my_rank==0) cout<<"all done!\n";
	
	MPI_Finalize();
	return 0;
}