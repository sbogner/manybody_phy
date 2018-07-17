#include "system.h"
#include "solver.h"


int main(int argc, char *argv[]){

	CSolver PairingModel(4, 4, 1.0, -1.0, 1.0, "pairingmodel");
	PairingModel.MBPT2();
	PairingModel.SRG();

	/*

	double d = 1;
	double g = atof(argv[1]);
	double smax = atof(argv[2]);
	double ds = atof(argv[3]);

	// set up Hamiltonian
	mat H0 = zeros<mat>(6,6);
	for(int i = 0; i < 6; ++i){
		for(int j = 0; j < 6; ++j){
			if(i+j != 5) H0(i,j) = -0.5*g; 
		}
	}
	H0(0,0) = 2*d-g;
	H0(1,1) = 4*d-g;
	H0(2,2) = 6*d-g;
	H0(3,3) = 6*d-g;
	H0(4,4) = 8*d-g;
	H0(5,5) = 10*d-g;

	SRG pairing_model(H0, 6, 170, smax, ds);

	vec snapshots = linspace<vec>(0.0, smax);
	pairing_model.srg(snapshots, "pairing_srg.dat");
	pairing_model.magnus(snapshots, "pairing_magnus.dat");

	*/


	return 0;
}