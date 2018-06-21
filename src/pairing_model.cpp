#include "srg.h"


int main(int argc, char *argv[]){

	double smax, ds;
	int n, count = 0;

	if(argc < 4){ cout << "Input maximum s-value 'smax', s-spacing 'ds', and number of mesh points 'n'." << endl; }
	else{
		smax = atof(argv[1]);
		ds = atof(argv[2]);
		n = atoi(argv[3]);
	}


	double dg = 2.0/n;           // spacing between 'g' values
	double d = 1.0;              // energy level spacing

	for(double g = -1.0; g <= 1.0; g += dg){

		cout << "g = " << g << endl;
		
		// set up Hamiltonian
		mat H = zeros<mat>(6,6);
		for(int i = 0; i < 6; ++i){
			for(int j = 0; j < 6; ++j){
				if(i+j != 5) H(i,j) = -0.5*g; 
			}
		}
		H(0,0) = 2*d-g;
		H(1,1) = 4*d-g;
		H(2,2) = 6*d-g;
		H(3,3) = 6*d-g;
		H(4,4) = 8*d-g;
		H(5,5) = 10*d-g;


		// diagonalize Hamiltonian
		srg(H, 6, smax, ds);

		count += 1;
	}

	return 0;
}