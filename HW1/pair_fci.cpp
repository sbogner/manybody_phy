// HW 1

#include "jacobi.h"

int main(int argc, char *argv[]){

	double d, g;

	if(argc < 2){ cout << "Input energy spacing 'd' and interaction strength 'g'." << endl; }
	else{
		d = atof(argv[1]);
		g = atof(argv[2]);
	}

	double gamma = d/g;
	cout << gamma << endl;

	// set up Hamiltonian
	mat H = zeros<mat>(6,6);
	for(int i = 0; i < 6; ++i){
		for(int j = 0; j < 6; ++j){
			if(i+j != 5) H(i,j) = -g; 
		}
	}

	print_matrix(H,6);


	// set-up matrix of eigenvectors
	mat U = eye<mat>(N,N);


}