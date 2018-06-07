
#include "jacobi.h"

void write_eigs(mat& D, mat& U, double d, double g, string filename){


	ofstream ofile;

	// open file
	cout << "Writing to '" << filename << "'... " << endl;
	ofile.open(filename);

	// print eigenvalues and eigenvectors
	ofile << "# energy spacing d = " << d << endl;
	ofile << "# interaction strength g = " << g << endl;
	ofile << "# eigenvalues E_i, eigenvectors (C_0i, C_1i, C_2i, C_3i, C_4i, C_5i)\n" << endl;
	for(int j = 0; j < 6; ++j){

		ofile << left << setw(15) << setprecision(5) << D(j,j);

		for(int i = 0; i < 6; ++i){
			ofile << setw(15) << setprecision(5) << U(i,j);
		}
		ofile << endl;
	}

	ofile.close();
}

int main(int argc, char *argv[]){

	int n, count = 0;

	if(argc < 2){ cout << "Input number of interaction strength 'g' mesh points." << endl; }
	else{
		n = atoi(argv[1]);
	}

	double dg = 2.0/n;      // spacing between 'g' values
	double d = 1.0;         // energy level spacing

	for(double g = -1.0; g <= 1.0; g += dg){

		// set up Hamiltonian
		mat H = zeros<mat>(6,6);
		for(int i = 0; i < 6; ++i){
			for(int j = 0; j < 6; ++j){
				if(i+j != 5) H(i,j) = -g; 
			}
		}
		H(0,0) = 2*d-2*g;
		H(1,1) = 4*d-2*g;
		H(2,2) = 6*d-2*g;
		H(3,3) = 6*d-2*g;
		H(4,4) = 8*d-2*g;
		H(5,5) = 10*d-2*g;

		// set-up matrix of eigenvectors
		mat U = eye<mat>(6,6);

		// solve for eigenvalues and eigenvectors
		jacobi(H, U, 6);

		// write to file
		write_eigs(H, U, d, g, "eigs"+to_string(count)+".dat");

		count += 1;
	}

	return 0;
}