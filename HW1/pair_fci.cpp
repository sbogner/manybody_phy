// HW 1

#include "jacobi.h"

void write_eigs(mat& D, mat& U, double d, double g, string filename){


	ofstream ofile;
	double lambda, u, uu;
	vec eigvals(6);
	ivec index(6);

	// sort eigenvalues by ascending order
	for(int i = 0; i < 6; i++){
		eigvals(i) = D(i,i);
	}
	sort(eigvals.begin(),eigvals.end());

	// get indices of ordered eigenvalues
	for(int i = 0; i < 6; i++){
		for(int j = 0; j < 6; j++){
			if(eigvals(i) == D(j,j)) index(i) = j;
		}
	}

	// open file
	cout << "Writing to '" << filename << "'... " << endl;
	ofile.open(filename);

	// print eigenvalues and eigenvectors
	ofile << "# energy spacing d = " << d << endl;
	ofile << "# interaction strength g = " << g << endl;
	ofile << "# eigenvalues E_i, eigenvectors (C_0i, C_1i, C_2i, C_3i, C_4i, C_5i)\n" << endl;
	for(int j = 0; j < 6; ++j){

		ofile << left << scientific << setw(15) << setprecision(5) << eigvals(j);

		for(int i = 0; i < 6; ++i){
			ofile << setw(15) << setprecision(5) << U(i,index(j));
		}
		ofile << endl;
	}

	ofile.close();
}

int main(int argc, char *argv[]){

	double d, g;

	if(argc < 2){ cout << "Input energy spacing 'd' and interaction strength 'g'." << endl; }
	else{
		d = atof(argv[1]);
		g = atof(argv[2]);
	}

	double gamma = d/g;

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
	write_eigs(H, U, d, g, "eigs.dat");

	return 0;
}