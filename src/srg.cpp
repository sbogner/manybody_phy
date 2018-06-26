// These are definitions of functions used to diagonalize an NxN 
// matrix H using the similarity renormalization group

#include "srg.h"

mat commutator(mat& A, mat& B){
	return A*B-B*A;
}

void split(mat& H, mat& Hd, mat& Hod, int N){

	for(int i = 0; i < N; ++i){
		for(int j = 0; j < N; ++j){
			if(i == j) Hd(i,j) = H(i,j);
			else Hod(i,j) = H(i,j);
		}
	}
}

void srg(mat& H, int N, double smax, double ds){

	for(double s = 0.0; s <= smax; s += ds){

		// split H into diag and off-diag parts
		mat Hd = zeros<mat>(N,N);
		mat Hod = zeros<mat>(N,N);
		split(H, Hd, Hod, N);

		// calculate generator eta and derivative dH
		mat eta = commutator(Hd, Hod);
		mat dH = commutator(eta, H);

		// construct ds*I where I is the identity matrix
		mat dsI = zeros<mat>(N,N);
		for(int i = 0; i < N; ++i) dsI(i,i) = ds;

		// forward euler
		H = H + dsI*dH;
	}

	// print
	//display_matrix(H, N);

}

void display_matrix(mat& H, int N){

	for(int i = 0; i < N; i++){
		cout << "[";
		for(int j = 0; j < N; j++){
			cout << setw(10) << setprecision(3) << H(i,j);
		}
		cout << "]" << endl;
	}
	cout << endl;
}
