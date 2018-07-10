#include <iostream>
#include <cmath>
#include <string> 
#include <fstream>
#include <iomanip>
#include <armadillo>
#include "time.h"

using namespace std;
using namespace arma;



static const int Bmax = 22;
static const vec B = {1.0, -1/2.0, 1/6.0, 0.0, -1/30.0, 
	                  0.0, 1/42.0, 0.0, -1/30.0, 0.0, 5/66.0, 
	                  0.0, -691/2730.0, 0.0, 7/6.0, 0.0, -3617/510.0, 
	                  0.0, 43867/798.0, 0.0, -174611/330, 0.0};

// displays NxN matrix
void display(mat& A, int N){
	
	for(int i = 0; i < N; i++){
		cout << "[";
		for(int j = 0; j < N; j++){
			cout << setw(10) << setprecision(3) << A(i,j);
		}
		cout << "]" << endl;
	}
	cout << endl;	
}

int factorial(int n){

	int factorial = 1;

	if(n > 1){
		for(int i = 2; i <= n; ++i) factorial *= i;
	}
	
	return factorial;
}

mat commutator(mat A, mat B){
	
	return A*B-B*A;
}

mat nested_commutator(mat A, mat B, int n){

	if(n == 0) return B;
	else{
		mat ad = B;
		for(int i = 0; i < n; ++i) ad = commutator(A,ad);
		return ad;
	}
}

// returns A after solving dA/ds=derivative(A(s),s)
// later add display statements for various s
mat RK4(mat A0, int N, double smax, double ds, mat (*derivative)(mat,mat)){

	mat k1, k2, k3, k4;

	// this is a dummy variable if srg_derivative() is used
	mat Omega  = zeros<mat>(N,N);

	// algorithm
	mat Atemp = A0;
	for(double s = 0; s <= smax; s += ds){

		k1 = ds*(*derivative)(Atemp,Omega);
		k2 = ds*(*derivative)(Atemp+0.5*k1,Omega);
		k3 = ds*(*derivative)(Atemp+0.5*k2,Omega);
		k4 = ds*(*derivative)(Atemp+k3,Omega);

		Atemp += (k1+2.0*k2+2.0*k3+k4)/6.0;
	}

	return Atemp;
}

// generator
mat get_eta(mat H){

	// split H into diag and off-diag parts
	mat Hd = diagmat(H);
	mat Hod = H-Hd;

	// calculate generator
	return commutator(Hd,Hod);	
}

// dH/ds 
// need dummy to pass function as argument in RK4
mat srg_derivative(mat H, mat dummy){

	return commutator(get_eta(H),H);
}

mat srg(mat& H0, int N, double smax, double ds){

	return RK4(H0, N, smax, ds, srg_derivative);

}

// dOmega/ds
mat magnus_derivative(mat H, mat Omega){

	double tolerance = 1.0E-5;
	mat eta = get_eta(H);

	// only one non-zero Bernoulli number
	mat sum = B[1]*nested_commutator(Omega,eta,1);

	// add for all even k
	for(int k = 0; k < Bmax; k += 2){
		sum = B[k]*nested_commutator(Omega,eta,k)/factorial(k);
	}

	return sum;

}

mat magnus(mat& H0, int N, double smax, double ds){

	return RK4(H0, N, smax, ds, magnus_derivative);

}


int main(int argc, char *argv[]){

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

	cout << "\nOriginal Hamiltonian:\n" << endl;
	display(H0, 6);	

	mat H = srg(H0, 6, smax, ds);

	cout << "\nDiagonalized Hamiltonian (SRG)\n" << endl;
	display(H, 6);

	H = magnus(H0, 6, smax, ds);

	cout << "\nDiagonalized Hamiltonian (SRG w/ Magnus Expansion)\n" << endl;
	display(H, 6);

	return 0;
}
