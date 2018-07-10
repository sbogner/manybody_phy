#include <iostream>
#include <cmath>
#include <string> 
#include <fstream>
#include <iomanip>
#include <armadillo>
#include "time.h"

using namespace std;
using namespace arma;



static const int Bmax = 32;
static const vec B = {1.0, -1/2.0, 1/6.0, 0.0, -1/30.0, 
	                  0.0, 1/42.0, 0.0, -1/30.0, 0.0, 
	                  5/66.0, 0.0, -691/2730.0, 0.0, 7/6.0, 
	                  0.0, -3617/510.0, 0.0, 43867/798.0, 0.0, 
	                  -174611/330, 0.0, 854513/138.0, 0.0, -236364091/2730.0, 
	                  0.0, 8553103/6.0, 0.0, -23749461029/870, 0.0,
	                  8615841276005/14322, 0.0};


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

	if(n > 1){
		return n*factorial(n-1);
	}
	else return 1;
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

// returns A(smax) after solving dA/ds = derivative
// if A = H, B0 is a dummy 
// if A = Omega, B0 is H0 
mat RK4(mat A0, mat B0, int N, double smax, double ds, mat (*derivative)(mat,mat)){

	mat k1, k2, k3, k4;

	// algorithm
	mat A = A0, B = B0;

	for(double s = 0; s <= smax; s += ds){

		k1 = ds*(*derivative)(A,B);
		k2 = ds*(*derivative)(A+0.5*k1,B);
		k3 = ds*(*derivative)(A+0.5*k2,B);
		k4 = ds*(*derivative)(A+k3,B);

		A += (k1+2.0*k2+2.0*k3+k4)/6.0;
		B = expmat(A)*B0*expmat(-A);
	}

	return A;
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
mat get_dH(mat H, mat dummy){

	return commutator(get_eta(H),H);
}

// dOmega/ds
mat get_dOmega(mat Omega, mat H){

	double tolerance = 1.0E-5;
	mat eta = get_eta(H);

	// only one non-zero Bernoulli number
	mat sum = 0.5*nested_commutator(Omega,eta,1);

	// add for all even k
	for(int k = 0; k < Bmax; k += 2){
		sum += B[k]*nested_commutator(Omega,eta,k)/factorial(k);
	}

	return sum;

}

mat srg(mat& H0, int N, double smax, double ds){

	mat dummy = zeros<mat>(6,6);
	return RK4(H0, dummy, N, smax, ds, get_dH);

}

mat magnus(mat& H0, int N, double smax, double ds){

	mat Omega0 = zeros<mat>(6,6);
	mat Omega = RK4(Omega0, H0, N, smax, ds, get_dOmega);
	return expmat(Omega)*H0*expmat(-Omega);
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
