#include "diag.h"

SRG::SRG(mat &H, int N, double smax, double ds){

	H0_ = H;
	N_ = N;
	smax_ = smax;
	ds_ = ds;
	dsI_ = zeros<mat>(N_,N_);
	for(int i = 0; i < N_; ++i) dsI(i,i) = ds;
	for(int n = 0; n < max_; ++n) prefactor_[n] = B_[n]/factorial(n);
}

void SRG::display_H(){
	
	for(int i = 0; i < N_; i++){
		cout << "[";
		for(int j = 0; j < N_; j++){
			cout << setw(10) << setprecision(3) << H_(i,j);
		}
		cout << "]" << endl;
	}
	cout << endl;	
}

mat SRG::commutator(mat& A, mat& B){
	
	return A*B-B*A;
}

void SRG::split(){

	Hd_ = diagmat(H_);
	Hod_ = H_-Hd_;
}

void SRG::direct(){

	H_ = H0_;

	for(s_ = 0.0; s_ <= smax_; s_ += ds_){

		// split H into diag and off-diag parts
		split();

		// calculate generator eta and derivative dH/ds
		eta_ = commutator(Hd_,Hod_);
		derivative_ = commutator(eta_,H_);

		// forward euler
		H_ = H_+ds_*derivative_;
	}

	display_H();
}

int SRG::factorial(int n){

	int factorial = 1;

	if(n > 1){
		for(int i = 2; i <= n; ++i) factorial *= i;
	}
	else return factorial;
}

mat SRG::nested_commutator(int n){

	if(n == 0) return eta_;
	else{
		mat ad = eta_;
		for(int i = 0; i < n; ++i) ad = commutator(omega_,ad);
		return ad;
	}
}

mat SRG::derivative(){

	int n = 0;
	double tolerance = 1.0E-5;
	mat sum, summand; 

	do{
 		summand = prefactor_[n]*nested_commutator(n);
 		sum += summand;
	}while(summand > tolerance);

	derivative_ = sum;
}

void SRG::magnus(){
	

}