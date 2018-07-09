#ifndef DIAG_H
#define DIAG_H

#include <iostream>
#include <cmath>
#include <string> 
#include <fstream>
#include <iomanip>
#include <armadillo>
#include <vector>
#include "time.h"

using namespace std;
using namespace arma;

const double pi = 4.0*atan(1.0);

class SRG{

public:

	static const int max_ = 22;
	static const vector<double> B_{1.0, -1/2.0, 1/6.0, 0.0, -1/30.0, 
		                  0.0, 1/42.0, 0.0, -1/30.0, 0.0, 5/66.0, 
		                  0.0, -691/2730.0, 0.0, 7/6.0, 0.0, -3617/510.0, 
		                  0.0, 43867/798.0, 0.0, -174611/330, 0.0}
    
    vector<double> prefactor_;
    int N_;
    double s_, smax_, ds_;

    mat H0_, H_, Hd_, Hod_;
    mat eta_, omega_, derivative_;


	SRG(mat &H, int N, double smax, double ds);
	~SRG(){}

	void display_H();
	mat commutator(mat& A, mat& B);

	void split();
	void direct();
	int factorial(int n);
	mat nested_commutator(int n);
	void derivative();
	void magnus();

};


#endif