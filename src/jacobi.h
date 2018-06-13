#ifndef JACOBI_H
#define JACOBI_H

#include <iostream>
#include <cmath>
#include <string> 
#include <fstream>
#include <iomanip>
#include <armadillo>
#include "time.h"

using namespace std;
using namespace arma;

const double pi = 4.0*atan(1.0);

double offdiag_sq(mat& A, int N);
double norm_sq(mat& A, int N);
void get_pivot(mat& A, int N, int& k, int& l);
void rotate(mat& A, mat& V, int k, int l, int N);
void jacobi(mat& A, mat& V, int N);
void print_matrix(mat& A, int N);
void write_eigs(mat& D, mat& U, double d, double g, string filename);

#endif