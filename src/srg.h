#ifndef SRG_H
#define SRG_H

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

mat commutator(mat& A, mat& B);
void split(mat& H, mat& Hd, mat& Hod, int N);
void srg(mat &H, int N, double smax, double ds);
void display_matrix(mat& A, int N);


#endif