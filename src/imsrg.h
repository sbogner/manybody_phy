#ifndef IMSRG_H
#define IMSRG_H

#include <iostream>
#include <cmath>
#include <string> 
#include <fstream>
#include <iomanip>
#include <armadillo>
#include "time.h"
#include "system.h"

using namespace std;
using namespace arma;

class CIMSRG{

public:

	CIMSRG(double smax, double ds, vec holes, vec parts);
	~IMCSRG(){}

	int nholes_, nparts_, nbasis_;
	double smax_, ds_;
	vec holes_, parts_;
	mat basis_, f_, Gamma_;

	void generate_basis();

	void RK4_imsrg();

	void imsrg(vec snapshots, string filename);

	mat imsrg();
};

#endif
