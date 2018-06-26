#ifndef SOLVER_H
#define SOLVER_H

#include <iostream>
#include <cmath>
#include <string> 
#include <fstream>
#include <iomanip>
#include <vector>
#include "time.h"
#include "system.h"
#include "srg.h"

using namespace std;

class CSolver{

public:

	int particles_, sp_states_;
	double d_, gmin_, gmax_, gstep_;
	string filename_;	

	CSolver(int particles, int sp_states, double d, double gmin, double gmax, string filename);
	~CSolver(){}

	void mbpt2();

};

#endif