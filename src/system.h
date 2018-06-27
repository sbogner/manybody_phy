#ifndef SYSTEM_H
#define SYSTEM_H

#include <iostream>
#include <cmath>
#include <string> 
#include <fstream>
#include <iomanip>
#include <vector>
#include "time.h"

using namespace std;

class CFermionSystem{

public:

	int N_, F_;
	int particles_, sp_states_;
	double d_, g_;
	

	vector<vector<int>> states_;
	vector<vector<double>> H_;

	CFermionSystem();
	CFermionSystem(int particles, int sp_states, double d, double g);
	~CFermionSystem(){}

	void generate_sp_basis();
	void generate_H();

	double h0(int q, int r);
	double v(int q, int r, int s, int t);
	double f(int q, int r);
	double eps(int q, int r);
	double eps(vector<int>& holes, vector<int>& parts);

};

#endif