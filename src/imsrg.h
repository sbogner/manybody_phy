// In-medium SRG for the pairing model

#ifndef IMSRG_H
#define IMSRG_H

#include <iostream>
#include <cmath>
#include <string> 
#include <fstream>
#include <iomanip>
#include <armadillo>
#include <map>

using namespace std;
using namespace arma;


// for ordering index map
struct StateComparator{
	bool operator()(const irowvec& left, const irowvec& right) const{
		if(left(0) > right(0)) return true;
		if(left(0) < right(0)) return false;
		if(left(0) == right(0)){
			if(left(1) > right(1)) return true;
			if(left(1) < right(1)) return false;
		}
	}
};


class CIMSRG{

public:

	CIMSRG(int dim1B, double d, double g, double smax, double ds, ivec holes, ivec parts);
	~CIMSRG(){}

	int dim1B_, dim2B_;
	int nholes_, nparts_;

	double d_, g_, E_;
	double smax_, ds_;

	ivec holes_, parts_;
	ivec basis1B_, occ1B_;

	imat basis2B_, ph_basis2B_;
	imat occ2B_1_, occ2B_2_, occ2B_3_;

	mat H1B_, H2B_, f_, Gamma_, eta;


	map<irowvec,int,StateComparator> index2B_;
	map<irowvec,int,StateComparator> ph_index2B_;


	void build_basis2B();
	void build_ph_basis2B();
	void build_occ1B();
	void build_occ2B_1();		// n_a - n_b
	void build_occ2B_2();		// 1 - n_a - n_b
	void build_occ2B_3();		// n_a * n_b
	void build_hamiltonian();	// pairing model hamiltonian

	void normal_order();

	double fod_norm();
	double Gammaod_norm();

	void ph_transform2B(); 
	void inverse_ph_transform2B();

	void calc_eta_imtime();
	void calc_eta_white();
	void calc_eta_wegner();

	void RK4_imsrg();

	void imsrg(vec snapshots, string filename);

	mat imsrg();
};

#endif
