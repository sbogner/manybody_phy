#include "system.h"

CFermionSystem::CFermionSystem(){

	// default: 4-fermion pairing model
	N_ = 8;
	F_ = 4;
	particles_ = 4;
	sp_states_ = 4;
	d_ = 1.0;
	g_ = 0.5;

	generate_sp_basis();
}

CFermionSystem::CFermionSystem(int particles, int sp_states, double d, double g){

	N_ = 2*sp_states;
	F_ = particles;
	particles_ = particles;
	sp_states_ = sp_states;
	d_ = d;
	g_ = g;

	generate_sp_basis();

	/*
	for(int q = 0; q < N_; ++q){
		for(int r = 0; r < N_; ++r){
			cout << q << "\t" << r << "\t" << h0(q,r) << endl;
		}
	}
	*/

	for(int q = 0; q < 8; ++q){
		for(int r = 0; r < 8; ++r){
			for(int s = 0; s < 8; ++s){
				for(int t = 0; t < 8; ++t){

					if(v(q,r,s,t) != 0){
						cout << q << "\t" <<r << "\t" <<s << "\t" <<t <<"\t"<< v(q,r,s,t) << endl;
					}
					
				}
			}
		}
	}
	
	
}

void CFermionSystem::generate_sp_basis(){

	states_.resize(N_);

	for(int i = 0; i < N_; ++i){


		states_[i].resize(2);

		// assign quantum number p
		states_[i][0] = floor(i/2.0);

		// assign spins
		if(i%2 == 0) states_[i][1] = 1;
		else states_[i][1] = -1;

	}

}

double CFermionSystem::h0(int q, int r){

	// return kinetic energy p*d
	if(q == r) return states_[q][0]*d_;
	else return 0.0;
}

double CFermionSystem::v(int q, int r, int s, int t){

	int p_q, p_r, p_s, p_t;
	int spin_q, spin_r, spin_s, spin_t;

	// retrieve quantum numbers 
	p_q = states_[q][0];
	p_r = states_[r][0];
	p_s = states_[s][0];
	p_t = states_[t][0];
	spin_q = states_[q][1];
	spin_r = states_[r][1];
	spin_s = states_[s][1];
	spin_t = states_[t][1];

	// return interaction energy
	if( (p_q == p_r) && (p_s == p_t) ){
		if( (spin_q == spin_r) || (spin_s == spin_t) ) return 0.0;
		if( (spin_q == spin_s) && (spin_r == spin_t) ) return -0.5*g_;
		if( (spin_q == spin_t) && (spin_r == spin_s) ) return 0.5*g_;	
	}
	else return 0.0;

}

// fock operator
double CFermionSystem::f(int q, int r){

	double E = h0(q,r);
	for(int i = 0; i < F_; ++i) E += v(q,i,r,i);

}

// energy of particle in mb state
double CFermionSystem::eps(int q, int r){

	if(q == r){
		double E = h0(q,q);
		for(int i = 0; i < F_; ++i) E += v(q,i,q,i);
		return E;
	}
	else return 0.0;
}

double CFermionSystem::eps(vector<int>& holes, vector<int>& parts){

	double E = 0.0;

	for(int i = 0; i < holes.size(); ++i) E += eps(holes[i],holes[i]);
	for(int a = 0; a < parts.size(); ++a) E -= eps(parts[a],parts[a]);

	return E;
}



