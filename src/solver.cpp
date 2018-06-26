#include "solver.h"


CSolver::CSolver(int particles, int sp_states, double d, double gmin, double gmax, string filename){

	particles_ = particles;
	sp_states_ = sp_states;
	d_ = d;
	gmin_ = gmin;
	gmax_ = gmax;
	gstep_ = (gmax-gmin)/100.0;
	filename_ = filename;
}


void CSolver::mbpt2(){

	ofstream outfile;
	outfile.open(filename_+"_mbpt2.dat");
	outfile << "# d = " << d_ << endl;

	for(double g = gmin_; g <= gmax_; g += gstep_){

		cout << "g = " << g << endl;

		CFermionSystem system(particles_, sp_states_, d_, g);

		double Ecorr = 0.0;
		for(int a = system.F_; a < system.N_; ++a){
			for(int b = system.F_; b < system.N_; ++b){
				for(int i = 0; i < system.F_; ++i){
					for(int j = 0; j < system.F_; ++j){
						vector<int> holes = {i,j};
						vector<int> parts = {a,b};
						Ecorr += system.v(a,b,i,j)*system.v(a,b,i,j)/system.eps(holes,parts);
					}
				}
			}
		}

		double exact = -0.25*g*g*(2/(4+g)+1/(6+g)+1/(2+g));
		outfile << g << "\t" << Ecorr << "\t" << exact << endl;
	}

	outfile.close();



}