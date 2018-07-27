// In-medium SRG for the pairing model
#include "imsrg.h"

CIMSRG::CIMSRG(int dim1B, double d, double g, double smax, double ds, ivec holes, ivec parts){

	d_ = d;
	g_ = g;
	smax_ = smax;
	ds_ = ds;

	holes_ = holes;
	parts_ = parts;
	nholes_ = holes.n_elem;
	nparts_ = parts.n_elem;

	dim1B_ = dim1B;
	dim2B_ = dim1B*dim1B;
	basis1B_.set_size(dim1B);
	for(int i = 0; i < dim1B; ++i) basis1B_(i) = i;

	build_basis2B();
	build_ph_basis2B();
	build_occ1B();
	build_occ2B_1();
	build_ph_occ2B_1();
	build_occ2B_2();
	build_occ2B_3();
	build_hamiltonian();
	normal_order();
}

void CIMSRG::build_basis2B(){

	basis2B_.set_size(dim2B_,2);
	irowvec state;
	int index = 0;

	for(int i = 0; i < nholes_; ++i){
		for(int j = 0; j < nholes_; ++j){
			basis2B_(index,0) = holes_(i);
			basis2B_(index,1) = holes_(j);
			state = basis2B_.row(index);
			index2B_[state] = index;
			index++;
		}
	}

	for(int i = 0; i < nholes_; ++i){
		for(int a = 0; a < nparts_; ++a){
			basis2B_(index,0) = holes_(i);
			basis2B_(index,1) = parts_(a);
			state = basis2B_.row(index);
			index2B_[state] = index;
			index++;
		}
	}

	for(int a = 0; a < nparts_; ++a){
		for(int i = 0; i < nholes_; ++i){
			basis2B_(index,0) = parts_(a);
			basis2B_(index,1) = holes_(i);
			state = basis2B_.row(index);
			index2B_[state] = index;
			index++;			
		}
	}

	for(int a = 0; a < nparts_; ++a){
		for(int b = 0; b < nparts_; ++b){
			basis2B_(index,0) = parts_(a);
			basis2B_(index,1) = parts_(b);
			state = basis2B_.row(index);
			index2B_[state] = index;
			index++;			
		}
	}
}

void CIMSRG::build_ph_basis2B(){

	ph_basis2B_.set_size(dim2B_,2);
	irowvec state;
	int index = 0;

	for(int i = 0; i < nholes_; ++i){
		for(int j = 0; j < nholes_; ++j){
			ph_basis2B_(index,0) = holes_(i);
			ph_basis2B_(index,1) = holes_(j);
			state = ph_basis2B_.row(index);
			ph_index2B_[state] = index;
			index++;
		}
	}

	for(int i = 0; i < nholes_; ++i){
		for(int a = 0; a < nparts_; ++a){
			ph_basis2B_(index,0) = holes_(i);
			ph_basis2B_(index,1) = parts_(a);
			state = ph_basis2B_.row(index);
			ph_index2B_[state] = index;
			index++;
		}
	}

	for(int a = 0; a < nparts_; ++a){
		for(int i = 0; i < nholes_; ++i){
			ph_basis2B_(index,0) = parts_(a);
			ph_basis2B_(index,1) = holes_(i);
			state = ph_basis2B_.row(index);
			ph_index2B_[state] = index;
			index++;			
		}
	}

	for(int a = 0; a < nparts_; ++a){
		for(int b = 0; b < nparts_; ++b){
			ph_basis2B_(index,0) = parts_(a);
			ph_basis2B_(index,1) = parts_(b);
			state = ph_basis2B_.row(index);
			ph_index2B_[state] = index;
			index++;			
		}
	}
}

void CIMSRG::build_occ1B(){

	occ1B_.zeros(dim1B_);
	for(int i = 0; i < nholes_; ++i) occ1B_(holes_(i)) = 1;
}

// n_a - n_b
void CIMSRG::build_occ2B_1(){

	occ2B_1_.zeros(dim2B_,dim2B_);
	for(int p = 0; p < dim2B_; ++p){
		occ2B_1_(p,p) = occ1B_(basis2B_(p,0))-occ1B_(basis2B_(p,1));
	}
}

// 1 - n_a - n_b
void CIMSRG::build_occ2B_2(){

	occ2B_2_.zeros(dim2B_,dim2B_);
	for(int p = 0; p < dim2B_; ++p){
		occ2B_2_(p,p) = 1-occ1B_(basis2B_(p,0))-occ1B_(basis2B_(p,1));
	}	
}

// n_a * n_b
void CIMSRG::build_occ2B_3(){

	occ2B_3_.zeros(dim2B_,dim2B_);
	for(int p = 0; p < dim2B_; ++p){
		occ2B_3_(p,p) = occ1B_(basis2B_(p,0))*occ1B_(basis2B_(p,1));
	}	
}

void CIMSRG::build_ph_occ2B_1(){

	ph_occ2B_1_.zeros(dim2B_,dim2B_);
	for(int p = 0; p < dim2B_; ++p){
		ph_occ2B_1_(p,p) = occ1B_(ph_basis2B_(p,0))-occ1B_(ph_basis2B_(p,1));
	}	
}

void CIMSRG::build_hamiltonian(){

	int i, j, k, l;

	// one-body
	H1B_.zeros(dim1B_,dim1B_);
	for(int index = 0; index < dim1B_; ++index){
		i = basis1B_(index);
		H1B_(i,i) = d_*floor(i/2.0);
	}

	// two-body
	H2B_.zeros(dim2B_,dim2B_);
	for(int index1 = 0; index1 < dim2B_; ++index1){

		i = basis2B_(index1,0);
		j = basis2B_(index1,1);

		if((i%2 == 0) && (j == i+1)){

			for(int index2 = 0; index2 < dim2B_; ++index2){

				k = basis2B_(index2,0);
				l = basis2B_(index2,1);

				if((k%2 == 0) && (l == k+1)){

					H2B_(index2B_[{i,j}],index2B_[{k,l}]) = -0.5*g_;
					H2B_(index2B_[{i,j}],index2B_[{l,k}]) = 0.5*g_;
					H2B_(index2B_[{j,i}],index2B_[{k,l}]) = 0.5*g_;
					H2B_(index2B_[{j,i}],index2B_[{l,k}]) = -0.5*g_;	
				}
			}
		}
	}
}

void CIMSRG::normal_order(){

	int index1, index2;

	// zero-body
	E_ = 0.0;
	for(int i = 0; i < nholes_; ++i){
		E_ += H1B_(holes_(i),holes_(i));
		for(int j = 0; j < nholes_; ++j){
			index1 = index2B_[{holes_(i),holes_(j)}];
			E_ += 0.5*H2B_(index1,index1);
		}
	}

	// one-body
	f_ = H1B_;
	for(int p = 0; p < dim1B_; ++p){
		for(int q = 0; q < dim1B_; ++q){
			for(int i = 0; i < nholes_; ++i){
				index1 = index2B_[{p,holes_(i)}];
				index2 = index2B_[{q,holes_(i)}];
				f_(p,q) += H2B_(index1,index2);
			}
		}
	}

	// two-body
	Gamma_ = H2B_;
}

double CIMSRG::fod_norm(){

	double norm = 0.0;

	for(int a = 0; a < nparts_; ++a){
		for(int i = 0; i < nholes_; ++i){
			norm += f_(parts_(a),holes_(i))*f_(parts_(a),holes_(i));
			norm += f_(holes_(i),parts_(a))*f_(holes_(i),parts_(a));
		}
	}

	return sqrt(norm);
}

double CIMSRG::Gammaod_norm(){

	double norm = 0.0;
	int index1, index2;

	for(int a = 0; a < nparts_; ++a){
		for(int b = 0; b < nparts_; ++b){
			for(int i = 0; i < nholes_; ++i){
				for(int j = 0; j < nholes_; ++j){
					index1 = index2B_[{parts_(a),parts_(b)}];
					index2 = index2B_[{holes_(i),holes_(j)}];
					norm += Gamma_(index1,index2)*Gamma_(index1,index2);
					norm += Gamma_(index2,index1)*Gamma_(index2,index1);
				}
			}
		}
	}

	return sqrt(norm);
}

mat CIMSRG::ph_transform2B(mat matrix2B){

	int i, j, k, l;
	mat ph_matrix2B = zeros<mat>(dim2B_,dim2B_);

	for(int p = 0; p < dim2B_; ++p){

		i = basis2B_(p,0);
		j = basis2B_(p,1);

		for(int q = 0; q < dim2B_; ++q){

			k = basis2B_(q,0);
			l = basis2B_(q,1);
			
			ph_matrix2B(p,q) -= matrix2B(index2B_[{i,j}],index2B_[{k,l}]);
		}
	}

	return ph_matrix2B;
}

mat CIMSRG::inverse_ph_transform2B(mat ph_matrix2B){

	int i, j, k, l;
	mat matrix2B = zeros<mat>(dim2B_,dim2B_);

	for(int p = 0; p < dim2B_; ++p){

		i = basis2B_(p,0);
		j = basis2B_(p,1);

		for(int q = 0; q < dim2B_; ++q){

			k = basis2B_(q,0);
			l = basis2B_(q,1);
			
			matrix2B(p,q) -= ph_matrix2B(ph_index2B_[{i,j}],ph_index2B_[{k,l}]);
		}
	}

	return matrix2B;	
}

mat CIMSRG::commutator(mat A, mat B){

	return A*B-B*A;
}

void CIMSRG::calc_eta_wegner(){

	int index1, index2, index3;

	// split into diagonal and off-diagonal parts
	mat fd, fod, Gammad, Gammaod, GammaGamma;
	fd.zeros(size(f_));
	fod.zeros(size(f_));
	Gammad.zeros(size(Gamma_));
	Gammaod.zeros(size(Gamma_));

	for(int a = 0; a < nparts_; ++a){
		for(int i = 0; i < nholes_; ++i){
			fod(parts_(a),holes_(i)) = f_(parts_(a),holes_(i));
			fod(holes_(i),parts_(a)) = f_(holes_(i),parts_(a));
		}
	}
	fd = f_-fod;

	for(int a = 0; a < nparts_; ++a){
		for(int b = 0; b < nparts_; ++b){
			for(int i = 0; i < nholes_; ++i){
				for(int j = 0; j < nholes_; ++j){

					index1 = index2B_[{parts_(a),parts_(b)}];
					index2 = index2B_[{holes_(i),holes_(j)}];

					Gammaod(index1,index2) = Gamma_(index1,index2);
					Gammaod(index2,index1) = Gamma_(index2,index1);
				}
			}
		}
	}
	Gammad = Gamma_-Gammaod;

	// 1B-1B interaction
	eta1B_ = commutator(fd, fod);

	// 1B-2B interaction
	for(int p = 0; p < dim1B_; ++p){
		for(int q = 0; q < dim1B_; ++q){
			for(int i = 0; i < nholes_; ++i){
				for(int a = 0; a < nparts_; ++a){

					index1 = index2B_[{parts_(a),p}];
					index2 = index2B_[{holes_(i),q}];

					eta1B_(p,q) += fd(holes_(i),parts_(a))*Gammaod(index1,index2)-fod(holes_(i),parts_(a))*Gammad(index1,index2);

					index1 = index2B_[{holes_(i),p}];
					index2 = index2B_[{parts_(a),q}];

					eta1B_(p,q) += fod(parts_(a),holes_(i))*Gammad(index1,index2)-fd(parts_(a),holes_(i))*Gammaod(index1,index2);
				}
			}
		}
	}

	// 2B-2B interaction
	GammaGamma = Gammad*occ2B_2_*Gammaod;
	mat GammaGammaT = GammaGamma.t();
	for(int p = 0; p < dim1B_; ++p){
		for(int q = 0; q < dim1B_; ++q){
			for(int i = 0; i < nholes_; ++i){

				index1 = index2B_[{holes_(i),p}];
				index2 = index2B_[{holes_(i),q}];

				eta1B_(p,q) += 0.5*(GammaGamma(index1,index2)-GammaGammaT(index1,index2));
			}
		}
	}

	GammaGamma = Gammad*occ2B_3_*Gammaod;
	GammaGammaT = GammaGamma.t();
	for(int p = 0; p < dim1B_; ++p){
		for(int q = 0; q < dim1B_; ++q){
			for(int r = 0; r < dim1B_; ++r){

				index1 = index2B_[{r,p}];
				index2 = index2B_[{r,q}];

				eta1B_(p,q) += 0.5*(GammaGamma(index1,index2)-GammaGammaT(index1,index2));
			}
		}
	}

	// 1B-2B interaction
	eta2B_.zeros(size(Gamma_));
	for(int p = 0; p < dim1B_; ++p){
		for(int q = 0; q < dim1B_; ++q){
			for(int r = 0; r < dim1B_; ++r){
				for(int s = 0; s < dim1B_; ++s){

					index1 = index2B_[{p,q}];
					index2 = index2B_[{r,s}];

					for(int t = 0; t < dim1B_; ++t){

						index3 = index2B_[{t,q}];
						eta2B_(index1,index2) += (fd(p,t)*Gammaod(index3,index2)-fod(p,t)*Gammad(index3,index2));
						
						index3 = index2B_[{p,t}];
						eta2B_(index1,index2) += (fd(q,t)*Gammaod(index3,index2)-fod(p,t)*Gammad(index3,index2));					

						index3 = index2B_[{t,s}];
						eta2B_(index1,index2) += (fod(t,r)*Gammad(index1,index3)-fd(t,r)*Gammaod(index1,index3));	

						index3 = index2B_[{r,t}];
						eta2B_(index1,index2) += (fod(t,s)*Gammad(index1,index3)-fd(t,s)*Gammaod(index1,index3));	
					}
				}
			}
		}
	}

	// 2B-2B interaction
	GammaGamma = Gammad*occ2B_2_*Gammaod;
	GammaGammaT = GammaGamma.t();
	eta2B_ += 0.5*(GammaGamma-GammaGammaT);

	// transform to particle-hole representation
	mat ph_Gammad = ph_transform2B(Gammad);
	mat ph_Gammaod = ph_transform2B(Gammaod);
	mat ph_GammaGamma = ph_Gammad*ph_occ2B_1_*ph_Gammaod;

	// transform back
	GammaGamma = inverse_ph_transform2B(ph_GammaGamma);

	// commutator/antisymmetrization
	int i, j, k, l;
	mat work = zeros<mat>(size(GammaGamma));
	for(int p = 0; p < dim2B_; ++p){

		i = basis2B_(p,0);
		j = basis2B_(p,1);

		for(int q = 0; q < dim2B_; ++q){

			k = basis2B_(q,0);
			l = basis2B_(q,1);

			work(p,q) += (GammaGamma(index2B_[{j,i}],q)+GammaGamma(p,index2B_[{l,k}])
				          -GammaGamma(p,q)-GammaGamma(index2B_[{j,i}],index2B_[{l,k}]));
		}
	}
	GammaGamma = work;
	eta2B_ += GammaGamma;
}

void CIMSRG::calc_derivatives(){

	// calculate generator
	calc_eta_wegner();

	// zero-body
	dE_ = 0.0;
	int k = 0;
	for(int i = 0; i < nholes_; ++i){
		for(int a = 0; a < nparts_; ++a){
			dE_ += (eta1B_(holes_(i),parts_(a))*f_(parts_(a),holes_(i))
			       - eta1B_(parts_(a),holes_(i))*f_(holes_(i),parts_(a)));
		}
	}

	int index1, index2, index3;
	for(int i = 0; i < nholes_; ++i){
		for(int j = 0; j < nholes_; ++j){
			for(int a = 0; a < nparts_; ++a){
				for(int b = 0; b <nparts_; ++b){

					index1 = index2B_[{holes_(i),holes_(j)}];
					index2 = index2B_[{parts_(a),parts_(b)}];

					dE_ += 0.5*eta2B_(index1,index2)*Gamma_(index2,index1);
				}
			}
		}
	}

	// one-body
	df_ = commutator(eta1B_,f_);
	for(int p = 0; p < dim1B_; ++p){
		for(int q = 0; q < dim1B_; ++q){
			for(int i = 0; i < nholes_; ++i){
				for(int a = 0; a <nparts_; ++a){

					index1 = index2B_[{parts_(a),p}];
					index2 = index2B_[{holes_(i),q}];

					df_(p,q) += (eta1B_(holes_(i),parts_(a))*Gamma_(index1,index2)-f_(holes_(i),parts_(a))*eta2B_(index1,index2));

					index1 = index2B_[{holes_(i),p}];
					index2 = index2B_[{parts_(a),q}];

					df_(p,q) += (f_(parts_(a),holes_(i))*eta2B_(index1,index2)-eta1B_(parts_(a),holes_(i))*Gamma_(index1,index2));
				}
			}
		}
	}

	mat etaGamma = eta2B_*occ2B_2_*Gamma_;
	mat etaGammaT = etaGamma.t();
	for(int p = 0; p < dim1B_; ++p){
		for(int q = 0; q < dim1B_; ++q){
			for(int r = 0; r < dim1B_; ++r){

				state1 = {r,p};
				state2 = {r,q};
				index1 = index2B_[state1];
				index2 = index2B_[state2];

				df_(p,q) += 0.5*(etaGamma(index1,index2)+etaGammaT(index1,index2));
			}
		}
	}

	// two-body
	dGamma_ = zeros<mat>(size(Gamma_));
	for(int p = 0; p < dim1B_; ++p){
		for(int q = 0; q < dim1B_; ++q){
			for(int r = 0; r < dim1B_; ++r){
				for(int s = 0; s < dim1B_; ++s){

					state1 = {p,q};
					state2 = {r,s};
					index1 = index2B_[state1];
					index2 = index2B_[state2];

					for(int t = 0; t < dim1B_; ++t){

						state3 = {t,q};
						index3 = index2B_[state3];
						dGamma_(index1,index2) += (eta1B_(p,t)*Gamma_(index3,index2)-f_(p,t)*eta2B_(index3,index2));

						state3 = {p,t};
						index3 = index2B_[state3];
						dGamma_(index1,index2) += (eta1B_(q,t)*Gamma_(index3,index2)-f_(q,t)*eta2B_(index3,index2));

						state3 = {t,s};
						index3 = index2B_[state3];
						dGamma_(index1,index2) += (f_(t,r)*eta2B_(index1,index3)-eta1B_(t,r)*Gamma_(index1,index3));					

						state3 = {r,t};
						index3 = index2B_[state3];
						dGamma_(index1,index2) += (f_(t,s)*eta2B_(index1,index3)-eta1B_(t,s)*Gamma_(index1,index3));
					}
				}
			}
		}
	}

	etaGamma = eta2B_*occ2B_2_*Gamma_;
	dGamma_ += 0.5*(etaGamma+trans(etaGamma));

	// transform matrices to ph representation
	mat ph_eta2B = ph_transform2B(eta2B_);
	mat ph_Gamma = ph_transform2B(Gamma_);
	mat ph_etaGamma = ph_eta2B*ph_occ2B_1_*ph_Gamma;

	// transform back
	etaGamma = inverse_ph_transform2B(ph_etaGamma);

	// antisymmetrization
	int i, j, l;
	mat work = zeros<mat>(size(etaGamma));
	for(int index1 = 0; index1 < dim2B_; ++index1){

		i = basis2B_(index1,0);
		j = basis2B_(index1,1);
		irowvec state1 = {j,i};

		for(int index2 = 0; index2 < dim2B_; ++index2){

			k = basis2B_(index2,0);
			l = basis2B_(index2,1);
			irowvec state2 = {l,k};

			work(index1,index2) += (etaGamma(index2B_[state1],index2) + etaGamma(index1,index2B_[state2])
				                    - etaGamma(index1,index2) - etaGamma(index2B_[state1],index2B_[state2]));
		}
	}
	etaGamma = work;
	dGamma_ += etaGamma;
}

/*
void CIMSRG::RK4(){

	double E_k1, E_k2, E_k3, E_k4;
	mat f_k1, f_k2, f_k3, f_k4;
	mat Gamma_k1, Gamma_k2, Gamma_k3, Gamma_k4;

	// store old values
	double E = E_;
	mat f = f_;
	mat Gamma = Gamma_;

	// K1
	calc_derivatives();
	E_k1 = ds_*dE_;
	f_k1 = ds_*df_;
	Gamma_k1 = ds_*dGamma_;

	// K2
	E_ = E+0.5*E_k1;
	f_ = f+0.5*f_k1;
	Gamma_ = Gamma+0.5*Gamma_k1;
	calc_derivatives();
	E_k2 = ds_*dE_;
	f_k2 = ds_*df_;
	Gamma_k2 = ds_*dGamma_;
	
	// K3
	E_ = E+0.5*E_k2;
	f_ = f+0.5*f_k2;
	Gamma_ = Gamma+0.5*Gamma_k2;
	calc_derivatives();
	E_k3 = ds_*dE_;
	f_k3 = ds_*df_;
	Gamma_k3 = ds_*dGamma_;

	// K4
	E_ = E+E_k3;
	f_ = f+f_k3;
	Gamma_ = Gamma+Gamma_k3;
	calc_derivatives();
	E_k4 = ds_*dE_;
	f_k4 = ds_*df_;
	Gamma_k4 = ds_*dGamma_;

	// update
	E_ = E+(E_k1+2.0*E_k2+2.0*E_k3+E_k4)/6.0;
	f_ = f+(f_k1+2.0*f_k2+2.0*f_k3+f_k4)/6.0;
	Gamma_ = Gamma+(Gamma_k1+2.0*Gamma_k2+2.0*Gamma_k3+Gamma_k4)/6.0;
}
*/

void CIMSRG::euler(){

	calc_derivatives();
	E_ += ds_*dE_;
	f_ += ds_*df_;
	Gamma_ += ds_*dGamma_;
}

void CIMSRG::imsrg(vec snapshots, string filename){

	int k = 0;

	ofstream outfile;
	outfile.open(filename);
	outfile << "# flow parameter s, zero-body component E" << endl;

	for(double s = 0; s <= smax_; s += ds_){

		printf("s = %4.3f   dE = %4.3f   E = %4.3f\n", s, dE_, E_);

		//RK4();
		euler();

		// write snapshots to file
		if(fabs(s-snapshots(k)) < 0.5*ds_){
		outfile << s << "\t" << E_ << endl;
			++k;
		}
	}
}

double CIMSRG::imsrg(){

	for(double s = 0; s <= smax_; s += ds_) euler();
	return E_;
}

int main(int argc, char *argv[]){

	double d = atof(argv[1]);
	double g = atof(argv[2]);
	double smax = atof(argv[3]);
	double ds = atof(argv[4]);

	ivec holes = {0,1,2,3};
	ivec parts = {4,5,6,7};
	CIMSRG PairingModel(8, d, g, smax, ds, holes, parts);

	vec snapshots = {0.001, 0.002, 0.003, 0.005, 0.007, 0.01, 0.02, 0.03, 0.05, 0.07, 0.1, 0.2, 0.3, 0.5, 0.7, 1.0, 2.0, 3.0, 5.0, 7.0, 10.0};
	PairingModel.imsrg(snapshots, "pairing_imsrg.dat");

	return 0;
}