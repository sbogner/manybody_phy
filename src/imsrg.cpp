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

void CIMSRG::build_hamiltonian(){

	int i, j, k, l;
	irowvec state1, state2;

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

					state1 = {i,j};
					state2 = {k,l};
					H2B_(index2B_[state1],index2B_[state2]) = -0.5*g_;

					state2 = {l,k};
					H2B_(index2B_[state1],index2B_[state2]) = 0.5*g_;

					state1 = {j,i};
					H2B_(index2B_[state1],index2B_[state2]) = -0.5*g_;	

					state2 = {k,l};
					H2B_(index2B_[state1],index2B_[state2]) = 0.5*g_;

				}

			}
		}
	}
}

void CIMSRG::normal_order(){

	int index1, index2;
	irowvec state1, state2;

	// zero-body
	E_ = 0.0;
	for(int i = 0; i < nholes_; ++i){
		E_ += H1B_(holes_(i),holes_(i));
		for(int j = 0; j < nholes_; ++j){
			state1 = {i,j};
			index1 = index2B_[state1];
			E_ += 0.5*H2B_(index1,index1);
		}
	}

	// one-body
	f_ = H1B_;
	for(int i = 0; i < dim1B_; ++i){
		for(int j = 0; j < dim2B_; ++j){
			for(int h = 0; h < nholes_; ++h){
				state1 = {i,h};
				state2 = {j,h};
				index1 = index2B_[state1];
				index2 = index2B_[state2];
				f_(i,j) += H2B_(index1,index2);
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
	irowvec state1, state2;

	for(int a = 0; a < nparts_; ++a){
		for(int b = 0; b < nparts_; ++b){
			for(int i = 0; i < nholes_; ++i){
				for(int j = 0; j < nholes_; ++j){
					state1 = {a,b};
					state2 = {i,j};
					index1 = index2B_[state1];
					index2 = index2B_[state2];
					norm += Gamma_(index1,index2)*Gamma_(index1,index2);
					norm += Gamma_(index2,index1)*Gamma_(index2,index1);
				}
			}
		}
	}
}

mat CIMSRG::ph_transform2B(mat matrix2B){

	int i, j, k, l;
	irowvec state1, state2;

	mat ph_matrix2B = zeros<mat>(dim2B_,dim2B_);
	for(int index1 = 0; index1 < dim2_; ++index1){

		i = basis2B_(index1,0);
		j = basis2B_(index1,1);


		for(int index2 = 0; index2 < dim2B_; ++index2){

			k = basis2B_(index2,0);
			l = basis2B_(index2,1);

			state1 = {i,j};
			state2 = {k,l};
			
			ph_matrix2B(p,q) -= matrix2B(index2B_[state1],index2B_[state2]);
		}
	}

	return ph_matrix2B;
}

mat CIMSRG::inverse_ph_transform2B(mat ph_matrix2B){

	int i, j, k, l;
	irowvec state1, state2;

	mat matrix2B = zeros<mat>(dim2B_,dim2B_);
	for(int index1 = 0; index1 < dim2_; ++index1){

		i = basis2B_(index1,0);
		j = basis2B_(index1,1);


		for(int index2 = 0; index2 < dim2B_; ++index2){

			k = basis2B_(index2,0);
			l = basis2B_(index2,1);

			state1 = {i,j};
			state2 = {k,l};
			
			matrix2B(p,q) -= ph_matrix2B(index2B_[state1],index2B_[state2]);
		}
	}

	return matrix2B;	
}

mat CIMSRG::commutator(mat A, mat B){

	return A*B-B*A;
}

void CIMSRG::calc_eta_wegner(){

	int index1, index2, index3;
	irowvec state1, state2, state3;

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

					state1 = {a,b};
					state2 = {i,j};
					index1 = index2B_[state1];
					index2 = index2B_[state2];

					Gammaod(index1,index2) = Gamma_(index1,index2);
					Gammaod(index2,index1) = Gamma_(index2,index1);
				}
			}
		}
	}
	Gammad = Gamma_-Gammaod;

	// ONE-BODY ETA
	// 1B-1B interaction
	eta1B_ = commutator(fd, fod);

	// 1B-2B interaction
	for(int p = 0; p < dim1B_; ++p){
		for(int q = 0; q < dim1B_; ++q){
			for(int i = 0; i < nholes_; ++i){
				for(int a = 0; a < nparts_; ++a){

					state1 = {a,p};
					state2 = {i,q};
					index1 = index2B_[state1];
					index2 = index2B_[state2];

					eta1B_(p,q) += fd(i,a)*Gammaod(index1,index2)-fod(i,a)*Gammad(index1,index2);

					state1 = {i,p};
					state2 = {a,q};
					index1 = index2B_[state1];
					index2 = index2B_[state2];

					eta1B_(p,q) += fod(a,i)*Gammad(index1,index2)-fd(a,i)*Gammaod(index1,index2);
				}
			}
		}
	}

	// 2B-2B interaction
	GammaGamma = Gammad*occ2B_2_*Gammaod;
	for(int p = 0; p < dim1B_; ++p){
		for(int q = 0; q < dim1B_; ++q){
			for(int i = 0; i < nholes_; ++i){

				state1 = {i,p};
				state2 = {i,q};
				index1 = index2B_[state1];
				index2 = index2B_[state2];

				eta1B_(p,q) += 0.5*(GammaGamma(index1,index2)-trans(GammaGamma)(index1,index2));
			}
		}
	}

	GammaGamma = Gammad*occ2B_3_*Gammaod;
	for(int p = 0; p < dim1B_; ++p){
		for(int q = 0; q < dim1B_; ++q){
			for(int r = 0; r < dim1B_; ++r){

				state1 = {r,p};
				state2 = {r,q};
				index1 = index2B_[state1];
				index2 = index2B_[state2];

				eta1B_(p,q) += 0.5*(GammaGamma(index1,index2)-trans(GammaGamma)(index1,index2));
			}
		}
	}

	// TWO-BODY ETA
	// 1B-2B interaction
	eta2B_.zeros(size(Gamma_));
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
						eta2B_(index1,index2) += (fd(p,t)*Gammaod(index3,index2)-fod(p,t)*Gammad(index3,index2));
						
						state3 = {p,t};
						index3 = index2B_[state3];
						eta2B_(index1,index2) += (fd(q,t)*Gammaod(index3,index2)-fod(p,t)*Gammad(index3,index2));					

						state3 = {t,s};
						index3 = index2B_[state3];
						eta2B_(index1,index2) += (fod(t,r)*Gammad(index1,index3)-fd(t,r)*Gammaod(index1,index3));	

						state3 = {r,t};
						index3 = index2B_[state3];
						eta2B_(index1,index2) += (fod(t,s)*Gammad(index1,index3)-fd(t,s)*Gammaod(index1,index3));	
					}
				}
			}
		}
	}

	// 2B-2B interaction
	GammaGamma = Gammad*occ2B_2_*Gammaod;
	eta2B_ += 0.5*(GammaGamma-trans(GammaGamma))

	// transform to particle-hole representation
	ph_Gammad = ph_transform2B(Gammad);
	ph_Gammaod = ph_transform2B(Gammaod);
	ph_GammaGamma = ph_Gammad*occ2B_1_*ph_Gammaod;

	// transform back
	GammaGamma = inverse_ph_transform2B(ph_GammaGamma);

}

void calc_eta_imtime();
void calc_eta_white();


int main(){

	ivec holes = {0,1,2,3};
	ivec parts = {4,5,6,7};
	CIMSRG Solver(8, 1.0, 0.5, 10.0, 0.001, holes, parts);


}