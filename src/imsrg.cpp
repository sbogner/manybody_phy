CIMSRG::CIMSRG(double smax, double ds, vec holes, vec parts){

	smax_ = smax;
	ds_ = ds;
	holes_ = holes;
	parts_ = parts;
	nholes_ = holes.n_elem;
	nparts_ = parts.n_elem;
	generate_basis();
}

void CIMSRG::generate_basis(){

	nbasis_ = pow(nholes_+nparts_,2.0);
	basis_.set_size(nbasis_,2);
	
	int index = 0;
	
	for(int i = 0; i < nholes_; ++i){
		for(int j = 0; j < nholes_; ++j){
			basis_(index,0) = holes_(i);
			basis_(index,1) = holes_(j);
			index++;
		}
	}

	for(int i = 0; i < nholes_; ++i){
		for(int a = 0; a < nparts_; ++a){
			basis_(index,0) = holes_(i);
			basis_(index,1) = parts_(a);
			index++;
		}
	}

	for(int a = 0; a < nparts_; ++a){
		for(int i = 0; i < nholes_; ++i){
			basis_(index,0) = parts_(a);
			basis_(index,1) = holes_(i);
			index++;			
		}
	}

	for(int a = 0; a < nparts_; ++a){
		for(int b = 0; b < nparts_; ++b){
			basis_(index,0) = parts_(a);
			basis_(index,1) = parts_(b);
			index++;			
		}
	}
}