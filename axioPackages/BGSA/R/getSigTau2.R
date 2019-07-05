getSigTau2 <-
function(tau2, nu, sc, s0, muSc, sigSc){

	w = 1
	m = 100
	
	z = getPostSigTau2(tau2, nu, sc, s0, muSc, sigSc) - rexp(1)

	# Stepping out to obtain the [L, R] range
	u = runif(1)
	L = sc - w*u
	R = L + w
	v = runif(1)
	J = floor(m*v)
	K = (m-1) - J
	
	L = max(0, L)
	while (J>0 && L>0 && z < getPostSigTau2(tau2, nu, L, s0, muSc, sigSc)) {
		L = L - w	
		L = max(0, L)		
		J = J - 1
	}

	while (K>0 && z < getPostSigTau2(tau2, nu, R, s0, muSc, sigSc)) {
		R = R+w
		K = K-1
	}


	# Shrinkage to obtain a sample
	u = runif(1)
	newParam = L + u*(R-L)
	
	while (z > getPostSigTau2(tau2, nu, newParam,  s0, muSc, sigSc)) {
		if (newParam < sc) {
			L = newParam
		}else{
			R = newParam
		}
    
		u = runif(1)
		newParam = L + u*(R-L)
	}

	return(newParam)
	}
