*!plssem_fimix version 0.3.0
*!Written 15Jan2019
*!Written by Sergio Venturini and Mehmet Mehmetoglu
*!The following code is distributed under GNU General Public License version 3 (GPL-3)

version 15.1

mata:

real scalar plssem_fimix_ll(real matrix eta, real matrix xi, real rowvector rho,
	struct plssem_struct_matrix colvector B,
	struct plssem_struct_matrix colvector Gamma,
	struct plssem_struct_matrix colvector Psi)
{
	/* Description:
		 ------------
		 Function that computes the incomplete log-likelihood for the FIMIX-PLS
		 method
	*/
	
	/* Arguments:
		 ----------
		 - eta				--> real matrix containing the endogenous latent variables
		 - xi			 		--> real matrix containing the exogenous latent variables
		 - rho				--> real colvector the mixing proportions
		 - B					--> struct plssem_struct_matrix colvector containing some
											inner model coefficients (endogenous latent variables)
		 - Gamma			--> struct plssem_struct_matrix colvector containing some
											inner model coefficients (exogenous latent variables)
		 - Psi				--> struct plssem_struct_matrix colvector containing some
											inner model coefficients (variances)
	*/
	
	/* Returned value:
		 ---------------
		 - ll					--> real scalar corresponding to the log-likelihood value
	*/

	real scalar i, n, k, K, P, ll, ll_i, f_ik
	real colvector zeta
	
	n = rows(eta)
	K = length(B)
	P = cols(eta)
	ll = 0
	
	for (i = 1; i <= n; i++) {
		ll_i = 0
		for (k = 1; k <= K; k++) {
			zeta = (I(P) - B[k, 1].mat)'*eta[i, .]' - Gamma[k, 1].mat*xi[i, .]'
			f_ik = lnmvnormalden(J(1, P, 0), Psi[k, 1].mat, zeta)
			ll_i = ll_i + rho[k]*exp(f_ik)
		}
		ll = ll + ln(ll_i)
	}
	
	return(ll)
}

real scalar plssem_fimix_ll_em(real matrix eta, real matrix xi,
	real matrix P_ik, real rowvector rho, 
	struct plssem_struct_matrix colvector B,
	struct plssem_struct_matrix colvector Gamma,
	struct plssem_struct_matrix colvector Psi)
{
	/* Description:
		 ------------
		 Function that computes the complete log-likelihood for the FIMIX-PLS method
	*/
	
	/* Arguments:
		 ----------
		 - eta				--> real matrix containing the endogenous latent variables
		 - xi			 		--> real matrix containing the exogenous latent variables
		 - P_ik			 	--> real matrix containing the posterior membership probability
		 - rho				--> real colvector the mixing proportions
		 - B					--> struct plssem_struct_matrix colvector containing some
											inner model coefficients (endogenous latent variables)
		 - Gamma			--> struct plssem_struct_matrix colvector containing some
											inner model coefficients (exogenous latent variables)
		 - Psi				--> struct plssem_struct_matrix colvector containing some
											inner model coefficients (variances)
	*/
	
	/* Returned value:
		 ---------------
		 - ll					--> real scalar corresponding to the log-likelihood value
	*/

	real scalar i, n, k, K, P, ll, f_ik
	real colvector zeta
	
	n = rows(eta)
	P = cols(eta)
	K = length(B)
	ll = 0
	
	for (i = 1; i <= n; i++) {
		for (k = 1; k <= K; k++) {
			zeta = (I(P) - B[k, 1].mat)'*eta[i, .]' - Gamma[k, 1].mat*xi[i, .]'
			f_ik = lnmvnormalden(J(1, P, 0), Psi[k, 1].mat, zeta)
			ll = ll + P_ik[i, k]*(f_ik + ln(rho[k]))
		}
	}
	
	return(ll)
}

real matrix plssem_fimix_estep(real rowvector rho, real matrix eta,
	real matrix xi, struct plssem_struct_matrix colvector B,
	struct plssem_struct_matrix colvector Gamma,
	struct plssem_struct_matrix colvector Psi)
{
	/* Description:
		 ------------
		 Function that computes the posterior mixture proportions in the FIMIX-PLS
		 method
	*/
	
	/* Arguments:
		 ----------
		 - rho				--> real colvector the mixing proportions
		 - eta				--> real matrix containing the endogenous latent variables
		 - xi			 		--> real matrix containing the exogenous latent variables
		 - B					--> struct plssem_struct_matrix colvector containing some
											inner model coefficients (endogenous latent variables)
		 - Gamma			--> struct plssem_struct_matrix colvector containing some
											inner model coefficients (exogenous latent variables)
		 - Psi				--> struct plssem_struct_matrix colvector containing some
											inner model coefficients (variances)
	*/
	
	/* Returned value:
		 ---------------
		 - P_ik				--> real matrix containing the posterior mixture proportions
	*/

	real scalar P_ik, i, n, k, K, P, f_ik, exponential, den
	real colvector zeta, a_, probabilities
	real matrix invert_
	
	n = rows(eta)
	P = cols(eta)
	K = length(B)
	P_ik = J(n, K, .)
	
	for (i = 1; i <= n; i++) {
		for (k = 1; k <= K; k++) {
			zeta = (I(P) - B[k, 1].mat)'*eta[i, .]' - Gamma[k, 1].mat*xi[i, .]'
			f_ik = lnmvnormalden(J(1, P, 0), Psi[k, 1].mat, zeta)
			P_ik[i, k] = rho[k]*exp(f_ik)

			a_ = B[k, 1].mat'*eta[i, .]' + Gamma[k, 1].mat*xi[i, .]'
			invert_ = invsym(Psi[k, 1].mat)
			exponential = exp(-.5 * a_' * invert_ * a_)
			den = ((2*pi())^(P/2))*sqrt(det(Psi[k, 1].mat))
			probabilities = exponential/den
// 			P_ik[i, k] = rho[k]*probabilities
		}
	}
	//P_ik = P_ik :/ rowsum(P_ik)
	P_ik = J(n, K, 1) :/ (J(n, K, 1) - (P_ik :- rowsum(P_ik)) :/ P_ik)
	
	return(P_ik)
}

real matrix plssem_fimix_ic(real scalar lnL, real scalar lnL_1, real scalar N,
	real scalar Q, real scalar R, real scalar K, real matrix P_ik)
{
	/* Description:
		 ------------
		 Function that computes the information criteria for a FIMIX-PLS analysis
	*/
	
	/* Arguments:
		 ----------
		 - lnL				--> real scalar providing the model's loglikelihood value
		 - lnL_1			--> real scalar providing the model's loglikelihood value
											when K = 1
		 - N					--> real scalar providing the sample size
		 - Q			 		--> real scalar providing the number of endogenous latent
											variables
		 - R					--> real scalar providing the number of coefficients in the
											inner model
		 - K					--> real scalar providing the number of classes/segments
		 - P_ik				--> real matrix containing the a posteriori probability of
											observation i belonging to segment k
	*/
	
	/* Returned value:
		 ---------------
		 - ic					--> real matrix containing the information criteria, that is
											- AIC		Akaike's information criterion
											- AIC3	Modified AIC with factor 3
											- AIC4	Modified AIC with factor 4
											- BIC		Bayesian information criterion
											- CAIC	Consistent AIC
											- HQ		Hannan-Quinn criterion
											- MDL5	Minimum description length with factor 5
											- LnL		Log-likelihood
											- EN		Entropy criterion (normed)
											- NFI		Non-fuzzy index
											- NEC		Normalized entropy criterion
	*/
	
	real scalar NK, ES
	real matrix ic
	
	NK = (K - 1) + K*R + K*Q
	ES = -sum(P_ik :* ln(P_ik))
	
	ic = J(11, 1, .)
	
	ic[1, 1] = -2*lnL + 2*NK
	ic[2, 1] = -2*lnL + 3*NK
	ic[3, 1] = -2*lnL + 4*NK
	ic[4, 1] = -2*lnL + ln(N)*NK
	ic[5, 1] = -2*lnL + (ln(N) + 1)*NK
	ic[6, 1] = -2*lnL + 2*ln(ln(N))*NK
	ic[7, 1] = -2*lnL + 5*ln(N)*NK
	ic[8, 1] = lnL
	ic[9, 1] = 1 - ES/(N*ln(K))
	ic[10, 1] = (K*sum(P_ik :^ 2) - N)/(N*(K - 1))
	ic[11, 1] = ES/(lnL - lnL_1)
	
	return(ic)
}

struct plssem_struct_fimix scalar plssem_fimix(real matrix Y, real matrix S,
	string scalar touse, real scalar K, real scalar maxit_fim, real scalar stop,
	|real scalar restart, real scalar seed, real scalar noisily,
	real scalar groups)
{
	/* Description:
		 ------------
		 Function that implements the FIMIX-PLS method using an EM algorithm
	*/
	
	/* Arguments:
		 ----------
		 - Y					--> real matrix containing the latent variable scores
		 - S			 		--> real matrix containing the structural model adjacency
											matrix
		 - touse	 		--> string scalar containing the name of the variable tracking
											the subset of the data to use
		 - K					--> real scalar providing the number of FIMIX classes
		 - maxit_fim	--> real scalar providing the maximum number of FIMIX
											iterations
		 - stop				--> real scalar providing the FIMIX stopping criterion
		 - restart		--> (optional) real scalar number of repetitions of the EM
											algorithm (default 10)
		 - seed				--> (optional) real scalar providing the seed
		 - noisily		--> (optional) real scalar; if not 0, it prints the EM
											iterations
		 - groups			--> (optional) real scalar; indicator of whether many
											FIMIX-PLS are required
	*/
	
	/* Returned value:
		 ---------------
		 - res_fimix	--> scalar plssem_struct_fimix structure containing the results
	*/

	struct plssem_struct_fimix scalar res_fimix
	struct plssem_struct_matrix scalar mat
	struct plssem_struct_matrix colvector B, Gamma, Psi, B_best, Gamma_best, ///
		Psi_best, path_coef, path_coef_best

	real matrix P_ik, P_ik_best, eta, xi, y, x, xx, xy, yhat, beta, beta_tmp, ///
		B_k, Gamma_k, Psi_k, ic
	real scalar iter, n, N, i, k, r, nendo, nexo, P, p, jj, sse, omega, ///
		rN0, rN_lte_5, delta_Q, old_Q, new_Q, new_ll, iter_best, ///
		lnL_best, lnL_c_best, R, Q, lnL_1, skip
	real colvector touse_vec, Sind, tau, iter_all, lnL_all, ///
		lnL_c_all, fimix_class, class_freq, lnL_iter, lnL_c_iter, lnL_iter_best, ///
		lnL_c_iter_best
	real rowvector rho, endo, exo, rho_best
	string scalar todisp, spaces
	
	if (restart == .) {
		restart = 10
	}
	if (seed != .) {
		rseed(seed)
	}
	if (noisily == .) {
		noisily = 1
	}
	
	touse_vec = st_data(., touse)
	P = cols(Y)
	n = rows(Y)
	N = sum(touse_vec)
	endo = (colsum(S)	:> 0)			 // indicators for the endogenous latent variables
	exo = (colsum(S) :== 0)			 // indicators for the exogenous latent variables
	eta = Y[selectindex(touse_vec), selectindex(endo)]
	xi = Y[selectindex(touse_vec), selectindex(exo)]
	nendo = cols(eta)
	nexo = cols(xi)
	
	B = J(K, 1, mat)
	Gamma = J(K, 1, mat)
	Psi = J(K, 1, mat)
	path_coef = J(K, 1, mat)
	iter_all = J(restart, 1, .)
	lnL_all = J(restart, 1, .)
	lnL_c_all = J(restart, 1, .)
	ic = J(11, 1, .)

	if ((noisily) & (restart > 1)) {
		printf("{txt}\n")
		printf("{txt}FIMIX-PLS --> EM algorithm runs (")
		printf("{res}%f", restart)
		printf("{txt})\n")
		printf("{hline 4}{c +}{hline 3} 1 ")
		printf("{hline 3}{c +}{hline 3} 2 ")
		printf("{hline 3}{c +}{hline 3} 3 ")
		printf("{hline 3}{c +}{hline 3} 4 ")
		printf("{hline 3}{c +}{hline 3} 5\n")
	}
	skip = 0
	for (r = 1; r <= restart; r++) {
		if ((noisily) & (restart > 1)) {
			if (mod(r, 50) == 0) {
				todisp = strtrim(strofreal(r, "%9.0f"))
				if (strlen(todisp) < strlen(strofreal(restart))) {
					skip = strlen(strofreal(restart)) - strlen(todisp)
				}
				spaces = ""
				for (i = 1; i <= (skip + 3); i++) {
					spaces = spaces + char(32)
				}
				printf("{txt}." + spaces + todisp + "\n")
				skip = 0
			}
			else {
				printf("{txt}.")
			}
			displayflush()
		}
		if (groups) {
			printf("{txt}.")
			displayflush()
		}
		
		// Reset likelihood vectors at each EM run
		lnL_iter = J(0, 1, .)
		lnL_c_iter = J(0, 1, .)
		
		// Random initialization of P_ik
		P_ik = rdirichlet(N, J(1, K, 1/K))

		// Initial M-step (compute B, Gamma and Psi)
		rho = mean(P_ik)
		for (k = 1; k <= K; k++) {
			beta = J(P, P, 0)
			B_k = J(nendo, nendo, 0)
			Gamma_k = J(nendo, nexo, 0)
			Psi_k = J(nendo, nendo, 0)
			jj = 1
			for (p = 1; p <= P; p++) {
				Sind = selectindex(S[., p])
				if (endo[1, p]) {
					//x = (J(N, 1, 1), Y[., Sind])
					x = Y[selectindex(touse_vec), Sind]
					y = Y[selectindex(touse_vec), p]
					xx = quadcross(x, P_ik[., k], x)
					xy = quadcross(x, P_ik[., k], y)
					tau = qrsolve(xx, xy)
					//beta[Sind, p] = tau[|2 \ .|]
					beta[Sind, p] = tau
					yhat = x * tau
					sse = quadcross((y - yhat), P_ik[., k], (y - yhat))
					omega = sse/(N*rho[k])
					Psi_k[jj, jj] = omega
					jj++
				}
			}
			B_k = beta[selectindex(endo), selectindex(endo)]
			Gamma_k = beta[selectindex(exo), selectindex(endo)]'
			B[k, 1].mat = B_k
			Gamma[k, 1].mat = Gamma_k
			Psi[k, 1].mat = Psi_k
		}
		
		// EM iterations
		iter = 1
		new_Q = plssem_fimix_ll_em(eta, xi, P_ik, rho, B, Gamma, Psi) //mindouble()/2
		delta_Q = maxdouble()/2
		while ((iter <= maxit_fim) & (delta_Q >= stop)) {
			/*
			if (noisily) {
				if (iter == 1) {
					printf("{txt}FIMIX iteration 1")
				}
				else if (mod(iter, 5) == 0) {
					todisp = strtrim(strofreal(iter, "%9.0f"))
					printf("{txt}" + todisp)
				}
				else {
					printf("{txt}.")
				}
				displayflush()
			}
			*/
			
			// E-step
			P_ik = plssem_fimix_estep(rho, eta, xi, B, Gamma, Psi)
			_editmissing(P_ik, 0)

			old_Q = new_Q
			
			// M-step (compute B, Gamma and Psi)
			rho = mean(P_ik)
			for (k = 1; k <= K; k++) {
				beta = J(P, P, 0)
				B_k = J(nendo, nendo, 0)
				Gamma_k = J(nendo, nexo, 0)
				Psi_k = J(nendo, nendo, 0)
				jj = 1
				for (p = 1; p <= P; p++) {
					Sind = selectindex(S[., p])
					if (endo[1, p]) {
						//x = (J(N, 1, 1), Y[., Sind])
						x = Y[selectindex(touse_vec), Sind]
						y = Y[selectindex(touse_vec), p]
						xx = quadcross(x, P_ik[., k], x)
						xy = quadcross(x, P_ik[., k], y)
						tau = qrsolve(xx, xy)
						//beta[Sind, p] = tau[|2 \ .|]
						beta[Sind, p] = tau
						yhat = x * tau
						sse = quadcross((y - yhat), P_ik[., k], (y - yhat))
						omega = sse/(N*rho[k])
						Psi_k[jj, jj] = omega
						jj++
					}
				}
				B_k = beta[selectindex(endo), selectindex(endo)]
				Gamma_k = beta[selectindex(exo), selectindex(endo)]'
				B[k, 1].mat = B_k
				Gamma[k, 1].mat = Gamma_k
				Psi[k, 1].mat = Psi_k
				beta_tmp = select(beta, rowsum(beta))
				beta_tmp = select(beta_tmp, colsum(beta_tmp))
				path_coef[k, 1].mat = beta_tmp
				_editvalue(path_coef[k, 1].mat, 0, .)
			}
			
			// Compute the complete and incomplete loglikelihood
			new_Q = plssem_fimix_ll_em(eta, xi, P_ik, rho, B, Gamma, Psi)
			delta_Q = abs(new_Q - old_Q)
			lnL_c_iter = (lnL_c_iter \ new_Q)
			new_ll = plssem_fimix_ll(eta, xi, rho, B, Gamma, Psi)
			lnL_iter = (lnL_iter \ new_ll)
			
			iter++
		}
		iter_all[r] = (iter - 1)
		lnL_all[r] = new_ll
		lnL_c_all[r] = new_Q
		
		if (r == 1) {
			iter_best = (iter - 1)
			lnL_best = new_ll
			lnL_c_best = new_Q
			B_best = B
			Gamma_best = Gamma
			path_coef_best = path_coef
			Psi_best = Psi
			rho_best = rho
			P_ik_best = P_ik
			lnL_iter_best = lnL_iter
			lnL_c_iter_best = lnL_c_iter
		}
		else {
			if (lnL_all[r] > lnL_best) {
				iter_best = (iter - 1)
				lnL_best = new_ll
				lnL_c_best = new_Q
				B_best = B
				Gamma_best = Gamma
				path_coef_best = path_coef
				Psi_best = Psi
				rho_best = rho
				P_ik_best = P_ik
				lnL_iter_best = lnL_iter
				lnL_c_iter_best = lnL_c_iter
			}
		}
	}
	if ((noisily) & (restart > 1)) {
		printf("{txt}\n")
		displayflush()
	}
	if (groups) {
		printf("{txt}\n")
		displayflush()
	}
	
	// Final classification
	fimix_class = which(P_ik_best, "max")
	class_freq = J(K, 1, 0)
	for (k = 1; k <= K; k++) {
		class_freq[k, 1] = sum(fimix_class :== k)
	}
	if (anyof(class_freq, 0)) {
		rN0 = 1
	}
	else {
		rN0 = 0
	}
	if (any(class_freq :<= 5) & all(class_freq :> 0)) {
		rN_lte_5 = 1
	}
	else {
		rN_lte_5 = 0
	}
	
	// Information criteria and fit indices
	B = J(1, 1, mat)
	Gamma = J(1, 1, mat)
	Psi = J(1, 1, mat)
	beta = J(P, P, 0)
	Psi_k = J(nendo, nendo, 0)
	jj = 1
	for (p = 1; p <= P; p++) {
		Sind = selectindex(S[., p])
		if (endo[1, p]) {
			//x = (J(N, 1, 1), Y[., Sind])
			x = Y[selectindex(touse_vec), Sind]
			y = Y[selectindex(touse_vec), p]
			xx = quadcross(x, x)
			xy = quadcross(x, y)
			tau = qrsolve(xx, xy)
			//beta[Sind, p] = tau[|2 \ .|]
			beta[Sind, p] = tau
			yhat = x * tau
			sse = quadcross((y - yhat), (y - yhat))
			omega = sse/N
			Psi_k[jj, jj] = omega
			jj++
		}
	}
	B[1, 1].mat = beta[selectindex(endo), selectindex(endo)]
	Gamma[1, 1].mat = beta[selectindex(exo), selectindex(endo)]'
	Psi[1, 1].mat = Psi_k
	lnL_1 = plssem_fimix_ll(eta, xi, J(1, 1, 1), B, Gamma, Psi)

	Q = length(selectindex(endo))
	R = sum(S[., selectindex(endo)])
	//R = R + Q  // to use if intercepts are included in the inner model regressions
	ic = plssem_fimix_ic(lnL_best, lnL_1, N, Q, R, K, P_ik_best)
	
	// Assigning results to structure's members
	res_fimix.niter = iter_best
	res_fimix.fimix_class = fimix_class
	res_fimix.touse_vec = touse_vec
	res_fimix.rN0 = rN0
	res_fimix.rN_lte_5 = rN_lte_5
	res_fimix.nclass = K
	res_fimix.maxiter = maxit_fim
	res_fimix.stop = stop
	res_fimix.restart = restart
	res_fimix.rho = rho_best
	res_fimix.path = path_coef_best
	res_fimix.classfreq = class_freq
	res_fimix.ll = lnL_iter_best
	res_fimix.ll_c = lnL_c_iter_best
	res_fimix.ll_restart = lnL_all
	res_fimix.ic = ic
	
	return(res_fimix)
}

end
