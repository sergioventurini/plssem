*!plssem_rebus version 0.3.0
*!Written 15Jan2019
*!Written by Sergio Venturini and Mehmet Mehmetoglu
*!The following code is distributed under GNU General Public License version 3 (GPL-3)

version 15.1

mata:

real colvector plssem_rebus_cm(real matrix e, real matrix f, real matrix loads, ///
	real rowvector r2)
{
	/* Description:
		 ------------
		 Function that computes the REBUS-PLS closeness measure
	*/
	
	/* Arguments:
		 ----------
		 - e 			--> real matrix containing the PLS-SEM measurement model residuals
		 - f			--> real matrix containing the PLS-SEM structural model residuals
		 - loads	--> real matrix containing the PLS-SEM loadings
		 - r2			--> real rowvector containing the PLS-SEM R-squares
	*/
	
	/* Returned value:
		 ---------------
		 - cm			--> real colvector containing the closeness measures for all
									observations
	*/
	
	real scalar N
	real colvector cm, left, right
	real matrix com, e2, f2, r2_mat
	
	N = rows(e)
	
	e2 = e :^ 2
	com = diag(rowsum(loads :^ 2))
	left = rowsum(e2 * invsym(com), .)
	left = left / colsum(left)
	
	f2 = f :^ 2
	r2_mat = diag(r2)
	right = rowsum(f2 * invsym(r2_mat), .)
	right = right / colsum(right)
	
	cm = (N - 2)*sqrt(left :* right)
	
	return(cm)
}

real scalar plssem_rebus_gqi(real matrix ind, real matrix lat, real matrix e, ///
	real matrix f, real colvector cl)
{
	/* Description:
		 ------------
		 Function that computes the Group Quality Index (GQI) for a REBUS-PLS
		 solution
	*/
	
	/* Arguments:
		 ----------
		 - ind		--> real matrix containing the indicator values stacked by
									REBUS-PLS group
		 - lat		--> real matrix containing the latent variable values stacked by
									REBUS-PLS group
		 - e 			--> real matrix containing the PLS-SEM measurement model residuals
		 - f			--> real matrix containing the PLS-SEM structural model residuals
		 - cl			--> real colvector containing the REBUS-PLS classification
	*/
	
	/* Returned value:
		 ---------------
		 - gqi		--> real scalar providing the Group Quality Index (GQI) value
	*/
	
	real scalar N, Ncl, k, gqi
	real colvector Nk, left, right
	real matrix e2, f2, e2_k, f2_k, ind_k, lat_k, ind_k_den, lat_k_den
	
	N = rows(ind)
	Ncl = max(cl)
	e2 = e:^2
	f2 = f:^2
	Nk = J(0, 1, .)
	left = J(0, 1, .)
	right = J(0, 1, .)
	for (k = 1; k <= Ncl; k++) {
		Nk = (Nk \ rows(cl[selectindex(cl :== k), 1]))
		e2_k = e2[selectindex(cl :== k), .]
		f2_k = f2[selectindex(cl :== k), .]
		ind_k = ind[selectindex(cl :== k), .]
		lat_k = lat[selectindex(cl :== k), .]
		ind_k_den = diagonal(crossdev(ind_k, mean(ind_k), ind_k, mean(ind_k)))
		lat_k_den = diagonal(crossdev(lat_k, mean(lat_k), lat_k, mean(lat_k)))
		left = (left \ (1 - mean(colsum(e2_k)' :/ ind_k_den)))
		right = (right \ (1 - mean(colsum(f2_k)' :/ lat_k_den)))
	}
	// gqi = sqrt((Nk' * (left :* right))/N) // this is the calculation according
																					 // to Sanchez's plspm R package

	gqi = sqrt((Nk' * left)*(Nk' * right))/N	// this is the calculation according
																						// to Trinchera's PhD thesis (2007)

	return(gqi)
}

struct plssem_struct_rebus scalar plssem_rebus(real matrix X, real matrix M,
	real matrix S, string scalar ind, string scalar stdind, string scalar latents,
	string scalar binary, real scalar tol, real scalar maxit, string scalar touse,
	string scalar scheme, string scalar crit, string scalar init,
	string scalar scale, real scalar structural, real scalar rawsum,
	real scalar consistent, string scalar rebus_cl, real scalar numclass,
	real scalar maxit_reb, real scalar stop, |real scalar noisily)
{
	/* Description:
		 ------------
		 Function that implements the REBUS-PLS algorithm
	*/
	
	/* Arguments:
		 ----------
		 - X					--> real matrix containing the observed manifest variables
		 - M					--> real matrix containing the measurement model adjacency
											matrix
		 - S			 		--> real matrix containing the structural model adjacency
											matrix
		 - ind				--> string scalar providing the list of all the unstandardized
											indicators
		 - stdind			--> string scalar providing the list of all the standardized
											indicators
		 - latents 		--> string scalar providing the list of all the latent
											variables
		 - binary	 		--> string scalar providing the list of the latent binary
											variables
		 - tol		 		--> real scalar providing the tolerance to use in the
											algorithm
		 - maxit	 		--> real scalar providing the maximum number of PLS iterations
		 - touse	 		--> string scalar containing the name of the variable tracking
											the subset of the data to use
		 - scheme	 		--> string scalar containing the weighting scheme to use (one
											of "centroid", "factorial" or "path")
		 - crit		 		--> string scalar containing the convergence criterion (either
											"relative" or "square")
		 - init	 			--> string scalar containing the initialization type (either
											"indsum" or "eigen")
		 - scale	 		--> string scalar indicating if the indicators must be scaled
		 - structural	--> real scalar equal to 1 if the model has a structural part,
											and 0 otherwise
		 - rawsum			--> real scalar equal to 1 if the 'rawsum' option has been
											chosen, and 0 otherwise
		 - consistent	--> real scalar equal to 1 if the consistent version (PLSc) is
											required, and 0 otherwise
		 - rebus_cl		--> string scalar providing the variable with the REBUS-PLS
											classes
		 - numclass		--> real scalar providing the number of REBUS-PLS classes
		 - maxit_reb	--> real scalar providing the maximum number of REBUS-PLS
											iterations
		 - stop				--> real scalar providing the REBUS-PLS stopping criterion
		 - noisily		--> (optional) real scalar; if not 0, it prints the bootstrap
											iterations
	*/
	
	/* Returned value:
		 ---------------
		 - res_rebus	--> scalar plssem_struct_rebus structure containing the results
	*/

	struct plssem_struct_rebus scalar res_rebus
	struct plssem_struct scalar res
	struct plssem_struct colvector localmodels

	real matrix cm, Xsc, Yinit, Xreb, Xrebsc, Xreb_all, Xrebsc_all, ow, path, ///
		loads, y, y_loc, block, x_hat, y_hat, out_res, inn_res
	real scalar iter, N, nchanged, k, P, p, rN0, rN_lte_5
	string scalar touse_loc_name, todisp
	real colvector modes, touse_vec, rebus_class, touse_loc, old_class, ///
		new_class, class_freq
	real rowvector r2, endo
	
	iter = 1
	P = cols(M)
	touse_vec = st_data(., touse)
	N = sum(touse_vec)
	nchanged = N
	
	rebus_class = st_data(., rebus_cl)
	localmodels = J(numclass, 1, res)
	modes = J(P, 1, 0)
	endo = (colsum(S)	:> 0)			 // indicators for the endogenous latent variables
	
	if (noisily == .) {
		noisily = 1
	}
	
	(void) st_addvar("byte", touse_loc_name = st_tempname())
	
	old_class = rebus_class
	
	while ((iter <= maxit_reb) & (nchanged > N*stop)) {
		/*
		if (noisily) {
			if (iter == 1) {
				printf("{txt}REBUS-PLS iteration 1")
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
		cm = J(N, 0, .)
		for (k = 1; k <= numclass; k++) {
			touse_loc = touse_vec :* (old_class :== k)
			st_store(., touse_loc_name, touse_loc)
			
			// Check that there are no zero-variance indicators
			if (anyof(sd(X[selectindex(touse_loc), .]), 0)) {
				printf("{err}some indicators during the REBUS-PLS calculations turned out to have zero variance\n")
				printf("{err}try reducing the number of classes\n")
				_error(1)
			}
			
			// Standardize the MVs (if required)
			if (scale == "") {
				Xsc = scale(X[selectindex(touse_loc), .])
			}
			else {
				Xsc = X[selectindex(touse_loc), .]
			}
			
			// Check that there are no zero-variance indicators
			if (anyof(sd(Xsc), 0)) {
				printf("{err}some indicators during the REBUS-PLS calculations turned out to have zero variance\n")
				printf("{err}try reducing the number of classes\n")
				_error(1)
			}
			
			// Initialize the LVs
			plssem_scale(Xsc, stdind, touse_loc_name, scale)
			Yinit = plssem_init_mat(Xsc, M, ind, stdind, latents, touse_loc_name, ///
				rawsum, init)
			
			// Run the PLS algorithm
			localmodels[k, 1] = plssem_base(Xsc, Yinit, M, S, modes, latents, ///
				binary, tol, maxit, touse_loc_name, scheme, crit, structural, ///
				rawsum, consistent)
			
			Xreb = X[selectindex(touse_loc), .]
			Xrebsc = scale(Xreb)
			ow = localmodels[k, 1].outer_weights
			path = localmodels[k, 1].path
			path = editmissing(path, 0)
			loads = localmodels[k, 1].loadings
			r2 = localmodels[k, 1].r2
			r2 = r2[., selectindex(r2 :!= .)]
			//y = Xrebsc * ow ???
			
			Xreb_all = X[selectindex(touse_vec), .]
			Xrebsc_all = scale(Xreb_all, 0, sd(Xreb), mean(Xreb))
			y_loc = Xrebsc_all * ow
			out_res = J(N, 0, .)
			for (p = 1; p <= P; p++) {
				block = selectindex(loads[., p] :!= .)
				x_hat = y_loc[., p] * loads[block, p]'
				out_res = (out_res, (Xrebsc_all[., block] - x_hat))
			}
			y_hat = y_loc * path[., selectindex(endo)]
			inn_res = (y_loc[., selectindex(endo)] - y_hat)
			cm = (cm, plssem_rebus_cm(out_res, inn_res, loads, r2))
		}
		new_class = J(length(old_class), 1, .)
		new_class[selectindex(touse_vec), 1] = which(cm, "min")
		nchanged = sum(old_class :!= new_class)
		old_class = new_class
		
		class_freq = uniqrows(new_class[selectindex(touse_vec), .], 1)[., 2]
		if (anyof(class_freq, 0)) {
			rN0 = 1
			break
		}
		else {
			rN0 = 0
		}
		if (any(class_freq :<= 5) & all(class_freq :> 0)) {
			rN_lte_5 = 1
			break
		}
		else {
			rN_lte_5 = 0
		}
		
		iter++
	}
	
	// Assigning results to structure's members
	res_rebus.niter = (iter - 1)
	res_rebus.rebus_class = new_class
	res_rebus.touse_vec = touse_vec
	res_rebus.rN0 = rN0
	res_rebus.rN_lte_5 = rN_lte_5
	res_rebus.nclass = numclass
	res_rebus.maxiter = maxit_reb
	res_rebus.stop = stop
	
	return(res_rebus)
}

real colvector plssem_rebus_ptest(real matrix X, real matrix M, real matrix S,
	string scalar ind, string scalar stdind, string scalar latents,
	string scalar binary, real scalar tol, real scalar maxit, string scalar touse,
	string scalar scheme, string scalar crit, string scalar init,
	string scalar scale, real scalar structural, real scalar rawsum,
	real scalar consistent, string scalar rebus_cl, real scalar numclass,
	|real scalar B, real scalar seed, real scalar noisily)
{
	/* Description:
		 ------------
		 Function that implements the permutation test for a REBUS-PLS solution
	*/
	
	/* Arguments:
		 ----------
		 - X					--> real matrix containing the observed manifest variables
		 - M					--> real matrix containing the measurement model adjacency
											matrix
		 - S			 		--> real matrix containing the structural model adjacency
											matrix
		 - ind				--> string scalar providing the list of all the unstandardized
											indicators
		 - stdind			--> string scalar providing the list of all the standardized
											indicators
		 - latents 		--> string scalar providing the list of all the latent
											variables
		 - binary	 		--> string scalar providing the list of the latent binary
											variables
		 - tol		 		--> real scalar providing the tolerance to use in the
											algorithm
		 - maxit	 		--> real scalar providing the maximum number of PLS iterations
		 - touse	 		--> string scalar containing the name of the variable tracking
											the subset of the data to use
		 - scheme	 		--> string scalar containing the weighting scheme to use (one
											of "centroid", "factorial" or "path")
		 - crit		 		--> string scalar containing the convergence criterion (either
											"relative" or "square")
		 - scale	 		--> string scalar indicating if the indicators must be scaled
		 - init	 			--> string scalar containing the initialization type (either
											"indsum" or "eigen")
		 - structural	--> real scalar equal to 1 if the model has a structural part,
											and 0 otherwise
		 - rawsum			--> real scalar equal to 1 if the 'rawsum' option has been
											chosen, and 0 otherwise
		 - consistent	--> real scalar equal to 1 if the consistent version (PLSc) is
											required, and 0 otherwise
		 - rebus_cl		--> string scalar providing the variable with the REBUS-PLS
											classes
		 - numclass		--> real scalar providing the number of REBUS-PLS classes
		 - B					--> (optional) real scalar number of permutation replications
											(default 100)
		 - seed				--> (optional) real scalar providing the permutation seed
		 - noisily		--> (optional) real scalar; if not 0, it prints the
											permutation iterations
	*/
	
	/* Returned value:
		 ---------------
		 - res --> scalar plssem_struct_rebus structure containing the results
	*/
	
	struct plssem_struct scalar res
	struct plssem_struct colvector localmodels
	
	real matrix Xsc, Yinit, Xreb, Xrebsc, Yendo, ow, path, ///
		loads, y_loc, block, x_hat, y_hat, out_res, inn_res, indstd_st, ///
		lat_st, out_res_st, inn_res_st
	real scalar P, Q, b, p, k, i, skip
	string scalar touse_loc_name, todisp, spaces
	real colvector gqi_perm, modes, touse_vec, rebus_class, touse_loc, ///
		rebclass_p, class_st, class_reb
	real rowvector endo
	
	P = cols(M)
	Q = rows(M)
	touse_vec = st_data(., touse)
	
	if (B == .) {
		B = 100
	}
	if (seed != .) {
		rseed(seed)
	}
	if (noisily == .) {
		noisily = 1
	}
	
	rebus_class = st_data(., rebus_cl)
	rebclass_p = rebus_class
	localmodels = J(numclass, 1, res)
	modes = J(P, 1, 0)
	endo = colsum(S)						 // indicators for the endogenous latent variables
	gqi_perm = J(B, 1, .)
	
	(void) st_addvar("byte", touse_loc_name = st_tempname())
	
	if (noisily) {
		printf("{txt}\n")
		printf("{txt}REBUS-PLS --> Permutation replications (")
		printf("{res}%f", B)
		printf("{txt})\n")
		printf("{hline 4}{c +}{hline 3} 1 ")
		printf("{hline 3}{c +}{hline 3} 2 ")
		printf("{hline 3}{c +}{hline 3} 3 ")
		printf("{hline 3}{c +}{hline 3} 4 ")
		printf("{hline 3}{c +}{hline 3} 5\n")
	}
	skip = 0
	for (b = 1; b <= B; b++) {
		if (noisily) {
			if (mod(b, 50) == 0) {
				todisp = strtrim(strofreal(b, "%9.0f"))
				if (strlen(todisp) < strlen(strofreal(B))) {
					skip = strlen(strofreal(B)) - strlen(todisp)
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
		
		rebclass_p[selectindex(touse_vec)] = jumble(rebus_class[selectindex(touse_vec)])
		
		indstd_st = J(0, Q, .)
		lat_st = J(0, sum(endo :> 0), .)
		out_res_st = J(0, Q, .)
		inn_res_st = J(0, sum(endo :> 0), .)
		class_st = J(0, 1, .)
		for (k = 1; k <= numclass; k++) {
			touse_loc = touse_vec :* (rebclass_p :== k)
			st_store(., touse_loc_name, touse_loc)
			
			// Check that there are no zero-variance indicators
			if (anyof(sd(X[selectindex(touse_loc), .]), 0)) {
				printf("{err}some indicators during the REBUS-PLS calculations turned out to have zero variance\n")
				printf("{err}try reducing the number of classes\n")
				_error(1)
			}
			
			// Standardize the MVs (if required)
			Xsc = plssem_scale_mat(X, touse_loc, scale)
			
			// Check that there are no zero-variance indicators
			if (anyof(sd(Xsc), 0)) {
				printf("{err}some indicators during the REBUS-PLS calculations turned out to have zero variance\n")
				printf("{err}try reducing the number of classes\n")
				_error(1)
			}
			
			// Initialize the LVs
			plssem_scale(Xsc, stdind, touse_loc_name, scale)
			Yinit = plssem_init_mat(Xsc, M, ind, stdind, latents, touse_loc_name, ///
				rawsum, init)
			
			// Run the PLS algorithm
			localmodels[k, 1] = plssem_base(Xsc, Yinit, M, S, modes, latents, ///
				binary, tol, maxit, touse_loc_name, scheme, crit, structural, ///
				rawsum, consistent)

			Xreb = X[selectindex(touse_loc), .]
			Xrebsc = scale(Xreb)
			Yendo = st_data(., latents, touse_loc_name)
			Yendo = Yendo[., selectindex(endo)]
			
			ow = localmodels[k, 1].outer_weights
			path = localmodels[k, 1].path
			path = editmissing(path, 0)
			loads = localmodels[k, 1].loadings
			
			y_loc = Xrebsc * ow
			out_res = J(rows(Xreb), 0, .)
			for (p = 1; p <= P; p++) {
				block = selectindex(loads[., p] :!= .)
				x_hat = y_loc[., p] * loads[block, p]'
				out_res = (out_res, (Xrebsc[., block] - x_hat))
			}
			y_hat = y_loc * path[., selectindex(endo)]
			inn_res = y_loc[., selectindex(endo)] - y_hat
			class_reb = J(rows(Xreb), 1, k)
			
			indstd_st = (indstd_st \ Xrebsc)
			lat_st = (lat_st \ Yendo)
			out_res_st = (out_res_st \ out_res)
			inn_res_st = (inn_res_st \ inn_res)
			class_st = (class_st \ class_reb)
		}
		gqi_perm[b, 1] = plssem_rebus_gqi(indstd_st, lat_st, out_res_st, inn_res_st, class_st)
	}
	
	return(gqi_perm)
}

end
