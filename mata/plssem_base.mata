*!plssem_base version 0.3.1
*!Written 16Apr2019
*!Written by Sergio Venturini and Mehmet Mehmetoglu
*!The following code is distributed under GNU General Public License version 3 (GPL-3)

version 15.1

mata:

void plssem_scale(real matrix X, string scalar stdind, string scalar touse,
	string scalar scale)
{
	/* Description:
		 ------------
		 Function that implements standardization of indicator variables
	*/
	
	/* Arguments:
		 ----------
		 - X					--> real matrix containing the observed manifest variables
		 - stdind	 		--> string scalar providing the list of the standardized
											indicators
		 - touse			--> string scalar containing the name of the variable tracking
											the subset of the data to use
		 - scale	 		--> string scalar indicating if the indicators are scaled
	*/
	
	/* Returned value:
		 ---------------
		 - void				--> the matrix of standardized indicators is saved directly
											in the data table
	*/
	
	real scalar V, v
	string vector std_ind
	
	std_ind = tokens(stdind)
	V = length(std_ind)
	
	for (v = 1; v <= V; v++) {
		stata("capture quietly generate " + std_ind[v] + " = .", 1)
	}
	if (scale == "") {
		st_store(., std_ind, touse, scale(X))
	}
	else {
		//st_store(., std_ind, touse, X :- mean(X))
		st_store(., std_ind, touse, X)
	}
}

real matrix plssem_scale_mat(real matrix X, real colvector touse,
	string scalar scale)
{
	/* Description:
		 ------------
		 Function that implements standardization of the indicators (returning the
		 result as a matrix)
	*/
	
	/* Arguments:
		 ----------
		 - X					--> real matrix containing the observed manifest variables
		 - touse			--> real colvector containing the data subset to use
		 - scale	 		--> string scalar indicating if the indicators are scaled
	*/
	
	/* Returned value:
		 ---------------
		 - Xsc				--> real matrix of the standardized indicators
	*/
	
	real matrix X_touse, Xsc
	
	X_touse = X[selectindex(touse), .]
	
	if (scale == "") {
		Xsc = scale(X_touse)
	}
	else {
		//Xsc = X_touse :- mean(X_touse)
		Xsc = X_touse
	}
	
	return(Xsc)
}

void plssem_init(real matrix X,  real matrix M, string scalar ind,
	string scalar stdind, string scalar latents, string scalar touse,
	real scalar rawsum, string scalar init)
{
	/* Description:
		 ------------
		 Function that initializes the LVs in a PLS-SEM model
	*/
	
	/* Arguments:
		 ----------
		 - X					--> real matrix containing the observed manifest variables
		 - M			 		--> real matrix containing the measurement model adjacency
											matrix
		 - ind				--> string scalar providing the list of all the unstandardized
											indicators
		 - stdind			--> string scalar providing the list of all the standardized
											indicators
		 - latents 		--> string scalar providing the list of all the latent
											variables
		 - touse	 		--> string scalar containing the name of the variable tracking
											the subset of the data to use
		 - rawsum			--> real scalar equal to 1 if the 'rawsum' option has been
											chosen, and 0 otherwise
		 - init	 			--> string scalar containing the initialization type (either
											"indsum" or "eigen")
	*/
	
	/* Returned value:
		 ---------------
		 - void				--> the matrix of initial LV scores is saved directly in the
											data table
	*/
	
	real matrix Y, Xunscaled
	real scalar V, v
	string scalar tmpscores
	real colvector scores
	real rowvector monoitem
	string rowvector lvs, sinds
	
	lvs = tokens(latents)
	V = length(lvs)
	
	monoitem = (colsum(M) :== 1)
	sinds = tokens(stdind)
	
	for (v = 1; v <= V; v++) {
		stata("capture quietly generate " + lvs[v] + " = .", 1)
		stata("capture quietly replace " + lvs[v] + " = .", 1)
	}
	
	if (!rawsum) {
		if (init == "indsum") {
			Y = X * M
		}
		else if (init == "eigen") {
			for (v = 1; v <= V; v++) {
				tmpscores = st_tempname()
				if (!monoitem[v]) {
					stata("quietly factor " + invtokens(sinds[selectindex(M[., v])]) + ///
						" if " + touse + ", factors(1) pcf", 1)
					stata("quietly predict " + tmpscores + " if " + touse, 1)
					st_store(., lvs[v], touse, st_data(., tmpscores, touse))
				}
				else {
					st_store(., lvs[v], touse, X[., selectindex(M[., v])])
				}
			}
			Y = st_data(., latents, touse)
		}
		st_store(., lvs, touse, scale(Y))
	}
	else {
		Xunscaled = st_data(., ind, touse)
		for (v = 1; v <= V; v++) {
			scores = rowsum(Xunscaled[., selectindex(M[., v])])
			st_store(., lvs[v], touse, scores)
			stata("capture quietly generate rs_" + lvs[v] + " = .", 1)
			st_store(., "rs_" + lvs[v], touse, scores)
		}
		Y = st_data(., latents, touse)
		st_store(., lvs, touse, scale(Y))
	}
}

real matrix plssem_init_mat(real matrix X,  real matrix M, string scalar ind,
	string scalar stdind, string scalar latents, string scalar touse,
	real scalar rawsum, string scalar init)
{
	/* Description:
		 ------------
		 Function that initializes the LVs in a PLS-SEM model (returning the result
		 as a matrix)
	*/
	
	/* Arguments:
		 ----------
		 - X					--> real matrix containing the observed manifest variables
		 - M			 		--> real matrix containing the measurement model adjacency
											matrix
		 - ind				--> string scalar providing the list of all the unstandardized
											indicators
		 - stdind			--> string scalar providing the list of all the standardized
											indicators
		 - latents 		--> string scalar providing the list of all the latent
											variables
		 - touse	 		--> string scalar containing the name of the variable tracking
											the subset of the data to use
		 - rawsum			--> real scalar equal to 1 if the 'rawsum' option has been
											chosen, and 0 otherwise
		 - init	 			--> string scalar containing the initialization type (either
											"indsum" or "eigen")
	*/
	
	/* Returned value:
		 ---------------
		 - Y					--> real matrix of initial LV scores
	*/
	
	real matrix Y, Xunscaled
	real scalar V, v
	string scalar tmpscores
	real rowvector monoitem
	string rowvector lvs, sinds, globalmodel
	
	lvs = tokens(latents)
	V = length(lvs)
	
	monoitem = (colsum(M) :== 1)
	sinds = tokens(stdind)
	
	if (!rawsum) {
		if (init == "indsum") {
			Y = X * M
		}
		else if (init == "eigen") {
			Y = J(rows(X), 0, .)
			globalmodel = st_tempname()
			stata("_estimates hold " + globalmodel, 1)
			for (v = 1; v <= V; v++) {
				tmpscores = st_tempname()
				if (!monoitem[v]) {
					stata("quietly factor " + invtokens(sinds[selectindex(M[., v])]) + ///
						" if " + touse + ", factors(1) pcf", 1)
					stata("quietly predict " + tmpscores + " if " + touse, 1)
					Y = (Y, st_data(., tmpscores, touse))
				}
				else {
					Y = X[., selectindex(M[., v])]
				}
			}
			stata("_estimates unhold " + globalmodel, 1)
		}
	}
	else {
		Xunscaled = st_data(., ind, touse)
		Y = J(rows(Xunscaled), 0, .)
		for (v = 1; v <= V; v++) {
			Y = (Y, rowsum(Xunscaled[., selectindex(M[., v])]))
		}
	}
	Y = scale(Y)
	
	return(Y)
}

struct plssem_struct scalar plssem_base(real matrix X, real matrix Yinit,
	real matrix M, real matrix S, real colvector mode, string scalar latents,
	string scalar binary, real scalar tol, real scalar maxit, string scalar touse,
	string scalar scheme, string scalar crit, real scalar structural,
	real scalar rawsum, real scalar consistent)
{
	/* Description:
		 ------------
		 Function that implements the basic PLS-SEM algorithm
	*/
	
	/* Arguments:
		 ----------
		 - X					--> real matrix containing the observed manifest variables
		 - Yinit			--> real matrix containing the initialized latent variables
		 - M					--> real matrix containing the measurement model adjacency
											matrix
		 - S			 		--> real matrix containing the structural model adjacency
											matrix
		 - mode		 		--> real colvector containing the modes of the latent
											variables (each element must be either 0 ["reflective"]
											or 1 ["formative"])
		 - latents 		--> string scalar providing the list of all the latent
											variables
		 - binary	 		--> string scalar providing the list of the latent binary
											variables
		 - tol		 		--> real scalar providing the tolerance to use in the
											algorithm
		 - maxit	 		--> real scalar providing the maximum number of iterations
		 - touse	 		--> string scalar containing the name of the variable tracking
											the subset of the data to use
		 - scheme	 		--> string scalar containing the weighting scheme to use (one
											of "centroid", "factorial" or "path")
		 - crit		 		--> string scalar containing the convergence criterion (either
											"relative" or "square")
		 - structural	--> real scalar equal to 1 if the model has a structural part,
											and 0 otherwise
		 - rawsum			--> real scalar equal to 1 if the 'rawsum' option has been
											chosen, and 0 otherwise
		 - consistent	--> real scalar equal to 1 if the consistent version (PLSc) is
											required, and 0 otherwise
	*/
	
	/* Returned value:
		 ---------------
		 - res --> scalar plssem_struct structure containing the estimated
							 quantities (loadings, path coefficients, etc.)
	*/
	
	struct plssem_struct scalar res
	
	real matrix W, w_old, w_new, C, E, Yhat, Ycor, Ytilde, Ydep, Yindep, mf, ///
		fscores, T, beta, beta_v, Yexo, Yendo, Ypred, Whist, xlambda, lambda, ///
		Snew, s2_mat, Xsc, beta0, y_tmp, outer_res, inner_res, x_hat, y_hat, ///
		block, lambda_c, xlambda_c, beta_c
	real scalar iter, delta, P, Q, p, n, minval, maxval, converged
	real rowvector endo, diff, isnotbinary, r2, r2_a
	real colvector Mind, Sind, predec
	string rowvector lvs_vec
	string scalar logitcmd
	
	iter = 1
	delta = 1e6
	P = cols(M)   			// number of latent variables
	Q = rows(M)   			// number of manifest variables
	n = rows(X)					// sample size
	W = J(Q, P, 0)
	w_old = J(Q, P, 0)
	C = S + S'
	endo = colsum(S)		// indicators for the endogenous latent variables
	predec = rowsum(S)	// indicators for the predecessor latent variables
	beta = J(P, P, .)
	beta_v = J(P, P, .)
	Whist = J(0, Q, .)
	diff = J(1, 0, .)
	xlambda = J(Q, P, .)
	lambda = J(Q, P, .)
	r2 = J(1, P, .)
	r2_a = J(1, P, .)
	lambda_c = J(Q, P, .)
	xlambda_c = J(Q, P, .)
	beta_c = J(P, P, .)
	
	// Indicators for the binary latent variables
	lvs_vec = tokens(latents)
	isnotbinary = J(1, P, 1)
	for (p = 1; p <= P; p++) {
		if (strpos(binary, lvs_vec[p])) {
			isnotbinary[p] = 0
		}
	}
	
	if (!rawsum) {
		if (structural) {
			// Step 0: latent scores initialization as the raw sum of indicators
			// 				 (the initial latent scores are passed directly)
			Yhat = Yinit
			W = M * diag(1 :/ sd(X * M))
			w_old = rowsum(W)'
			Whist = w_old
			
			while ((delta >= tol) & (iter <= maxit)) {
				// Step 1: inner weights estimation
				E = J(P, P, 0)
				Ycor = correlation(Yhat)
				
				if (scheme == "centroid") {
					E = sign(C :* Ycor)
				}
				else if (scheme == "factorial") {
					E = C :* Ycor
				}
				else if (scheme == "path") {
					for (p = 1; p <= P; p++) {
						if (endo[1, p]) {
							Sind = selectindex(S[., p])
							Yexo = (J(n, 1, 1), Yhat[., Sind])
							Yendo = Yhat[., p]
							E[Sind, p] = qrsolve(cross(Yexo, Yexo), cross(Yexo, Yendo))[|2 \ .|]
						}
						if (predec[p, 1]) {
							Sind = selectindex(S[p, .])
							Ypred = (Yhat[., p], Yhat[., Sind])
							E[Sind, p] = correlation(Ypred)[|2, 1 \ ., 1|]
						}
					}
				}
				
				// Step 2: inner approximation
				Ytilde = Yhat * E
				Ytilde = scale(Ytilde)
				
				// Step 3: outer weights estimation
				W = M
				for (p = 1; p <= P; p++) {
					Mind = selectindex(M[., p])
					if (rows(Mind) == 1) {
						continue
					}
					mf = X[., Mind]
					fscores = Ytilde[., p]
					/* (this is the code according to the algortihm) */
					if (!mode[p, 1]) {	// mode A
						W[Mind, p] = correlation((fscores, mf))[|2, 1 \ ., 1|]
					}
					else {							// mode B
						mf = (J(n, 1, 1), mf)
						W[Mind, p] = qrsolve(cross(mf, mf), cross(mf, fscores))[|2 \ .|]
					}
					/* (the following code is as in the semPLS package) */
					/*
					if (mode[p, 1]) {	// mode B
						T = luinv(cholesky(correlation(mf)))'
						mf = mf * T
					}
					W[Mind, p] = correlation((fscores, mf))[|2, 1 \ ., 1|]
					if (mode[p, 1]) {	// mode B
						W[Mind, p] = T * W[Mind, p]
					}
					*/
				}
				
				// Step 4: outer approximation
				Ytilde = X * W
				Yhat = scale(Ytilde)
				W = W * diag(1 :/ sd(Ytilde))
				w_new = rowsum(W)'
				Whist = (Whist \ w_new)

				// Step 5: convergence check
				if (crit == "relative") {
					delta = max(abs((w_old - w_new) :/ w_new))
					// delta = mreldif(w_old, w_new)
				}
				else if (crit == "square") {
					delta = max((w_old - w_new) :^ 2)
				}
				diff = (diff, delta)
				w_old = w_new
				iter++
			}
			converged = ((delta < tol) & (iter - 1 <= maxit))
		}
		else {
			W = M
			for (p = 1; p <= P; p++) {
				Mind = selectindex(M[., p])
				if (rows(Mind) == 1) {
					continue
				}
				if (!mode[p, 1]) {	// mode A
					W[Mind, p] = correlation((Yinit[., p], X[., Mind]))[|2, 1 \ ., 1|]
				}
				else {							// mode B
					Yindep = (J(n, 1, 1), X[., Mind])
					Ydep = Yinit[., p]
					W[Mind, p] = qrsolve(cross(Yindep, Yindep), cross(Yindep, Ydep))[|2 \ .|]
				}
			}
			Yhat = scale(X * W)
			iter = .
			converged = 0
			diff = .
			E = .
		}
	}
	else {
		Yhat = Yinit
		iter = .
		converged = 0
		diff = .
		E = .
		W = M
	}
	
	// Make the binary latent variables equal to 0 or 1
	for (p = 1; p <= P; p++) {
		if (!isnotbinary[p]) {
			minval = min(Yhat[., p])
			maxval = max(Yhat[., p])
			if (minval < maxval) {
				Yhat[., p] = editvalue(Yhat[., p], minval, 0)
				Yhat[., p] = editvalue(Yhat[., p], maxval, 1)
			}
			else {
				Yhat[., p] = editvalue(Yhat[., p], minval, 1)
				Yhat[., p] = editvalue(Yhat[., p], maxval, 0)
			}
		}
	}
	
	// Store latent scores in the data set
	st_store(., tokens(latents), touse, Yhat)
	
	// Path coefficients estimation
	if (structural) {
		real matrix xx, xxinv, xy, H, V
		real scalar nendo, sse, sst, s2, jj
		
		nendo = sum(endo :> 0)
		s2_mat = J(nendo, nendo, 0)
		jj = 1
		
		for (p = 1; p <= P; p++) {
			if (endo[1, p]) {
				Sind = selectindex(S[., p])
				Yexo = (J(n, 1, 1), Yhat[., Sind])
				Yendo = Yhat[., p]
				if (isnotbinary[p]) {
					xx = quadcross(Yexo, Yexo)
					xxinv = invsym(xx)
					xy = quadcross(Yexo, Yendo)
					beta[Sind, p] = qrsolve(xx, xy)[|2 \ .|]
					H = Yexo * xxinv * Yexo'
					sse = Yendo' * (I(n) - H) * Yendo
					sst = quadcrossdev(Yendo, mean(Yendo), Yendo, mean(Yendo))
					r2[p] = 1 - sse/sst
					r2_a[p] = 1 - (sse/sst)*(n - 1)/(n - cols(Yexo))
					s2 = sse/(n - cols(Yexo))
					V = s2 * xxinv
					beta_v[Sind, p] = diagonal(V)[|2 \ .|]
					s2_mat[jj, jj] = s2
					jj++
				}
				else {
					logitcmd = "quietly logit " + lvs_vec[p] + " " + ///
						invtokens(lvs_vec[selectindex(S[., p])]) + ///
						" if " + "`" + "touse" + "'"
					stata(logitcmd, 0)
					beta[Sind, p] = st_matrix("e(b)")[|1 \ rows(Sind)|]
					r2_a[p] = st_numscalar("e(r2_p)")
					beta_v[Sind, p] = diagonal(st_matrix("e(V)"))[|1 \ rows(Sind)|]
				}
			}
		}

		// Total effects estimation
		real matrix beta_tmp, ret, step
		
		beta_tmp = editmissing(beta, 0)
		ret = beta_tmp
		step = beta_tmp
		for (p = 2; p <= P; p++) {
			step = step * beta_tmp
			ret = step + ret
		}
		
		// Resizing matrices
		/*
		beta = beta[selectindex(predec), selectindex(endo)]
		beta_v = beta_v[selectindex(predec), selectindex(endo)]
		r2 = r2[selectindex(endo)]
		r2_a = r2_a[selectindex(endo)]
		*/
		Snew = S
	}
	else {
		beta = .
		beta_v = .
		s2_mat = .
		ret = .
		Snew = .
	}
	
	// Loadings estimation
	for (p = 1; p <= P; p++) {
		T = (Yhat[., p], X)
		xlambda[., p] = correlation(T)[|2, 1 \ ., 1|]
		Mind = selectindex(M[., p])
		lambda[Mind, p] = xlambda[Mind, p]
	}
	
	// Residuals calculation
	Xsc = scale(X)
	y_tmp = Xsc * W
	outer_res = J(n, 0, .)
	for (p = 1; p <= P; p++) {
		if (!mode[p, 1]) {	// mode A
			block = selectindex(lambda[., p] :!= .)
			x_hat = y_tmp[., p] * lambda[block, p]'
			outer_res = (outer_res, (Xsc[., block] - x_hat))
		}
	}
	if (structural) {
		beta0 = editmissing(beta, 0)
		y_hat = y_tmp * beta0[., selectindex(endo)]
		inner_res = (y_tmp[., selectindex(endo)] - y_hat)
	}
	else {
		inner_res = J(n, 0, .)
	}
	
	// Partial least-squares consistent estimation (PLSc)
	real matrix R, rhoA_mat, Ratt, beta0_c, Yhat_c, inner_res_c
	real colvector rhoA
	real rowvector r2_c, r2_a_c
	real scalar nmiss
	
	if (consistent) {
		y_tmp = scale(X) * W
		if (structural) {
			R = correlation(Yhat)
			rhoA = consistentrho(X, M, W, mode)
			
			r2_c = J(1, P, .)
			r2_a_c = J(1, P, .)
	
			rhoA_mat = rhoA' * rhoA
			for (p = 1; p <= P; p++) {
				rhoA_mat[p, p] = 1
			}
			Ratt = R :/ sqrt(rhoA_mat)
			for (p = 1; p <= P; p++) {
				if (endo[1, p]) {
					if (isnotbinary[p]) {
						Sind = selectindex(S[., p])
						beta_c[Sind, p] = invsym(Ratt[Sind, Sind]) * Ratt[Sind, p]
						r2_c[p] = Ratt[Sind, p]' * invsym(Ratt[Sind, Sind]) * Ratt[Sind, p]
						r2_a_c[p] = 1 - (1 - r2_c[p])*(n - 1)/(n - rows(Sind) - 1)
					}
					else {
						printf("[TODO]\n")
					}
				}
			}

			// Inner model's residuals calculation
			beta0_c = editmissing(beta_c, 0)
			Yhat_c = y_tmp * beta0_c
			inner_res_c = (y_tmp - Yhat_c)

			// Total effects estimation
			real matrix beta_tmp_c, ret_c
			
			beta_tmp_c = editmissing(beta_c, 0)
			ret_c = beta_tmp_c
			step = beta_tmp_c
			for (p = 2; p <= P; p++) {
				step = step * beta_tmp_c
				ret_c = step + ret_c
			}
		}
		else {
			Ratt = .
			beta_c = .
			ret_c = .
			r2_c = .
			r2_a_c = .
			inner_res_c = .
		}
		
		lambda_c = W * diag(sqrt(rhoA) :/ (diagonal(W' * W))')
		xlambda_c = xlambda * diag(1 :/ sqrt(rhoA))
		xlambda_c = xlambda_c :* (lambda :== .) + lambda_c
		for (p = 1; p <= P; p++) {
			Mind = selectindex(M[., p])
			if (length(Mind) == 1) {
				lambda_c[Mind, p] = 1
				xlambda_c[Mind, p] = 1
			}
			if (mode[p, 1]) {	// mode B (no correction)
				lambda_c[Mind, p] = lambda[Mind, p]
				xlambda_c[Mind, p] = xlambda[Mind, p]
			}
			nmiss = missing(lambda[., p])
			lambda_c[selectindex(lambda[., p] :== .), p] = J(nmiss, 1, .)
		}

		// Outer model's residuals calculation
		real matrix x_hat_c, outer_res_c

		outer_res_c = J(n, 0, .)
		for (p = 1; p <= P; p++) {
			if (!mode[p, 1]) {	// mode A
				block = selectindex(lambda_c[., p] :!= .)
				x_hat_c = y_tmp[., p] * lambda_c[block, p]'
				outer_res_c = (outer_res_c, (Xsc[., block] - x_hat_c))
			}
		}

		// Store latent scores in the data set
		//st_store(., tokens(latents), touse, Yhat_c)		// not clear what to save!
	}
	else {
		Ratt = .
		beta_c = .
		lambda_c = .
		xlambda_c = .
		ret_c = .
		r2_c = .
		r2_a_c = .
		outer_res_c = .
		inner_res_c = .
	}
	
	// Assigning results to structure's members
	res.n = n
	res.niter = (iter - 1)
	res.converged = converged
	res.diff = diff
	res.loadings = lambda
	res.xloadings = xlambda
	res.path = beta
	res.path_v = beta_v
	res.inner_v = s2_mat
	res.scores = Yhat
	res.total_effects = ret
	res.inner_weights = E
	res.outer_weights = W
	res.evo = Whist
	res.r2 = r2
	res.r2_a = r2_a
	res.M = M
	res.S = Snew
	res.e = outer_res
	res.f = inner_res
	res.R_c = Ratt
	res.loadings_c = lambda_c
	res.xloadings_c = xlambda_c
	res.path_c = beta_c
	res.scores_c = Yhat_c
	res.total_effects_c = ret_c
	res.r2_c = r2_c
	res.r2_a_c = r2_a_c
	res.e_c = outer_res_c
	res.f_c = inner_res_c
	
	return(res)
}

real matrix plssem_lv(real matrix X, real matrix Y, string scalar latents,
	string scalar binary)
{
	/* Description:
		 ------------
		 Function that computes the PLS-SEM (cross) loading variances assuming
		 normally distributed errors
	*/
	
	/* Arguments:
		 ----------
		 - X			 --> real matrix containing the observed manifest variables
		 - Y			 --> real matrix containing the latent variables scores
		 - latents --> string scalar providing the list of all the latent variables
		 - binary	 --> string scalar providing the list of the latent binary
									 variables
	*/
	
	/* Returned value:
		 ---------------
		 - xload_v --> real matrix containing the (cross) loading variances
	*/

	real matrix xload_v, y, x, xx, xxinv, H, V
	real scalar P, Q, p, q, n, sse, s2
	real rowvector isnotbinary
	string rowvector lvs_vec

	P = cols(Y)   			// number of latent variables
	Q = cols(X)   			// number of manifest variables
	n = rows(X)					// sample size
	xload_v = J(Q, P, 0)

	// Indicators for the binary latent variables
	lvs_vec = tokens(latents)
	isnotbinary = J(1, P, 1)
	for (p = 1; p <= P; p++) {
		if (strpos(binary, lvs_vec[p])) {
			isnotbinary[p] = 0
		}
	}

	// Variance computation
	for (p = 1; p <= P; p++) {
		if (isnotbinary[p]) {
			y = scale(Y[., p])
			for (q = 1; q <= Q; q++) {
				x = (J(n, 1, 1), scale(X[., q]))
				xx = quadcross(x, x)
				xxinv = invsym(xx)
				H = x * xxinv * x'
				sse = y' * (I(n) - H) * y
				s2 = sse/(n - 2)
				V = s2 * xxinv
				xload_v[q, p] = V[2, 2]
			}
		}
	}
	
	return(xload_v)
}

real matrix plssem_pval(real matrix path, real matrix path_se,
	string scalar latents, string scalar binary, real scalar n, real scalar boot)
{
	/* Description:
		 ------------
		 Function that computes the path coefficients' p-values for a PLS-SEM model
	*/
	
	/* Arguments:
		 ----------
		 - path				--> real matrix containing the path coefficients
		 - path_se		--> real matrix containing the path coefficients' standard
											errors
		 - latents		--> string scalar providing the list of all the latent
											variables
		 - binary			--> string scalar providing the list of the latent binary
											variables
		 - n					--> real scalar that provides the sample size
		 - boot				--> real scalar equal to 1 if the bootstrap variances are
											provided and 0 otherwise
	*/
	
	/* Returned value:
		 ---------------
		 - pval --> real matrix containing the path coefficients' variances
	*/

	real matrix pval
	real scalar P, p, q, stat, df
	real rowvector isnotbinary
	string rowvector lvs_vec

	P = cols(path)   			// number of latent variables
	pval = J(P, P, .)

	// Indicators for the binary latent variables
	lvs_vec = tokens(latents)
	isnotbinary = J(1, P, 1)
	for (p = 1; p <= P; p++) {
		if (strpos(binary, lvs_vec[p])) {
			isnotbinary[p] = 0
		}
	}

	// P-value computation
	for (p = 1; p <= P; p++) {
		for (q = 1; q <= P; q++) {
			if (!missing(path[q, p])) {
				stat = abs(path[q, p]) / path_se[q, p]
				pval[q, p] = 2*(1 - normal(stat))
				if ((isnotbinary[p]) & (!boot)) {
					df = n - colnonmissing(path[., p]) - 1
					pval[q, p] = 2*ttail(df, stat)
				}
			}
		}
	}

	return(pval)
}

real matrix plssem_pathtab(real matrix path, real matrix path_pval,
	real rowvector r2_a)
{
	/* Description:
		 ------------
		 Function for setting up the path coefficients' table for a PLS-SEM model
	*/
	
	/* Arguments:
		 ----------
		 - path				--> real matrix containing the path coefficients
		 - path_pval	--> real matrix containing the path coefficients' standard
											errors
		 - r2_a				--> real rowvector containing the structural model's adjusted
											R2's
	*/
	
	/* Returned value:
		 ---------------
		 - pathtab --> real matrix containing the path coefficients' table
	*/

	real matrix pathtab
	real colvector tokeep
	real scalar P, p

	P = rows(path)   						// number of latent variables
	pathtab = J(2*P + 1, P, .)
	tokeep = J(0, 1, .)

	for (p = 1; p <= P; p++) {
		pathtab[2*p - 1, .] = path[p, .]
		pathtab[2*p, .] = path_pval[p, .]
		if (missing(path[p, .]) != P) {
			tokeep = (tokeep \ (2*p - 1))
			tokeep = (tokeep \ 2*p)
		}
	}
	pathtab[2*P + 1, .] = r2_a
	tokeep = (tokeep \ (2*P + 1))
	pathtab = pathtab[tokeep, selectindex(colnonmissing(pathtab))]

	return(pathtab)
}

struct plssem_struct_boot scalar plssem_boot(real matrix X, real matrix Yinit,
	real matrix M, real matrix S, real colvector mode, string scalar latents,
	string scalar binary, real scalar tol, real scalar maxit, string scalar touse,
	string scalar scheme, string scalar crit, real scalar structural,
	real scalar rawsum, real scalar consistent, |real scalar B, real scalar seed,
	real scalar noisily)
{
	/* Description:
		 ------------
		 Function that implements the bootstrap version of the PLS-SEM algorithm
	*/
	
	/* Arguments:
		 ----------
		 - X					--> real matrix containing the observed manifest variables
		 - Yinit			--> real matrix containing the initialized latent variables
		 - M					--> real matrix containing the measurement model adjacency
											matrix
		 - S					--> real matrix containing the structural model adjacency
											matrix
		 - mode				--> real colvector containing the modes of the latent
											variables (each element must be either 0 ["reflective"]
											or 1 ["formative"])
		 - latents		--> string scalar providing the list of all the latent
											variables
		 - binary			--> string scalar providing the list of the latent binary
											variables
		 - tol				--> real scalar providing the tolerance to use in the
											algorithm
		 - maxit			--> real scalar providing the maximum number of iterations
		 - touse			--> string scalar containing the name of the variable tracking
											the subset of the data to use
		 - scheme			--> string scalar containing the weighting scheme to use (one
											of "centroid", "factorial" or "path")
		 - crit		 		--> string scalar containing the convergence criterion (either
											"relative" or "square")
		 - structural	--> real scalar equal to 1 if the model has a structural part,
											and 0 otherwise
		 - rawsum			--> real scalar equal to 1 if the 'rawsum' option has been
											chosen, and 0 otherwise
		 - consistent	--> real scalar equal to 1 if the consistent version (PLSc) is
											required, and 0 otherwise
		 - B					--> (optional) real scalar number of bootstrap replications
											(default 100)
		 - seed				--> (optional) real scalar providing the bootstrap seed
		 - noisily		--> (optional) real scalar; if not 0, it prints the bootstrap
											iterations
	*/
	
	/* Returned value:
		 ---------------
		 - res --> scalar plssem_struct_boot structure containing the bootstrap
							 replications for all unknown quantities
	*/

	struct plssem_struct scalar res
	struct plssem_struct_boot scalar res_bs
	
	real matrix X_bs, Yinit_bs, lambda, xlambda, beta, load_v, xload_v
	real vector bs_ind
	real scalar P, Q, n, b, skip, i, p, skip_b
	string scalar todisp, spaces
	real rowvector isnotbinary
	string rowvector lvs_vec

	P = cols(M)   			// number of latent variables
	Q = rows(M)   			// number of manifest variables
	n = rows(X)
	
	// Boostrap resampling settings
	if (B == .) {
		B = 100
	}
	if (seed != .) {
		rseed(seed)
	}
	
	lambda = J(B, Q, .)
	xlambda = J(B, Q*P, .)
	beta = J(B, P*P, .)

	// Indicators for the binary latent variables
	lvs_vec = tokens(latents)
	isnotbinary = J(1, P, 1)
	for (p = 1; p <= P; p++) {
		if (strpos(binary, lvs_vec[p])) {
			isnotbinary[p] = 0
		}
	}
	
	if (noisily == .) {
		noisily = 1
	}
	
	if (noisily) {
		printf("{txt}\n")
		printf("{txt}PLS-SEM --> Bootstrap replications (")
		printf("{res}%f", B)
		printf("{txt})\n")
		printf("{hline 4}{c +}{hline 3} 1 ")
		printf("{hline 3}{c +}{hline 3} 2 ")
		printf("{hline 3}{c +}{hline 3} 3 ")
		printf("{hline 3}{c +}{hline 3} 4 ")
		printf("{hline 3}{c +}{hline 3} 5\n")
	}
	skip = 0
	skip_b = 0
	for (b = 1; b <= B; b++) {
		if (noisily & (skip_b == 0)) {
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
		skip_b = 0
		bs_ind = runiformint(n, 1, 1, n)
		X_bs = X[bs_ind, .]
		Yinit_bs = Yinit[bs_ind, .]
		if (anyof(sd(X_bs), 0)) {
			skip_b = 1
			b--
			continue
		}
		res = plssem_base(X_bs, Yinit_bs, M, S, mode, latents, binary, tol, ///
			maxit, touse, scheme, crit, structural, rawsum, consistent)
		if (structural) {
			if (consistent) {
				beta[b, .] = vec(res.path_c)'
			}
			else {
				beta[b, .] = vec(res.path)'
			}
		}
		if (consistent) {
			lambda[b, .] = rowsum(res.loadings_c)'
			xlambda[b, .] = vec(res.xloadings_c)'
		} else {
			lambda[b, .] = rowsum(res.loadings)'
			xlambda[b, .] = vec(res.xloadings)'
		}
	}
	
	// Storing the latent scores using the original full sample
	res = plssem_base(X, Yinit, M, S, mode, latents, binary, tol, ///
			maxit, touse, scheme, crit, structural, rawsum, consistent)
	
	// Assigning results to structure's members
	res_bs.reps = B
	res_bs.loadings_reps = lambda
	res_bs.loadings_bs = plssem_boot_lm(lambda, M, "outer")
	load_v = plssem_boot_lv(lambda, M, "outer")
	for (p = 1; p <= P; p++) {
		if (!isnotbinary[p]) {
			load_v[., p] = J(Q, 1, 0)
		}
	}
	res_bs.loadings_v = load_v
	res_bs.xloadings_reps = xlambda
	res_bs.xloadings_bs = plssem_boot_lv(xlambda, M, "cross")
	xload_v = plssem_boot_lv(xlambda, M, "cross")
	for (p = 1; p <= P; p++) {
		if (!isnotbinary[p]) {
			xload_v[., p] = J(Q, 1, 0)
		}
	}
	res_bs.xloadings_v = xload_v
	if (structural) {
		res_bs.path_reps = beta
		res_bs.path_bs = plssem_boot_pm(beta)
		res_bs.path_v = plssem_boot_pv(beta)
	}
	else {
		res_bs.path_reps = .
		res_bs.path_bs = .
		res_bs.path_v = .
	}
	
	return(res_bs)
}

real matrix plssem_boot_lm(real matrix lambda, real matrix M,
	string scalar type)
{
	/* Description:
		 ------------
		 Function that computes the bootstrap estimates of the outer loadings
	*/
	
	/* Arguments:
		 ----------
		 - lambda		--> real matrix whose rows contain the bootstrap replications
										of the (stacked) loading matrix
		 - M				--> real matrix containing the measurement model adjacency
										matrix
		 - type			--> string scalar indicating whether outer or cross loadings
										are passed
	*/
	
	/* Returned value:
		 ---------------
		 - loading 	--> real matrix containing the loadings bootstrap estimates
	*/

	real matrix loading
	real colvector lambda_mean, Mind
	real scalar P, Q, p
	
	P = cols(M)   			// number of latent variables
	Q = rows(M)   			// number of manifest variables
	
	lambda_mean = mean(lambda)'
	loading = J(Q, P, 0)

	if (type == "outer") {
		for (p = 1; p <= P; p++) {
			Mind = selectindex(M[., p])
			loading[Mind, p] = lambda_mean[Mind]
		}
	}
	else {
		loading = rowshape(lambda_mean, P)'
	}
	
	return(loading)
}

real matrix plssem_boot_pm(real matrix path)
{
	/* Description:
		 ------------
		 Function that computes the bootstrap estimates of the path coefficients
	*/
	
	/* Arguments:
		 ----------
		 - path				--> real matrix whose rows contain the bootstrap replications
											of the (stacked) path coefficients matrix
	*/
	
	/* Returned value:
		 ---------------
		 - pathcoef 	--> real matrix containing the path coefficients bootstrap
											estimates
	*/

	real matrix pathcoef
	real vector path_mean
	real scalar P, j
	
	P = sqrt(cols(path))
	
	path_mean = J(1, P*P, .)
	pathcoef = J(P, P, 0)
	
	for (j = 1; j <= P*P; j++) {
		if (!hasmissing(path[., j])) {
			path_mean[j] = mean(path[., j])
		}
	}
	
	pathcoef = rowshape(path_mean, P)'
	
	return(pathcoef)
}

real matrix plssem_boot_lv(real matrix lambda, real matrix M,
	string scalar type)
{
	/* Description:
		 ------------
		 Function that computes the variance of the outer loadings bootstrap
		 estimates
	*/
	
	/* Arguments:
		 ----------
		 - lambda			--> real matrix whose rows contain the bootstrap replications
											of the (stacked) loading matrix
		 - M					--> real matrix containing the measurement model adjacency
											matrix
		 - type				--> string scalar indicating whether outer or cross loadings
											are passed
	*/
	
	/* Returned value:
		 ---------------
		 - loading_v 	--> real matrix containing the loadings bootstrap variances
	*/

	real matrix loading_v
	real colvector lambda_var, Mind
	real scalar P, Q, p
	
	P = cols(M)   			// number of latent variables
	Q = rows(M)   			// number of manifest variables
	
	lambda_var = diagonal(variance(lambda))
	loading_v = J(Q, P, .)
	
	if (type == "outer") {
		for (p = 1; p <= P; p++) {
			Mind = selectindex(M[., p])
			loading_v[Mind, p] = lambda_var[Mind]
		}
	}
	else {
		loading_v = rowshape(lambda_var, P)'
	}
	
	return(loading_v)
}

real matrix plssem_boot_pv(real matrix path)
{
	/* Description:
		 ------------
		 Function that computes the variance of the path coefficients bootstrap
		 estimates
	*/
	
	/* Arguments:
		 ----------
		 - path					--> real matrix whose rows contain the bootstrap
												replications of the (stacked) path coefficients matrix
	*/
	
	/* Returned value:
		 ---------------
		 - pathcoef_v 	--> real matrix containing the path coefficients bootstrap
												estimates
	*/

	real matrix pathcoef_v
	real vector path_var
	real scalar P, j
	
	P = sqrt(cols(path))
	
	path_var = J(1, P*P, .)
	pathcoef_v = J(P, P, 0)
	
	for (j = 1; j <= P*P; j++) {
		if (!hasmissing(path[., j])) {
			path_var[j] = variance(path[., j])
		}
	}
	
	pathcoef_v = rowshape(path_var, P)'
	
	return(pathcoef_v)
}

real matrix plssem_reliability(real matrix X, real matrix load,
	real colvector mode, real matrix outerw)
{
	/* Description:
		 ------------
		 Function that computes the reliability coefficients for a PLS-SEM model
	*/
	
	/* Arguments:
		 ----------
		 - X					--> real matrix containing the observed manifest variables
		 - load				--> real matrix containing the estimated outer loadings
		 - mode				--> real colvector containing the modes of the latent
											variables (each element must be either 0 ["reflective"]
											or 1 ["formative"])
		 - outerw			--> real matrix containing the outer weigths
	*/
	
	/* Returned value:
		 ---------------
		 - relcoef --> real matrix containing the latent variable reliability
									 coefficients
	*/

	real matrix M, relcoef
	real scalar P

	P = cols(load)   			// number of latent variables
	M = (load :!= .)
	relcoef = J(3, P, 0)

	relcoef[1, .] = cronbach(X, M, mode)
	relcoef[2, .] = dillongoldstein(X, load, mode)
	if (outerw == J(0, 0, .)) {
		relcoef[3, .] = J(1, P, .)
	}
	else {
		relcoef[3, .] = consistentrho(X, M, outerw, mode)
	}
	
	return(relcoef)
}

real matrix plssem_vif(real matrix R, real matrix S)
{
	/* Description:
		 ------------
		 Function that computes the VIFs for the structural part of a PLS-SEM model
	*/
	
	/* Arguments:
		 ----------
		 - R			--> real matrix containing the latent variable correlations
		 - S			--> real matrix containing the structural model adjacency
									matrix
	*/
	
	/* Returned value:
		 ---------------
		 - vif		--> real matrix containing the VIF values
	*/

	real matrix vif, Sind
	real scalar P, p
	real rowvector endo

	P = cols(S)   			// number of manifest variables
	endo = colsum(S)		// indicators for the endogenous latent variables
	vif = J(P, P, .)
	
	for (p = 1; p <= P; p++) {
		if (endo[1, p]) {
			Sind = selectindex(S[., p])
			vif[Sind, p] = diagonal(invsym(R[Sind, Sind]))
		}
	}
	
	return(vif)
}

real rowvector cronbach(real matrix X, real matrix M, real colvector mode)
{
	/* Description:
		 ------------
		 Function that computes the standardizes Cronbach's reliability coefficients
	*/
	
	/* Arguments:
		 ----------
		 - X			--> real matrix containing the latent variable scores
		 - M			--> measurement model's adjacency matrix
		 - mode		--> real colvector containing the modes of the latent variables
									(each element must be either 0 for reflective or 1 for
									formative)
	*/
	
	/* Returned value:
		 ---------------
		 - alpha	--> real matrix containing the Cronbach's reliability coefficients
	*/
	
	real matrix Mind
	real scalar P, n, p, k, l, r_bar
	real rowvector alpha
	
	P = cols(M)
	n = rows(X)
	alpha = J(1, P, 1)
	
	for (p = 1; p <= P; p++) {
		Mind = selectindex(M[., p])
		k = rows(Mind)
		if (k == 1) {
			continue
		}
		l = k*(k - 1)/2
		r_bar = sum(lowertriangle(correlation(X[., Mind]), 0))/l
		alpha[1, p] = k/(1/r_bar + (k - 1))
		if (alpha[p] < 0) {
			alpha[1, p] = 0
		}
		if (mode[p]) {
			alpha[1, p] = .
		}
	}
	
	return(alpha)
}

real rowvector dillongoldstein(real matrix X, real matrix load,
	real colvector mode)
{
	/* Description:
		 ------------
		 Function that computes the Dillon-Goldstein's reliability coefficients
	*/
	
	/* Arguments:
		 ----------
		 - X			--> real matrix containing the latent variable scores
		 - load		--> real matrix containing the estimated outer loadings
		 - mode		--> real colvector containing the modes of the latent variables
									(each element must be either 0 for reflective or 1 for
									formative)
	*/
	
	/* Returned value:
		 ---------------
		 - rho		--> real matrix containing the Dillon-Goldstein's reliability
									coefficients
	*/
	
	real matrix M, Mind, V, scores
	real scalar P, p, k, num, den
	real rowvector S, lambdas, rho

	M = (load :!= .)
	P = cols(M)   			// number of latent variables
	/*
	V = .
	S = .
	*/
	rho = J(1, P, 1)

	for (p = 1; p <= P; p++) {
		Mind = selectindex(M[., p])
		k = rows(Mind)
		if (k == 1) {
			continue
		}
		/* // this is based on equation (9) from Tenenhaus et al. (2005)
		symeigensystem(correlation(X[., Mind]), V, S)
		scores = X[., Mind] * V
		lambdas = correlation((X[., Mind], scores[, 1]))[|(k + 1), 1 \ (k + 1), k|]
		*/
		lambdas = load[Mind, p]
		num = sum(lambdas)^2
		den = num + (k - sum(lambdas :^ 2))
		rho[1, p] = num/den
		if (mode[p]) {
			rho[1, p] = .
		}
	}
	
	return(rho)
}

real rowvector consistentrho(real matrix X, real matrix M, real matrix outerw,
	real colvector mode)
{
	/* Description:
		 ------------
		 Function that computes the Dillon-Goldstein's reliability coefficients
	*/
	
	/* Arguments:
		 ----------
		 - X			--> real matrix containing the latent variable scores
		 - M			--> measurement model's adjacency matrix
		 - outerw	--> real matrix containing the model's outer weights
		 - mode		--> real colvector containing the modes of the latent variables
									(each element must be either 0 for reflective or 1 for
									formative)
	*/
	
	/* Returned value:
		 ---------------
		 - rhoA		--> real matrix containing the consistent reliability coefficients
	*/
	
	real matrix Mind, S, w, w_wp
	real scalar P, p, k, wp_w, num, den
	real rowvector rhoA

	P = cols(M)   			// number of latent variables
	rhoA = J(1, P, 1)

	for (p = 1; p <= P; p++) {
		Mind = selectindex(M[., p])
		k = rows(Mind)
		if (k == 1) {
			continue
		}
		w = outerw[Mind, p]
		S = correlation(X[., Mind])
		wp_w = w' * w
		w_wp = w * w'
		num = wp_w^2*(w' * (S - diag(S)) * w)
		den = w' * (w_wp - diag(w_wp)) * w
		rhoA[1, p] = num/den
		if (mode[p]) {
			rhoA[1, p] = 1
		}
	}
	
	return(rhoA)
}

void meanimp(real matrix X, string scalar vars, |string scalar categ,
	string scalar touse)
{
	/* Description:
		 ------------
		 Function that performs mean imputation for a set of variables
	*/
	
	/* Arguments:
		 ----------
		 - X					--> real matrix containing the variable to impute
		 - vars	 			--> string scalar providing the list of the variables to
											impute
		 - categ	 		--> (optional) string scalar providing the list of the
											categorical variables
		 - touse	 		--> (optional) string scalar containing the name of the
											variable tracking the subset of the data to use
	*/
	
	/* Returned value:
		 ---------------
		 - void				--> the matrix of imputed variables is saved directly in the
											data table
	*/

	real matrix Ximp, ur, xvals, xfreqs
	real scalar Q, q, mostfreq
	string rowvector vars_vec
	
	Q = cols(X)
	
	Ximp = X
	vars_vec = tokens(vars)

	for (q = 1; q <= Q; q++) {
		if (hasmissing(Ximp[., q])) {
			if (strpos(categ, vars_vec[q])) {
				ur = uniqrows(Ximp[., q], 1)
				xvals = ur[., 1]
				xfreqs = ur[., 2]
				mostfreq = xvals[which(xfreqs, "max")]
				Ximp[., q] = editmissing(Ximp[., q], mostfreq)
			}
			else {
				if (st_vartype(vars_vec[q]) != "float") {
					stata("quietly recast float " + vars_vec[q])
				}
				Ximp[., q] = editmissing(Ximp[., q], mean(Ximp[., q]))
			}
		}
	}
	
	st_store(., vars_vec, touse, Ximp)
}

void knnimp(real matrix X, string scalar vars, real scalar k,
	|string scalar categ, string scalar touse, real scalar noisily,
	real scalar ties)
{
	/* Description:
		 ------------
		 Function that performs mean imputation for a set of variables
	*/
	
	/* Arguments:
		 ----------
		 - X					--> real matrix containing the variable to impute
		 - vars	 			--> string scalar providing the list of variables to impute
		 - k	 				--> real scalar providing the neighbor size
		 - categ	 		--> (optional) string scalar providing the list of the
											categorical variables
		 - touse	 		--> (optional) string scalar containing the name of the
											variable tracking the subset of the data to use
		 - noisily		--> (optional) real scalar; if not 0, it prints the bootstrap
											iterations
		 - ties				--> (optional) real scalar; if 1, it breaks ties randomly,
											otherwise ties are ordered according to observation
											number
	*/
	
	/* Returned value:
		 ---------------
		 - void				--> the matrix of imputed variables is saved directly in the
											data table
	*/

	real matrix Ximp, Xsub, ur, xvals, xfreqs
	real scalar n, Q, N, R, i, skip, skip_b, j, q, mostfreq, D_idx, knew,
		disti_num
	real colvector complete_idx, incomplete_idx, disti, tmp_idx, disti_uniq,
		disti_j
	real rowvector vars_nonmiss
	string rowvector vars_vec, vars_touse
	string scalar todisp, spaces
	
	n = rows(X)												// number of observations
	Q = cols(X)
	
	if (noisily == .) {
		noisily = 1
	}
	if (ties == .) {
		ties = 1
	}
	
	Ximp = X
	vars_vec = tokens(vars)
	complete_idx = selectindex(rowmissing(X) :== 0)
	N = rows(complete_idx)						// number of complete observations
	incomplete_idx = selectindex(rowmissing(X) :> 0)
	R = n - N													// number of incomplete observations
	
	if (R > 0) {
		if (noisily) {
			printf("{txt}\n")
			printf("{txt}PLS-SEM --> Imputing ")
			printf("{res}%f", R)
			printf("{txt} missing values using k nearest neighbors ")
			printf("{txt}(with k = " + strtrim(strofreal(k, "%9.0f")))
			printf("{txt})\n")
			printf("{hline 4}{c +}{hline 3} 1 ")
			printf("{hline 3}{c +}{hline 3} 2 ")
			printf("{hline 3}{c +}{hline 3} 3 ")
			printf("{hline 3}{c +}{hline 3} 4 ")
			printf("{hline 3}{c +}{hline 3} 5\n")
		}
		skip = 0
		skip_b = 0
		
		for (i = 1; i <= R; i++) {
			if (noisily & (skip_b == 0)) {
				if (mod(i, 50) == 0) {
					todisp = strtrim(strofreal(i, "%9.0f"))
					if (strlen(todisp) < strlen(strofreal(R))) {
						skip = strlen(strofreal(R)) - strlen(todisp)
					}
					spaces = ""
					for (j = 1; j <= (skip + 3); j++) {
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
			skip_b = 0
			
			vars_nonmiss = selectindex(colnonmissing(X[incomplete_idx[i], .]))
			vars_touse = vars_vec[vars_nonmiss]
			tmp_idx = selectindex(rowmissing(X[., vars_nonmiss]) :== 0)
			tmp_idx = (tmp_idx, range(1, rows(tmp_idx), 1))
			D_idx = tmp_idx[selectindex(incomplete_idx[i] :== tmp_idx[., 1]), 2]
			disti = euclidean_dist(X[tmp_idx[., 1], vars_nonmiss], D_idx)
			disti[D_idx] = 1e+100
			disti = (disti, tmp_idx)
			disti = sort(disti, (1, 2))				// ties ordered according to obs. number
			if (ties) {
				disti_uniq = uniqrows(disti[., 1], 0)
				disti_num = rows(disti_uniq)
				for (j = 1; j <= disti_num; j++) {
					disti_j = selectindex(disti[., 1] :== disti_uniq[j])
					disti[disti_j, .] = jumble(disti[disti_j, .])	 // ties randomly broken
				}
			}
			if (k > rows(disti)) {
				printf("\n")
				printf("{err}value of k too large\n")
			}
			Xsub = X[disti[1..k, 2], .]
			knew = k
			while (anyof(colmissing(Xsub), knew)) {
				knew++
				Xsub = X[disti[1..knew, 2], .]
			}
			for (q = 1; q <= Q; q++) {
				if (hasmissing(X[incomplete_idx[i], q])) {
					if (strpos(categ, vars_vec[q])) {
						ur = uniqrows(Xsub[., q], 1)
						xvals = ur[., 1]
						xfreqs = ur[., 2]
						mostfreq = xvals[which(xfreqs, "max")]
						Ximp[incomplete_idx[i], q] = ///
							editmissing(Ximp[incomplete_idx[i], q], mostfreq)
					}
					else {
						if (st_vartype(vars_vec[q]) != "float") {
							stata("quietly recast float " + vars_vec[q])
						}
						Ximp[incomplete_idx[i], q] = ///
							editmissing(Ximp[incomplete_idx[i], q], mean(Xsub[., q]))
					}
				}
			}
		}
		if (noisily) {
			printf("\n")
			displayflush()
		}
		
		stata("capture quietly matrix drop D", 1)
		st_store(., vars_vec, touse, Ximp)
	}
}

end
