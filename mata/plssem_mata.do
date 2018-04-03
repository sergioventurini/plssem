*!plssem_mata version 0.3.0
*!Written 03Apr2018
*!Written by Sergio Venturini and Mehmet Mehmetoglu
*!The following code is distributed under GNU General Public License version 3 (GPL-3)

capture mata: mata drop ///
	plssem_struct() ///
	plssem_scale() ///
	plssem_scale_mat() ///
	plssem_init() ///
	plssem_init_mat() ///
	plssem_base() ///
	plssem_lv() ///
	plssem_pval() ///
	plssem_pathtab() ///
	plssem_struct_boot() ///
	plssem_boot() ///
	plssem_boot_lm() ///
	plssem_boot_pm() ///
	plssem_boot_lv() ///
	plssem_boot_pv() ///
	plssem_reliability() ///
	plssem_vif() ///
	plssem_mga_perm() ///
	plssem_mga_perm_diff() ///
	plssem_mga_boot() ///
	plssem_mga_boot_diff() ///
	scale() ///
	unscale() ///
	sd() ///
	which() ///
	cronbach() ///
	meanimp() ///
	knnimp() ///
	euclidean_dist() ///
	rdirichlet() ///
	plssem_rebus_cm() ///
	plssem_rebus_gqi() ///
	plssem_struct_rebus() ///
	plssem_rebus() ///
	plssem_rebus_ptest() ///
	plssem_struct_matrix() ///
	plssem_struct_fimix() ///
	plssem_fimix_ll() ///
	plssem_fimix_ll_em() ///
	plssem_fimix_estep() ///
	plssem_fimix_ic() ///
	plssem_fimix() ///
	print_struct_matrix() ///
	prod() ///
	cleanup()

version 14.2
mata:
mata set matastrict on
mata set matafavor speed

struct plssem_struct {
	real scalar n								// number of observations used in the analysis
	real scalar niter						// number of iterations to reach convergence
	real scalar converged				// indicator that the algorithm did converge
	real rowvector diff					// row vector of the max outer weights differences
	real matrix loadings				// matrix of estimated outer loadings
	real matrix xloadings				// matrix of estimated cross loadings
	real matrix path						// matrix of estimated path coefficients
	real matrix path_v					// matrix of estimated path coefficient variances
	real matrix inner_v					// diagonal matrix of inner relation variances
	real matrix scores					// matrix of estimated latent scores
	real matrix total_effects		// matrix of estimated total effects
	real matrix inner_weights		// matrix of estimated inner weights
	real matrix outer_weights		// matrix of estimated outer weights
	real matrix evo							// matrix of the outer weights evolution
	real rowvector r2						// row vector of the structural model's R2's
	real rowvector r2_a					// row vector of the structural model's adjusted
															// R2's
	real matrix M								// measurement model's adjacency matrix
	real matrix S								// structural model's adjacency matrix (if any)
}

struct plssem_struct_boot {
	real scalar reps						// number of bootstrap replications
	real matrix loadings_reps		// matrix of outer loadings bootstrap replications
	real matrix loadings_bs			// matrix of outer loadings bootstrap estimates
	real matrix loadings_v			// matrix of outer loadings bootstrap variances
	real matrix xloadings_reps	// matrix of cross loadings bootstrap replications
	real matrix xloadings_bs		// matrix of cross loadings bootstrap estimates
	real matrix xloadings_v			// matrix of cross loadings bootstrap variances
	real matrix path_reps				// matrix of path coefficients bootstrap
															// replications
	real matrix path_bs					// matrix of path coefficients bootstrap estimates
	real matrix path_v					// matrix of path coefficients bootstrap variances
}

struct plssem_struct_rebus {
	real scalar niter						// number of iterations to reach convergence
	real colvector rebus_class	// column vector of the final class memberships
	real colvector touse_vec		// column vector for the observations used
	real scalar rN0							// indicator that any of the classes is empty
	real scalar rN_lte_5				// indicator that any of the classes has 5
															// observations or less
	real scalar nclass					// number of classes
	real scalar maxiter					// maximum number of REBUS iterations
	real scalar stop						// REBUS stopping criterion
}

struct plssem_struct_matrix {
	real matrix mat							// generic real matrix
}

struct plssem_struct_fimix {
	real scalar niter						// number of iterations to reach convergence
	real colvector fimix_class	// column vector of the final class memberships
	real colvector touse_vec		// column vector of the observations used
	real scalar rN0							// indicator that any of the classes is empty
	real scalar rN_lte_5				// indicator that any of the classes has 5
															// observations or less
	real scalar nclass					// number of classes
	real scalar maxiter					// maximum number of FIMIX iterations
	real scalar stop						// FIMIX stopping criterion
	real scalar restart					// number of restarts of the EM algorithm
	real rowvector rho					// row vector of mixture proportions
	struct plssem_struct_matrix ///
		colvector path						// column vector structure with path coefficients
	real colvector classfreq		// column vector of class frequencies
	real colvector ll						// column vector of incomplete loglikelihood values
	real colvector ll_c					// column vector of complete loglikelihood values
	real colvector ll_restart		// column vector of incomplete loglikelihood
															//   values across the different restarts
	real matrix ic							// real matrix of information criteria
}

void plssem_scale(real matrix X, string scalar stdind, string scalar touse,
	string scalar scale)
{
	/* Description:
		 ------------
		 Function that implements standardization of the indicators
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
		 - touse			--> real colvector containing the name of the variable
											tracking the subset of the data to use
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
											chosen
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
											chosen
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
	string rowvector lvs, sinds
	
	lvs = tokens(latents)
	V = length(lvs)
	
	monoitem = (colsum(M) :== 1)
	sinds = tokens(stdind)
	
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
					Y = st_data(., tmpscores, touse)
				}
				else {
					Y = X[., selectindex(M[., v])]
				}
			}
		}
	}
	else {
		Xunscaled = st_data(., ind, touse)
		for (v = 1; v <= V; v++) {
			Y = rowsum(Xunscaled[., selectindex(M[., v])])
		}
	}
	Y = scale(Y)
	
	return(Y)
}

struct plssem_struct scalar plssem_base(real matrix X, real matrix Yinit,
	real matrix M, real matrix S, real colvector mode, string scalar latents,
	string scalar binary, real scalar tol, real scalar maxit, string scalar touse,
	string scalar scheme, string scalar crit, real scalar structural,
	real scalar rawsum)
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
											chosen
	*/
	
	/* Returned value:
		 ---------------
		 - res --> scalar plssem_struct structure containing the estimated
							 quantities (loadings, path coefficients, etc.)
	*/
	
	struct plssem_struct scalar res
	
	real matrix W, w_old, w_new, C, E, Yhat, Ycor, Ytilde, Ydep, Yindep, mf, ///
		fscores, T, beta, beta_v, Yexo, Yendo, Ypred, Whist, xlambda, lambda, ///
		Snew, s2_mat
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
			
			while (delta >= tol & iter <= maxit) {
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
	
	if (structural) {
		// Path coefficients estimation
		real matrix xx, xxm1, xy, H, V
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
					xxm1 = invsym(xx)
					xy = quadcross(Yexo, Yendo)
					beta[Sind, p] = qrsolve(xx, xy)[|2 \ .|]
					H = Yexo * xxm1 * Yexo'
					sse = Yendo' * (I(n) - H) * Yendo
					sst = quadcrossdev(Yendo, mean(Yendo), Yendo, mean(Yendo))
					r2[p] = 1 - sse/sst
					r2_a[p] = 1 - (sse/sst)*(n - 1)/(n - cols(Yexo))
					s2 = sse/(n - cols(Yexo))
					V = s2 * xxm1
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

	real matrix xload_v, y, x, xx, xxm1, H, V
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
				xxm1 = invsym(xx)
				H = x * xxm1 * x'
				sse = y' * (I(n) - H) * y
				s2 = sse/(n - 2)
				V = s2 * xxm1
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
				if (isnotbinary[p] & !boot) {
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
		 Function that sets up the path coefficients' table for a PLS-SEM model
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
	real scalar P, p

	P = rows(path)   						// number of latent variables
	pathtab = J(2*P + 1, P, .)

	for (p = 1; p <= P; p++) {
		pathtab[2*p - 1, .] = path[p, .]
		pathtab[2*p, .] = path_pval[p, .]
	}
	pathtab[2*P + 1, .] = r2_a
	pathtab = pathtab[selectindex(rownonmissing(pathtab)), ///
		selectindex(colnonmissing(pathtab))]

	return(pathtab)
}

struct plssem_struct_boot scalar plssem_boot(real matrix X, real matrix Yinit,
	real matrix M, real matrix S, real colvector mode, string scalar latents,
	string scalar binary, real scalar tol, real scalar maxit, string scalar touse,
	string scalar scheme, string scalar crit, real scalar structural,
	real scalar rawsum, |real scalar B, real scalar seed, real scalar noisily)
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
											chosen
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
			maxit, touse, scheme, crit, structural, rawsum)
		if (structural) {
			beta[b, .] = vec(res.path)'
		}
		lambda[b, .] = rowsum(res.loadings)'
		xlambda[b, .] = vec(res.xloadings)'
	}
	
	// Storing the latent scores using the original full sample
	res = plssem_base(X, Yinit, M, S, mode, latents, binary, tol, ///
			maxit, touse, scheme, crit, structural, rawsum)
	
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
	real colvector mode)
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
	*/
	
	/* Returned value:
		 ---------------
		 - relcoef --> real matrix containing the (cross) loading variances
	*/

	real matrix relcoef, M, Mind, V, scores
	real scalar P, p, k, num, den
	real rowvector S, lambdas

	P = cols(load)   			// number of latent variables
	relcoef = J(2, P, 0)
	V = .
	S = .

	M = (load :!= .)
	relcoef[1, .] = cronbach(X, M)
	
	for (p = 1; p <= P; p++) {
		Mind = selectindex(M[., p])
		k = rows(Mind)
		if (k == 1) {
			relcoef[2, p] = 1
			continue
		}
		symeigensystem(correlation(X[., Mind]), V, S)
		scores = scale(X[., Mind], 1) * V
		lambdas = correlation((X[., Mind], scores[, 1]))[|(k + 1), 1 \ (k + 1), k|]
		num = sum(lambdas)^2
		den = num + (k - sum(lambdas :^ 2))
		relcoef[2, p] = num/den
		if (mode[p]) {
			relcoef[1, p] = .
			relcoef[2, p] = .
		}
	}
	
	return(relcoef)
}

real matrix plssem_vif(real matrix Y, real matrix S)
{
	/* Description:
		 ------------
		 Function that computes the VIFs for the structural part of a PLS-SEM model
	*/
	
	/* Arguments:
		 ----------
		 - Y			--> real matrix containing the initialized latent variables
		 - S			--> real matrix containing the structural model adjacency
									matrix
	*/
	
	/* Returned value:
		 ---------------
		 - vif		--> real matrix containing the VIF values
	*/

	real matrix vif, Sind, X, y, x, xx, xxm1, H
	real scalar P, p, q, n, lv_i, sse, sst
	real rowvector endo

	P = cols(S)   			// number of manifest variables
	n = rows(Y)					// sample size
	endo = colsum(S)		// indicators for the endogenous latent variables
	vif = J(P, P, .)
	
	for (p = 1; p <= P; p++) {
		if (endo[1, p]) {
			Sind = selectindex(S[., p])
			lv_i = 1
			for (q = 1; q <= P; q++) {
				if (S[q, p]) {
					X = Y[., Sind]
					y = X[., lv_i]
					X[., lv_i] = J(n, 1, 1)
					x = X
					xx = quadcross(x, x)
					xxm1 = invsym(xx)
					H = x * xxm1 * x'
					sse = y' * (I(n) - H) * y
					sst = quadcrossdev(y, mean(y), y, mean(y))
					vif[q, p] = sst/sse
					lv_i++
				}
			}
		}
	}
	
	return(vif)
}

class AssociativeArray scalar plssem_mga_perm(real matrix X, real matrix Yinit,
	real matrix M, real matrix S, real colvector mode, string scalar latents,
	string scalar binary, real scalar tol, real scalar maxit, string scalar touse,
	string scalar scheme, string scalar crit, real scalar structural,
	real scalar rawsum, real colvector group, |real scalar B, ///
	real scalar seed, real scalar noisily)
{
	/* Description:
		 ------------
		 Function that performs the calculations for the permutation version of the
		 PLS-MGA algorithm
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
											chosen
		 - group			--> real colvector containing the group variable
		 - B					--> (optional) real scalar number of permutation replications
											(default 100)
		 - seed				--> (optional) real scalar providing the permutation seed
		 - noisily		--> (optional) real scalar; if not 0, it prints the
											permutation iterations
	*/
	
	/* Returned value:
		 ---------------
		 - res_mga	--> scalar plssem_struct_boot structure containing the
										permutation replications for all unknown quantities
	*/

	struct plssem_struct scalar res
	
	class AssociativeArray scalar res_mga
	res_mga.reinit("string", 2)
	
	real matrix bs_mat, alldata, group_tmp, X_bs, Yinit_bs, res_p, res_l
	real vector bs_ind
	real scalar P, Q, N, n, G, b, skip, i, g, groupi
	string scalar todisp, spaces, touse_new_name
	real colvector touse_vec, groupvals, groupsizes, __id__, __group__, ///
		touse_new, touse_tmp, res_tmp

	P = cols(M)   									// number of latent variables
	Q = rows(M)   									// number of manifest variables
	N = rows(X)											// overall sample size
	touse_vec = st_data(., touse)
	n = sum(touse_vec)							// sample size
	
	if (B == .) {
		B = 100
	}
	if (seed != .) {
		rseed(seed)
	}
	if (noisily == .) {
		noisily = 1
	}
		
	__id__ = range(1, N, 1)
	__id__[selectindex(!touse_vec)] = J(N - n, 1, .)
	bs_mat = J(N, B, .)							// used to track how the permutations evolve
	touse_new = J(N, 1, .)
	
	group_tmp = uniqrows(group[selectindex(touse_vec)], 1)
	groupvals = group_tmp[., 1]
	groupsizes = group_tmp[., 2]
	G = rows(groupvals)
	for (g = 1; g <= G; g++) {
		res_mga.put((strofreal(g), "path"), J(B, sum(S), .))
		res_mga.put((strofreal(g), "loadings"), J(B, sum(M), .))
	}
	
	groupi = 1
	__group__ = J(N, 1, .)
	for (g = 1; g <= G; g++) {
		__group__[|groupi, 1 \ (groupi + groupsizes[g] - 1), 1|] = ///
			J(groupsizes[g], 1, groupvals[g])
		groupi = groupi + groupsizes[g]
	}

	(void) st_addvar("byte", touse_new_name = st_tempname())
	
	if (noisily) {
		printf("{txt}\n")
		printf("{txt}Multigroup analysis --> Permutation test replications (")
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
		bs_ind = jumble(__id__[selectindex(touse_vec)])
		__id__[selectindex(touse_vec)] = bs_ind
		alldata = (__id__, X, Yinit, group, touse_vec)
		alldata = sort(alldata, 1)
		alldata = (alldata, __group__)
		bs_mat[., b] = alldata[., 1]
		X_bs = alldata[range(1, sum(groupsizes), 1), range(2, (Q + 1), 1)']
		Yinit_bs = alldata[range(1, sum(groupsizes), 1), ///
			range((Q + 2), (Q + P + 1), 1)']
		touse_tmp = alldata[., (Q + P + 3)]
		for (g = 1; g <= G; g++) {
			touse_new = touse_tmp :* (__group__ :== groupvals[g])
			st_store(., touse_new_name, touse_new)
			res = plssem_base(X_bs[selectindex(touse_new), .], ///
				Yinit_bs[selectindex(touse_new), .], M, S, mode, ///
				latents, binary, tol, maxit, touse_new_name, scheme, crit, ///
				structural, rawsum)
			if (structural) {
				res_tmp = res.path
				res_tmp = res_tmp[selectindex(rownonmissing(res_tmp)), ///
					selectindex(colnonmissing(res_tmp))]
				res_tmp = vec(res_tmp)
				res_tmp = res_tmp[selectindex(rownonmissing(res_tmp))]
				res_p = res_mga.get((strofreal(g), "path"))
				res_p[b, .] = res_tmp'
				res_mga.put((strofreal(g), "path"), res_p)
			}
			res_tmp = res.loadings
			res_tmp = vec(res_tmp)
			res_tmp = res_tmp[selectindex(rownonmissing(res_tmp))]
			res_l = res_mga.get((strofreal(g), "loadings"))
			res_l[b, .] = res_tmp'
			res_mga.put((strofreal(g), "loadings"), res_l)
		}
	}
	
	return(res_mga)
}

real matrix plssem_mga_perm_diff(class AssociativeArray res_mga, 
	real matrix obstest, string scalar what)
{
	/* Description:
		 ------------
		 Utility function used for returning results after MGA through permutations
	*/
	
	/* Arguments:
		 ----------
		 - res_mga		--> class AssociativeArray containing the resulsts of a MGA
											using permutations
		 - obstest		--> real matrix containing the observed test statistic values
											for all groups
		 - what				--> scalar string providing whether the path coefficients or
											the loadings must be used
	*/
	
	/* Returned value:
		 ---------------
		 - dtest			--> real matrix containing the results
	*/
	
	real matrix dtest, coef_1, coef_ng, diff_ng, dtest_ng
	real scalar neff, G, g
	
	neff = rows(obstest)
	G = cols(obstest) + 1
	
	dtest = J(G - 1, neff, .)
	
	coef_1 = res_mga.get((strofreal(1), what))
	for (g = 2; g <= G; g++) {
		coef_ng = res_mga.get((strofreal(g), what))
		diff_ng = abs(coef_1 - coef_ng)
		dtest_ng = (diff_ng :> (obstest[., g - 1]'))
		dtest[g - 1, .] = colsum(dtest_ng)
	}
	
	return(dtest)
}

class AssociativeArray scalar plssem_mga_boot(real matrix X, real matrix Yinit,
	real matrix M, real matrix S, real colvector mode, string scalar latents,
	string scalar binary, real scalar tol, real scalar maxit, string scalar touse,
	string scalar scheme, string scalar crit, real scalar structural,
	real scalar rawsum, real colvector group, |real scalar B, ///
	real scalar seed, real scalar noisily)
{
	/* Description:
		 ------------
		 Function that performs the calculations for the bootstrap version of the
		 PLS-MGA algorithm
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
											chosen
		 - group			--> real colvector containing the group variable
		 - B					--> (optional) real scalar number of bootstrap replications
											(default 100)
		 - seed				--> (optional) real scalar providing the bootstrap seed
		 - noisily		--> (optional) real scalar; if not 0, it prints the
											bootstrap iterations
	*/
	
	/* Returned value:
		 ---------------
		 - res_mga	--> scalar plssem_struct_boot structure containing the
										bootstrap replications for all unknown quantities
	*/
	
	struct plssem_struct scalar res
	
	class AssociativeArray scalar res_mga
	res_mga.reinit("string", 2)
	
	real matrix bs_mat, alldata, group_tmp, X_bs, Yinit_bs, alldata_bs, ///
		alldata_tmp, res_p, res_l
	real vector bs_ind
	real scalar P, Q, N, n, G, b, skip, i, g, skip_b
	string scalar todisp, spaces, touse_new_name
	real colvector touse_vec, groupvals, groupsizes, __id__, touse_new, ///
		touse_tmp, group_idx, res_tmp

	P = cols(M)   									// number of latent variables
	Q = rows(M)   									// number of manifest variables
	N = rows(X)											// overall sample size
	touse_vec = st_data(., touse)
	n = sum(touse_vec)							// sample size
	
	if (B == .) {
		B = 100
	}
	if (seed != .) {
		rseed(seed)
	}
	if (noisily == .) {
		noisily = 1
	}
	
	__id__ = range(1, N, 1)
	__id__[selectindex(!touse_vec)] = J(N - n, 1, .)
	bs_mat = J(N, B, .)							// used to track how the bootstrap evolves
	touse_new = J(N, 1, .)
	
	group_tmp = uniqrows(group[selectindex(touse_vec)], 1)
	groupvals = group_tmp[., 1]
	groupsizes = group_tmp[., 2]
	G = rows(groupvals)
	for (g = 1; g <= G; g++) {
		res_mga.put((strofreal(g), "path"), J(B, sum(S), .))
		res_mga.put((strofreal(g), "loadings"), J(B, sum(M), .))
	}
	
	(void) st_addvar("byte", touse_new_name = st_tempname())
	
	if (noisily) {
		printf("{txt}\n")
		printf("{txt}Multigroup analysis --> Bootstrap replications (")
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
		alldata = (__id__, X, Yinit, group, touse_vec)
		alldata_bs = alldata
		for (g = 1; g <= G; g++) {
			group_idx = (group :== groupvals[g]) :* touse_vec
			alldata_tmp = alldata[selectindex(group_idx), .]
			bs_ind = runiformint(groupsizes[g], 1, 1, groupsizes[g])
			alldata_bs[selectindex(group_idx), .] = alldata_tmp[bs_ind, .]
		}
		bs_mat[., b] = alldata_bs[., 1]
		touse_tmp = alldata_bs[., (Q + P + 3)]
		for (g = 1; g <= G; g++) {
			group_idx = (alldata_bs[., (Q + P + 2)] :== groupvals[g])
			touse_new = group_idx :* touse_tmp
			st_store(., touse_new_name, touse_new)
			X_bs = alldata_bs[., range(2, (Q + 1), 1)']
			Yinit_bs = alldata_bs[., range((Q + 2), (Q + P + 1), 1)']
			if (anyof(sd(X_bs[selectindex(touse_new), .]), 0)) {
				skip_b = 1
				break
			}
			res = plssem_base(X_bs[selectindex(touse_new), .], ///
				Yinit_bs[selectindex(touse_new), .], M, S, mode, ///
				latents, binary, tol, maxit, touse_new_name, scheme, crit, ///
				structural, rawsum)
			if (structural) {
				res_tmp = res.path
				res_tmp = res_tmp[selectindex(rownonmissing(res_tmp)), ///
					selectindex(colnonmissing(res_tmp))]
				res_tmp = vec(res_tmp)
				res_tmp = res_tmp[selectindex(rownonmissing(res_tmp))]
				res_p = res_mga.get((strofreal(g), "path"))
				res_p[b, .] = res_tmp'
				res_mga.put((strofreal(g), "path"), res_p)
			}
			res_tmp = res.loadings
			res_tmp = vec(res_tmp)
			res_tmp = res_tmp[selectindex(rownonmissing(res_tmp))]
			res_l = res_mga.get((strofreal(g), "loadings"))
			res_l[b, .] = res_tmp'
			res_mga.put((strofreal(g), "loadings"), res_l)
		}
		if (skip_b) {
			b--
		}
	}
	
	return(res_mga)
}

real matrix plssem_mga_boot_diff(class AssociativeArray res_mga, 
	real matrix groupsizes, real scalar neff, string scalar what)
{
	/* Description:
		 ------------
		 Utility function used for returning results after MGA through bootstrap
	*/
	
	/* Arguments:
		 ----------
		 - res_mga		--> class AssociativeArray containing the resulsts of a MGA
											using bootstrap
		 - groupsizes	--> real matrix containing the observed test statistic values
											for all groups
		 - neff				--> real scalar containing the number of effects to consider
		 - what				--> scalar string providing whether the path coefficients or
											the loadings must be used
	*/
	
	/* Returned value:
		 ---------------
		 - dtest			--> real matrix containing the results
	*/
	
	real matrix dtest, coef_1, coef_ng, coef_var_1, coef_var_ng, diff_ng
	real scalar k0, k1, k2, G, g, k3
	
	G = rows(groupsizes)
	k0 = sum(groupsizes) - 2
	k1 = ((groupsizes[1, 1] - 1)^2)/k0
	
	dtest = J(G - 1, neff, .)
	
	coef_1 = res_mga.get((strofreal(1), what))
	coef_var_1 = diagonal(variance(coef_1))'
	for (g = 2; g <= G; g++) {
		coef_ng = res_mga.get((strofreal(g), what))
		coef_var_ng = diagonal(variance(coef_ng))'
		diff_ng = abs(mean(coef_1) - mean(coef_ng))
		k2 = ((groupsizes[g, 1] - 1)^2)/k0
		k3 = sqrt(1/groupsizes[1, 1] + 1/groupsizes[g, 1])
		dtest[g - 1, .] = diff_ng :/ (sqrt(k1*coef_var_1 + k2*coef_var_ng)*k3)
	}
	
	return(dtest)
}

real matrix scale(real matrix X, |real scalar biased, real rowvector scale, ///
	real rowvector center)
{
	/* Description:
		 ------------
		 Function that standardizes the columns of a real matrix X
	*/
	
	/* Arguments:
		 ----------
		 - X			--> real matrix whose columns must be scaled
		 - biased	--> (optional) real scalar; if not zero the biased sample
									variances are computed
		 - scale	--> (optional) real rowvector containing the column scales
									(after centering)
		 - center	--> (optional) real rowvector containing the column centers
	*/
	
	/* Returned value:
		 ---------------
		 - Xs 		--> real matrix with columns that have been centered and/or scaled
	*/

	real matrix cen, sc, sc_inv, Xc, Xcs
	real scalar n, factor
	
	n = rows(X)
	if (center == J(1, 0, .)) {
		cen = J(n, 1, mean(X))
	}
	else {
		cen = J(n, 1, center)
	}
	if (scale == J(1, 0, .)) {
		sc = sqrt(diag(variance(X)))
	}
	else {
		sc = diag(scale)
	}
	if (biased == .) {
		biased = 0
	}
	if (biased) {
		factor = sqrt((n - 1)/n)
		sc = sc :* factor
	}
	Xc = X - cen
	sc_inv = invsym(sc)
	Xcs = Xc * sc_inv

	return(Xcs)
}

real matrix unscale(real matrix Xs, real rowvector scale, ///
	real rowvector center)
{
	/* Description:
		 ------------
		 Function that unstandardizes the columns of a real matrix Xs
	*/
	
	/* Arguments:
		 ----------
		 - Xs			--> real matrix whose columns are standardized
		 - scale	--> real rowvector containing the column scales
									(after centering)
		 - center	--> real rowvector containing the column centers
	*/
	
	/* Returned value:
		 ---------------
		 - X 			--> real matrix with columns that have been unstandardized
	*/

	real matrix cen, sc, Xunscaled, X
	real scalar n
	
	n = rows(Xs)
	cen = J(n, 1, center)
	sc = diag(scale)
	Xunscaled = Xs * sc
	X = Xunscaled + cen

	return(X)
}

real rowvector sd(real matrix X, |real scalar biased)
{
	/* Description:
		 ------------
		 Function that computes the standard deviation of columns of X
	*/
	
	/* Arguments:
		 ----------
		 - X			--> real matrix of data
		 - biased	--> (optional) real scalar; if not zero the biased sample variance
									is computed
	*/
	
	/* Returned value:
		 ---------------
		 - sc --> real rowvector containing the X column standard deviations
	*/
	
	real rowvector sc
	real scalar n, bias, factor
	
	n = rows(X)

	sc = sqrt(diagonal(variance(X)))'
	
	if (biased == . | biased == 0) {
		bias = 0
	}
	if (bias) {
		factor = sqrt((n - 1)/n)
		sc = sc:*factor
	}

	return(sc)
}

real colvector which(real matrix X, |string scalar minmax)
{
	/* Description:
		 ------------
		 Function that determines the index of the (first) minimum or maximum of a
		 numeric matrix
	*/
	
	/* Arguments:
		 ----------
		 - X			--> real matrix of data
		 - minmax	--> (optional) string scalar for minimum ("min") or maximum
									("max"); deafult is "min"
	*/
	
	/* Returned value:
		 ---------------
		 - index --> real vector containing the indexes at which min or max occur
	*/
	
	real scalar i, j, N, w
	real colvector index
	
	pragma unset i		// this is needed to avoid useless warning messages
	pragma unset w		// this is needed to avoid useless warning messages

	N = rows(X)
	
	if (args() == 1) {
		minmax = "min"
	}

	index = J(0, 1, .)
	for (j = 1; j <= N; j++) {
		if (minmax == "min") {
			minindex(X[j, .], 1, i, w)
		}
		else if (minmax == "max") {
			maxindex(X[j, .], 1, i, w)
		}
		else {
			_error(3000, "the 'minmax' argument value is unrecognized")
		}
		index = (index \ i)
	}
	
	return(index)
}

real rowvector cronbach(real matrix X, real matrix M)
{
	/* Description:
		 ------------
		 Function that computes the standardizes Cronbach's reliability coefficients
	*/
	
	/* Arguments:
		 ----------
		 - X			--> real matrix containing the latent variable scores
		 - M			--> measurement model's adjacency matrix
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
		alpha[p] = k/(1/r_bar + (k - 1))
		if (alpha[p] < 0) {
			alpha[p] = 0
		}
	}
	
	return(alpha)
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
											imputes
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
				mostfreq = xvals[which(xfreqs, "max")[1]]
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
	|string scalar categ, string scalar touse, real scalar noisily)
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
	*/
	
	/* Returned value:
		 ---------------
		 - void				--> the matrix of imputed variables is saved directly in the
											data table
	*/

	real matrix Ximp, Xsub, ur, xvals, xfreqs
	real scalar n, Q, N, R, V, v, i, skip, skip_b, j, q, mostfreq, D_idx, knew
	real colvector complete_idx, incomplete_idx, disti, tmp_idx, D
	real rowvector vars_nonmiss
	string rowvector vars_vec, vars_touse
	string scalar meas, todisp, spaces
	
	n = rows(X)												// number of observations
	Q = cols(X)
	
	if (noisily == .) {
		noisily = 1
	}
	meas = "L2"
	
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
			disti[D_idx] = max(disti)
			disti = (disti, tmp_idx)
			disti = sort(disti, 1)					// ties are randomly broken
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
						mostfreq = xvals[which(xfreqs, "max")[1]]
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

real matrix rdirichlet(real scalar n, real rowvector alpha)
{
	/* Description:
		 ------------
		 Function that generates random deviates from a Dirichlet distribution
	*/
	
	/* Arguments:
		 ----------
		 - n 			--> real scalar providing the number of deviates to draw
		 - alpha	--> real rowvector containing the Dirichlet parameters
	*/
	
	/* Returned value:
		 ---------------
		 - out		--> real colvector containing the generated Dirichlet deviates
	*/
	
	real scalar l
	real matrix out_mat
	
	l = length(alpha)
	
	out_mat = rgamma(n, 1, alpha, 1)
	
	return(out_mat :/ rowsum(out_mat))
}

real colvector euclidean_dist(real matrix X, real scalar idx)
{
	/* Description:
		 ------------
		 Function that computes the dissimilarities among a set of observations
	*/
	
	/* Arguments:
		 ----------
		 - X					--> real matrix containing the data
		 - idx 				--> real scalar providing the row index of X for which to
											calculate the dissimilarities
	*/
	
	/* Returned value:
		 ---------------
		 - D					--> vector of dissimilarities
	*/
	
	real colvector D
	real scalar n
	
	n = rows(X)
	D = sqrt(rowsum((X - J(n, 1, 1)*X[idx, .]):^2))
	
	return(D)
}

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
	left = rowsum(e2 * luinv(com), .)
	left = left / colsum(left)
	
	f2 = f :^ 2
	r2_mat = diag(r2)
	right = rowsum(f2 * luinv(r2_mat), .)
	right = right / colsum(right)
	
	cm = (N - 2)*sqrt(left :* right)
	
	return(cm)
}

real scalar plssem_rebus_gqi(real matrix ind, real matrix lat, real matrix e, ///
	real matrix f, real colvector cl)
{
	/* Description:
		 ------------
		 Function that computes the Group Quality Index (GQI) for a REBUS solution
	*/
	
	/* Arguments:
		 ----------
		 - ind		--> real matrix containing the indicator values stacked by REBUS
									group
		 - lat		--> real matrix containing the latent variable values stacked by
									REBUS group
		 - e 			--> real matrix containing the PLS-SEM measurement model residuals
		 - f			--> real matrix containing the PLS-SEM structural model residuals
		 - cl			--> real colvector containing the REBUS classification
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
	string scalar rebus_cl, real scalar numclass, real scalar maxit_reb,
	real scalar stop, |real scalar noisily)
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
											chosen
		 - rebus_cl		--> string scalar providing the variable with the REBUS
											classes
		 - numclass		--> real scalar providing the number of REBUS classes
		 - maxit_reb	--> real scalar providing the maximum number of REBUS
											iterations
		 - stop				--> real scalar providing the REBUS stopping criterion
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
		loads, y, y_local, block, x_hat, y_hat, out_res, inn_res
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
				printf("{txt}REBUS iteration 1")
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
				printf("{err}some variables during the REBUS calculations turned out to have zero variance\n")
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
				printf("{err}some variables during the REBUS calculations turned out to have zero variance\n")
				printf("{err}try reducing the number of classes\n")
				_error(1)
			}
			
			// Initialize the LVs
			Yinit = plssem_init_mat(Xsc,  M, ind, stdind, latents, touse_loc_name, ///
				rawsum, init)
			
			// Run the PLS algorithm
			localmodels[k, 1] = plssem_base(Xsc, Yinit, M, S, modes, latents, ///
				binary, tol, maxit, touse_loc_name, scheme, crit, structural, rawsum)
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
			y_local = Xrebsc_all * ow
			out_res = J(N, 0, .)
			for (p = 1; p <= P; p++) {
				block = selectindex(loads[., p] :!= .)
				x_hat = y_local[., p] * loads[block, p]'
				out_res = (out_res, (Xrebsc_all[., block] - x_hat))
			}
			y_hat = y_local * path[., selectindex(endo)]
			inn_res = (y_local[., selectindex(endo)] - y_hat)
			cm = (cm, plssem_rebus_cm(out_res, inn_res, loads, r2))
		}
		new_class = J(length(old_class), 1, .)
		new_class[selectindex(touse_vec), 1] = which(cm, "min")
		nchanged = sum(old_class :!= new_class)
		old_class = new_class
		
		class_freq = uniqrows(new_class, 1)[., 2]
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
	string scalar rebus_cl, real scalar numclass, |real scalar B,
	real scalar seed, real scalar noisily)
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
											chosen
		 - rebus_cl		--> string scalar providing the variable with the REBUS
											classes
		 - numclass		--> real scalar providing the number of REBUS classes
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
		loads, y_local, block, x_hat, y_hat, out_res, inn_res, indstd_st, ///
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
				printf("{err}some variables during the REBUS calculations turned out to have zero variance\n")
				printf("{err}try reducing the number of classes\n")
				_error(1)
			}
			
			// Standardize the MVs (if required)
			Xsc = plssem_scale_mat(X, touse_loc, scale)
			
			// Check that there are no zero-variance indicators
			if (anyof(sd(Xsc), 0)) {
				printf("{err}some variables during the REBUS calculations turned out to have zero variance\n")
				printf("{err}try reducing the number of classes\n")
				_error(1)
			}
			
			// Initialize the LVs
			Yinit = plssem_init_mat(Xsc,  M, ind, stdind, latents, touse_loc_name, ///
				rawsum, init)
			
			// Run the PLS algorithm
			localmodels[k, 1] = plssem_base(Xsc, Yinit, M, S, modes, latents, ///
				binary, tol, maxit, touse_loc_name, scheme, crit, structural, rawsum)

			Xreb = X[selectindex(touse_loc), .]
			Xrebsc = scale(Xreb)
			Yendo = st_data(., latents, touse_loc_name)
			Yendo = Yendo[., selectindex(endo)]
			
			ow = localmodels[k, 1].outer_weights
			path = localmodels[k, 1].path
			path = editmissing(path, 0)
			loads = localmodels[k, 1].loadings
			
			y_local = Xrebsc * ow
			out_res = J(rows(Xreb), 0, .)
			for (p = 1; p <= P; p++) {
				block = selectindex(loads[., p] :!= .)
				x_hat = y_local[., p] * loads[block, p]'
				out_res = (out_res, (Xrebsc[., block] - x_hat))
			}
			y_hat = y_local * path[., selectindex(endo)]
			inn_res = y_local[., selectindex(endo)] - y_hat
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

	real scalar P_ik, i, n, k, K, P, f_ik
	real colvector zeta
	
	n = rows(eta)
	P = cols(eta)
	K = length(B)
	P_ik = J(n, K, .)
	
	for (i = 1; i <= n; i++) {
		for (k = 1; k <= K; k++) {
			zeta = (I(P) - B[k, 1].mat)'*eta[i, .]' - Gamma[k, 1].mat*xi[i, .]'
			f_ik = lnmvnormalden(J(1, P, 0), Psi[k, 1].mat, zeta)
			P_ik[i, k] = rho[k]*exp(f_ik)
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

struct plssem_struct_fimix scalar plssem_fimix(real matrix Y, real matrix M,
	real matrix S, string scalar ind, string scalar stdind, string scalar latents,
	string scalar binary, real scalar tol, real scalar maxit, string scalar touse,
	string scalar scheme, string scalar crit, string scalar init,
	string scalar scale, real scalar structural, real scalar rawsum,
	real scalar K, real scalar maxit_fim, real scalar stop, |real scalar restart,
	real scalar seed, real scalar noisily, real scalar groups)
{
	/* Description:
		 ------------
		 Function that implements the FIMIX-PLS method using an EM algorithm
	*/
	
	/* Arguments:
		 ----------
		 - Y					--> real matrix containing the latent variables
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
											chosen
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
		rN0, rN_lte_5, delta_Q, old_Q, new_Q, new_ll, old_ll, iter_best, ///
		lnL_best, lnL_c_best, ll_restart, R, Q, lnL_1, skip
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
		new_Q = mindouble()/2 //plssem_fimix_ll_em(eta, xi, P_ik, rho, B, Gamma, Psi)
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
			old_ll = new_ll
			
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

void print_struct_matrix(struct plssem_struct_matrix colvector A)
{
	/* Description:
		 ------------
		 Function that prints the elements of a struct of plssem_struct_matrix
		 objects
	*/
	
	/* Arguments:
		 ----------
		 - A					--> struct of plssem_struct_matrix objects in a colvector
											layout
	*/
	
	/* Returned value:
		 ---------------
		 - void 			--> no object returned
	*/

	real scalar i, j, nr, nc
	
	nr = rows(A)
	nc = cols(A)
	
	for (i = 1; i <= nr; i++) {
		for (j = 1; j <= nc; j++) {
			A[i, j].mat
		}
	}
}

real scalar prod(real vector v)
{
	/* Description:
		 ------------
		 Function that computes the product of a vector's elements
	*/
	
	/* Arguments:
		 ----------
		 - v					--> vector of real elements
	*/
	
	/* Returned value:
		 ---------------
		 - res 				--> product of the v vector's elements
	*/

	real scalar i, nr, res
	
	nr = length(v)
	
	res = 1
	for (i = 1; i <= nr; i++) {
		res = res*v[i]
	}
	
	return(res)
}

void cleanup()
{
	/* Description:
		 ------------
		 Function that cleans up Mata from temporary objects 
	*/
	
	/* Arguments:
		 ----------
		 - void				--> no argument passed
	*/
	
	/* Returned value:
		 ---------------
		 - void 			--> no object returned
	*/

	string colvector names
	real scalar i
	
	names = direxternal("__*")

	for (i = 1; i <= rows(names); i++) {
		rmexternal(names[i])
	}
}

end

********************************************************************************

mata: mata mlib create lplssem, dir(PERSONAL) replace
// mata: mata describe *()
// mata: mata mlib add lplssem *()
mata: mata mlib add lplssem ///
	plssem_struct() ///
	plssem_scale() ///
	plssem_scale_mat() ///
	plssem_init() ///
	plssem_init_mat() ///
	plssem_base() ///
	plssem_lv() ///
	plssem_pval() ///
	plssem_pathtab() ///
	plssem_struct_boot() ///
	plssem_boot() ///
	plssem_boot_lm() ///
	plssem_boot_pm() ///
	plssem_boot_lv() ///
	plssem_boot_pv() ///
	plssem_reliability() ///
	plssem_vif() ///
	plssem_mga_perm() ///
	plssem_mga_perm_diff() ///
	plssem_mga_boot() ///
	plssem_mga_boot_diff() ///
	scale() ///
	unscale() ///
	sd() ///
	which() ///
	cronbach() ///
	meanimp() ///
	knnimp() ///
	euclidean_dist() ///
	rdirichlet() ///
	plssem_rebus_cm() ///
	plssem_rebus_gqi() ///
	plssem_struct_rebus() ///
	plssem_rebus() ///
	plssem_rebus_ptest() ///
	plssem_struct_matrix() ///
	plssem_struct_fimix() ///
	plssem_fimix_ll() ///
	plssem_fimix_ll_em() ///
	plssem_fimix_estep() ///
	plssem_fimix_ic() ///
	plssem_fimix() ///
	print_struct_matrix() ///
	prod() ///
	cleanup()
mata: mata mlib index
mata: mata describe using lplssem
