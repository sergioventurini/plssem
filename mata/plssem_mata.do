* set trace on

	/* TO DO:
		 1. 
	*/

	/* ISSUES:
		 1. 
	*/

*!plssem_mata version 0.3.0
*!Written 19Aug2017
*!Written by Sergio Venturini and Mehmet Mehmetoglu
*!The following code is distributed under GNU General Public License version 3 (GPL-3)

capture mata: mata drop plssem_struct() plssem_base() plssem_lv() ///
	plssem_pval() plssem_pathtab() plssem_struct_boot() plssem_boot() ///
	plssem_boot_lm() plssem_boot_pm() plssem_boot_lv() plssem_boot_pv() ///
	plssem_reliability() plssem_vif() scale() sd() which() cronbach() ///
	rebus_cm() rebus_gqi()

version 10
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
	real scalar reps						// number of bootstrap replicates
	real matrix loadings_reps		// matrix of outer loadings bootstrap replicates
	real matrix loadings_bs			// matrix of outer loadings bootstrap estimates
	real matrix loadings_v			// matrix of outer loadings bootstrap variances
	real matrix xloadings_reps	// matrix of cross loadings bootstrap replicates
	real matrix xloadings_bs		// matrix of cross loadings bootstrap estimates
	real matrix xloadings_v			// matrix of cross loadings bootstrap variances
	real matrix path_reps				// matrix of path coefficients bootstrap replicates
	real matrix path_bs					// matrix of path coefficients bootstrap estimates
	real matrix path_v					// matrix of path coefficients bootstrap variances
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
		fscores, T, beta, beta_v, Yexo, Yendo, Ypred, Whist, xlambda, lambda, Snew
	real scalar iter, delta, P, Q, i, j, p, q, n, minval, maxval, converged
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
					/* (this is the code as in the semPLS package) */
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
				Yhat = Ytilde
				Yhat = scale(Yhat)
				W = W * diag(1 :/ sd(Ytilde))
				w_new = rowsum(W)'
				Whist = (Whist \ w_new)

				// Step 5: convergence check
				if (crit == "relative") {
					delta = max(abs(w_old - w_new) :/ w_new)
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
		real scalar sse, sst, s2

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
	real scalar P, Q, p, q, stat, df
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
		 - seed				--> (optional) real scalar bootstrap seed
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
	
	real matrix X_b, Yinit_b, lambda, xlambda, beta, load_v, xload_v
	real vector bs_ind
	real scalar P, Q, n, b, skip, i, p
	string scalar todisp, spaces
	real rowvector isnotbinary
	string rowvector lvs_vec

	P = cols(M)   			// number of latent variables
	Q = rows(M)   			// number of manifest variables
	n = rows(X)
	
	// Boostrap resampling
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
		printf("{txt}Bootstrap replications (")
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
		bs_ind = runiformint(n, 1, 1, n)
		X_b = X[bs_ind, .]
		Yinit_b = Yinit[bs_ind, .]
		res = plssem_base(X_b, Yinit_b, M, S, mode, latents, binary, tol, ///
			maxit, touse, scheme, crit, structural, rawsum)
		if (structural) {
			beta[b, .] = vec(res.path)'
		}
		lambda[b, .] = rowsum(res.loadings)'
		xlambda[b, .] = vec(res.xloadings)'
	}
	
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
		 - lambda		--> real matrix whose rows contain the bootstrap replicates
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
		 - path				--> real matrix whose rows contain the bootstrap replicates
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
		if (!any(path[., j] :== .)) {
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
		 - lambda			--> real matrix whose rows contain the bootstrap replicates
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
		 - path					--> real matrix whose rows contain the bootstrap replicates
												of the (stacked) path coefficients matrix
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
		if (!any(path[., j] :== .)) {
			path_var[j] = variance(path[., j])
		}
	}
	
	pathcoef_v = rowshape(path_var, P)'
	
	return(pathcoef_v)
}

real matrix plssem_reliability(real matrix X, real matrix load, real matrix ave,
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
		 - ave				--> real matrix containing the AVE
		 - mode				--> real colvector containing the modes of the latent
											variables (each element must be either 0 ["reflective"]
											or 1 ["formative"])
	*/
	
	/* Returned value:
		 ---------------
		 - relcoef --> real matrix containing the (cross) loading variances
	*/

	real matrix relcoef, M
	real scalar P, p
	real rowvector nind, num, den

	P = cols(load)   			// number of latent variables
	relcoef = J(2, P, 0)

	M = (load :!= .)
	relcoef[1, .] = cronbach(X, M)
	nind = colsum(M)
	num = colsum(load) :^ 2
	den = num + nind :* (J(1, P, 1) - ave)
	relcoef[2, .] = num :/ den
	for (p = 1; p <= P; p++) {
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
		 - Xs --> real matrix with columns that have been centered and/or scaled
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
		sc = sc:*factor
	}
	Xc = X - cen
	sc_inv = luinv(sc)
	Xcs = Xc * sc_inv

	return(Xcs)
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
	real scalar n, factor
	
	n = rows(X)

	sc = sqrt(diagonal(variance(X)))'
	
	if (biased == .) {
		biased = 0
	}
	if (biased) {
		factor = sqrt((n - 1)/n)
		sc = sc:*factor
	}

	return(sc)
}

real vector which(real matrix X, |string scalar minmax)
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
	real vector index
	
	pragma unset i // this is needed to avoid useless warning messages
	pragma unset w // this is needed to avoid useless warning messages

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

	real matrix alpha, Mind
	real scalar P, n, p, k, l, r_bar
	
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
		alpha[p] = k*r_bar/(1 + (k - 1)*r_bar)
	}
	
	return(alpha)
}

real colvector rebus_cm(real matrix e, real matrix f, real matrix loads, ///
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
		 - loads	--> real matrix containing the PLS-SEM communalities
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
	
	e2 = e:^2
	com = diag(rowsum(loads:^2))
	left = rowsum(e2 * luinv(com), .)
	left = left / colsum(left)

	f2 = f:^2
	r2_mat = diag(r2)
	right = rowsum(f2 * luinv(r2_mat), .)
	right = right / colsum(right)

	cm = (N - 2)*sqrt(left :* right)
	
	return(cm)
}

real scalar rebus_gqi(real matrix ind, real matrix lat, real matrix e, ///
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
	// sqrt((Nk' * (left :* right))/N) // this is the calculation according
																		 // to Sanchez's plspm R package

	gqi = sqrt((Nk' * left)*(Nk' * right))/N	// this is the calculation according
																						// to Trinchera's PhD thesis (2007)

	return(gqi)
}
end

********************************************************************************

mata: mata mosave plssem_struct(), dir(PERSONAL) replace
mata: mata mosave plssem_base(), dir(PERSONAL) replace
mata: mata mosave plssem_lv(), dir(PERSONAL) replace
mata: mata mosave plssem_pval(), dir(PERSONAL) replace
mata: mata mosave plssem_pathtab(), dir(PERSONAL) replace
mata: mata mosave plssem_struct_boot(), dir(PERSONAL) replace
mata: mata mosave plssem_boot(), dir(PERSONAL) replace
mata: mata mosave plssem_boot_lm(), dir(PERSONAL) replace
mata: mata mosave plssem_boot_pm(), dir(PERSONAL) replace
mata: mata mosave plssem_boot_lv(), dir(PERSONAL) replace
mata: mata mosave plssem_boot_pv(), dir(PERSONAL) replace
mata: mata mosave plssem_reliability(), dir(PERSONAL) replace
mata: mata mosave plssem_vif(), dir(PERSONAL) replace
mata: mata mosave scale(), dir(PERSONAL) replace
mata: mata mosave sd(), dir(PERSONAL) replace
mata: mata mosave which(), dir(PERSONAL) replace
mata: mata mosave cronbach(), dir(PERSONAL) replace
mata: mata mosave rebus_cm(), dir(PERSONAL) replace
mata: mata mosave rebus_gqi(), dir(PERSONAL) replace
