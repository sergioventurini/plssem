*!plssem_gas version 0.3.0
*!Written 15Jan2019
*!Written by Sergio Venturini and Mehmet Mehmetoglu
*!The following code is distributed under GNU General Public License version 3 (GPL-3)

version 15.1

mata:

real scalar plssem_gas_fitness(struct plssem_struct_matrix colvector f,
	|struct plssem_struct_matrix colvector e)
{
	/* Description:
		 ------------
		 Function that computes the fitness function for the PLS-GAS method
	*/
	
	/* Arguments:
		 ----------
		 - f			--> struct plssem_struct_matrix colvector containing structural
									model's residuals
		 - e 			--> (optional) struct plssem_struct_matrix colvector containing
									measurement model's residuals (only for reflective constructs)
	*/
	
	/* Returned value:
		 ---------------
		 - F			--> real scalar providing the value of the fitness function
	*/
	
	real scalar K, k, F
	K = length(e)
	F = 0
	
	for (k = 1; k <= K; k++) {
		if (args() == 1) {
			// only the structural model's residuals are used
			F = F + sum(f[k, 1].mat :^ 2)
		}
		else {
			// both the measurement and structural model's residuals are used
			F = F + sum(e[k, 1].mat :^ 2) + sum(f[k, 1].mat :^ 2)
		}
	}
	
	return(F)
}

real scalar plssem_gas_validity_check(real colvector child, 
	real matrix chromo, real scalar ind_count, real scalar K, real matrix X)
{
	/* Description:
		 ------------
		 Function that checks whether an individual generated during a PLS-GAS
		 analysis is valid or not
	*/
	
	/* Arguments:
		 ----------
		 - child			--> real colvector containing the individual to test for
											validity
		 - chromo			--> real matrix containing the generated individuals so far
		 - ind_count	--> real scalar providing the number of individuals generated
											so far
		 - K					--> real scalar providing the number of clusters
		 - X					--> real matrix containing the observed manifest variables
	*/
	
	/* Returned value:
		 ---------------
		 - isvalid	--> real scalar equal to 1 or 0 depending on whether the
										individual is valid or not
	*/

	real scalar isvalid, k, h, l
	real matrix tab
	real colvector ref_ind

	isvalid = 1
	
	// Check that child includes exactly K clusters
	if (length(uniqrows(child)) < K) {
		isvalid = 0
		return(isvalid)
	}

	// Check that all clusters include at least 10 observations
	if (any(uniqrows(child, 1)[., 2] :< 10)) {
		isvalid = 0
		return(isvalid)
	}
	
	// Check that child is not a duplicate of any of the existing individuals
	for (l = 1; l <= ind_count; l++) {
		tab = J(K, K, .)
		ref_ind = chromo[., l]
		for (k = 1; k <= K; k++) {
			for (h = 1; h <= K; h++) {
				tab[h, k] = sum(child[selectindex(ref_ind :== h), 1] :== k)
			}
		}
		if (nonmissing(tab) == K) {
			isvalid = 0
			return(isvalid)
		}
	}

	// Check that in each cluster there are no zero-variance indicators
	for (k = 1; k <= K; k++) {
		if (anyof(sd(X[selectindex(child :== k), .]), 0)) {
			isvalid = 0
			return(isvalid)
		}
	}
	
	return(isvalid)
}

struct plssem_struct_gas scalar plssem_gas(real matrix X, real matrix M, 
	real matrix S, string scalar ind, string scalar stdind,
	string scalar latents, string scalar binary, real scalar tol,
	real scalar maxit, string scalar touse, string scalar scheme,
	string scalar crit, string scalar init, string scalar scale,
	real scalar structural, real scalar rawsum, real scalar consistent,
	real scalar K, real scalar I, real scalar G, real scalar p_m, real scalar p_t,
	real scalar maxit_gas, |real scalar seed, real scalar noisily)
{
	/* Description:
		 ------------
		 Function that implements the PLS-GAS method based on a genetic algortihm
		 approach
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
		 - K					--> real scalar providing the number of classes
		 - I					--> real scalar providing the number of individuals in the
											current population
		 - G					--> real scalar providing the number of generations
		 - p_m				--> real scalar providing the probability of mutation
		 - p_t				--> real scalar providing the probability of transformation
		 - maxit_gas	--> real scalar providing the maximum number of iterations in
											the PLS-GAS second stage
		 - seed				--> (optional) real scalar providing the seed
		 - noisily		--> (optional) real scalar; if not 0, it prints the EM
											iterations
	*/
	
	/* Returned value:
		 ---------------
		 - res_gas	--> scalar plssem_struct_gas structure containing the results
	*/

	struct plssem_struct_gas scalar res_gas
	struct plssem_struct_matrix scalar mat
	struct plssem_struct_matrix colvector e, f
  struct plssem_struct scalar res
  struct plssem_struct colvector localmodels

	real matrix Xsc, Yinit, chromo, chromo_new, Xgas, ow, path, loads, ///
		Xgas_all, Xgassc_all, y_loc, out_res, block, x_hat, y_hat, inn_res
	real scalar n, N, i, k, j, P, p, skip, a, b, gen_count, ind_count, ///
		mut_count, isvalid, t, best_k, nchanged, F_tmp, F_best, j_nomiss
	real colvector touse_vec, modes, touse_loc, best, parent, child, best_tmp, ///
		alloc_best
	real rowvector F, endo, rank, rankF, p_i, F_k
	string scalar todisp, spaces, touse_loc_name
	
	if (seed != .) {
		rseed(seed)
	}
	if (noisily == .) {
		noisily = 1
	}
	
	touse_vec = st_data(., touse)
	P = cols(S)
	n = rows(X)
	N = sum(touse_vec)
	endo = (colsum(S)	:> 0)			 // indicators for the endogenous latent variables

  localmodels = J(K, 1, res)
  e = J(K, 1, mat)
  f = J(K, 1, mat)
	F = J(1, I, .)
  modes = J(P, 1, 0)

	(void) st_addvar("byte", touse_loc_name = st_tempname())

	if (noisily) {
		printf("{txt}\n")
		printf("{txt}Computing PLS-GAS solution...\n")
	}

	// 1st STAGE
	// ---------
	if (noisily) {
		printf("{txt}\n")
		printf("{txt}First stage --> Creating generations (")
		printf("{res}%f", G)
		printf("{txt})\n")
		printf("{hline 4}{c +}{hline 3} 1 ")
		printf("{hline 3}{c +}{hline 3} 2 ")
		printf("{hline 3}{c +}{hline 3} 3 ")
		printf("{hline 3}{c +}{hline 3} 4 ")
		printf("{hline 3}{c +}{hline 3} 5\n")
	}
	skip = 0

	// Step 1: Randomly select a starting population of size I individuals
	chromo = J(n, I, .)
	chromo[selectindex(touse_vec), .] = runiformint(N, I, 1, K)

	gen_count = 0
	while (gen_count < G) {
		if (noisily) {
			if (mod(gen_count + 1, 50) == 0) {
				todisp = strtrim(strofreal(gen_count + 1, "%9.0f"))
				if (strlen(todisp) < strlen(strofreal(G))) {
					skip = strlen(strofreal(G)) - strlen(todisp)
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
		
		// Step 2: Determine the rank fitness of each individual
		for (i = 1; i <= I; i++) {
			for (k = 1; k <= K; k++) {
				touse_loc = touse_vec :* (chromo[., i] :== k)
				st_store(., touse_loc_name, touse_loc)
				
				// Check that there are no zero-variance indicators
				if (anyof(sd(X[selectindex(touse_loc), .]), 0)) {
					printf("{err}some indicators during the PLS-GAS calculations turned out to have zero variance\n")
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
					printf("{err}some indicators during the PLS-GAS calculations turned out to have zero variance\n")
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
				e[k, 1].mat = localmodels[k, 1].e
				f[k, 1].mat = localmodels[k, 1].f
			}
			
			// Fitness calculation for individual i
			F[1, i] = plssem_gas_fitness(f) //, e) // only structural model's resid
		}
		rank = invorder(order(F', 1))'
		rankF = J(1, I, I + 1) - rank
		p_i = rankF/(I*(I + 1)/2)
		
		// Step 3: Create a new generation
		chromo_new = J(n, I, .)
		ind_count = 1
		//  - Step 3.1: Select the fittest individual as a member of the next
		//						  generation (elitist selection)
		best = chromo[., which(rankF, "max")]
		chromo_new[., ind_count] = best
		while (ind_count < I) {
			//  - Step 3.2: Randomly select (using rank fitness) a parent from the old
			//						  population
			parent = chromo[., rdiscrete(1, 1, p_i)]
			//  - Step 3.3: Draw a random number
			a = runiform(1, 1)
			child = parent
			//  - Step 3.4: Generate individuals for the new population
			if (a <= p_m) {
				// mutation
				mut_count = 0
				for (j = 1; j <= n; j++) {
					if (touse_vec[j]) {   // skip missing values
						b = runiform(1, 1)
						if (b <= p_t) {
							child[j, 1] = runiformint(1, 1, 1, K)
							mut_count++
						}
					}
				}
				if (mut_count == 0) {
					j_nomiss = runiformint(1, 1, 1, n)
					while (!touse_vec[j_nomiss]) {   // skip missing values
						j_nomiss = runiformint(1, 1, 1, n)
					}
					child[j_nomiss, 1] = runiformint(1, 1, 1, K)
				}
			}
			else {
				// reproduction (do nothing)
			}

			//  - Steps 3.5-3.6: Test the individual for validity
			if (ind_count > 1) {
				isvalid = plssem_gas_validity_check(child, chromo_new, ind_count, K, X)
			}
			else {
				isvalid = 1
			}

			// if the individual is valid, it is included in the new generation
			if (isvalid) {
				ind_count++
				chromo_new[., ind_count] = child
			}
		}

		chromo = chromo_new
		gen_count++
	}
	if (noisily) printf("{txt}\n")

	// 2nd STAGE
	// ---------
	// Select the fittest individual from 1st stage
	for (i = 1; i <= I; i++) {
		for (k = 1; k <= K; k++) {
			touse_loc = touse_vec :* (chromo[., i] :== k)
			st_store(., touse_loc_name, touse_loc)
			
			// Standardize the MVs (if required)
			if (scale == "") {
				Xsc = scale(X[selectindex(touse_loc), .])
			}
			else {
				Xsc = X[selectindex(touse_loc), .]
			}
			
			// Initialize the LVs
			plssem_scale(Xsc, stdind, touse_loc_name, scale)
			Yinit = plssem_init_mat(Xsc, M, ind, stdind, latents, touse_loc_name, ///
				rawsum, init)
			
			// Run the PLS algorithm
			localmodels[k, 1] = plssem_base(Xsc, Yinit, M, S, modes, latents, ///
				binary, tol, maxit, touse_loc_name, scheme, crit, structural, ///
				rawsum, consistent)
			e[k, 1].mat = localmodels[k, 1].e
			f[k, 1].mat = localmodels[k, 1].f
		}
		
		// Fitness calculation for individual i
		F[1, i] = plssem_gas_fitness(f) //, e) // only structural model's resid
	}
	rank = invorder(order(F', 1))'
	rankF = J(1, I, I + 1) - rank
	best = chromo[., which(rankF, "max")]
	
	// Check if reassigning each observation to a different group does improve
	// the fitness
	if (noisily) {
		printf("{txt}\n")
		printf("{txt}Second stage --> Local improvement (hill-climbing) (")
		printf("{res}%f", maxit_gas)
		printf("{txt})\n")
		printf("{hline 4}{c +}{hline 3} 1 ")
		printf("{hline 3}{c +}{hline 3} 2 ")
		printf("{hline 3}{c +}{hline 3} 3 ")
		printf("{hline 3}{c +}{hline 3} 4 ")
		printf("{hline 3}{c +}{hline 3} 5\n")
	}
	skip = 0

	best_tmp = best
	t = 1
	F_best = plssem_gas_fitness(f) //, e) // only structural model's resid
	alloc_best = best
	
	while ((t <= maxit_gas) & (nchanged > 0)) {
		if (noisily) {
			if (mod(t, 50) == 0) {
				todisp = strtrim(strofreal(t, "%9.0f"))
				if (strlen(todisp) < strlen(strofreal(maxit_gas))) {
					skip = strlen(strofreal(maxit_gas)) - strlen(todisp)
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
		
		nchanged = 0
		for (j = 1; j <= n; j++) {
			if (touse_vec[j]) {   // skip missing values
				F_k = J(1, K, .)
				
				// the computation for each observation uses the PLS-SEM results of the
				// segment to which the observation has been assigned
				ow = localmodels[best[j, 1], 1].outer_weights
				path = localmodels[best[j, 1], 1].path
				path = editmissing(path, 0)
				loads = localmodels[best[j, 1], 1].loadings

				for (k = 1; k <= K; k++) {
					best_tmp[j, 1] = k
					touse_loc = touse_vec :* (best_tmp :== k)
					Xgas = X[selectindex(touse_loc), .]
					
					Xgas_all = X[selectindex(touse_vec), .]
					Xgassc_all = scale(Xgas_all, 0, sd(Xgas), mean(Xgas))
					y_loc = Xgassc_all * ow
					out_res = J(N, 0, .)
					for (p = 1; p <= P; p++) {
						block = selectindex(loads[., p] :!= .)
						x_hat = y_loc[., p] * loads[block, p]'
						out_res = (out_res, (Xgassc_all[., block] - x_hat))
					}
					y_hat = y_loc * path[., selectindex(endo)]
					inn_res = (y_loc[., selectindex(endo)] - y_hat)
					F_k[1, k] = sum(out_res :^ 2) + sum(inn_res :^ 2)
				}
				best_k = which(F_k, "min")
				if (best_k != best[j, 1]) {
					nchanged++
					best[j, 1] = best_k
					best_tmp = best
					
					for (k = 1; k <= K; k++) {
						touse_loc = touse_vec :* (best[., 1] :== k)
						st_store(., touse_loc_name, touse_loc)
						
						// Check that there are no zero-variance indicators
						if (anyof(sd(X[selectindex(touse_loc), .]), 0)) {
							if (noisily) {
								printf("{txt}\n\n")
							}
							printf("{err}some indicators during the PLS-GAS calculations turned out to have zero variance\n")
							printf("{err}try reducing the number of classes\n")
							//_error(1)

							printf("{err}best results found so far returned\n")
							
							// Assigning results to structure's members
							res_gas.niter = (t - 1)
							res_gas.gas_class = alloc_best
							res_gas.touse_vec = touse_vec
							res_gas.nclass = K
							res_gas.popsize = I
							res_gas.ngen = G
							
							return(res_gas)
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
							if (noisily) {
								printf("{txt}\n\n")
							}
							printf("{err}some indicators during the PLS-GAS calculations turned out to have zero variance\n")
							printf("{err}try reducing the number of classes\n")
							//_error(1)

							printf("{err}best results found so far returned\n")
							
							// Assigning results to structure's members
							res_gas.niter = (t - 1)
							res_gas.gas_class = alloc_best
							res_gas.touse_vec = touse_vec
							res_gas.nclass = K
							res_gas.popsize = I
							res_gas.ngen = G
							
							return(res_gas)
						}
						
						// Initialize the LVs
						plssem_scale(Xsc, stdind, touse_loc_name, scale)
						Yinit = plssem_init_mat(Xsc, M, ind, stdind, latents, touse_loc_name, ///
							rawsum, init)
						
						// Run the PLS algorithm
						localmodels[k, 1] = plssem_base(Xsc, Yinit, M, S, modes, latents, ///
							binary, tol, maxit, touse_loc_name, scheme, crit, structural, ///
							rawsum, consistent)
						e[k, 1].mat = localmodels[k, 1].e
						f[k, 1].mat = localmodels[k, 1].f
					}
					
					F_tmp = plssem_gas_fitness(f) //, e) // only structural model's resid
					if (F_tmp < F_best) {
						F_best = F_tmp
						alloc_best = best
					}
				}
			}
		}
		
		t++
	}
	
	// Assigning results to structure's members
	res_gas.niter = (t - 1)
	res_gas.gas_class = alloc_best
	res_gas.touse_vec = touse_vec
	res_gas.nclass = K
	res_gas.popsize = I
	res_gas.ngen = G
	
	return(res_gas)
}

end
