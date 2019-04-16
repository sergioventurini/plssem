*!plssem_mga version 0.3.0
*!Written 15Jan2019
*!Written by Sergio Venturini and Mehmet Mehmetoglu
*!The following code is distributed under GNU General Public License version 3 (GPL-3)

version 15.1

mata:

class AssociativeArray scalar plssem_mga_perm(real matrix X, real matrix Yinit,
	real matrix M, real matrix S, real colvector mode, string scalar latents,
	string scalar binary, real scalar tol, real scalar maxit, string scalar touse,
	string scalar scheme, string scalar crit, real scalar structural,
	real scalar rawsum, real scalar consistent, real colvector group,
	|real scalar B, real scalar seed, real scalar noisily)
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
											chosen, and 0 otherwise
		 - consistent	--> real scalar equal to 1 if the consistent version (PLSc) is
											required, and 0 otherwise
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
				structural, rawsum, consistent)
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
	real scalar rawsum, real scalar consistent, real colvector group,
	|real scalar B, real scalar seed, real scalar noisily)
{
	/* Description:
		 ------------
		 Function that performs the bootstrap version of the PLS-MGA algorithm
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
				structural, rawsum, consistent)
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

end
