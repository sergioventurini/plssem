*!plssem_estat version 0.3.1
*!Written 03Nov2018
*!Written by Sergio Venturini and Mehmet Mehmetoglu
*!The following code is distributed under GNU General Public License version 3 (GPL-3)

program plssem_estat, rclass
	version 15.1
	gettoken subcmd rest : 0 , parse(", ")
	local lsubcmd = length("`subcmd'")
	
	if ("`subcmd'" == substr("indirect", 1, max(2, `lsubcmd'))) {
		indirect `rest'
	}
	else if ("`subcmd'" == substr("total", 1, max(2, `lsubcmd'))) {
		total `rest'
	}
	else if ("`subcmd'" == substr("vif", 1, max(2, `lsubcmd'))) {
		plssem_vif `rest'
	}
	else if ("`subcmd'" == substr("unobshet", 1, max(2, `lsubcmd'))) {
		unobshet `rest'
	}
	else if ("`subcmd'" == substr("mediate", 1, max(2, `lsubcmd'))) {
		mediate `rest'
	}
	else {
		// estat_default `0'
		display as error "the `subcmd' postestimation command is not implemented for plssem"
		exit
	}

	return add
end

program indirect, rclass
	version 15.1
	syntax , Effects(string) [ Boot(numlist integer >0 min=1 max=1) ///
		Seed(numlist max=1) Level(real 0.95) DIGits(integer 3) ]
	
	/* Options:
	   --------
		 effects(string)						--> list of indirect effects
		 boot(numlist integer >0 min=1 max=1)
																--> bootstrap estimation (# of repetions;
																		default 50)
		 seed(numlist max=1)				--> bootstrap seed number
		 level(real 0.95)						--> confidence level (default 0.95)
		 digits(integer 3)					--> number of digits to display (default 3)
	 */
	 
	 /* Description:
			------------
			This postestimation command provides the estimates for the indirect
			effects mediated by only one LV.
	 */
	
	if ("`effects'" == "") {
		display as error "effects() option must be provided"
		exit
	}
	local struct "structural"
	local props = e(properties)
	local isstruct : list struct in props
	if (!`isstruct') {
		display as error "the fitted plssem model includes only the measurement part"
		exit
	}
	
	tempvar __touse__
	quietly generate `__touse__' = e(sample)
	local reg3eqs = e(struct_eqs)
	if (`level' <= 0 | `level' >= 1) {
		display as error "confidence level must be in the range (0, 1)"
		exit
	}
	tempname normperc alpha_cl
	scalar `alpha_cl' = 1 - ((1 - `level')/2)
	scalar `normperc' = invnormal(`alpha_cl')
	
	/* Parse indirect statements (dep1 med1 indep1, dep2 med2 indep2, etc.) */
	local effects : list clean effects
	tokenize `"`effects'"', parse(",")
	
	local num_effects = 0
	local tok_i = 1
	while ("``tok_i''" != "") {
		local 0 ``tok_i''
		syntax varlist(min=3 max=3)
		
		local ind_dvar `ind_dvar' `: word 1 of `varlist''
		local ind_mvar `ind_mvar' `: word 2 of `varlist''
		local ind_ivar `ind_ivar' `: word 3 of `varlist''

		local ++num_effects

		local tok_i = `tok_i' + 2
	}
	if (`num_effects' > 5) {
		display as error "a maximum of 5 indirect effects are allowed"
		error 103
	}
	/* End parsing indirect statements */
	
	tempname sobel_se sobel_z sobel_pv indirect sobel_lci sobel_uci
	tempname indmeas reg3coef reg3var
	tempname dommat dommat2 moimat moimat2
	if ("`boot'" == "") {
		matrix `indmeas' = J(6, `num_effects', .)
	}
	else {
		matrix `indmeas' = J(10, `num_effects', .)
	}
	
	tempname ehold
	_estimates hold `ehold'

	if ("`boot'" != "") {
		display as text "Bootstrapping indirect effects..."
	}
	
	capture {
		forvalues ll = 1/`num_effects' {
			local dep `: word `ll' of `ind_dvar''
			local med `: word `ll' of `ind_mvar''
			local indep `: word `ll' of `ind_ivar''
			local indmeas_colnames "`indmeas_colnames' `dep'_`med'_`indep'"
			
			if ("`dep'" != "" & "`med'" != "" & "`indep'" != "") {
				local doi `dep':`indep'
				local moi `med':`indep'
				local dom `dep':`med'

				if ("`boot'" == "") {
					quietly reg3 `reg3eqs' if `__touse__', mvreg corr(independent)
					matrix `reg3coef' = e(b)
					matrix `reg3var' = e(V)

					tempname moicoef moivar domcoef domvar
					matrix `moimat' = `reg3coef'[1, "`moi'"]
					scalar `moicoef' = `moimat'[1, 1]
					matrix `moimat2' = `reg3var'["`moi'", "`moi'"]
					scalar `moivar' = `moimat2'[1, 1]
					matrix `dommat' = `reg3coef'[1, "`dom'"]
					scalar `domcoef' = `dommat'[1, 1]
					matrix `dommat2' = `reg3var'["`dom'", "`dom'"]
					scalar `domvar' = `dommat2'[1, 1]

					scalar `indirect' = `moicoef'*`domcoef'
					scalar `sobel_se' = sqrt(((`domcoef')^2)*`moivar' + ((`moicoef')^2)*`domvar')
					scalar `sobel_z' = `indirect'/`sobel_se'
					scalar `sobel_pv' =  2*(1 - normal(abs(`sobel_z')))	
					scalar `sobel_lci' = `indirect' - `normperc'*`sobel_se'
					scalar `sobel_uci' = `indirect' + `normperc'*`sobel_se'
					
					matrix `indmeas'[1, `ll'] = `indirect'
					matrix `indmeas'[2, `ll'] = `sobel_se'
					matrix `indmeas'[3, `ll'] = `sobel_z'
					matrix `indmeas'[4, `ll'] = `sobel_pv'
					matrix `indmeas'[5, `ll'] = `sobel_lci'
					matrix `indmeas'[6, `ll'] = `sobel_uci'
				}
				else {
					tempname reg3ciP reg3ciBC
					quietly bootstrap indeff=(_b[`dom']*_b[`moi']), reps(`boot') ///
						seed(`seed'): reg3 `reg3eqs' if `__touse__', mvreg corr(independent)
					quietly estat bootstrap, all
					matrix `reg3coef' = e(b)
					matrix `reg3var' = e(V)
					matrix `reg3ciP' = e(ci_percentile)
					matrix `reg3ciBC' = e(ci_bc)

					tempname iecoef ievar perc_lci perc_uci bc_lci bc_uci
					matrix `moimat' = `reg3coef'[1, 1]
					scalar `iecoef' = `moimat'[1, 1]
					matrix `moimat2' = `reg3var'[1, 1]
					scalar `ievar' = `moimat2'[1, 1]

					scalar `indirect' = `iecoef'
					scalar `sobel_se' = sqrt(`ievar')
					scalar `sobel_z' = `indirect'/`sobel_se'
					scalar `sobel_pv' =  2*(1 - normal(abs(`sobel_z')))	
					scalar `sobel_lci' = `indirect' - `normperc'*`sobel_se'
					scalar `sobel_uci' = `indirect' + `normperc'*`sobel_se'
					scalar `perc_lci' = `reg3ciP'[1, 1]
					scalar `perc_uci' = `reg3ciP'[2, 1]
					scalar `bc_lci' = `reg3ciBC'[1, 1]
					scalar `bc_uci' = `reg3ciBC'[2, 1]
					
					matrix `indmeas'[1, `ll'] = `indirect'
					matrix `indmeas'[2, `ll'] = `sobel_se'
					matrix `indmeas'[3, `ll'] = `sobel_z'
					matrix `indmeas'[4, `ll'] = `sobel_pv'
					matrix `indmeas'[5, `ll'] = `sobel_lci'
					matrix `indmeas'[6, `ll'] = `sobel_uci'
					matrix `indmeas'[7, `ll'] = `perc_lci'
					matrix `indmeas'[8, `ll'] = `perc_uci'
					matrix `indmeas'[9, `ll'] = `bc_lci'
					matrix `indmeas'[10, `ll'] = `bc_uci'
				}
			}
			else {
				display as error "provide a single dependent, independent and mediator variable for each indirect effect"
				exit
			}
		}
		if ("`boot'" == "") {
			matrix rownames `indmeas' = "Indirect effect" "Standard error" ///
				"Z statistic" "P-value" "Lower CI" "Upper CI"
		}
		else {
			matrix rownames `indmeas' = "Indirect effect" "Standard error" ///
				"Z statistic" "P-value" "Lower CI (N)" "Upper CI (N)" ///
				"Lower CI (P)" "Upper CI (P)" "Lower CI (BC)" "Upper CI (BC)"
		}
		matrix colnames `indmeas' = `indmeas_colnames'
	} // end of -capture-
	local rc = _rc
	_estimates unhold `ehold'
	error `rc'
	
	/* Display results */
	local ind_line = 5
	if ("`boot'" != "") {
		local ind_line = `ind_line' + 2
	}
	mktable_indirect, matrix(`indmeas') digits(`digits') boot(`boot') ///
		title("Significance testing of (standardized) indirect effects") ///
		firstcolname("Statistics") firstcolwidth(24) colwidth(18) ///
		hlines(`ind_line') level(`level') //novlines

	/* Return values */
	return matrix indirect `indmeas'
end

program total, rclass
	version 15.1
	syntax [ , DIGits(integer 3) Plot ]
	
	/* Options:
	   --------
		 digits(integer 3)		--> number of digits to display (default 3)
		 plot									--> bar plot of the effects
	 */
	
	local struct "structural"
	local props = e(properties)
	local isstruct : list struct in props
	if (!`isstruct') {
		display as error "the fitted plssem model includes only the measurement part"
		exit
	}
	if (`digits' < 0) {
		display as error "number of digits must be a nonnegative integer"
		exit
	}
		
	tempname B path direct indirect total
	matrix `B' = e(pathcoef)
	local nlv = colsof(`B')
	local nlv_m1 = `nlv' - 1
	matrix `indirect' = J(`nlv', `nlv', 0)
	if (`nlv' <= 2) {
		matrix `total' = `B'
	}
	else {
		matrix `path' = `B'
		forvalues k = 2/`nlv_m1' {
			matrix `path' = `path' * `B'
			matrix `indirect' = `indirect' + `path'
		}
	}
	matrix `total' = `B' + `indirect'
	
	tempname alleffects alllen maxlen alleffects_new
	local nall = `nlv'*(`nlv' - 1)
	matrix `alleffects' = J(`nall', 3, .)
	matrix `alllen' = J(`nall', 1, .)
	local lvnames : colnames `B'
	local i = 1
	foreach var_f in `lvnames' {
		local lvnames_rest : list lvnames - var_f
		foreach var_t in `lvnames_rest' {
			local allnames "`allnames' `var_f'->`var_t'"
			local rowi : list posof "`var_f'" in lvnames
			local coli : list posof "`var_t'" in lvnames
			if (`B'[`rowi', `coli'] == 0) {
				matrix `alleffects'[`i', 1] = .
			}
			else {
				matrix `alleffects'[`i', 1] = `B'[`rowi', `coli']
			}
			if (`indirect'[`rowi', `coli'] == 0) {
				matrix `alleffects'[`i', 2] = .
			}
			else {
				matrix `alleffects'[`i', 2] = `indirect'[`rowi', `coli']
			}
			if (`total'[`rowi', `coli'] == 0) {
				matrix `alleffects'[`i', 3] = .
			}
			else {
				matrix `alleffects'[`i', 3] = `total'[`rowi', `coli']
			}
			local ++i
		}
	}
	local firstnonmiss = 0
	forvalues i = 1/`nall' {
		if (!missing(`alleffects'[`i', 3])) {
			if (`firstnonmiss' == 0) {
				matrix `alleffects_new' = `alleffects'[`i', 1..3]
				local ++firstnonmiss
			}
			else {
				matrix `alleffects_new' = (`alleffects_new' \ `alleffects'[`i', 1..3])
			}
			local efftoleave : word `i' of `allnames'
			local allnames_new "`allnames_new' `efftoleave'"
		}
	}
	matrix rownames `alleffects_new' = `allnames_new'
	matrix colnames `alleffects_new' = "Direct" "Indirect" "Total"

	local nall = rowsof(`alleffects_new')
	forvalues i = 1/`nall' {
		local nm : word `i' of `allnames_new'
		local nm_len : strlen local nm
		matrix `alllen'[`i', 1] = `nm_len' + 2
		local nm_new = subinstr("`nm'", "->", " -> ", 1)
		local lblbar "`lblbar' `i' "`nm_new'""
	}
	mata: st_numscalar("`maxlen'", max(st_matrix("`alllen'")))
	local `maxlen' = `maxlen' + 2
	mktable, matrix(`alleffects_new') digits(`digits') firstcolname("Effect") ///
		title("Direct, Indirect (overall) and Total Effects") ///
		firstcolwidth(``maxlen'') colwidth(12) hlines(`nall') novlines total
	
	if ("`plot'" != "") {
		tempvar __id__ // __touse__
		// quietly generate `__touse__' = e(sample)
		quietly generate `__id__' = _n // if `__touse__'
		quietly svmat `alleffects_new', names(col)
		graph bar (asis) Direct Indirect, over(`__id__', relabel(`lblbar') ///
			label(angle(90) labsize(small))) stack nofill legend(position(12) ///
			region(style(none))) ylabel(, nogrid) yline(0, lwidth(thin)) ///
			title("Direct and Indirect (overall) Effects") scheme(sj)
		quietly drop Direct Indirect Total
	}

	/* Return values */
	return matrix total `alleffects_new'
end

program plssem_vif, rclass
	version 15.1
	syntax [ , DIGits(integer 3) ]
	
	/* Options:
	   --------
		 digits(integer 3)					--> number of digits to display (default 3)
	 */
	
	local struct "structural"
	local boot "bootstrap"
	local props = e(properties)
	local isstruct : list struct in props
	local isboot : list boot in props
	if (!`isstruct') {
		display as error "the fitted plssem model includes only the measurement part"
		exit
	}
	else if (`isboot') {
		display as error "estat vif is not allowed after bootstrap"
		exit
	}
	else {
		tempname strvif
		matrix `strvif' = e(struct_vif)
		local hline_path = rowsof(`strvif')
		mktable, matrix(`strvif') digits(`digits') firstcolname("Variable") ///
			title("Structural model - Multicollinearity check (VIFs)") ///
			firstcolwidth(14) colwidth(12) hlines(`hline_path') novlines
	}
	
	/* Return values */
	return matrix strvif `strvif'
end

program unobshet
	version 15.1
	syntax [ , Method(string) Numclass(numlist integer >=1 max=1) ///
		POPsize(numlist integer >=1 max=1) MAXCLass(integer 20) Dendrogram ///
		MAXITer(integer 50) NUMGen(integer 1000) PMut(real 0.3) ///
		PTransf(real 0.1) Stop(numlist >0 max=1) Test ///
		Reps(numlist integer >1 max=1) REStart(numlist integer >=1 max=1) ///
		SEed(numlist max=1) Plot name(string) DIGits(integer 3) ///
		GRoups(numlist integer >=1 <=20 min=1 sort) MAXITGas(integer 30) ]

	/* Options:
	   --------
		 method(string)									--> method to use for assessing unobserved
																				heterogeneity; available methods are
																				'rebus', 'fimix' and 'gas'. default is
																				'rebus'.
		 numclass(integer integer >=1
							max=1))								--> number of classes to use; if empty, it
																				is chosen automatically using a Ward
																				hierarchical algorithm
		 popsize(integer integer >=1
						 max=1))								--> number of individual in each population;
																				this is an option to use with
																				method(gas)
		 maxclass(integer 20)						--> maximum number of classes to test in
																				choosing automatically the number of
																				classes (default 20)
		 dendrogram											--> display the dendrogram for the Ward
																				cluster analysis
		 maxiter(integer 50)						--> maximum number of iterations (default
																				50)
		 numgen(integer 1000)						--> number of generations (default 1000)
		 pmut(real 0.3)									--> probability of mutation (default 0.3)
		 ptransf(real 0.1)							--> probability of transformation (default
																				0.1)
		 stop(real >0)									--> stopping criterion
		 test														--> permutation test
		 reps(numlist integer >1 max=1)	--> number of permutation test replications
																				(default is 50)
		 restart(numlist integer >=1
						 max=1)									--> number of repetitions of the EM
																				algorithm in FIMIX-PLS (default is 10)
		 seed(numlist max=1)						--> permutation test seed number
		 plot														--> plot of the empirical distribution
																				for the permutation test statistic
		 name(string)										--> variable name where to store the final
																				rebus classification
		 digits(integer 3)							--> number of digits to display (default 3)
		 groups(numlist integer >=1
						<=20 min=1 sort)				--> list of group numbers for which to
																				compare the corresponding fit indices
																				using FIMIX-PLS
		 maxitgas(integer 30)						--> maximum number of PLS-GAS second stage
																				iterations (default 30)
	 */
	 
	 /* Description:
			------------
			This postestimation command provides various methods to assess the
			presence of unobserved heterogeneity.
			
			Currently it implements the REBUS-PLS and FIMIX-PLS approaches.
	 */
	
	if ("`method'" == "") {
		local method "rebus"
	}
	
	if ("`method'" == "rebus") {
		REBUS `0'
	}
	else if ("`method'" == "fimix") {
		FIMIX `0'
	}
	else if ("`method'" == "gas") {
		GAS `0'
	}
	else {
		display as error "currently 'estat unobshet' implements only the " _continue
		display as error "REBUS-PLS and FIMIX-PLS methods"
		error 198
	}
	
	/* Clean up */
	capture mata: cleanup()
end

program REBUS, rclass
	version 15.1
	syntax [ , Method(string) Numclass(numlist integer >=1 max=1) ///
		MAXCLass(integer 20) Dendrogram MAXITer(integer 50) Stop(real 0.005) ///
		Test Reps(numlist integer >1 max=1) SEed(numlist max=1) Plot ///
		name(string) DIGits(integer 3) ]

	/* Options:
	   --------
		 method(string)									--> method to use for assessing unobserved
																				heterogeneity; available methods are
																				'rebus', 'fimix' and 'gas'. default is
																				'rebus'.
		 numclass(integer integer >=1
							max=1))								--> number of classes to use; if empty, it
																				is chosen automatically using a Ward
																				hierarchical algorithm
		 maxclass(integer 20)						--> maximum number of classes to test in
																				choosing automatically the number of
																				classes (default 20)
		 dendrogram											--> display the dendrogram for the Ward
																				cluster analysis
		 maxiter(integer 50)						--> maximum number of iterations (default
																				50)
		 stop(real 0.005)								--> stopping criterion (default 0.005)
		 test														--> permutation test
		 reps(numlist integer >1 max=1)	--> number of permutation test replications
																				(default is 50)
		 seed(numlist max=1)						--> permutation test seed number
		 plot														--> plot of the empirical distribution
																				for the permutation test statistic
		 name(string)										--> variable name where to store the final
																				REBUS-PLS classification
		 digits(integer 3)							--> number of digits to display (default 3)
	 */
	 
	 /* Description:
			------------
			This postestimation command implements the REBUS-PLS approach to assess
			the presence of unobserved heterogeneity.
	 */
	
	if ("`e(formative)'" != "") {
		display as error "the REBUS-PLS approach can't be used with formative blocks"
		error 198
	}
	local props = e(properties)
	local struct "structural"
	local isstruct : list struct in props
	if (!`isstruct') {
		display as error "the fitted plssem model includes only the measurement part"
		error 198
	}
	local rawsum "rawsum"
	local israwsum : list rawsum in props
	local cc "relative"
	local isrelative : list cc in props
	if (`isrelative') {
		local convcrit "relative"
	}
	else {
		local convcrit "square"
	}
	local initmet "indsum"
	local isindsum : list initmet in props
	if (`isindsum') {
		local init "indsum"
	}
	else {
		local init "eigen"
	}
	local indscale "scaled"
	local isscaled : list indscale in props
	if (`isscaled') {
		local scale ""
	}
	else {
		local scale "noscale"
	}
	local boot "bootstrap"
	local isboot : list boot in props
	if (`isboot') {
		display as error "the global model used the 'boot()' option, which slows down excessively the REBUS-PLS calculations"
		display as error "try refitting the global model without bootstrap"
		exit
	}
	if ("`test'" != "") {
		if ("`seed'" != "") {
			if (`seed' < 0 | `seed' >  2^31-1) {
				display as error "'seed()' option requires a value between 0 and 2^31-1"
				exit
			}
		}
	}
	if ("`name'" != "") & (`: word count `name'' != 1) {
		display as error "the 'name' option must include a single word"
		exit
	}
	if ("`name'" == "") {
		local name "rebus_class"
	}
	if ("`e(binarylvs)'" != "") {
		display as error "currently the REBUS-PLS implementation doesn't allow for binary latent variables"
		exit
	}
	local knnimp "knn"
	local isknnimp : list knnimp in props
	local meanimp "mean"
	local ismeanimp : list meanimp in props
	if (`isknnimp') {
		local missing "knn"
	}
	else if (`ismeanimp') {
		local missing "mean"
	}
	else {
		local missing ""
	}
	
	/* Global model */
	tempname globalmodel gof_global
	_estimates hold `globalmodel', copy
	preserve
	capture quietly predict, residuals
	if (_rc == 908) {
		display as error "matsize too small"
    display as error "    The command has attempted to create a matrix "  _continue
		display as error "with more than 400 rows or columns."
    display as error "    To run the command increase matsize by using the " _continue
    display as result "set matsize" _continue
		display as error " command; see help " _continue
		display as smcl "{help matsize}" _continue
		display as error "."
		restore
		exit
	}
	local todrop "`: colnames r(meas_res)' `: colnames r(struct_res)'"
	
	/* Collecting results to display later */
	matrix `gof_global' = e(assessment)
	matrix `gof_global' = `gof_global'[1, 3]
	local gof_loc = strofreal(`gof_global'[1, 1], "%9.4f")
	
	tempvar __touse__
	quietly generate byte `__touse__' = e(sample)
	quietly count if `__touse__'
	local N = r(N)
	
	if ("`missing'" != "") {
		/* Save original data set */
		local allindicators = e(mvs)
		tempname original_data
		mata: `original_data' = st_data(., "`: list uniq allindicators'", "`__touse__'")
		
		/* Recovery of missing values */
		mata: st_store(., tokens("`: list uniq allindicators'"), "`__touse__'", ///
			st_matrix("e(imputed_data)"))
	}
	
	/* Clustering residuals using Ward hierarchical linkage */
	tempvar rebus_clus
	capture cluster drop `rebus_clus'
	cluster wardslinkage `todrop' if `__touse__', name(`rebus_clus') ///
		measure(Euclidean)
	if ("`dendrogram'" != "") {
		tempname dendplot
		capture {
			quietly cluster dendrogram `rebus_clus' if `__touse__', ///
				xlabel(, angle(90) labsize(*.5)) ytitle("Euclidean distance") ///
				subtitle("(based on residuals from the global model)") ///
				title("Dendrogram for REBUS-PLS analysis") scheme(sj) ///
				name(`dendplot', replace)
		}
		if (_rc != 0) {
			display as error "dendrogram is not available because of too many leaves"
		}
	}
	if ("`numclass'" == "") {
		if (`N' < 5*`maxclass') {
			local maxclass = max(2, floor(min(`maxclass', `N'/5)))
		}
		tempname rebus_stop i w
		/*
		quietly cluster stop `rebus_clus', rule(duda) groups(2/`maxclass') ///
			matrix(`rebus_stop')
		*/
		quietly cluster stop `rebus_clus', rule(calinski) groups(2/`maxclass') ///
			matrix(`rebus_stop')
		mata: `i' = 0; `w' = 0
		mata: maxindex(st_matrix("`rebus_stop'")[., 2], 1, `i', `w')
		mata: st_local("numclass", strofreal(`i' + 1))
	}
	if (`N' < 5*`numclass') {
		display as text "warning: the number of classes chosen seems to be too large"
		display as text "calculations may abort; " _continue
		display as text "in this case, consider reducing the number of classes"
	}
	tempvar rebus_class
	capture drop `rebus_class'
	cluster generate `rebus_class' = groups(`numclass'), name(`rebus_nm')
	quietly drop `todrop'
	
	/* Parse global model e(cmdline) */
	local cmdline = e(cmdline)
	local trash
	gettoken cmdline options : cmdline, parse(",")
	gettoken trash options : options, parse(",")

	local mm_start = strpos(`"`cmdline'"', "(")
	if (strpos(`"`cmdline'"', "if")) {
		local mm_end = strpos(`"`cmdline'"', "if") - `mm_start' - 1
	}
	else {
		local mm_end = strlen(`"`cmdline'"')
	}
	local mm = substr(`"`cmdline'"', `mm_start', `mm_end')
	local mm : list clean mm

	tokenize `"`options'"', parse(")")
	local tok_i = 1
	while (substr(`"``tok_i''"', 1, 3) != "str") & (`"``tok_i''"' != "") {
		local ++tok_i
	}
	local sm_full = `"``tok_i''"' + ")"
	tokenize `"``tok_i''"', parse("(")
	local sm = `"`3'"'
	local sm : list clean sm
	local options : list options - sm_full
	local options : list clean options

	// the following code doesn't work
	/*
	tokenize `"`options'"', parse(")")
	local tok_i = 1
	while (substr(`"``tok_i''"', 1, 3) != "dig") & (`"``tok_i''"' != "") {
		local ++tok_i
	}
	local options_digits = `"``tok_i''"' + ")"
	local options : list options - options_digits
	local options : list clean options
	*/
	if (`israwsum') {
		local options "tol(`e(tolerance)') wscheme("
	}
	else {
		local options "tol(`e(tolerance)') maxiter(`e(maxiter)') wscheme("
	}
	local ws_centroid "centroid"
	local ws_factor "factor"
	local ws_path "path"
	if (`: list ws_centroid in props') {
		local scheme "centroid"
	}
	else if (`: list ws_factor in props') {
		local scheme "factor"
	}
	else if (`: list ws_path in props') {
		local scheme "path"
	}
	local options "`options'`scheme')"
	if ("`e(binarylvs)'" != "") {
		local options "`options' binary(`e(binarylvs)')"
	}
	
	/* Set temporary variables */
	local allindicators = e(mvs)
	local alllatents = e(lvs)
	local allreflective = e(reflective)
	local num_ind : word count `allindicators'
	local num_lv : word count `alllatents'
	tempname cm ow path ind indstd y_local loads r2 block x_hat out_res endo ///
		y_hat inn_res
	tempvar __touseloc__
	quietly generate byte `__touseloc__' = .
	mata: `endo' = colsum(st_matrix("e(adj_struct)"))
	mata: `endo' = (`endo' :> 0)
	foreach var in `allindicators' {
		local allstdindicators "`allstdindicators' std`var'"
	}
	local allstdindicators : list clean allstdindicators
	
	/* Run the REBUS-PLS algorithm */
	tempname res_rebus
	capture noisily {
		mata: `res_rebus' = ///
			plssem_rebus( ///
				st_data(., "`allindicators'"), ///			 note: `__touse__' not used here
				st_matrix("e(adj_meas)"), ///
				st_matrix("e(adj_struct)"), ///
				"`allindicators'", ///
				"`allstdindicators'", ///
				"`alllatents'", ///
				"`e(binarylvs)'", ///
				st_numscalar("e(tolerance)"), ///
				st_numscalar("e(maxiter)"), ///
				"`__touse__'", ///
				"`scheme'", ///
				"`convcrit'", ///
				"`init'", ///
				"`scale'", ///
				strtoreal("`isstruct'"), ///
				strtoreal("`israwsum'"), ///
				0, ///
				"`rebus_class'", ///
				strtoreal("`numclass'"), ///
				strtoreal("`maxiter'"), ///
				strtoreal("`stop'"), ///
				1)
		
		mata: st_local("iter", strofreal(`res_rebus'.niter))
		mata: st_local("rN0", strofreal(`res_rebus'.rN0))
		mata: st_local("rN_lte_5", strofreal(`res_rebus'.rN_lte_5))
		mata: st_store(., "`rebus_class'", `res_rebus'.rebus_class)
	}
	if (_rc != 0) {
		if (mod(`iter', 5) == 0) {
			display as error " aborting"
		}
		else {
			display as error "aborting"
		}
		if (real("`rN0'")) {
			display as error "at least one class is empty"
		}
		else if (real("`rN_lte_5'")) {
			display as error "too few observations (5 or less) assigned to a single class"
		}
		else if (_rc == 409) {
			display as error "at least one indicator has zero variance in one iteration"
		}
		else {
			display as error "something went wrong in the REBUS-PLS calculations"
		}
		display as error "try reducing the number of classes " _continue
		display as error "or relaxing the stopping criteria"
		restore
		_estimates unhold `globalmodel'
		exit
	}
	/*
	else {
		if (mod(`iter', 5) == 0) {
			display as text " done!"
		}
		else {
			display as text "done!"
		}
	}
	*/
	/* End of REBUS-PLS algorithm */
	
	/* Checking that the established classes have enough observations */
	forvalues k = 1/`numclass' {
		quietly count if (`rebus_class' == `k' & `__touse__')
		if (r(N) < 10) {
			display as error "less than 10 observations in class " + `k'
			display as error "at least 10 complete observations required " _continue
			display as error "in each class to proceed with the calculations"
			display as error "try reducing the number of classes with " _continue
			display as error "the 'numclass()' option"
			error 2001
		}
	}
	
	/* Once stability is attained, final local models are estimated */
	local allendogenous "`: colnames e(struct_b)'"
	tempname lat class indstd_st lat_st out_res_st inn_res_st class_st gqi
	mata: `indstd_st' = J(0, strtoreal("`num_ind'"), .)
	mata: `lat_st' = J(0, sum(`endo'), .)
	mata: `out_res_st' = J(0, strtoreal("`num_ind'"), .)
	mata: `inn_res_st' = J(0, sum(`endo'), .)
	mata: `class_st' = J(0, 1, .)
	local donotcleanup "nocleanup"
	forvalues k = 1/`numclass' {
		tempname localmodel_`k'
		quietly plssem `mm' if (`rebus_class' == `k' & `__touse__'), ///
			structural(`sm') `options' `donotcleanup'
		_estimates hold `localmodel_`k'', copy
		quietly drop `__touseloc__'
		quietly generate byte `__touseloc__' = e(sample)
		
		mata: st_view(`lat' = ., ., "`allendogenous'", "`__touseloc__'")
		mata: st_view(`ind' = ., ., "`allindicators'", "`__touseloc__'")
		mata: `indstd' = scale(`ind')
		
		mata: `ow' = st_matrix("e(outerweights)")
		mata: `path' = st_matrix("e(pathcoef)")
		mata: `loads' = st_matrix("e(loadings)")
		mata: `r2' = st_matrix("e(rsquared)")
		mata: `r2' = `r2'[., selectindex(`r2' :!= .)]
		
		mata: `y_local' = `indstd' * `ow'
		mata: `out_res' = J(rows(`ind'), 0, .)
		forvalues j = 1/`num_lv' {
			mata: `block' = selectindex(`loads'[., `j'] :!= .)
			mata: `x_hat' = `y_local'[., `j'] * `loads'[`block', `j']'
			mata: `out_res' = (`out_res', (`indstd'[., `block'] - `x_hat'))
		}
		mata: `y_hat' = `y_local' * `path'[., selectindex(`endo')]
		mata: `inn_res' = `y_local'[., selectindex(`endo')] - `y_hat'
		mata: `class' = J(rows(`ind'), 1, `k')

		mata: `indstd_st' = (`indstd_st' \ `indstd')
		mata: `lat_st' = (`lat_st' \ `lat')
		mata: `out_res_st' = (`out_res_st' \ `out_res')
		mata: `inn_res_st' = (`inn_res_st' \ `inn_res')
		mata: `class_st' = (`class_st' \ `class')
	}
	mata: `gqi' = ///
		plssem_rebus_gqi(`indstd_st', `lat_st', `out_res_st', `inn_res_st', `class_st')
	tempname gqi_final
	mata: st_numscalar("`gqi_final'", `gqi')
	/* End final REBUS-PLS calculations */

	/* Permutation test */
	if ("`test'" != "") {
		tempname res_ptest
		if ("`reps'" == "") {
			local reps = 100
		}
	
	capture noisily {
		mata: `res_ptest' = ///
			plssem_rebus_ptest( ///
				st_data(., "`allindicators'"), ///					 note: `touse' not used here
				st_matrix("e(adj_meas)"), ///
				st_matrix("e(adj_struct)"), ///
				"`allindicators'", ///
				"`allstdindicators'", ///
				"`alllatents'", ///
				"`e(binarylvs)'", ///
				st_numscalar("e(tolerance)"), ///
				st_numscalar("e(maxiter)"), ///
				"`__touse__'", ///
				"`scheme'", ///
				"`convcrit'", ///
				"`init'", ///
				"`scale'", ///
				strtoreal("`isstruct'"), ///
				strtoreal("`israwsum'"), ///
				0, ///
				"`rebus_class'", ///
				strtoreal("`numclass'"), ///
				strtoreal("`reps'"), ///
				strtoreal("`seed'"), ///
				1)
		}
		if (_rc != 0) {
			if (_rc == 409) {
				display as error "at least one indicator has zero variance " _continue
				display as error "in one iteration of the permutation test"
			}
			else {
				display as error "something went wrong in the REBUS-PLS permutation test"
			}
			restore
			_estimates unhold `globalmodel'
			exit
		}
		
		local oldN = _N
		if (`oldN' < `reps') {
			quietly set obs `reps'
		}
		tempname pvalue_sc permdist
		tempvar gqi_dist
		quietly generate `gqi_dist' = .
		mata: st_store(range(1, `reps', 1), "`gqi_dist'", `res_ptest')
		quietly count if (`gqi_dist' > `gqi_final') & !missing(`gqi_dist')
		local pvalue = strofreal(r(N)/`reps', "%9.4f")
		scalar `pvalue_sc' = r(N)/`reps'
		local gqi_final_loc = strofreal(`gqi_final', "%9.4f")
		if ("`plot'" != "") {
			quietly twoway ///
				histogram `gqi_dist', fraction || ///
				scatteri 0 `gqi_final_loc' (12) "GQI", msymbol(D) msize(medlarge) ///
				mfcolor(gs12) mlcolor(black) mlabsize(medsmall) mlabgap(*2) || ///
				scatteri 0 `gof_loc' (12) "GoF", msymbol(O) msize(medlarge) ///
				mfcolor(gs12) mlcolor(black) mlabsize(medsmall) mlabgap(*2) || , ///
				xtitle("Statistic value") ytitle("") ///
				title("Empirical distribution of the Group Quality Index (GQI)") ///
				subtitle("(based on `reps' replications)") legend(off) ///
				note("GoF: `gof_loc'" "GQI: `gqi_final_loc' - p-value = `pvalue'") ///
				scheme(sj) name(`permdist', replace)
		}
		if (`oldN' < `reps') {
			local firsttodelete = `oldN' + 1
			quietly drop in `firsttodelete'/l
		}
	}
	/* End permutation test */

	/* Display results */
	local skip1 = 1
	local skip3 = 3

	tempname tmp strb nonmiss alllen numeff
	matrix `tmp' = e(struct_b)
	matrix `strb' = vec(`tmp')
	local nrows = rowsof(`tmp')
	local ncols = colsof(`tmp')
	local totrows = `nrows'*`ncols'
	matrix `nonmiss' = J(`totrows', 1, 0)
	local strb_rn : rowfullnames `strb'
	forvalues i = 1/`totrows' {
		if (!missing(`strb'[`i', 1])) {
			matrix `nonmiss'[`i', 1] = 1
			local nm_tmp `: word `i' of `strb_rn''
			local tok_i = 1
			tokenize `"`nm_tmp'"', parse(":")
			while ("``tok_i''" != "") {
				if (`tok_i' == 1) {
					local nm_Y "``tok_i''"
				}
				else if (`tok_i' == 3) {
					local nm_X "``tok_i''"
				}
				local ++tok_i
			}
			local strbok_rn "`strbok_rn' `nm_X':`nm_Y'"
		}
	}	
	mata: st_numscalar("`numeff'", colsum(st_matrix("`nonmiss'")))
	local neff = `numeff'
	matrix `alllen' = J(`neff', 1, .)
	local resnm : subinstr local strbok_rn ":" "->", all
	forvalues j = 1/`neff' {
		local nm : word `j' of `resnm'
		local nm_len : strlen local nm
		matrix `alllen'[`j', 1] = `nm_len' + 2
	}

	tempname results_n results_p results_l tmp_mat
	local nind : word count `allindicators'
	forvalues k = 1/`numclass' {
		local grp_cn `grp_cn' "Class_`k'"
	}
	matrix `results_n' = J(3, 1 + `numclass', .)
	matrix rownames `results_n' = "Observations" "Percentage" "GoF"
	matrix colnames `results_n' = "Global" `grp_cn'
	matrix `results_p' = J(`neff', 1 + `numclass', .)
	matrix rownames `results_p' = `resnm'
	matrix colnames `results_p' = "Global" `grp_cn'
	matrix `results_l' = J(`nind', 1 + `numclass', .)
	matrix rownames `results_l' = `allindicators'
	matrix colnames `results_l' = "Global" `grp_cn'

	_estimates unhold `globalmodel'
	_estimates hold `globalmodel', copy  // needed to keep it in memory
	matrix `tmp' = e(struct_b)
	matrix `strb' = vec(`tmp')

	tempname global_N
	scalar `global_N' = e(N)
	matrix `results_n'[1, 1] = e(N)
	matrix `results_n'[2, 1] = 100
	matrix `results_n'[3, 1] = `gof_global'[1, 1]

	local b_i = 1
	forvalues i = 1/`totrows' {
		if (`nonmiss'[`i', 1]) {
			matrix `results_p'[`b_i', 1] = `strb'[`i', 1]
			local ++b_i
		}
	}
	mata: st_matrix("`tmp_mat'", rowsum(st_matrix("e(loadings)")))
	matrix `results_l'[1, 1] = `tmp_mat'
	
	tempname gof_tmp
	forvalues k = 1/`numclass' {
		_estimates unhold `localmodel_`k''
		matrix `tmp' = e(struct_b)
		matrix `strb' = vec(`tmp')

		matrix `gof_tmp' = e(assessment)
		matrix `gof_tmp' = `gof_tmp'[1, 3]
		matrix `results_n'[1, 1 + `k'] = e(N)
		matrix `results_n'[2, 1 + `k'] = e(N)/`global_N'*100
		matrix `results_n'[3, 1 + `k'] = `gof_tmp'[1, 1]
		
		local b_i = 1
		forvalues i = 1/`totrows' {
			if (`nonmiss'[`i', 1]) {
				matrix `results_p'[`b_i', 1 + `k'] = `strb'[`i', 1]
				local ++b_i
			}
		}

		mata: st_matrix("`tmp_mat'", rowsum(st_matrix("e(loadings)")))
		matrix `results_l'[1, 1 + `k'] = `tmp_mat'
	}
	
	local gqi_final_loc = `gqi_final'
	mkheader, digits(5) rebus_it(`iter') rebus_gqi(`gqi_final_loc')

	tempname maxlen
	mata: st_numscalar("`maxlen'", max(st_matrix("`alllen'")))
  local firstcollbl ""	
	local `maxlen' = max(max(strlen("`firstcollbl'"), `maxlen') + 2, 14)
	local colw = max(9, `digits' + 5)

	local title "REBUS-PLS classes"
	local firstcollbl ""
	mktable, matrix(`results_n') digits(`digits') firstcolname(`firstcollbl') ///
		title(`title') firstcolwidth(``maxlen'') colwidth(`colw') hlines(3) ///
		novlines total rebus

	local title "Loadings"
	local firstcollbl ""
	mktable, matrix(`results_l') digits(`digits') firstcolname(`firstcollbl') ///
		title(`title') firstcolwidth(``maxlen'') colwidth(`colw') hlines(`nind') ///
		novlines total

	local title "Path coefficients"
	local firstcollbl ""
	mktable, matrix(`results_p') digits(`digits') firstcolname(`firstcollbl') ///
		title(`title') firstcolwidth(``maxlen'') colwidth(`colw') hlines(`neff') ///
		novlines total

	if ("`test'" != "") {
		tempname results_t
		matrix `results_t' = J(2, 1, .)
		matrix rownames `results_t' = "Replications" "P-value"
		matrix colnames `results_t' = "Value"
		matrix `results_t'[1, 1] = `reps'
		matrix `results_t'[2, 1] = `pvalue_sc'
		
		local title "Permutation test"
		local firstcollbl ""
		mktable, matrix(`results_t') digits(`digits') firstcolname(`firstcollbl') ///
			title(`title') firstcolwidth(``maxlen'') colwidth(`colw') hlines(2) ///
			novlines total rebus
	}
	/* End of display */

	/* Restore global model results */
	tempname rebus_c
	mata: `rebus_c' = st_data(., "`rebus_class'", "`__touse__'")
	restore
	_estimates unhold `globalmodel'
	quietly generate byte `__touse__' = e(sample)

	/* Save final classification */
	capture drop `name'
	quietly generate int `name' = .
	mata: st_store(., "`name'", "`__touse__'", `rebus_c')
	local now "`c(current_date)', `c(current_time)'"
	local now : list clean now
	label variable `name' "REBUS-PLS classification [`now']"
	
	/* Clean up */
	if ("`missing'" != "") {
		mata: st_store(., tokens("`: list uniq allindicators'"), "`__touse__'", ///
			`original_data')
	}
	
	/* Return values */
	return scalar nclasses = `numclass'
	return scalar GQI = `gqi_final'
	return scalar GoF = `gof_global'[1, 1]
	
	/* Maximum number of iterations reached */
	if (`iter' > `maxiter') {
		display as error "warning: REBUS-PLS algorithm did not converge"
		display as error "the solution provided may not be acceptable; " _continue
		display as error "try relaxing the stopping criteria"
	}
end

program FIMIX, rclass
	version 15.1
	syntax [ , Method(string) Numclass(numlist integer >=1 max=1) ///
		MAXITer(integer 30000) Stop(real 1.e-5) ///
		REStart(numlist integer >=1 max=1) SEed(numlist max=1) ///
		name(string) DIGits(integer 3) GRoups(numlist integer >=1 <=20 min=1 sort) ]

	/* Options:
	   --------
		 method(string)									--> method to use for assessing unobserved
																				heterogeneity; available methods are
																				'rebus', 'fimix' and 'gas'. default is
																				'rebus'.
		 numclass(integer integer >=1
							max=1))								--> number of classes to use; if empty, it
																				is chosen automatically using a Ward
																				hierarchical algorithm
		 maxiter(integer 30000)					--> maximum number of iterations (default
																				30000)
		 stop(real 1e-5)								--> stopping criterion (default 1e-5)
		 restart(numlist integer >=1
						 max=1)									--> number of repetitions of the EM
																				algorithm in FIMIX-PLS (default is 10)
		 seed(numlist max=1)						--> seed number
		 name(string)										--> variable name where to store the final
																				FIMIX-PLS classification
		 digits(integer 3)							--> number of digits to display (default 3)
		 groups(numlist integer >=1
						<=20 min=1 sort)				--> list of group numbers for which to
																				compare the corresponding fit indices
																				using FIMIX-PLS
	 */
	 
	 /* Description:
			------------
			This postestimation command implements the FIMIX-PLS approach to assess
			the presence of unobserved heterogeneity.
	 */
	
	local props = e(properties)
	local struct "structural"
	local isstruct : list struct in props
	if (!`isstruct') {
		display as error "the fitted plssem model includes only the measurement part"
		error 198
	}
	local rawsum "rawsum"
	local israwsum : list rawsum in props
	local cc "relative"
	local isrelative : list cc in props
	if (`isrelative') {
		local convcrit "relative"
	}
	else {
		local convcrit "square"
	}
	local initmet "indsum"
	local isindsum : list initmet in props
	if (`isindsum') {
		local init "indsum"
	}
	else {
		local init "eigen"
	}
	local indscale "scaled"
	local isscaled : list indscale in props
	if (`isscaled') {
		local scale ""
	}
	else {
		local scale "noscale"
	}
	local boot "bootstrap"
	local isboot : list boot in props
	if (`isboot') {
		display as error "the global model used the 'boot()' option, which slows down excessively the FIMIX-PLS calculations"
		display as error "try refitting the global model without bootstrap"
		exit
	}
	if ("`seed'" != "") {
		if (`seed' < 0 | `seed' >  2^31-1) {
			display as error "'seed()' option requires a value between 0 and 2^31-1"
			exit
		}
	}
	if ("`name'" != "") & (`: word count `name'' != 1) {
		display as error "the 'name' option must include a single word"
		exit
	}
	if ("`name'" == "") {
		local name "fimix_class"
	}
	if ("`e(binarylvs)'" != "") {
		display as error "currently the FIMIX-PLS implementation doesn't allow for binary latent variables"
		exit
	}
	local knnimp "knn"
	local isknnimp : list knnimp in props
	local meanimp "mean"
	local ismeanimp : list meanimp in props
	if (`isknnimp') {
		local missing "knn"
	}
	else if (`ismeanimp') {
		local missing "mean"
	}
	else {
		local missing ""
	}
	if ("`restart'" == "") {
		local restart = 10
	}
	if (("`numclass'" == "") & ("`groups'" == "")) {
		display as error "the 'numclass' option must include one value for FIMIX-PLS"
		exit
	}
	
	tempvar __touse__
	quietly generate byte `__touse__' = e(sample)
	quietly count if `__touse__'
	local N = r(N)

	if ("`missing'" != "") {
		/* Save original data set */
		local allindicators = e(mvs)
		tempname original_data
		mata: `original_data' = st_data(., "`: list uniq allindicators'", "`__touse__'")
		
		/* Recovery of missing values */
		mata: st_store(., tokens("`: list uniq allindicators'"), "`__touse__'", ///
			st_matrix("e(imputed_data)"))
	}

  /* Parse global model e(cmdline) */
  local cmdline = e(cmdline)
  local trash
  gettoken cmdline options : cmdline, parse(",")
  gettoken trash options : options, parse(",")

  local mm_start = strpos(`"`cmdline'"', "(")
  if (strpos(`"`cmdline'"', "if")) {
    local mm_end = strpos(`"`cmdline'"', "if") - `mm_start' - 1
  }
  else {
    local mm_end = strlen(`"`cmdline'"')
  }
  local mm = substr(`"`cmdline'"', `mm_start', `mm_end')
  local mm : list clean mm

  tokenize `"`options'"', parse(")")
  local tok_i = 1
  while (substr(`"``tok_i''"', 1, 3) != "str") & (`"``tok_i''"' != "") {
    local ++tok_i
  }
  local sm_full = `"``tok_i''"' + ")"
  tokenize `"``tok_i''"', parse("(")
  local sm = `"`3'"'
  local sm : list clean sm
  local options : list options - sm_full
  local options : list clean options

  // the following code doesn't work
  /*
  tokenize `"`options'"', parse(")")
  local tok_i = 1
  while (substr(`"``tok_i''"', 1, 3) != "dig") & (`"``tok_i''"' != "") {
    local ++tok_i
  }
  local options_digits = `"``tok_i''"' + ")"
  local options : list options - options_digits
  local options : list clean options
  */
	if (`israwsum') {
		local options "tol(`e(tolerance)') wscheme("
	}
	else {
		local options "tol(`e(tolerance)') maxiter(`e(maxiter)') wscheme("
	}
  local ws_centroid "centroid"
  local ws_factor "factor"
  local ws_path "path"
  if (`: list ws_centroid in props') {
    local scheme "centroid"
  }
  else if (`: list ws_factor in props') {
    local scheme "factor"
  }
  else if (`: list ws_path in props') {
    local scheme "path"
  }
  local options "`options'`scheme')"
  if ("`e(binarylvs)'" != "") {
    local options "`options' binary(`e(binarylvs)')"
  }
  
	/* Set temporary variables */
	local allindicators = e(mvs)
	local alllatents = e(lvs)
	local allreflective = e(reflective)
	local num_ind : word count `allindicators'
	local num_lv : word count `alllatents'
	tempname endo
	mata: `endo' = colsum(st_matrix("e(adj_struct)"))
	mata: `endo' = (`endo' :> 0)
	foreach var in `allindicators' {
		local allstdindicators "`allstdindicators' std`var'"
	}
	local allstdindicators : list clean allstdindicators

	/* Run the FIMIX-PLS algorithm */
	tempname res_fimix
	if ("`groups'" == "") {
		capture noisily {
			mata: `res_fimix' = ///
				plssem_fimix( ///
					st_data(., "`alllatents'"), ///			 	 note: `__touse__' not used here
					st_matrix("e(adj_struct)"), ///
					"`__touse__'", ///
					strtoreal("`numclass'"), ///
					strtoreal("`maxiter'"), ///
					strtoreal("`stop'"), ///
					strtoreal("`restart'"), ///
					strtoreal("`seed'"), ///
					1, 0)
			
			mata: st_local("iter", strofreal(`res_fimix'.niter))
			mata: st_local("rN0", strofreal(`res_fimix'.rN0))
			mata: st_local("rN_lte_5", strofreal(`res_fimix'.rN_lte_5))
			mata: st_local("fimix_ll", strofreal(max(`res_fimix'.ll)))
		}
		if (_rc != 0) {
			if (mod(`iter', 5) == 0) {
				display as error " aborting"
			}
			else {
				display as error "aborting"
			}
			if (_rc == 409) {
				display as error "at least one indicator has zero variance in one iteration"
			}
			else {
				display as error "something went wrong in the FIMIX-PLS calculations"
			}
			display as error "try reducing the number of classes " _continue
			display as error "or relaxing the stopping criteria"
			exit
		}
	}
	else {
		display
		display as text "Computing FIMIX-PLS solution with..."
		
		tempname ic_groups ic_tmp
		foreach g of numlist `groups' {
			if (`g' == 1) {
				display as text "  --> " as result `g' as text " group  " _continue
			}
			else {
				display as text "  --> " as result `g' as text " groups " _continue
			}
			capture noisily {
				mata: `res_fimix' = ///
					plssem_fimix( ///
						st_data(., "`alllatents'"), ///			 note: `__touse__' not used here
						st_matrix("e(adj_struct)"), ///
						"`__touse__'", ///
						`g', ///
						strtoreal("`maxiter'"), ///
						strtoreal("`stop'"), ///
						strtoreal("`restart'"), ///
						strtoreal("`seed'"), ///
						0, 1)
				
				mata: st_local("iter", strofreal(`res_fimix'.niter))
				mata: st_local("rN0", strofreal(`res_fimix'.rN0))
				mata: st_local("rN_lte_5", strofreal(`res_fimix'.rN_lte_5))
				mata: st_local("fimix_ll", strofreal(max(`res_fimix'.ll)))
			}
			if (_rc != 0) {
				if (mod(`iter', 5) == 0) {
					display as error " aborting"
				}
				else {
					display as error "aborting"
				}
				if (_rc == 409) {
					display as error "at least one indicator has zero variance in one iteration"
				}
				else {
					display as error "something went wrong in the FIMIX-PLS calculations"
				}
				display as error "try relaxing the stopping criteria"
				exit
			}
			
			mata: st_matrix("`ic_tmp'", `res_fimix'.ic)
			matrix `ic_groups' = (nullmat(`ic_groups'), `ic_tmp')
			local ic_cn "`ic_cn' `g'"
		}
	}
	/* End of FIMIX-PLS algorithm */
	
  /* Display results */
  local skip1 = 1
  local skip3 = 3
	
	if ("`groups'" == "") {
		tempname tmp strb nonmiss alllen numeff
		matrix `tmp' = e(struct_b)
		matrix `strb' = vec(`tmp')
		local nrows = rowsof(`tmp')
		local ncols = colsof(`tmp')
		local totrows = `nrows'*`ncols'
		matrix `nonmiss' = J(`totrows', 1, 0)
		local strb_rn : rowfullnames `strb'
		forvalues i = 1/`totrows' {
			if (!missing(`strb'[`i', 1])) {
				matrix `nonmiss'[`i', 1] = 1
				local nm_tmp `: word `i' of `strb_rn''
				local tok_i = 1
				tokenize `"`nm_tmp'"', parse(":")
				while ("``tok_i''" != "") {
					if (`tok_i' == 1) {
						local nm_Y "``tok_i''"
					}
					else if (`tok_i' == 3) {
						local nm_X "``tok_i''"
					}
					local ++tok_i
				}
				local strbok_rn "`strbok_rn' `nm_X':`nm_Y'"
			}
		}
		mata: st_numscalar("`numeff'", colsum(st_matrix("`nonmiss'")))
		local neff = `numeff'
		matrix `alllen' = J(`neff', 1, .)
		local resnm : subinstr local strbok_rn ":" "->", all
		forvalues j = 1/`neff' {
			local nm : word `j' of `resnm'
			local nm_len : strlen local nm
			matrix `alllen'[`j', 1] = `nm_len' + 2
		}

		tempname results_n results_p results_ic tmp_mat
		local nind : word count `allindicators'
		forvalues k = 1/`numclass' {
			local grp_cn `grp_cn' "Class_`k'"
		}
		matrix `results_n' = J(3, 1 + `numclass', .)
		matrix rownames `results_n' = "Observations" "Percentage" "Mixture%%weight"
		matrix colnames `results_n' = "Global" `grp_cn'
		matrix `results_p' = J(`neff', 1 + `numclass', .)
		matrix rownames `results_p' = `resnm'
		matrix colnames `results_p' = "Global" `grp_cn'

		tempname global_N
		scalar `global_N' = e(N)
		matrix `results_n'[1, 1] = e(N)
		matrix `results_n'[2, 1] = 100

		local b_i = 1
		forvalues i = 1/`totrows' {
			if (`nonmiss'[`i', 1]) {
				matrix `results_p'[`b_i', 1] = `strb'[`i', 1]
				local ++b_i
			}
		}
		
		tempname fimixclassfreq fimixrho
		mata: st_matrix("`fimixclassfreq'", `res_fimix'.classfreq)
		mata: st_matrix("`fimixrho'", `res_fimix'.rho)
		forvalues k = 1/`numclass' {
			mata: st_matrix("`tmp'", `res_fimix'.path[`k', 1].mat)
			matrix `strb' = vec(`tmp')

			matrix `results_n'[1, 1 + `k'] = `fimixclassfreq'[`k', 1]
			matrix `results_n'[2, 1 + `k'] = `fimixclassfreq'[`k', 1]/`global_N'*100
			matrix `results_n'[3, 1 + `k'] = `fimixrho'[1, `k']
			
			local b_i = 1
			forvalues i = 1/`totrows' {
				if (`nonmiss'[`i', 1]) {
					matrix `results_p'[`b_i', 1 + `k'] = `strb'[`i', 1]
					local ++b_i
				}
			}
		}
		
		mata: st_matrix("`results_ic'", `res_fimix'.ic)
		matrix rownames `results_ic' = "AIC" "AIC3" "AIC4" "BIC" "CAIC" "HQ" ///
																	 "MDL5" "LnL" "EN" "NFI" "NEC"
		matrix colnames `results_ic' = "Value"
		
		mkheader, digits(5) fimix_it(`iter') fimix_ll(`fimix_ll')
		
		tempname maxlen
		mata: st_numscalar("`maxlen'", max(st_matrix("`alllen'")))
		local firstcollbl ""
		local `maxlen' = max(max(strlen("`firstcollbl'"), `maxlen') + 2, 16)
		local colw = max(10, `digits' + 5)

		local title "FIMIX-PLS classes"
		mktable, matrix(`results_n') digits(`digits') firstcolname(`firstcollbl') ///
			title(`title') firstcolwidth(``maxlen'') colwidth(`colw') hlines(3) ///
			novlines total rebus fimix

		local title "Path coefficients"
		local firstcollbl ""
		mktable, matrix(`results_p') digits(`digits') firstcolname(`firstcollbl') ///
			title(`title') firstcolwidth(``maxlen'') colwidth(`colw') hlines(`neff') ///
			novlines total
		
		local title "Fit indices"
		local firstcollbl ""
		mktable, matrix(`results_ic') digits(`digits') firstcolname(`firstcollbl') ///
			title(`title') firstcolwidth(``maxlen'') colwidth(`colw') hlines(11) ///
			novlines total
	}
	else {
		matrix rownames `ic_groups' = "AIC" "AIC3" "AIC4" "BIC" "CAIC" "HQ" ///
																	"MDL5" "LnL" "EN" "NFI" "NEC"
		matrix colnames `ic_groups' = `ic_cn'

		local colw = max(11, `digits' + 5)
		local title "FIMIX-PLS comparison"
		local firstcollbl "Groups"
		mktable, matrix(`ic_groups') digits(`digits') firstcolname(`firstcollbl') ///
			title(`title') firstcolwidth(8) colwidth(`colw') hlines(11) ///
			novlines total
	}
  /* End of display */

	/* Save final classification */
	if ("`groups'" == "") {
		capture drop `name'
		quietly generate int `name' = .
		mata: st_store(., "`name'", "`__touse__'", `res_fimix'.fimix_class)
		local now "`c(current_date)', `c(current_time)'"
		local now : list clean now
		label variable `name' "FIMIX-PLS classification [`now']"
	}
	
	/* Clean up */
	if ("`missing'" != "") {
		mata: st_store(., tokens("`: list uniq allindicators'"), "`__touse__'", ///
			`original_data')
	}
	
	/* Return values */
	if ("`groups'" == "") {
		return scalar nclasses = `numclass'
		tempname ll ll_c ll_restart
		mata: st_matrix("`ll'", `res_fimix'.ll)
		mata: st_matrix("`ll_c'", `res_fimix'.ll_c)
		mata: st_matrix("`ll_restart'", `res_fimix'.ll_restart)
		return matrix loglik_restart `ll_restart'
		return matrix loglik_c `ll_c'
		return matrix loglik `ll'
		return matrix fitindices `results_ic'
	}
	else {
		return matrix fitindices `ic_groups'
	}
	
	/* Maximum number of iterations reached */
	if ("`groups'" == "") {
		if (real("`rN0'")) {
			display as error "warning: at least one class is empty"
		}
		if (real("`rN_lte_5'")) {
			display as error "warning: too few observations (5 or less) assigned to a single class"
		}
	}
	if (`iter' > `maxiter') {
		display as error "warning: FIMIX-PLS algorithm did not converge"
		display as error "the solution provided may not be acceptable; " _continue
		display as error "try relaxing the stopping criteria"
	}
end

program GAS, rclass
	version 15.1
	syntax [ , Method(string) Numclass(numlist integer >=1 max=1) ///
		POPsize(numlist integer >=1 max=1) NUMGen(integer 1000) ///
		PMut(real 0.3) PTransf(real 0.1) SEed(numlist max=1) name(string) ///
		DIGits(integer 3) MAXITGas(integer 30) ]

	/* Options:
	   --------
		 method(string)									--> method to use for assessing unobserved
																				heterogeneity; available methods are
																				'rebus', 'fimix' and 'gas'. default is
																				'rebus'.
		 numclass(integer integer >=1
							max=1))								--> number of classes to use; if empty, it
																				is chosen automatically using a Ward
																				hierarchical algorithm
		 popsize(integer integer >=1
						 max=1))								--> number of individual in each population;
																				this is an option to use with
																				method(gas)
		 numgen(integer 1000)						--> number of generations (default 1000)
		 pmut(real 0.3)									--> probability of mutation (default 0.3)
		 ptransf(real 0.1)							--> probability of transformation (default
																				0.1)
		 stop(real 0.005)								--> stopping criterion (default 0.005)
		 seed(numlist max=1)						--> permutation test seed number
		 name(string)										--> variable name where to store the final
																				rebus classification
		 digits(integer 3)							--> number of digits to display (default 3)
		 maxitgas(integer 30)						--> maximum number of PLS-GAS second stage
																				iterations (default 30)
	 */
	 
	 /* Description:
			------------
			This postestimation command implements the PLS-GAS approach to assess
			the presence of unobserved heterogeneity.
	 */
	
	local props = e(properties)
	local struct "structural"
	local isstruct : list struct in props
	if (!`isstruct') {
		display as error "the fitted plssem model includes only the measurement part"
		error 198
	}
	local rawsum "rawsum"
	local israwsum : list rawsum in props
	local cc "relative"
	local isrelative : list cc in props
	if (`isrelative') {
		local convcrit "relative"
	}
	else {
		local convcrit "square"
	}
	local initmet "indsum"
	local isindsum : list initmet in props
	if (`isindsum') {
		local init "indsum"
	}
	else {
		local init "eigen"
	}
	local indscale "scaled"
	local isscaled : list indscale in props
	if (`isscaled') {
		local scale ""
	}
	else {
		local scale "noscale"
	}
	local boot "bootstrap"
	local isboot : list boot in props
	if (`isboot') {
		display as error "the global model used the 'boot()' option, which slows down excessively the PLS-GAS calculations"
		display as error "try refitting the global model without bootstrap"
		exit
	}
	if ("`seed'" != "") {
		if (`seed' < 0 | `seed' >  2^31-1) {
			display as error "'seed()' option requires a value between 0 and 2^31-1"
			exit
		}
	}
	if ("`name'" != "") & (`: word count `name'' != 1) {
		display as error "the 'name' option must include a single word"
		exit
	}
	if ("`name'" == "") {
		local name "gas_class"
	}
	if ("`e(binarylvs)'" != "") {
		display as error "currently the PLS-GAS implementation doesn't allow for binary latent variables"
		exit
	}
	local knnimp "knn"
	local isknnimp : list knnimp in props
	local meanimp "mean"
	local ismeanimp : list meanimp in props
	if (`isknnimp') {
		local missing "knn"
	}
	else if (`ismeanimp') {
		local missing "mean"
	}
	else {
		local missing ""
	}
	if ("`numclass'" == "") {
		display as error "the 'numclass' option must include one value for PLS-GAS"
		exit
	}
	if ("`popsize'" == "") {
		display as error "the 'popsize' option must include one value for PLS-GAS"
		exit
	}
	if (`numgen' < 1) {
		display as error "the 'numgen' option for PLS-GAS must be an integer at least equal to 1"
		exit
	}
	if ((`pmut' > 1) | (`pmut' < 0)) {
		display as error "the 'pmut' option for PLS-GAS must be in between 0 and 1"
		exit
	}
	if ((`ptransf' > 1) | (`ptransf' < 0)) {
		display as error "the 'ptransf' option for PLS-GAS must be in between 0 and 1"
		exit
	}
	
	tempvar __touse__
	quietly generate byte `__touse__' = e(sample)
	quietly count if `__touse__'
	local N = r(N)

	/* Global model */
	tempname globalmodel gof_global
	_estimates hold `globalmodel', copy
	
	/* Collecting results to display later */
	matrix `gof_global' = e(assessment)
	matrix `gof_global' = `gof_global'[1, 3]

	if ("`missing'" != "") {
		/* Save original data set */
		local allindicators = e(mvs)
		tempname original_data
		mata: `original_data' = st_data(., "`: list uniq allindicators'", "`__touse__'")
		
		/* Recovery of missing values */
		mata: st_store(., tokens("`: list uniq allindicators'"), "`__touse__'", ///
			st_matrix("e(imputed_data)"))
	}

  /* Parse global model e(cmdline) */
  local cmdline = e(cmdline)
  local trash
  gettoken cmdline options : cmdline, parse(",")
  gettoken trash options : options, parse(",")

  local mm_start = strpos(`"`cmdline'"', "(")
  if (strpos(`"`cmdline'"', "if")) {
    local mm_end = strpos(`"`cmdline'"', "if") - `mm_start' - 1
  }
  else {
    local mm_end = strlen(`"`cmdline'"')
  }
  local mm = substr(`"`cmdline'"', `mm_start', `mm_end')
  local mm : list clean mm

  tokenize `"`options'"', parse(")")
  local tok_i = 1
  while (substr(`"``tok_i''"', 1, 3) != "str") & (`"``tok_i''"' != "") {
    local ++tok_i
  }
  local sm_full = `"``tok_i''"' + ")"
  tokenize `"``tok_i''"', parse("(")
  local sm = `"`3'"'
  local sm : list clean sm
  local options : list options - sm_full
  local options : list clean options

  // the following code doesn't work
  /*
  tokenize `"`options'"', parse(")")
  local tok_i = 1
  while (substr(`"``tok_i''"', 1, 3) != "dig") & (`"``tok_i''"' != "") {
    local ++tok_i
  }
  local options_digits = `"``tok_i''"' + ")"
  local options : list options - options_digits
  local options : list clean options
  */
	if (`israwsum') {
		local options "tol(`e(tolerance)') wscheme("
	}
	else {
		local options "tol(`e(tolerance)') maxiter(`e(maxiter)') wscheme("
	}
  local ws_centroid "centroid"
  local ws_factor "factor"
  local ws_path "path"
  if (`: list ws_centroid in props') {
    local scheme "centroid"
  }
  else if (`: list ws_factor in props') {
    local scheme "factor"
  }
  else if (`: list ws_path in props') {
    local scheme "path"
  }
  local options "`options'`scheme')"
  if ("`e(binarylvs)'" != "") {
    local options "`options' binary(`e(binarylvs)')"
  }
  
	/* Set temporary variables */
	local allindicators = e(mvs)
	local alllatents = e(lvs)
	local allreflective = e(reflective)
	local num_ind : word count `allindicators'
	local num_lv : word count `alllatents'
	tempname endo
	mata: `endo' = colsum(st_matrix("e(adj_struct)"))
	mata: `endo' = (`endo' :> 0)
	foreach var in `allindicators' {
		local allstdindicators "`allstdindicators' std`var'"
	}
	local allstdindicators : list clean allstdindicators

	/* Run the PLS-GAS algorithm */
	tempname res_gas
	capture noisily {
		mata: `res_gas' = ///
			plssem_gas( ///
				st_data(., "`allindicators'"), ///			 note: `__touse__' not used here
				st_matrix("e(adj_meas)"), ///
				st_matrix("e(adj_struct)"), ///
				"`allindicators'", ///
				"`allstdindicators'", ///
				"`alllatents'", ///
				"`e(binarylvs)'", ///
				st_numscalar("e(tolerance)"), ///
				st_numscalar("e(maxiter)"), ///
				"`__touse__'", ///
				"`scheme'", ///
				"`convcrit'", ///
				"`init'", ///
				"`scale'", ///
				strtoreal("`isstruct'"), ///
				strtoreal("`israwsum'"), ///
				0, ///
				strtoreal("`numclass'"), ///
				strtoreal("`popsize'"), ///
				strtoreal("`numgen'"), ///
				strtoreal("`pmut'"), ///
				strtoreal("`ptransf'"), ///
				strtoreal("`maxitgas'"), ///
				strtoreal("`seed'"), ///
				1)
		
		mata: st_local("iter", strofreal(`res_gas'.niter))
	}
	if (_rc != 0) {
		if (mod(`iter', 5) == 0) {
			display as error " aborting"
		}
		else {
			display as error "aborting"
		}
		if (_rc == 409) {
			display as error "at least one indicator has zero variance in one iteration"
		}
		else {
			display as error "something went wrong in the PLS-GAS calculations"
		}
		display as error "try reducing the number of classes " _continue
		display as error "or relaxing the stopping criteria"
		exit
	}
	/* End of PLS-GAS algorithm */

	/* Save final classification */
	tempvar tmp_class
	quietly generate int `tmp_class' = .
	mata: st_store(., "`tmp_class'", `res_gas'.gas_class)

	/* Once stability is attained, final local models are estimated */
	local donotcleanup "nocleanup"
	forvalues k = 1/`numclass' {
		tempname localmodel_`k'
		quietly plssem `mm' if (`tmp_class' == `k' & `__touse__'), ///
			structural(`sm') `options' `donotcleanup'
		_estimates hold `localmodel_`k''
	}
	
	/* Restore global model results */
	_estimates unhold `globalmodel'

  /* Display results */
	display
	
	local skip1 = 1
	local skip3 = 3

	tempname tmp strb nonmiss alllen numeff
	matrix `tmp' = e(struct_b)
	matrix `strb' = vec(`tmp')
	local nrows = rowsof(`tmp')
	local ncols = colsof(`tmp')
	local totrows = `nrows'*`ncols'
	matrix `nonmiss' = J(`totrows', 1, 0)
	local strb_rn : rowfullnames `strb'
	forvalues i = 1/`totrows' {
		if (!missing(`strb'[`i', 1])) {
			matrix `nonmiss'[`i', 1] = 1
			local nm_tmp `: word `i' of `strb_rn''
			local tok_i = 1
			tokenize `"`nm_tmp'"', parse(":")
			while ("``tok_i''" != "") {
				if (`tok_i' == 1) {
					local nm_Y "``tok_i''"
				}
				else if (`tok_i' == 3) {
					local nm_X "``tok_i''"
				}
				local ++tok_i
			}
			local strbok_rn "`strbok_rn' `nm_X':`nm_Y'"
		}
	}	
	mata: st_numscalar("`numeff'", colsum(st_matrix("`nonmiss'")))
	local neff = `numeff'
	matrix `alllen' = J(`neff', 1, .)
	local resnm : subinstr local strbok_rn ":" "->", all
	forvalues j = 1/`neff' {
		local nm : word `j' of `resnm'
		local nm_len : strlen local nm
		matrix `alllen'[`j', 1] = `nm_len' + 2
	}

	tempname results_n results_p results_l tmp_mat
	local nind : word count `allindicators'
	forvalues k = 1/`numclass' {
		local grp_cn `grp_cn' "Class_`k'"
	}
	matrix `results_n' = J(3, 1 + `numclass', .)
	matrix rownames `results_n' = "Observations" "Percentage" "GoF"
	matrix colnames `results_n' = "Global" `grp_cn'
	matrix `results_p' = J(`neff', 1 + `numclass', .)
	matrix rownames `results_p' = `resnm'
	matrix colnames `results_p' = "Global" `grp_cn'
	matrix `results_l' = J(`nind', 1 + `numclass', .)
	matrix rownames `results_l' = `allindicators'
	matrix colnames `results_l' = "Global" `grp_cn'

	_estimates hold `globalmodel', copy  // needed to keep it in memory
	matrix `tmp' = e(struct_b)
	matrix `strb' = vec(`tmp')

	tempname global_N
	scalar `global_N' = e(N)
	matrix `results_n'[1, 1] = e(N)
	matrix `results_n'[2, 1] = 100
	matrix `results_n'[3, 1] = `gof_global'[1, 1]

	local b_i = 1
	forvalues i = 1/`totrows' {
		if (`nonmiss'[`i', 1]) {
			matrix `results_p'[`b_i', 1] = `strb'[`i', 1]
			local ++b_i
		}
	}
	mata: st_matrix("`tmp_mat'", rowsum(st_matrix("e(loadings)")))
	matrix `results_l'[1, 1] = `tmp_mat'
	
	tempname gof_tmp
	forvalues k = 1/`numclass' {
		_estimates unhold `localmodel_`k''
		matrix `tmp' = e(struct_b)
		matrix `strb' = vec(`tmp')

		matrix `gof_tmp' = e(assessment)
		matrix `gof_tmp' = `gof_tmp'[1, 3]
		matrix `results_n'[1, 1 + `k'] = e(N)
		matrix `results_n'[2, 1 + `k'] = e(N)/`global_N'*100
		matrix `results_n'[3, 1 + `k'] = `gof_tmp'[1, 1]
		
		local b_i = 1
		forvalues i = 1/`totrows' {
			if (`nonmiss'[`i', 1]) {
				matrix `results_p'[`b_i', 1 + `k'] = `strb'[`i', 1]
				local ++b_i
			}
		}

		mata: st_matrix("`tmp_mat'", rowsum(st_matrix("e(loadings)")))
		matrix `results_l'[1, 1 + `k'] = `tmp_mat'
	}
	
	mkheader, digits(5) gas

	tempname maxlen
	mata: st_numscalar("`maxlen'", max(st_matrix("`alllen'")))
  local firstcollbl ""	
	local `maxlen' = max(max(strlen("`firstcollbl'"), `maxlen') + 2, 14)
	local colw = max(9, `digits' + 5)

	local title "PLS-GAS classes"
	local firstcollbl ""
	mktable, matrix(`results_n') digits(`digits') firstcolname(`firstcollbl') ///
		title(`title') firstcolwidth(``maxlen'') colwidth(`colw') hlines(3) ///
		novlines total rebus

	local title "Loadings"
	local firstcollbl ""
	mktable, matrix(`results_l') digits(`digits') firstcolname(`firstcollbl') ///
		title(`title') firstcolwidth(``maxlen'') colwidth(`colw') hlines(`nind') ///
		novlines total

	local title "Path coefficients"
	local firstcollbl ""
	mktable, matrix(`results_p') digits(`digits') firstcolname(`firstcollbl') ///
		title(`title') firstcolwidth(``maxlen'') colwidth(`colw') hlines(`neff') ///
		novlines total
	/* End of display */

	/* Restore global model results */
	_estimates unhold `globalmodel'

	/* Save final classification */
	capture drop `name'
	quietly generate int `name' = .
	mata: st_store(., "`name'", `res_gas'.gas_class)
	local now "`c(current_date)', `c(current_time)'"
	local now : list clean now
	label variable `name' "PLS-GAS classification [`now']"

	/* Clean up */
	if ("`missing'" != "") {
		mata: st_store(., tokens("`: list uniq allindicators'"), "`__touse__'", ///
			`original_data')
	}
	
	/* Return values */
	return scalar nclasses = `numclass'
	return scalar niter = `iter'
end

program mediate, rclass
	version 15.1
	syntax , INDep(string) MED(string) DEP(string) ///
		[ BReps(numlist integer >0 min=1 max=1) Seed(numlist max=1) ZLC RIT RID ///
		BCa Level(real 0.95) DIGits(integer 3) ]
	
	/* Options:
	   --------
		 indep(string)							--> independent variable
		 med(string)								--> mediator variable
		 dep(string)								--> dependent variable
		 breps(numlist integer >0 min=1 max=1)
																--> number of bootstrap replications (default
																		is 50)
		 seed(numlist max=1)				--> bootstrap seed number
		 zlc												--> mediation procedures described by Zhao
																		et al. (2010)
		 rit												--> ratio of the indirect effect to the total
																		effect
		 rid												--> ratio of the indirect effect to the direct
																		effect
		 bca												--> returns bias-corrected accelerated bootstrap
																		confidence intervals instead of percentile
																		confidence intervals, which is the default
		 level(real 0.95)						--> confidence level (default 0.95)
		 digits(integer 3)					--> number of digits to display (default 3)
	 */
	 
	 /* Description:
			------------
			This postestimation command provides the details for a single indirect
			effect mediated by only one LV.
	 */

	tempname coef coef_el
	matrix `coef' = e(struct_b)
	
	local exo_nm : rownames `coef'
	local endo_nm : colnames `coef'
	
	if (!`: list dep in endo_nm') {
		display as error "`dep' is not an endogenous variable"
		exit
	}
	if (!`: list med in endo_nm') {
		display as error "`med' is not an endogenous variable"
		exit
	}
	if (!`: list indep in exo_nm') {
		display as error "`indep' is not an exogenous variable"
		exit
	}

	local coef_el = `coef'[rownumb(`coef', "`indep'"), colnumb(`coef', "`dep'")]
	if (missing(`coef_el')) {
		display as error "the `indep' variable has no direct effect on `dep'"
		exit
	}
	local coef_el = `coef'[rownumb(`coef', "`med'"), colnumb(`coef', "`dep'")]
	if (missing(`coef_el')) {
		display as error "the `med' variable has no direct effect on `dep'"
		exit
	}
	local coef_el = `coef'[rownumb(`coef', "`indep'"), colnumb(`coef', "`med'")]
	if (missing(`coef_el')) {
		display as error "the `indep' variable has no direct effect on `med'"
		exit
	}

	local struct "structural"
	local props = e(properties)
	local isstruct : list struct in props
	if (!`isstruct') {
		display as error "the fitted plssem model includes only the measurement part"
		exit
	}

	tempvar __touse__
	quietly generate `__touse__' = e(sample)
	local reg3eqs = e(struct_eqs)
	if (`level' <= 0 | `level' >= 1) {
		display as error "confidence level must be in the range (0, 1)"
		exit
	}

	if ("`breps'" == "") {
		local breps = 50
	}

	tempname normperc alpha_cl
	scalar `alpha_cl' = 1 - ((1 - `level')/2)
	scalar `normperc' = invnormal(`alpha_cl')

	tempname ehold
	_estimates hold `ehold'

	display as text "Bootstrapping mediation effect..."

	capture {
		if ("`dep'" != "" & "`med'" != "" & "`indep'" != "") {
			local doi `dep':`indep'
			local moi `med':`indep'
			local dom `dep':`med'

			tempname reg3coef_b reg3var_b reg3ciP reg3ciBC reg3coef reg3var
			quietly bootstrap indeff=(_b[`dom']*_b[`moi']), reps(`breps') ///
				seed(`seed'): reg3 `reg3eqs' if `__touse__', mvreg corr(independent)
			quietly estat bootstrap, all
			matrix `reg3coef_b' = e(b)
			matrix `reg3var_b' = e(V)
			matrix `reg3ciP' = e(ci_percentile)
			matrix `reg3ciBC' = e(ci_bc)

			tempname perc_lci perc_uci bc_lci bc_uci
			scalar `perc_lci' = `reg3ciP'[1, 1]
			scalar `perc_uci' = `reg3ciP'[2, 1]
			scalar `bc_lci' = `reg3ciBC'[1, 1]
			scalar `bc_uci' = `reg3ciBC'[2, 1]

			quietly reg3 `reg3eqs' if `__touse__', mvreg corr(independent)
			matrix `reg3coef' = e(b)
			matrix `reg3var' = e(V)

			tempname moimat moicoef moimat2 moivar dommat domcoef dommat2 domvar
			matrix `moimat' = `reg3coef'[1, "`moi'"]
			scalar `moicoef' = `moimat'[1, 1]
			matrix `moimat2' = `reg3var'["`moi'", "`moi'"]
			scalar `moivar' = `moimat2'[1, 1]
			matrix `dommat' = `reg3coef'[1, "`dom'"]
			scalar `domcoef' = `dommat'[1, 1]
			matrix `dommat2' = `reg3var'["`dom'", "`dom'"]
			scalar `domvar' = `dommat2'[1, 1]
		}
		else {
			display as error "provide a single dependent, independent and mediator variable for each indirect effect"
			exit
		}
	} // end of -capture-
	local rc = _rc
	error `rc'
	
	/* Description of mediating effect */
	tempname indmeas mat_bk mat_zlc mat_rit mat_rid
	matrix `indmeas' = J(6, 3, .)
	matrix `mat_bk' = J(7, 1, .)
	matrix `mat_zlc' = J(4, 1, .)
	matrix `mat_rit' = J(3, 1, .)
	matrix `mat_rid' = J(3, 1, .)
	
	tempname coef_moi se_moi var_moi coef_dom se_dom var_dom prodterm
	scalar `coef_moi' = `moicoef'
	scalar `var_moi' = `moivar'
	scalar `se_moi' = sqrt(`var_moi')
	scalar `coef_dom' = `domcoef'
	scalar `var_dom' = `domvar'
	scalar `se_dom' = sqrt(`var_dom')
	scalar `prodterm' = `coef_moi'*`coef_dom'

	// Sobel's method
	tempname sobel_se sobel_z sobel_pv sobel_lci sobel_uci
	scalar `sobel_se' = sqrt(((`coef_dom')^2)*`var_moi' + ((`coef_moi')^2)*`var_dom')
	scalar `sobel_z' = `prodterm'/`sobel_se'
	scalar `sobel_pv' =  2*(1 - normal(abs(`sobel_z')))
	scalar `sobel_lci' = `prodterm' - `normperc'*`sobel_se'
	scalar `sobel_uci' = `prodterm' + `normperc'*`sobel_se'

	matrix `indmeas'[1, 1] = `prodterm'
	matrix `indmeas'[2, 1] = `sobel_se'
	matrix `indmeas'[3, 1] = `sobel_z'
	matrix `indmeas'[4, 1] = `sobel_pv'
	matrix `indmeas'[5, 1] = `sobel_lci'
	matrix `indmeas'[6, 1] = `sobel_uci'
	
	// Delta method
	tempname v delta_var delta_se delta_z delta_pv delta_lci delta_uci
	quietly nlcom _b[`moi']*_b[`dom']
	matrix `v' = r(V)
	scalar `delta_var' = `v'[1, 1]
	scalar `delta_se' = sqrt(`delta_var')
	scalar `delta_z' = `prodterm'/`delta_se'
	scalar `delta_pv' =  2*(1 - normal(abs(`delta_z')))
	scalar `delta_lci' = `prodterm' - `normperc'*`delta_se'
	scalar `delta_uci' = `prodterm' + `normperc'*`delta_se'

	matrix `indmeas'[1, 2] = `prodterm'
	matrix `indmeas'[2, 2] = `delta_se'
	matrix `indmeas'[3, 2] = `delta_z'
	matrix `indmeas'[4, 2] = `delta_pv'
	matrix `indmeas'[5, 2] = `delta_lci'
	matrix `indmeas'[6, 2] = `delta_uci'
	
	// Bootstrap method
	tempname boot_prod boot_se boot_z boot_pv
	scalar `boot_prod' = `reg3coef_b'[1, 1]
	scalar `boot_se' = sqrt(`reg3var_b'[1, 1])
	scalar `boot_z' = `boot_prod'/`boot_se'
	scalar `boot_pv' =  2*(1 - normal(abs(`boot_z'))) 						// CHECK THIS!!!
	
	matrix `indmeas'[1, 3] = `boot_prod'
	matrix `indmeas'[2, 3] = `boot_se'
	matrix `indmeas'[3, 3] = `boot_z'
	matrix `indmeas'[4, 3] = `boot_pv'
	if ("`bca'" == "") {
		matrix `indmeas'[5, 3] = `perc_lci'
		matrix `indmeas'[6, 3] = `perc_uci'
	}
	else {
		matrix `indmeas'[5, 3] = `bc_lci'
		matrix `indmeas'[6, 3] = `bc_uci'
	}
	
	matrix rownames `indmeas' = "Indirect effect" "Standard error" ///
		"Z statistic" "P-value" "Lower CI" "Upper CI"
	matrix colnames `indmeas' = "Sobel" "Delta" "Bootstrap"

	/* Clean up */
	_estimates unhold `ehold'

	/* Baron-Kenny mediation testing adjusted by Iacobucci et al. */
	tempname coef st_tab pval stderr
	matrix `coef' = e(struct_b)
	local struct_b_rn : rownames `coef'
	local struct_b_cn : colnames `coef'
	matrix `st_tab' = e(struct_table)
	local st_nrow : rowsof(`st_tab')
	local st_nrow = (`st_nrow' - 1)/2
	local st_ncol : colsof(`st_tab')
	matrix `pval' = J(`st_nrow', `st_ncol', .)
	forvalues i = 1/`st_nrow' {
		matrix `pval'[`i', 1] = `st_tab'[`i'*2, 1...]
	}
	matrix rownames `pval' = `struct_b_rn'
	matrix colnames `pval' = `struct_b_cn'
	mat `stderr' = e(struct_se)

	// X -> M
	tempname coef_moi_m coef_moi moi_pval_m moi_pval
	matrix `coef_moi_m' = `coef'["`indep'", "`med'"]
	scalar `coef_moi' = `coef_moi_m'[1, 1]
	matrix `moi_pval_m' = `pval'["`indep'", "`med'"]
	scalar `moi_pval' = `moi_pval_m'[1, 1]
	
	// M -> Y
	tempname coef_dom_m coef_dom dom_pval_m dom_pval
	matrix `coef_dom_m' = `coef'["`med'", "`dep'"]
	scalar `coef_dom' = `coef_dom_m'[1, 1]
	matrix `dom_pval_m' = `pval'["`med'", "`dep'"]
	scalar `dom_pval' = `dom_pval_m'[1, 1]
	
	// X -> Y
	tempname coef_doi_m coef_doi doi_pval_m doi_pval
	matrix `coef_doi_m' = `coef'["`indep'", "`dep'"]
	scalar `coef_doi' = `coef_doi_m'[1, 1]
	matrix `doi_pval_m' = `pval'["`indep'", "`dep'"]
	scalar `doi_pval' = `doi_pval_m'[1, 1]
	
	matrix `mat_bk'[1, 1] = `coef_doi'
	matrix `mat_bk'[2, 1] = `coef_moi'
	matrix `mat_bk'[3, 1] = `coef_dom'
	matrix `mat_bk'[4, 1] = `doi_pval'
	matrix `mat_bk'[5, 1] = `moi_pval'
	matrix `mat_bk'[6, 1] = `dom_pval'
	matrix `mat_bk'[7, 1] = `sobel_pv'	

	/* Zhao et al. mediation testing */
	if ("`zlc'" != "") {
		tempname axbxc
		scalar `axbxc' = `coef_moi'*`coef_dom'*`coef_doi'
		
		matrix `mat_zlc'[1, 1] = `boot_pv'
		matrix `mat_zlc'[2, 1] = `doi_pval'
		matrix `mat_zlc'[3, 1] = `coef_doi'
		matrix `mat_zlc'[4, 1] = `axbxc'
	}
	
	if ("`rit'" != "") {
		tempname prod toteff
		scalar `prod' = abs(`prodterm')
		scalar `toteff' = abs(`prodterm' + `coef_doi')
		
		matrix `mat_rit'[1, 1] = `prod'
		matrix `mat_rit'[2, 1] = `toteff'
		matrix `mat_rit'[3, 1] = `prod'/`toteff'
	}
	
	if ("`rid'" != "") {
		tempname prod absprod
		scalar `prod' = abs(`prodterm')
		scalar `absprod' = abs(`coef_doi')
		
		matrix `mat_rid'[1, 1] = `prod'
		matrix `mat_rid'[2, 1] = `absprod'
		matrix `mat_rid'[3, 1] = `prod'/`absprod'
	}
	
	/* Display results */
	local med_line = 5
	mktable_mediate, matrix(`indmeas') matrix_bk(`mat_bk') ///
		matrix_zlc(`mat_zlc') matrix_rit(`mat_rit') matrix_rid(`mat_rid') ///
		depv("`dep'") medv("`med'") indepv("`indep'") ///
		title("Significance testing of (standardized) indirect effect") ///
		firstcolname("Statistics") firstcolwidth(24) colwidth(20) ///
		hlines(`med_line') reps(`breps') digits(`digits') level(`level') novlines
	
	/* Return values */
	return matrix mediate `indmeas'
end
