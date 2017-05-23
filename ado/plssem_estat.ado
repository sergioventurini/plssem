*!plssem_estat version 0.2.0
*!Written 23May2017
*!Written by Sergio Venturini and Mehmet Mehmetoglu
*!The following code is distributed under GNU General Public License version 3 (GPL-3)

program plssem_estat, rclass
	version 10
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
	else {
		// estat_default `0'
		display as error "the `subcmd' postestimation command is not implemented for plssem"
		exit
	}

	return add
end

program indirect, rclass
	version 10
	syntax , Effects(string) [ Boot(numlist min=1 max=1) Seed(numlist max=1) ///
		Level(real 0.95) DIGits(integer 3) ]
	
	/* Options:
	   --------
		 effects(string)						--> list of indirect effects
		 boot(numlist min=1 max=1)	--> bootstrap estimation (# of repetions;
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

	tempname sobel_se sobel_z sobel_pv indirect sobel_lci sobel_uci
	tempname indmeas reg3coef reg3coef_bs reg3var
	tempname dommat dommat2 moimat moimat2 //doimat doimat2
	if ("`boot'" == "") {
		matrix `indmeas' = J(6, `num_effects', .)
	}
	else {
		matrix `indmeas' = J(10, `num_effects', .)
	}

	tempname ehold
	_estimates hold `ehold'
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
					local reg3eq_nm : colfullnames `reg3coef'

					tempname moicoef moivar domcoef domvar
					matrix `moimat' = `reg3coef'[1, "`moi'"]
					scalar `moicoef' = `moimat'[1, 1]
					matrix `moimat2' = `reg3var'["`moi'","`moi'"]
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
					matrix `reg3coef_bs' = e(b_bs) // this is not used for now
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
	return matrix indirect = `indmeas'
end

program total, rclass
	version 10
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
	return matrix total = `alleffects_new'
end

program plssem_vif, rclass
	version 10
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
	return matrix strvif = `strvif'
end

program unobshet, rclass
	version 10
	syntax [ , Method(string) Numclass(numlist integer >1 max=1) ///
		MAXCLass(integer 20) Dendrogram MAXITer(integer 50) Stop(real 0.005) ///
		Test Reps(numlist integer >1 max=1) SEed(numlist max=1) Plot ///
		name(string) DIGits(integer 3) ]
	
	/* Options:
	   --------
		 method(string)									--> method to use for assessing unobserved
																				heterogeneity (default is 'rebus')
		 numclass(integer)							--> number of classes to use; if empty, it
																				is chosen automatically using a Ward
																				hierarchical algorithm
		 maxclass(integer 20)						--> maximum number of classes to test in
																				choosing automatically the number of
																				classes (default 20)
		 dendrogram											--> display the dendrogram for the Ward
																				cluster analysis
		 maxiter(integer 50)						--> maximum number of iterations (default
																				50)
		 stop(real 0.05)								--> stopping criterion (default 0.005)
		 test														--> permutation test
		 reps(numlist integer >1 max=1)	--> number of permutation test replications
																				(default is 50)
		 seed(numlist max=1)						--> permutation test seed number
		 plot														--> plot of the empirical distribution
																				for the permutation test statistic
		 name(string)										--> variable name where to store the final
																				rebus classification
		 digits(integer 3)							--> number of digits to display (default 3)
	 */
	 
	 /* Description:
			------------
			This postestimation command provides various methods to assess the
			presence of unobserved heterogeneity.
			
			Currently it implements only the REBUS approach.
	 */
	
	if ("`method'" == "") {
		local method "rebus"
	}
	if ("`method'" != "rebus") {
		display as error "currently 'estat unobshet' implements only the REBUS method"
		exit
	}
	if ("`e(formative)'" != "") {
		display as error "the REBUS approach can't be used with formative blocks; " _continue
		display as error "only reflective blocks are allowed"
		exit
	}
	local props = e(properties)
	local struct "structural"
	local isstruct : list struct in props
	if (!`isstruct') {
		display as error "the fitted plssem model includes only the measurement part"
		exit
	}
	local boot "bootstrap"
	local isboot : list boot in props
	if (`isboot') {
		display as error "the global model used the 'boot()' option, which slows down excessively the REBUS calculations"
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
		display as error "currently the REBUS approach implementation doesn't allow for binary latent variables"
		exit
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

	// collect results for later to display
	matrix `gof_global' = e(assessment)
	matrix `gof_global' = `gof_global'[1, 3]
	local gof_loc = strofreal(`gof_global'[1, 1], "%9.4f")

	tempvar __touse__
	quietly generate byte `__touse__' = e(sample)
	quietly count if `__touse__'
	if (r(N) == 0) {
		error 2000
	}
	local N = r(N)

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
				title("Dendrogram for PLS-SEM REBUS analysis") scheme(sj) ///
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
	local options "tol(`e(tolerance)') maxiter(`e(maxiter)') wscheme("
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
	local num_lv : word count `e(lvs)'
	tempname cm ow path ind indstd y y_local loads r2 ind_all indstd_all block ///
		x_hat out_res endo y_hat inn_res old_class old_class_i new_class counts
	tempvar __touseloc__
	quietly generate byte `__touseloc__' = .
	mata: `endo' = colsum(st_matrix("e(adj_struct)"))
	mata: `endo' = (`endo' :> 0)
	local iter = 1
	local nchanged = `N'

	/* REBUS iterations */
	display
	capture {
		while (`iter' <= `maxiter') & (`nchanged' > `N'*`stop') {
			if (`iter' == 1) {
				noisily display as text "REBUS iteration 1" _continue
			}
			else if (mod(`iter', 5) == 0) {
				noisily display as text "`iter'" _continue
			}
			else if (`iter') {
				noisily display as text "." _continue
			}

			mata: `cm' = J(strtoreal(st_local("N")), 0, .)
			forvalues k = 1/`numclass' {
				tempname localmodel_`k'
				quietly plssem `mm' if (`rebus_class' == `k' & `__touse__'), ///
					structural(`sm') `options'
				_estimates hold `localmodel_`k'', copy
				quietly drop `__touseloc__'
				quietly generate byte `__touseloc__' = e(sample)

				mata: st_view(`ind' = ., ., "`allindicators'", "`__touseloc__'")					// split.DM[[k]]

				mata: `ow' = st_matrix("e(outerweights)")																	// w.locals[[k]]
				mata: `path' = st_matrix("e(pathcoef)")																		// path.locals[[k]]
				mata: `loads' = st_matrix("e(loadings)")																	// loads.locals[[k]]
				mata: `r2' = st_matrix("e(rsquared)")	
				mata: `r2' = `r2'[., selectindex(`r2' :!= .)]															// R2.locals[[k]]

				// mata: `indstd' = scale(`ind') 																					// split.X[[k]]
				// mata: `y' = `indstd' * `ow'																						// Y.k

				mata: st_view(`ind_all' = ., ., "`allindicators'", "`__touse__'")					// DM
				mata: `indstd_all' = scale(`ind_all', mean(`ind'), sd(`ind'))							// X.k
				mata: `y_local' = `indstd_all' * `ow'																			// Y.locals[[k]]
				mata: `out_res' = J(strtoreal(st_local("N")), 0, .)
				forvalues j = 1/`num_lv' {
					mata: `block' = selectindex(`loads'[., `j'] :!= .)											// blocklist == j
					mata: `x_hat' = `y_local'[., `j'] * `loads'[`block', `j']'							// X.hat
					mata: `out_res' = (`out_res', (`indstd_all'[., `block'] - `x_hat'))			// out.res[, blocklist == j]
				}
				mata: `y_hat' = `y_local' * `path'[., selectindex(`endo')]								// Y.hat
				mata: `inn_res' = `y_local'[., selectindex(`endo')] - `y_hat'							// innres.locals[[k]]
				mata: `cm' = (`cm', rebus_cm(`out_res', `inn_res', `loads', `r2'))				// CM.locals[, k]
				
				_estimates drop `localmodel_`k''
			}

			mata: st_view(`old_class' = ., ., "`rebus_class'", "`__touse__'")
			mata: `old_class_i' = st_viewobs(`old_class')
			mata: `new_class' = which(`cm', "min")
			mata: st_local("nchanged", strofreal(sum(`old_class' :!= `new_class')))
			mata: st_store(`old_class_i', "`rebus_class'", ., `new_class')
			forvalues k = 1/`numclass' {
				quietly count if `rebus_class' == `k'
				if (r(N) == 0) {
					local rN0 = 1
					error 148
				}
				else {
					local rN0 = 0
				}
				if (r(N) <= 5) {
					local rN_lte_5 = 1
					error 148
				}
				else {
					local rN_lte_5 = 0
				}
			}

			local ++iter
		}
	}
	if (_rc != 0) {
		if (mod(`iter', 5) == 0) {
			display as error " aborting"
		}
		else {
			display as error "aborting"
		}
		if (`rN0') {
			display as error "one of the classes is empty"
		}
		else if (`rN_lte_5') {
			display as error "too few observations (5 or less) assigned to a single class"
		}
		else {
			display as error "something went wrong in the REBUS calculations"
		}
		display as error "try reducing the number of classes " _continue
		display as error "or relaxing any of the stopping criteria"
		restore
		_estimates unhold `globalmodel'
		exit
	}
	else {
		if (mod(`iter' - 1, 5) == 0) {
			display as text " done!"
		}
		else {
			display as text "done!"
		}
	}
	display
	/* End REBUS iterations */

	/* Once stability is attained, final local models are estimated */
	local allendogenous "`: colnames e(struct_b)'"
	tempname lat class indstd_st lat_st out_res_st inn_res_st class_st gqi
	mata: `indstd_st' = J(0, cols(`ind'), .)
	mata: `lat_st' = J(0, sum(`endo'), .)
	mata: `out_res_st' = J(0, cols(`ind'), .)
	mata: `inn_res_st' = J(0, sum(`endo'), .)
	mata: `class_st' = J(0, 1, .)
	forvalues k = 1/`numclass' {
		quietly plssem `mm' if (`rebus_class' == `k' & `__touse__'), ///
			structural(`sm') `options'
		_estimates hold `localmodel_`k'', copy
		quietly drop `__touseloc__'
		quietly generate byte `__touseloc__' = e(sample)

		mata: st_view(`lat' = ., ., "`allendogenous'", "`__touseloc__'")
		mata: st_view(`ind' = ., ., "`allindicators'", "`__touseloc__'")						// split.DM[[k]]
		mata: `indstd' = scale(`ind') 																							// split.X[[k]]

		mata: `ow' = st_matrix("e(outerweights)")																		// w.locals[[k]]
		mata: `path' = st_matrix("e(pathcoef)")																			// path.locals[[k]]
		mata: `loads' = st_matrix("e(loadings)")																		// loads.locals[[k]]
		mata: `r2' = st_matrix("e(rsquared)")	
		mata: `r2' = `r2'[., selectindex(`r2' :!= .)]																// R2.locals[[k]]

		mata: `y_local' = `indstd' * `ow'																						// LV.locals[[k]]
		mata: `out_res' = J(rows(`ind'), 0, .)
		forvalues j = 1/`num_lv' {
			mata: `block' = selectindex(`loads'[., `j'] :!= .)												// blocklist == j
			mata: `x_hat' = `y_local'[., `j'] * `loads'[`block', `j']'								// X.hat
			mata: `out_res' = (`out_res', (`indstd'[., `block'] - `x_hat'))						// out.res[, blocklist == j]
		}
		mata: `y_hat' = `y_local' * `path'[., selectindex(`endo')]									// Y.hat
		mata: `inn_res' = `y_local'[., selectindex(`endo')] - `y_hat'								// innres.locals[[k]]
		mata: `class' = J(rows(`ind'), 1, `k')

		mata: `indstd_st' = (`indstd_st' \ `indstd')
		mata: `lat_st' = (`lat_st' \ `lat')
		mata: `out_res_st' = (`out_res_st' \ `out_res')
		mata: `inn_res_st' = (`inn_res_st' \ `inn_res')
		mata: `class_st' = (`class_st' \ `class')
	}
	mata: `gqi' = rebus_gqi(`indstd_st', `lat_st', `out_res_st', `inn_res_st', `class_st')
	tempname gqi_final
	mata: st_numscalar("`gqi_final'", `gqi')
	/* End final REBUS calculations */

	/* Permutation test */
	if ("`test'" != "") {
		tempname rebclass rebclass_p gqi_perm
		if ("`reps'" == "") {
			local reps = 50
		}
		if ("`seed'" != "") {
			set seed `seed'
		}
		quietly generate `rebclass_p' = .
		mata: `rebclass' = st_data(., "`rebus_class'", "`__touse__'")
		mata: `gqi_perm' = J(`reps', 1, .)

		noisily display as text "Permutation replications (" _continue
		noisily display as result string(`reps') _continue
		noisily display as text ")"
		noisily display as text ///
			"{hline 4}{c +}{hline 3}" " 1 " ///
			"{hline 3}{c +}{hline 3}" " 2 " ///
			"{hline 3}{c +}{hline 3}" " 3 " ///
			"{hline 3}{c +}{hline 3}" " 4 " ///
			"{hline 3}{c +}{hline 3}" " 5"
		
		local tmp_skip = 0
		forvalues m = 1/`reps' {
			if (mod(`m', 50) == 0) {
				local todisp = strtrim(string(`m', "%9.0f"))
				if (strlen("`todisp'") < strlen(string(`reps'))) {
					local tmp_skip = strlen(string(`reps')) - strlen("`todisp'")
				}
				noisily display as text "." _skip(3) _skip(`tmp_skip') "`todisp'"
				local tmp_skip = 0
			}
			else {
				noisily display as text "." _continue
			}

			mata: st_store(., "`rebclass_p'", "`__touse__'", jumble(`rebclass'))

			mata: `indstd_st' = J(0, cols(`ind'), .)
			mata: `lat_st' = J(0, sum(`endo'), .)
			mata: `out_res_st' = J(0, cols(`ind'), .)
			mata: `inn_res_st' = J(0, sum(`endo'), .)
			mata: `class_st' = J(0, 1, .)
			forvalues k = 1/`numclass' {
				quietly plssem `mm' if (`rebclass_p' == `k' & `__touse__'), ///
					structural(`sm') `options'
				quietly drop `__touseloc__'
				quietly generate byte `__touseloc__' = e(sample)

				mata: st_view(`lat' = ., ., "`allendogenous'", "`__touseloc__'")
				mata: st_view(`ind' = ., ., "`allindicators'", "`__touseloc__'")				// split.DM[[k]]
				mata: `indstd' = scale(`ind') 																					// split.X[[k]]

				mata: `ow' = st_matrix("e(outerweights)")																// w.locals[[k]]
				mata: `path' = st_matrix("e(pathcoef)")																	// path.locals[[k]]
				mata: `loads' = st_matrix("e(loadings)")																// loads.locals[[k]]
				mata: `r2' = st_matrix("e(rsquared)")	
				mata: `r2' = `r2'[., selectindex(`r2' :!= .)]														// R2.locals[[k]]

				mata: `y_local' = `indstd' * `ow'																				// LV.locals[[k]]
				mata: `out_res' = J(rows(`ind'), 0, .)
				forvalues j = 1/`num_lv' {
					mata: `block' = selectindex(`loads'[., `j'] :!= .)										// blocklist == j
					mata: `x_hat' = `y_local'[., `j'] * `loads'[`block', `j']'						// X.hat
					mata: `out_res' = (`out_res', (`indstd'[., `block'] - `x_hat'))				// out.res[, blocklist == j]
				}
				mata: `y_hat' = `y_local' * `path'[., selectindex(`endo')]							// Y.hat
				mata: `inn_res' = `y_local'[., selectindex(`endo')] - `y_hat'						// innres.locals[[k]]
				mata: `class' = J(rows(`ind'), 1, `k')

				mata: `indstd_st' = (`indstd_st' \ `indstd')
				mata: `lat_st' = (`lat_st' \ `lat')
				mata: `out_res_st' = (`out_res_st' \ `out_res')
				mata: `inn_res_st' = (`inn_res_st' \ `inn_res')
				mata: `class_st' = (`class_st' \ `class')
			}
			mata: `gqi_perm'[`m', 1] = ///
				rebus_gqi(`indstd_st', `lat_st', `out_res_st', `inn_res_st', `class_st')
		}
		
		local oldN = _N
		if (`oldN' < `reps') {
			quietly set obs `reps'
		}
		tempname pvalue_sc permdist
		tempvar gqi_dist
		quietly generate `gqi_dist' = .
		mata: st_store(range(1, `reps', 1), "`gqi_dist'", `gqi_perm')
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
		display
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
	
	local iterm1 = `iter' - 1
	local gqi_final_loc = `gqi_final'
	mkheader, digits(5) rebus_it(`iterm1') rebus_gqi(`gqi_final_loc')

	tempname maxlen
	mata: st_numscalar("`maxlen'", max(st_matrix("`alllen'")))
	local `maxlen' = max(strlen("`firstcollbl'"), `maxlen') + 2
	if (`numclass' == 2) {
		local colw = 11
	}
	else {
		local colw = 9
	}

	local title "REBUS classes"
	local firstcollbl ""
	mktable, matrix(`results_n') digits(`digits') firstcolname(`firstcollbl') ///
		title(`title') firstcolwidth(``maxlen'') colwidth(`colw') hlines(3) ///
		novlines total rebus

	local title "Path coefficients"
	local firstcollbl ""
	mktable, matrix(`results_p') digits(`digits') firstcolname(`firstcollbl') ///
		title(`title') firstcolwidth(``maxlen'') colwidth(`colw') hlines(`neff') ///
		novlines total

	local title "Loadings"
	local firstcollbl ""
	mktable, matrix(`results_l') digits(`digits') firstcolname(`firstcollbl') ///
		title(`title') firstcolwidth(``maxlen'') colwidth(`colw') hlines(`nind') ///
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
	/* End display */

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
	label variable `name' "REBUS classification [`now']"

	/* Return values */
	return scalar nclasses = `numclass'
	return scalar GQI = `gqi_final'
	return scalar GoF = `gof_global'[1, 1]
	
	/* Maximum number of iterations reached */
	if (`iter' > `maxiter') {
		display as error "warning: REBUS algorithm did not converge"
		display as error "the solution provided may not be acceptable; " _continue
		display as error "try to relax any of the stopping criteria"
	}
end
