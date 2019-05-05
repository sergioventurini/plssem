*!plssemc_estat version 0.3.0
*!Written 03Nov2018
*!Written by Sergio Venturini and Mehmet Mehmetoglu
*!The following code is distributed under GNU General Public License version 3 (GPL-3)

program plssemc_estat, rclass
	version 15.1
	gettoken subcmd rest : 0 , parse(", ")
	local lsubcmd = length("`subcmd'")
	
	if ("`subcmd'" == substr("indirect", 1, max(2, `lsubcmd'))) {
		indirect_c `rest'
	}
	else if ("`subcmd'" == substr("total", 1, max(2, `lsubcmd'))) {
		total_c `rest'
	}
	else if ("`subcmd'" == substr("vif", 1, max(2, `lsubcmd'))) {
		plssem_vif_c `rest'
	}
	else if ("`subcmd'" == substr("mediate", 1, max(2, `lsubcmd'))) {
		mediate_c `rest'
	}
	else {
		// estat_default `0'
		display as error "the `subcmd' postestimation command is not implemented for plssemc"
		exit
	}

	return add
end

program indirect_c, rclass
	version 15.1
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
	
	display as error "estat indirect after plssemc will be available soon! :)"
	exit
	
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

	if ("`boot'" != "") {
		display
		display as text "Computing indirect effects bootstrap distribution..."
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
	return matrix indirect `indmeas'
end

program total_c, rclass
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

program plssem_vif_c, rclass
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

program mediate_c, rclass
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
