*!plssemc_estat version 0.3.0
*!Written 27Apr2018
*!Written by Sergio Venturini and Mehmet Mehmetoglu
*!The following code is distributed under GNU General Public License version 3 (GPL-3)

program plssemc_estat, rclass
	version 14.2
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
	else {
		// estat_default `0'
		display as error "the `subcmd' postestimation command is not implemented for plssemc"
		exit
	}

	return add
end

program indirect_c, rclass
	version 14.2
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
	version 14.2
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
	version 14.2
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
