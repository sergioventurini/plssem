*!plssem_estat version 0.1
*!Written 05Jan2017
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
	else {
		// estat_default `0'
		display as error "the `subcmd' postestimation command is not implemented for plssem"
		error 197
	}

	return add
end

program indirect, eclass
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
		error 198
	}
	local struct "structural"
	local props = e(properties)
	local isstruct : list struct in props
	if (!`isstruct') {
		display as error "the fitted plssem model includes only the measurement part"
		error 198
	}
	
	local touse = e(sample)
	local reg3eqs = e(struct_eqs)
	if (`level' <= 0 | `level' >= 1) {
		display as error "confidence level must be in the range (0, 1)"
		error 197
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
					quietly reg3 `reg3eqs' if `touse', mvreg corr(independent)
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
						seed(`seed'): reg3 `reg3eqs' if `touse', mvreg corr(independent)
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
				error 198
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
	// ereturn matrix indirect = `indmeas'
end

program total, eclass
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
		error 198
	}
	if (`digits' < 0) {
		display as error "number of digits must be a nonnegative integer"
		error 197
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
	// ereturn matrix total = `alleffects_new'
end

program mktable_indirect
	version 10
	syntax , matrix(string) [ FIRSTCOLName(string) FIRSTCOLWidth(integer 25) ///
		COLWidth(integer 15) Title(string) HLines(numlist >0 integer sort) ///
		NOVLines DIGits(integer 3) Boot(numlist min=1 max=1) Level(real 0.95) ]

	/* Options:
	   --------
		 matrix(string)							--> matrix containing the numbers to display
		 firstcolname(string)				--> first column's name
		 firstcolwidth(integer 25)	--> first column's width (default 25)
		 colwidth(integer 15)				--> other columns' width (default 15)
		 title(string)							--> table's main title
		 hlines											--> rows at which to insert an horizontal line
		 novlines										--> no vertical lines between columns
		 digits(integer 3)					--> number of digits to display (default 3)
		 boot												--> flag for bootstrap
		 level(real 0.95)						--> confidence level (default 0.95)
	 */
	
	local skip0 = 0
	local skip1 = 1
	local skip2 = 2
	local skip3 = 3
	local skip4 = 4

	if (strlen("`firstcolname'") > `firstcolwidth' - 2*`skip1') {
		local firstcolname = abbrev("`firstcolname'", `firstcolwidth' - 2*`skip1')
	}
	if (`colwidth' < 10) {
		display as error "colwidth must be larger than 9 to properly show the table"
		error 198
	}

	local firstcolwidth_p1 = `firstcolwidth' + 1
	local ncols = colsof(`matrix')
	local nrows = rowsof(`matrix')
	local tabrows = 5
	local matcolnames : colnames `matrix'
	if ("`boot'" == "") {
		local matrownames "Indirect effect" "_" "Standard error" "_" ///
			"Z statistic" "_" "P-value" "_" "Conf. interval"
	}
	else {
		local matrownames "Indirect effect" "_" "Standard error" "_" ///
			"Z statistic" "_" "P-value" "_" "Conf. interval (N)" "_" ///
			"Conf. interval (P)" "_" "Conf. interval (BC)"
	}
	if ("`boot'" != "") {
		local tabrows = `tabrows' + 2
	}
	local ncols_m1 = `ncols' - 1
	local usable = `colwidth' - 2*`skip1'
	local usable2 = floor(`usable'/2 - 1)
	if (`digits' >= `usable') {
		display as error "the number of digits chosen is too large"
		error 198
	}
	if ("`hlines'" != "") {
		numlist "`hlines'"
		local hlinestodisp = r(numlist)
	}
	if ("`novlines'" == "") {
		local tvlines: display "{c TT}"
		local bvlines: display "{c BT}"
		local mvlines: display "{c +}"
		local vlines: display _skip(`skip1') "{c |}"
	}
	else {
		local tvlines: display "{c -}"
		local bvlines: display "{c -}"
		local mvlines: display "{c -}"
		local vlines: display _skip(`skip2')
	}
	local firstline: display "{hline ""`firstcolwidth'""}{c TT}"
	local title1: display _col(`firstcolwidth_p1') "{c |}"
	local todisp "`firstcolname'"
	if (strlen("`todisp'") > `firstcolwidth' - 2*`skip1') {
		local todisp = abbrev("`todisp'", `firstcolwidth' - 2*`skip1')
	}
	local tmp_skip = `firstcolwidth' - strlen("`todisp'") - 2*`skip1'
	local title2: display _skip(`skip1') "`firstcolname'" _skip(`tmp_skip') _col(`firstcolwidth_p1') "{c |}"
	local title3: display _col(`firstcolwidth_p1') "{c |}"
	local secondline: display "{hline ""`firstcolwidth'""}{c +}"
	local lastline: display "{hline ""`firstcolwidth'""}{c BT}"

	forvalues j = 1/`ncols_m1' {
		local firstline "`firstline'" "`: display "{hline ""`colwidth'""}`tvlines'"'"
		local colname_j : word `j' of `matcolnames'
		tokenize `"`colname_j'"', parse("_")
		local dep `1'
		macro shift
		local mod `2'
		macro shift
		local indep `3'
		local todisp "`dep'"
		if (strlen("`todisp'") > `usable' - 3) {
			local todisp = abbrev("`todisp'", `usable' - 3)
		}
		local tmp_skip = `usable' - strlen("`todisp'") - 3
		local title1 "`title1'" "`: display _skip(`skip1') _skip(`tmp_skip')'" "`todisp' <-" "`vlines'"
		local todisp "`mod'"
		if (strlen("`todisp'") > `usable' - 3) {
			local todisp = abbrev("`todisp'", `usable' - 3)
		}
		local tmp_skip = `usable' - strlen("`todisp'") - 3
		local title2 "`title2'" "`: display _skip(`skip1') _skip(`tmp_skip')'" "`todisp' <-" "`vlines'"
		local todisp "`indep'"
		if (strlen("`todisp'") > `usable') {
			local todisp = abbrev("`todisp'", `usable')
		}
		local tmp_skip = `usable' - strlen("`todisp'")
		local title3 "`title3'" "`: display _skip(`skip1') _skip(`tmp_skip')'" "`todisp'" "`vlines'"
		local secondline "`secondline'" "`: display "{hline ""`colwidth'""}`mvlines'"'"
		local lastline "`lastline'" "`: display "{hline ""`colwidth'""}`bvlines'"'"
	}
	local firstline "`firstline'" "`: display "{hline ""`colwidth'""}"'"
	local colname_j : word `ncols' of `matcolnames'
	tokenize `"`colname_j'"', parse("_")
	local dep `1'
	macro shift
	local mod `2'
	macro shift
	local indep `3'
	local todisp "`dep'"
	if (strlen("`todisp'") > `usable' - 3) {
		local todisp = abbrev("`todisp'", `usable' - 3)
	}
	local tmp_skip = `usable' - strlen("`todisp'") - 3
	local title1 "`title1'" "`: display _skip(`skip1') _skip(`tmp_skip')'" "`todisp' <-"
	local todisp "`mod'"
	if (strlen("`todisp'") > `usable' - 3) {
		local todisp = abbrev("`todisp'", `usable' - 3)
	}
	local tmp_skip = `usable' - strlen("`todisp'") - 3
	local title2 "`title2'" "`: display _skip(`skip1') _skip(`tmp_skip')'" "`todisp' <-"
	local todisp "`indep'"
	if (strlen("`todisp'") > `usable') {
		local todisp = abbrev("`todisp'", `usable')
	}
	local tmp_skip = `usable' - strlen("`todisp'")
	local title3 "`title3'" "`: display _skip(`skip1') _skip(`tmp_skip')'" "`todisp'"
	local secondline "`secondline'" "`: display "{hline ""`colwidth'""}"'"
	local lastline "`lastline'" "`: display "{hline ""`colwidth'""}"'"
	
	display
	if ("`boot'" == "") {
		display as txt _skip(0) "`title'"
	}
	else {
		display as txt _skip(0) "`title'" " (Bootstrap)"
	}
	display as txt "`firstline'"
	display as txt "`title1'"
	display as txt "`title2'"
	display as txt "`title3'"
	display as txt "`secondline'"
	
	tokenize `"`matrownames'"', parse("_")
	local disp_ci = 1
	forvalues i = 1/`tabrows' {
		local rownametodisp "``i''"
		macro shift
		macro shift
		if (strlen("`rownametodisp'") >= (`firstcolwidth' - 2*`skip1')) {
			local rownametodisp = abbrev("`rownametodisp'", `firstcolwidth' - 2*`skip1')
		}
		local tmp_skip = `firstcolwidth' - strlen("`rownametodisp'") - 2*`skip1'
		local rownametodisp "`: display _skip(`skip1')'" "`rownametodisp'" "`: display _skip(`tmp_skip') _skip(`skip1') "{c |}"'"
		display as txt "`rownametodisp'" _continue

		local rowtodisp ""
		if (`i' <= 4) {
			forvalues j = 1/`ncols_m1' {
				if (!missing(`matrix'[`i', `j'])) {
					local todisp = strtrim(string(`matrix'[`i', `j'], "%`usable'.`digits'f"))
					if (strlen("`todisp'") > `usable') {
						local todisp = strtrim(string(`matrix'[`i', `j'], "%`usable'.3e"))
					}
					local tmp_skip = `usable' - strlen("`todisp'")
					local rowtodisp "`rowtodisp'" "`: display _skip(`skip1') _skip(`tmp_skip') "`todisp'`vlines'"'"
				}
				else {
					local rowtodisp "`rowtodisp'" "`: display _skip(`skip1') _skip(`usable') "`vlines'"'"
				}
			}
			if (!missing(`matrix'[`i', `ncols'])) {
				local todisp = strtrim(string(`matrix'[`i', `ncols'], "%`usable'.`digits'f"))
				if (strlen("`todisp'") > `usable') {
					local todisp = strtrim(string(`matrix'[`i', `ncols'], "%`usable'.3e"))
				}
				local tmp_skip = `usable' - strlen("`todisp'")
				local rowtodisp "`rowtodisp'" "`: display _skip(`skip1') _skip(`tmp_skip') "`todisp'"'"
			}
		}
		else {
			forvalues j = 1/`ncols_m1' {
				if (!missing(`matrix'[`i' + `disp_ci' - 1, `j']) & ///
						!missing(`matrix'[`i' + `disp_ci', `j'])) {
					local todisp1 = strtrim(string(`matrix'[`i' + `disp_ci' - 1, `j'], "%`usable'.`digits'f"))
					local todisp2 = strtrim(string(`matrix'[`i' + `disp_ci', `j'], "%`usable'.`digits'f"))
					if (strlen("`todisp1'") >= `usable2') {
						local todisp1 = strtrim(string(`todisp1', "%`usable'.2f"))
					}
					if (strlen("`todisp2'") >= `usable2') {
						local todisp2 = strtrim(string(`todisp2', "%`usable'.2f"))
					}
					local todisp = "(`todisp1', `todisp2')"
					local tmp_skip = `usable' - strlen("`todisp'")
					local rowtodisp "`rowtodisp'" "`: display _skip(`skip1') _skip(`tmp_skip') "`todisp'`vlines'"'"
				}
				else {
					local rowtodisp "`rowtodisp'" "`: display _skip(`skip1') _skip(`usable') "`vlines'"'"
				}
			}
			if (!missing(`matrix'[`i' + `disp_ci' - 1, `ncols']) & ///
					!missing(`matrix'[`i' + `disp_ci', `ncols'])) {
				local todisp1 = strtrim(string(`matrix'[`i' + `disp_ci' - 1, `ncols'], "%`usable'.`digits'f"))
				local todisp2 = strtrim(string(`matrix'[`i' + `disp_ci', `ncols'], "%`usable'.`digits'f"))
				if (strlen("`todisp1'") >= `usable2') {
					local todisp1 = strtrim(string(`todisp1', "%`usable'.2f"))
				}
				if (strlen("`todisp2'") >= `usable2') {
					local todisp2 = strtrim(string(`todisp2', "%`usable'.2f"))
				}
				local todisp = "(`todisp1', `todisp2')"
				local tmp_skip = `usable' - strlen("`todisp'")
				local rowtodisp "`rowtodisp'" "`: display _skip(`skip1') _skip(`tmp_skip') "`todisp'"'"
			}
			local ++disp_ci
		}
		display as result "`rowtodisp'"
		
		foreach num in `hlinestodisp' {
			if (`num' == `i') {
				if (`num' < `tabrows') {
					display as txt "`secondline'"
				}
				else if (`num' == `tabrows') {
					display as txt "`lastline'"
				}
			}
			continue, break
		}
	}
	
	if (`level' <= 0 | `level' >= 1) {
		display as error "confidence level must be in the range (0, 1)"
		error 197
	}
	local conflevel = strtrim(string(`level'*100, "%9.2g"))
	display as txt _skip(`skip1') "confidence level: `conflevel'%"
	if ("`boot'" != "") {
		display as txt _skip(`skip1') "(N)    normal confidence interval"
		display as txt _skip(`skip1') "(P)    percentile confidence interval"
		display as txt _skip(`skip1') "(BC)   bias-corrected confidence interval"
	}
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
		error 198
	}
	else if (`isboot') {
		display as error "estat vif is not allowed after bootstrap"
		error 198
	}
	else {
		tempname strvif
		matrix `strvif' = e(struct_vif)
		local hline_path = rowsof(`strvif')
		mktable, matrix(`strvif') digits(`digits') firstcolname("Variable") ///
			title("Structural model - Multicollinearity check (VIFs)") ///
			firstcolwidth(14) colwidth(12) hlines(`hline_path') novlines
	}
end
