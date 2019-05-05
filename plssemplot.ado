*!plssemplot version 0.3.0
*!Written 27Apr2018
*!Written by Sergio Venturini and Mehmet Mehmetoglu
*!The following code is distributed under GNU General Public License version 3 (GPL-3)

program plssemplot
	version 15.1

	if (("`e(cmd)'" != "plssem") & ("`e(cmd)'" != "plssemc")) {
		error 301
	}
	
	tempname results eresults
	_return hold `results'
	_return restore `results' , hold
	_estimates hold `eresults', copy restore
	capture nobreak noisily {
		_plssemplot `0'
	}
	local rc = _rc
	_return restore `results'
	_estimates unhold `eresults'
	if (`rc') {
		exit `rc'
	}
end

program _plssemplot
	syntax [, INnermodel OUTermodel STats(varlist min=1 max=1) SCores ///
		Crossloadings Loadings OUTERWeights ]
	local 0 `", `options'"'

	tempvar __touse__
	quietly generate `__touse__' = e(sample)

	// check that estimation sample has not changed
	// checkestimationsample
	
	local props = e(properties)
	local struct "structural"
	local isstruct : list struct in props
	local rawsum "rawsum"
	local israwsum : list rawsum in props
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

	if ("`innermodel'" != "") {
		if (!`isstruct') {
			display as error "the fitted plssem model includes only the measurement part"
			exit
		}

		tempname adjmat
		mata: `adjmat' = st_matrix("e(adj_struct)")
		
		preserve
		
		capture quietly nwset, clear
		if ((_rc != 0) & (_rc != 6001)) {
			display as error "to use the 'innermodel' option you need to install the nwcommands suite"
			display as error "type " _continue
			display as smcl "{stata net install nwcommands-ado.pkg}"
			exit
		}
		
		quietly nwset, mat(`adjmat') directed labs(`e(lvs)')
		quietly nwplot, lab layout(circle) title("Structural model" " ") ///
			labelopt(mlabposition(0)) nodefactor(4.5) arcstyle(automatic) ///
			scatteropt(msymbol(O)) arrowgap(1) scheme(sj)
		
		restore
		mata: mata drop `adjmat'
	}
	
	if ("`outermodel'" != "") {
		display as error "this feature will be available soon! :)"
		/*
		tempname adjmat adjmat_full
		local nlvs : word count `e(lvs)'
		local nmvs : word count `e(mvs)'
		local all_nm "`: colnames e(adj_meas)' `: rownames e(adj_meas)'"
		matrix `adjmat' = J(`nlvs' + `nmvs', `nlvs' + `nmvs', 0)
		matrix `adjmat'[`nlvs' + 1, 1] = e(adj_meas)
		matrix `adjmat'[1, 1] = e(adj_struct)
		matrix rownames `adjmat' = `all_nm'
		matrix colnames `adjmat' = `all_nm'
		mata: `adjmat_full' = st_matrix("`adjmat'")
		
		preserve
		
		capture quietly nwset, clear
		if ((_rc != 0) & (_rc != 6001)) {
			display as error "to use the 'outermodel' option you need to install the nwcommands suite"
			display as error "type: net install nwcommands-ado.pkg"
			exit
		}

		quietly nwset, mat(`adjmat_full') directed labs(`all_nm')
		quietly nwplot, lab layout(mds) title("Measurement model" " ") ///
			labelopt(mlabposition(0)) nodefactor(2.5) arcstyle(automatic) ///
			scatteropt(msymbol(O)) arrowgap(1) scheme(sj)
		
		restore
		mata: mata drop `adjmat_full'
		*/
	}
	
	if ("`stats'" != "") {
		tempname loads lv
		matrix `loads' = e(loadings)
		local num_ind = rowsof(`loads')
		local num_lv = colsof(`loads')
		local loads_cn : colnames `loads'
		local loads_rn : rownames `loads'
		local pos_lv : list posof "`stats'" in loads_cn
		matrix `lv' = `loads'[1..`num_ind', `pos_lv']
		forvalues i = 1/`num_ind' {
			if (!missing(`lv'[`i', 1])) {
				local ind_lv : word `i' of `loads_rn'
				local block_lv "`block_lv' `ind_lv'"
			}
		}
		
		if ("`missing'" != "") {
			/* Save original data set */
			local allindicators = e(mvs)
			tempname original_data
			mata: `original_data' = st_data(., "`: list uniq allindicators'", "`__touse__'")
			
			/* Recovery of missing values */
			mata: st_store(., tokens("`: list uniq allindicators'"), "`__touse__'", ///
				st_matrix("e(imputed_data)"))
		}
		
		quietly summarize `block_lv' if `__touse__'
		local jsize = (r(max) - r(min))*.15
		if ("`missing'" == "") {
			local stats_title "`stats'"
		}
		else {
			local stats_title "`stats' (imputed)"
		}
		graph matrix `block_lv' if `__touse__', half title(`stats_title') ///
			maxes(ylab(#6, grid) xlab(#6, grid)) msymbol(oh) jitter(`jsize') ///
			jitterseed(1406) scheme(sj)
		
		/*
		local num_block : word count `block_lv'
		forvalues i = 1/`num_block' {
			forvalues j = 1/`num_block' {
				if (`i' == `j') {
					tempname hist`j'
					local hist_lv : word `j' of `block_lv'
					quietly tabulate `hist_lv' if `__touse__'
					if (r(r) > 10) {
						local type ""
						local freq ""
					}
					else {
						local type "discrete"
						local freq "frequency"
					}
					if (`i' == 1) {
						local xlabel "xlabel(none)"
						local ylabel ""
						local xscale ""
						local yscale ""
					}
					else if (`i' == `num_block') {
						local xlabel ""
						local ylabel ""
						local xscale ""
						local yscale "yscale(alt)"
					}
					else {
						local xlabel "xlabel(none)"
						local ylabel "ylabel(none)"
						local xscale ""
						local yscale ""
					}
					quietly histogram `hist_lv' if `__touse__', `type' name(`hist`j'') ///
						ytitle("") xtitle("") `xlabel' `ylabel' `xscale' `yscale' `freq' ///
						plotregion(lcolor(black) lwidth(thin)) ///
						title(`hist_lv', ring(0)) nodraw scheme(sj)
					local toplot "`toplot' `hist`j''"
				}
				else if (`i' < `j') {
					tempname scatter`i'`j'
					tempvar empty
					quietly generate `empty' = . if `__touse__'
					local sc_x : word `j' of `block_lv'
					local sc_y : word `i' of `block_lv'
					quietly correlate `sc_y' `sc_x' if `__touse__'
					local corr_todisp = string(r(rho), "%9.4g")
					local corr_todisp "r = `corr_todisp'"
					local used = r(N)
					quietly count if `__touse__'
					local size = r(N)
					local N_todisp = string(`used'/`size'*100, "%9.2g")
					local N_todisp "N = `used' (`N_todisp'%)"
					twoway scatter `empty' `empty' if `__touse__', legend(off) ///
						ytitle("") xtitle("") xlabel(none) ylabel(none) ///
						xscale(range(0 2)) yscale(range(0 3)) ///
						text(1.75 1 "`corr_todisp'", size(medium)) ///
						text(1.25 1 "`N_todisp'", size(medium)) ///
						plotregion(lcolor(black) lwidth(thin)) ///
						name(`scatter`i'`j'') nodraw scheme(sj)
					local toplot "`toplot' `scatter`i'`j''"
				}
				else {
					tempname scatter`i'`j'
					local sc_x : word `j' of `block_lv'
					local sc_y : word `i' of `block_lv'
					if ((`i' > 1) & (`i' < `num_block')) {
						local xlabel "xlabel(none)"
						local ylabel ""
						local xscale ""
						local yscale ""
					}
					else if ((`j' == 1) & (`i' == `num_block')) {
						local xlabel ""
						local ylabel ""
						local xscale ""
						local yscale ""
					}
					else if ((`i' > 1) & (`i' == `num_block')) {
						local xlabel ""
						local ylabel "ylabel(none)"
						local xscale ""
						local yscale ""
					}
					else {
						local xlabel "xlabel(none)"
						local ylabel "ylabel(none)"
						local xscale ""
						local yscale ""
					}
					twoway (scatter `sc_y' `sc_x' if `__touse__', jitter(`jsize') ///
						ylab(#6, grid) xlab(#6, grid) msymbol(oh) jitterseed(1406)) ///
						(lfit `sc_y' `sc_x', lpattern(solid)), legend(off) ///
						ytitle("") xtitle("") `xlabel' `ylabel' `xscale' `yscale' ///
						plotregion(lcolor(black) lwidth(thin)) ///
						name(`scatter`i'`j'') nodraw scheme(sj)
					local toplot "`toplot' `scatter`i'`j''"
				}
			}
		}
		graph combine `toplot', imargin(vsmall) ///
			graphregion(margin(l=22 r=22)) title(`stats') scheme(sj)
		*/
		
		/* Clean up */
		if ("`missing'" != "") {
			mata: st_store(., tokens("`: list uniq allindicators'"), "`__touse__'", ///
				`original_data')
		}
	}
	
	if ("`scores'" != "") {
		if (!`isstruct') {
			display as error "the fitted plssem model includes only the measurement part"
			exit
		}
		
		graph matrix `e(lvs)' if `__touse__', half title("PLS-SEM scores") ///
			maxes(ylab(#6, grid) xlab(#6, grid)) msymbol(Oh) scheme(sj)
	}

	if ("`crossloadings'" != "") {
		tempname loads
		matrix `loads' = e(loadings)
		local num_ind = rowsof(`loads')
		local num_lv = colsof(`loads')
		tempname xload C
		matrix `xload' = J(`num_ind', `num_lv', .)
		matrix rownames `xload' = `e(mvs)'
		matrix colnames `xload' = `e(lvs)'
		local j = 1
		local last = `num_ind' + 1
		foreach var in `e(lvs)' {
			quietly correlate `var' `e(mvs)' if `__touse__'
			matrix `C' = r(C)
			matrix `xload'[1, `j'] = `C'[2..`last', 1]
			local ++j
		}
		mktable, matrix(`xload') digits(4) firstcolname("") ///
			title("Cross loadings") firstcolwidth(12) colwidth(9) ///
			hlines(`num_ind') novlines corr
		
		preserve
		quietly drop _all
		quietly svmat double `xload', names("lv")
		quietly generate indicator = .
		quietly generate block = ""
		quietly generate order = .
		local i 1
		local j 1
		foreach ind in `e(mvs)' {
			quietly replace indicator = `i' in `i'
			quietly replace order = `i' if indicator == `i'
			foreach var in `e(lvs)' {
				if (!missing(`loads'[`i', `j'])) {
					local lv_nm : word `j' of `e(lvs)'
					quietly replace block = "`var'" in `i'
				}
				local ++j
			}
			local lblbar "`lblbar' `i' `ind'"
			local j 1
			local ++i
		}
		quietly generate id = _n
		quietly reshape long lv, i(id) j(latent) string
		local i 1
		foreach var in `e(lvs)' {
			quietly replace latent = "`var'" if latent == "`i'"
			quietly clonevar lv`i' = lv
			quietly replace lv`i' = . if block != "`var'"
			quietly label variable lv`i' "`var'"
			local ++i
		}
		tempname minlv maxlv
		mata: st_numscalar("`minlv'", min(st_data(., "lv?")))
		mata: st_numscalar("`maxlv'", max(st_data(., "lv?")))
		if (`minlv' < -.5) {
			local yline_neg = "dash"
		}
		else {
			local yline_neg = "blank"
		}
		graph bar (asis) lv?, over(indicator, label(angle(90) labsize(small)) ///
			relabel(`lblbar') gap(0)) bargap(-90) by(block latent, imargin(tiny) ///
			note(""/*"Graphs by block and latent variable"*/) ///
			title("Cross loadings") legend(position(12)) iscale(*.885)) ///
			ylabel(``minlv'(.25)`maxlv'', labsize(small) nogrid) ///
			yline(0.5, lpattern(dash) lwidth(thin)) ///
			yline(-0.5, lpattern(`yline_neg') lwidth(thin)) ///
			yline(0, lpattern(solid) lwidth(thin)) legend(rows(1) ///
			region(style(none)) size(3) keygap(1) symxsize(5)) scheme(sj)
		restore
	}

	if ("`loadings'" != "") {
		tempname loads
		matrix `loads' = e(loadings)
		local num_ind = rowsof(`loads')
		local num_lv = colsof(`loads')
		mktable, matrix(`loads') digits(4) firstcolname("") ///
			title("Loadings") firstcolwidth(14) colwidth(14) ///
			hlines(`num_ind') novlines
		
		tempvar lv indicator
		quietly svmat double `loads', names(`lv')
		quietly generate `indicator' = .
		local i 1
		foreach ind in `e(mvs)' {
			quietly replace `indicator' = `i' in `i'
			local lblbar "`lblbar' `i' `ind'"
			local ++i
		}
		local i 1
		foreach var in `e(lvs)' {
			quietly label variable `lv'`i' "`var'"
			local ++i
		}
		tempname minlv maxlv
		mata: st_numscalar("`minlv'", min(st_data(., "`lv'?")))
		mata: st_numscalar("`maxlv'", max(st_data(., "`lv'?")))
		if (`minlv' < -.7) {
			local yline_neg = "dash"
		}
		else {
			local yline_neg = "blank"
		}
		graph bar (asis) `lv'?, over(`indicator', label(angle(90) labsize(small)) ///
			relabel(`lblbar') gap(*.3)) title("Loadings") nofill ///
			ylabel(``minlv'(.1)`maxlv'', labsize(small) nogrid) ///
			yline(0.7, lpattern(dash) lwidth(thin)) ///
			yline(-0.7, lpattern(`yline_neg') lwidth(thin)) ///
			yline(0, lpattern(solid) lwidth(thin)) legend(position(12) rows(1) ///
			region(style(none)) size(3) keygap(1) symxsize(5)) scheme(sj)
	}
	
	if ("`outerweights'" != "") {
		if (!`isstruct') {
			display as error "the fitted plssem model includes only the measurement part"
			exit
		}
		if (`israwsum') {
			display as error "the plssem model has been fitted using the rawsum option"
			exit
		}
	
		tempname ow
		matrix `ow' = e(ow_history)
		local niter = e(iterations)
		
		preserve
		
		tempvar iter lv
		quietly drop _all
		quietly svmat double `ow', names(eqcol)
		quietly generate `iter' = _n - 1
		capture quietly reshape long `e(mvs)', i(`iter') j(`lv') string
		if (_rc != 0) {
			display as error "some of the indicators have been used more than once; " _continue
			display as error "the outer weights diagram can't be displayed"
			restore
			exit
		}
		// here we use colors (s2gcolor) to emphasize differences
		twoway (line `e(mvs)' `iter', lpattern(solid)), ///
			by(`lv', legend(position(3)) title("Outer weights convergence") ///
			note("")) xtitle("Iteration") ytitle("Outer weights") xlabel(#`niter') ///
			legend(cols(1) size(vsmall)) scheme(s2gcolor)

		restore
	}	
end
