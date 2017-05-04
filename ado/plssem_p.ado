*!plssem_p version 0.1
*!Written 04May2017
*!Written by Sergio Venturini and Mehmet Mehmetoglu
*!The following code is distributed under GNU General Public License version 3 (GPL-3)

program plssem_p, rclass sortpreserve
	version 10

	syntax , [ xb RESiduals NOOUTer NOINner ]
	
	/* Options:
	   --------
		 xb								--> fitted values
		 residuals				--> residuals
		 noouter					--> no results returned for the outer model
		 noinner					--> no results returned for the inner model
	 */
	 
	 /* Description:
			------------
			This postestimation command computes fitted values and residuals for a
			PLSSEM model.
	 */
	
	local touse = e(sample)
	
	// check that estimation sample has not changed
	checkestimationsample

	if ("`e(formative)'" != "") {
		display as txt "(fitted values and residuals computed only " _continue
		display as txt "for reflective latent variables)"
	}
	if ("`noinner'" != "") & ("`noouter'" != "") {
		display as txt "noinner and noouter option chosen simultaneously; " _continue
		display as txt "nothing saved"
	}
	else {
		if ("`xb'" == "") & ("`residuals'" == "") {
			display as txt "(option " _continue
			display as result "xb" _continue
			display as txt " assumed; fitted values)"
			local xb "xb"
		}
	}

	local refl = e(reflective)
	local empty = .
	local noreflective : list refl == empty
	local struct "structural"
	local props = e(properties)
	local isstruct : list struct in props
	if (!`isstruct') {
		display as txt "(the model doesn't include any structural part; " _continue
		display as txt "quantities saved only for the measurement part)"
	}
	
	/* Fitted values */
	tempname latents
	local lvs `e(lvs)'
	mkmat `lvs' if `touse', matrix(`latents')

	/* 	- outer model */
	if (!`noreflective') {
		tempname adj_meas loadings b2use ohat modeA_scores indicators ///
			indicators_std tmp
		local modeA `e(reflective)'
		local lvs_adj : colnames e(adj_meas)
		matrix `adj_meas' = e(adj_meas)
		matrix `loadings' = e(loadings)
		local mvs : rownames `adj_meas'
		local n_lvs : word count `modeA'
		local n_mvs : word count `mvs'
		mkmat `mvs' if `touse', matrix(`indicators')
		mkmat `modeA' if `touse', matrix(`modeA_scores')
		local idx = 1
		local start = 1
		/*
		quietly misstable patterns `mvs' if `touse'
		local nobs = r(N_complete)
		*/
		local nobs = _N
		mata: st_matrix("`ohat'", J(`nobs', 0, .))
		mata: mvs_idx = J(0, 1, .)
		foreach var in `modeA' {
			local lv_col : list posof "`var'" in lvs_adj
			mata: load_var = st_matrix("`loadings'")[., `lv_col']
			mata: load_idx = selectindex(load_var :!= .)
			mata: mvs_idx = (mvs_idx \ load_idx)
			mata: load_nm = st_matrixrowstripe("`loadings'")[load_idx, .]
			mata: nblock = rows(load_idx)
			mata: st_local("nblock", strofreal(nblock))
			mata: st_matrix("`b2use'", transposeonly(load_var[load_idx, .]))
			mata: st_matrixcolstripe("`b2use'", load_nm)
			mata: st_matrix("`tmp'", ///
						st_matrix("`modeA_scores'")[., `idx'] * st_matrix("`b2use'"))
			mata: st_matrix("`ohat'", (st_matrix("`ohat'"), st_matrix("`tmp'")))
			local ++idx
			local end = `start' + `nblock' - 1
			matrix `ohat'[1, `start'] = `tmp'
			local start = `end' + 1
			local tmp_nm = ""
			forvalues j = 1/`nblock' {
				local tmp_nm "`tmp_nm' `: word `j' of `: colnames `b2use'''_hat"
			}
			local tmp_nm : list clean tmp_nm
			local ohat_nm "`ohat_nm' `tmp_nm'"
		}
		local ohat_nm : list clean ohat_nm
		matrix colnames `ohat' = `ohat_nm'
	}
	
	/*	- inner model */
	if (`isstruct') {
		tempname adj_struct path ihat
		matrix `adj_struct' = e(adj_struct)
		matrix `path' = e(pathcoef)
		mata: endo = colsum(st_matrix("`adj_struct'"))
		mata: endo = (endo :> 0)
		mata: st_local("n_endo", strofreal(sum(endo)))
		mata: path = st_matrix("`path'")[., selectindex(endo)]
		mata: st_matrix("`ihat'", st_matrix("`latents'") * path)
		mata: lvs_endo_nm = st_matrixcolstripe("`path'")[selectindex(endo), 2]
		forvalues j = 1/`n_endo' {
			mata: st_local("lvs_endo", lvs_endo_nm[`j'])
			local endovs "`endovs' `lvs_endo'"
			local ihat_nm "`ihat_nm' `lvs_endo'_hat"
		}
		local ihat_nm : list clean ihat_nm
		matrix colnames `ihat' = `ihat_nm'
	}

	/* Residuals */
	/*	- outer model */
	if (!`noreflective') {
		tempname ores
		mata: st_matrix("`indicators_std'", ///
			scale(st_matrix("`indicators'")[., mvs_idx]))
		matrix `ores' = `indicators_std' - `ohat'
		mata: mvs_nm = st_matrixcolstripe("`indicators'")[mvs_idx, 2]
		mata: st_local("n_ind", strofreal(rows(mvs_idx)))
		forvalues j = 1/`n_ind' {
			mata: st_local("ind_nm", mvs_nm[`j'])
			local indvs "`indvs' `ind_nm'"
			local ores_nm "`ores_nm' `ind_nm'_res"
		}
		matrix colnames `ores' = `ores_nm'
	}

	/*	- inner model */
	if (`isstruct') {
		tempname ires
		mata: st_matrix("`latents'", st_matrix("`latents'")[., selectindex(endo)])
		matrix `ires' = `latents' - `ihat'
		forvalues j = 1/`n_endo' {
			mata: st_local("lvs_endo", lvs_endo_nm[`j'])
			local ires_nm "`ires_nm' `lvs_endo'_res"
		}
		matrix colnames `ires' = `ires_nm'
	}

	/* All fitted values */
	tempname yhat
	if (!`noreflective') & (`isstruct') {
		matrix `yhat' = (`ohat', `ihat')
		matrix colnames `yhat' = `ohat_nm' `ihat_nm'
	}
	else if (!`noreflective') & (!`isstruct') {
		matrix `yhat' = `ohat'
		matrix colnames `yhat' = `ohat_nm'
	}
	else if (`noreflective') & (`isstruct') {
		matrix `yhat' = `ihat'
		matrix colnames `yhat' = `ihat_nm'
	}
	else {
		// nothing to save
	}
	
	/* All residuals */
	tempname res
	if (!`noreflective') & (`isstruct') {
		matrix `res' = (`ores', `ires')
		matrix colnames `res' = `ores_nm' `ires_nm'
	}
	else if (!`noreflective') & (!`isstruct') {
		matrix `res' = `ores'
		matrix colnames `res' = `ores_nm'
	}
	else if (`noreflective') & (`isstruct') {
		matrix `res' = `ires'
		matrix colnames `res' = `ires_nm'
	}
	else {
		// nothing to save
	}
	
	/* Copy quantities in data set */
	local now "`c(current_date)', `c(current_time)'"
	local now : list clean now
	if ("`xb'" != "") {
		if ("`noouter'" == "") & (!`noreflective') {
			forvalues j = 1/`n_ind' {
				local varname "`: word `j' of `ohat_nm''"
				capture confirm new variable `varname'
				if (_rc == 110) {
					quietly drop `varname'
				}
			}
			svmat double `ohat', names(col)
			forvalues j = 1/`n_ind' {
				local varname "`: word `j' of `ohat_nm''"
				label variable `varname' "Fitted values for outer model (`: word `j' of `indvs'') [`now']"
			}
		}
		if ("`noinner'" == "") & (`isstruct') {
			forvalues j = 1/`n_endo' {
				local varname "`: word `j' of `ihat_nm''"
				capture confirm new variable `varname'
				if (_rc == 110) {
					quietly drop `varname'
				}
			}
			svmat double `ihat', names(col)
			forvalues j = 1/`n_endo' {
				local endoname "`: word `j' of `ihat_nm''"
				label variable `endoname' "Fitted values for inner model (`: word `j' of `endovs'') [`now']"
			}
		}
	}
	if ("`residuals'" != "") {
		if ("`noouter'" == "") & (!`noreflective') {
			forvalues j = 1/`n_ind' {
				local varname "`: word `j' of `ores_nm''"
				capture confirm new variable `varname'
				if (_rc == 110) {
					quietly drop `varname'
				}
			}
			svmat double `ores', names(col)
			forvalues j = 1/`n_ind' {
				local varname "`: word `j' of `ores_nm''"
				label variable `varname' "Residuals for outer model (`: word `j' of `indvs'') [`now']"
			}
		}
		if ("`noinner'" == "") & (`isstruct') {
			forvalues j = 1/`n_endo' {
				local varname "`: word `j' of `ires_nm''"
				capture confirm new variable `varname'
				if (_rc == 110) {
					quietly drop `varname'
				}
			}
			svmat double `ires', names(col)
			forvalues j = 1/`n_endo' {
				local varname "`: word `j' of `ires_nm''"
				label variable `varname' "Residuals for inner model (`: word `j' of `endovs'') [`now']"
			}
		}
	}
	
	/* Clean up */
	if (!`noreflective') {
		mata: mata drop mvs_idx load_var load_idx load_nm nblock mvs_nm
	}
	if (`isstruct') {
		mata: mata drop endo path lvs_endo_nm
	}

	/* Return values */
	if (!`noreflective') | (`isstruct') {
		return matrix residuals = `res'
		return matrix fitted = `yhat'
	}
end

version 10
mata:
real matrix scale(real matrix X)
{
	/* Description:
		 ------------
		 Function that standardizes the columns of a real matrix X.
	*/
	
	/* Arguments:
		 ----------
		 - X --> real matrix whose columns must be standardized.
	*/
	
	/* Returned value:
		 ---------------
		 - Xstd --> real matrix with columns that have been standardized.
	*/
	
	real matrix m, Xc, sd_inv, Xstd
	real scalar n
	
	n = rows(X)
	m = J(n, 1, mean(X))
	Xc = X - m
	sd_inv = luinv(sqrt(diag(variance(Xc))))
	Xstd = Xc * sd_inv

	return(Xstd)
}
end
