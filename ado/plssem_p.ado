*!plssem_p version 0.3.0
*!Written 01Oct2017
*!Written by Sergio Venturini and Mehmet Mehmetoglu
*!The following code is distributed under GNU General Public License version 3 (GPL-3)

program plssem_p, rclass sortpreserve
	version 15.1

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
		 PLS-SEM model.
	*/
	
	tempvar __touse__
	quietly generate `__touse__' = e(sample)
	
	// check that estimation sample has not changed
	// checkestimationsample

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
	local props = e(properties)
	local struct "structural"
	local isstruct : list struct in props
	if (!`isstruct') {
		display as txt "(the model doesn't include any structural part; " _continue
		display as txt "quantities saved only for the measurement part)"
	}
	local noscale "unscaled"
	local isnoscale : list noscale in props
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
	
	/* Save original data set */
	local allindicators = e(mvs)
	tempname original_data original_data_exp
	mata: `original_data' = st_data(., "`: list uniq allindicators'", "`__touse__'")
	mata: `original_data_exp' = st_data(., "`allindicators'", "`__touse__'")
	
	/* Recovery of missing values */
	if ("`missing'" != "") {
		mata: st_store(., tokens("`: list uniq allindicators'"), "`__touse__'", ///
			st_matrix("e(imputed_data)"))
	}
	
	/* Save fitted values */
	tempname latents
	local lvs `e(lvs)'
	mkmat `lvs' if `__touse__', matrix(`latents')
	
	/* -- outer model -- */
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
		mkmat `mvs' if `__touse__', matrix(`indicators')
		mkmat `modeA' if `__touse__', matrix(`modeA_scores')
		local idx = 1
		local start = 1
		quietly count if `__touse__'
		local nobs = r(N)
		
		mata: st_matrix("`ohat'", J(`nobs', 0, .))
		tempname mvs_idx load_var load_idx load_nm num_block
		mata: `mvs_idx' = J(0, 1, .)
		foreach var in `modeA' {
			local lv_col : list posof "`var'" in lvs_adj
			mata: `load_var' = st_matrix("`loadings'")[., `lv_col']
			mata: `load_idx' = selectindex(`load_var' :!= .)
			mata: `mvs_idx' = (`mvs_idx' \ `load_idx')
			mata: `load_nm' = st_matrixrowstripe("`loadings'")[`load_idx', .]
			mata: `num_block' = rows(`load_idx')
			mata: st_local("nblock", strofreal(`num_block'))
			mata: st_matrix("`b2use'", transposeonly(`load_var'[`load_idx', .]))
			mata: st_matrixcolstripe("`b2use'", `load_nm')
			mata: st_matrix("`tmp'", ///
							st_matrix("`modeA_scores'")[., `idx'] * st_matrix("`b2use'"))
			mata: st_matrix("`ohat'", (st_matrix("`ohat'"), st_matrix("`tmp'")))
			local ++idx
			local end = `start' + `nblock' - 1
			matrix `ohat'[1, `start'] = `tmp'
			local start = `end' + 1
			local tmp_nm = ""
			forvalues j = 1/`nblock' {
				local tmp_nm "`tmp_nm' `var'`: word `j' of `: colnames `b2use'''_hat"
			}
			local tmp_nm : list clean tmp_nm
			local ohat_nm "`ohat_nm' `tmp_nm'"
		}
		local ohat_nm : list clean ohat_nm
		matrix colnames `ohat' = `ohat_nm'
	}

	/* -- inner model -- */
	if (`isstruct') {
		tempname adj_struct path endo path_mata ihat lvs_endo_nm
		matrix `adj_struct' = e(adj_struct)
		matrix `path' = e(pathcoef)
		mata: `endo' = colsum(st_matrix("`adj_struct'"))
		mata: `endo' = (`endo' :> 0)
		mata: st_local("n_endo", strofreal(sum(`endo')))
		mata: `path_mata' = st_matrix("`path'")[., selectindex(`endo')]
		mata: st_matrix("`ihat'", st_matrix("`latents'") * `path_mata')
		mata: `lvs_endo_nm' = st_matrixcolstripe("`path'")[selectindex(`endo'), 2]
		forvalues j = 1/`n_endo' {
			mata: st_local("lvs_endo", `lvs_endo_nm'[`j'])
			local endovs "`endovs' `lvs_endo'"
			local ihat_nm "`ihat_nm' `lvs_endo'_hat"
		}
		local ihat_nm : list clean ihat_nm
		matrix colnames `ihat' = `ihat_nm'
	}
	/* End of saving fitted values */

	/* Save residuals */
	/* -- outer model -- */
	if (!`noreflective') {
		tempname ores mvs_nm
		if (!`isnoscale') {
			mata: st_matrix("`indicators_std'", ///
				scale(st_matrix("`indicators'")[., `mvs_idx']))
			matrix `ores' = `indicators_std' - `ohat'
		}
		else {
			matrix `ores' = `indicators' - `ohat'
		}
		mata: `mvs_nm' = st_matrixcolstripe("`indicators'")[`mvs_idx', 2]
		mata: st_local("n_ind", strofreal(rows(`mvs_idx')))
		forvalues j = 1/`n_ind' {
			mata: st_local("ind_nm", `mvs_nm'[`j'])
			local indvs "`indvs' `ind_nm'"
		}
		local ores_nm : subinstr local ohat_nm "_hat" "_res", all
		matrix colnames `ores' = `ores_nm'
	}

	/* -- inner model -- */
	if (`isstruct') {
		tempname ires
		mata: st_matrix("`latents'", st_matrix("`latents'")[., selectindex(`endo')])
		matrix `ires' = `latents' - `ihat'
		forvalues j = 1/`n_endo' {
			mata: st_local("lvs_endo", `lvs_endo_nm'[`j'])
			local ires_nm "`ires_nm' `lvs_endo'_res"
		}
		matrix colnames `ires' = `ires_nm'
	}
	/* End of saving residuals */

	/*
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
	*/
	
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
				quietly generate `varname' = .
				label variable `varname' "Fitted values for outer model (`: word `j' of `indvs'') [`now']"
			}
			tempname ohat_mata orig_m orig_s
			mata: `ohat_mata' = st_matrix("`ohat'")
			mata: `orig_m' = mean(`original_data_exp')
			mata: `orig_s' = sd(`original_data_exp')
			mata: st_matrix("`ohat'", unscale(scale(`ohat_mata'), `orig_s', `orig_m'))
			mata: st_store(., tokens("`ohat_nm'"), "`__touse__'", st_matrix("`ohat'"))
		}
		if ("`noinner'" == "") & (`isstruct') {
			forvalues j = 1/`n_endo' {
				local varname "`: word `j' of `ihat_nm''"
				capture confirm new variable `varname'
				if (_rc == 110) {
					quietly drop `varname'
				}
				quietly generate `varname' = .
				label variable `varname' "Fitted values for inner model (`: word `j' of `endovs'') [`now']"
			}
			mata: st_store(., tokens("`ihat_nm'"), "`__touse__'", st_matrix("`ihat'"))
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
				quietly generate `varname' = .
				label variable `varname' "Residuals for outer model (`: word `j' of `indvs'') [`now']"
			}
			mata: st_store(., tokens("`ores_nm'"), "`__touse__'", st_matrix("`ores'"))
		}
		if ("`noinner'" == "") & (`isstruct') {
			forvalues j = 1/`n_endo' {
				local varname "`: word `j' of `ires_nm''"
				capture confirm new variable `varname'
				if (_rc == 110) {
					quietly drop `varname'
				}
				quietly generate `varname' = .
				label variable `varname' "Residuals for inner model (`: word `j' of `endovs'') [`now']"
			}
			mata: st_store(., tokens("`ires_nm'"), "`__touse__'", st_matrix("`ires'"))
		}
	}
	/* End of copying quantities in data set */
	
	/* Clean up */
	if ("`missing'" != "") {
		mata: st_store(., tokens("`: list uniq allindicators'"), "`__touse__'", ///
			`original_data')
	}
	capture mata: cleanup()
	
	/* Return values */
	if (!`noreflective') | (`isstruct') {
		return matrix struct_res = `ires'
		return matrix meas_res = `ores'
		return matrix struct_fit = `ihat'
		return matrix meas_fit = `ohat'
	}
end
