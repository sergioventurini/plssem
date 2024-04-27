*!plssemc_mi version 0.1.1
*!Written 03Aug2023
*!Written by Sergio Venturini and Mehmet Mehmetoglu
*!The following code is distributed under GNU General Public License version 3 (GPL-3)

program plssemc_mi
	version 15.1
	syntax [anything] [if] [in] [, * ]
	
	if replay() {
		if ("`e(cmd)'" != "mi estimate") {
			error 301
		}
		if (_by()) {
			error 190
		}
		
		// Extract all the options from 'cmdline'
		local cmdline = `"`e(cmdline)'"'
		local cmdlen = strlen(`"`cmdline'"')
		local cmdcomma = strpos(`"`cmdline'"', ",")
		local old_options = substr(`"`cmdline'"', `cmdcomma' + 1, `cmdlen' - `cmdcomma')
		local old_options : list clean old_options
		local hasstruct = strpos(`"`old_options'"', "str")
		if (!`hasstruct') {
			local nostructural "nostructural"
		}
		
		if (`"`options'"' != "") {
      if (`"`options'"' == "midisplay") {
        mi estimate
      }
      else {
  			Display_c_mi, `options' `nostructural'
      }
		}
		else {
      Display_c_mi, `old_options' `nostructural'
		}
		exit
	}
	
	if (_caller() < 8) {
		local version : display "version " string(_caller()) ", missing :"
	}
	else {
		local version : display "version " string(_caller()) " :"
	}
	
  `version' Estimate_c_mi `0'  // version is not necessary
end

program Estimate_c_mi, eclass
	version 15.1
	syntax anything [if] [in], [ MIOptions(string) MIDisplay * ]
	
	/* Options:
	   --------
		 mioptions(string)							--> options to pass to mi estimate
		 midisplay							        --> options to show results of mi estimate
	 */
	
  local tempnamelist

  /* Recover some plssem properties */
  quietly plssemc `anything', `options'
	local props = e(properties)
  local allindicators = e(mvs)
  local modeA = e(reflective)
  local binary `"`e(binarylvs)'"'
	local init : word 1 of `props'
	local wscheme : word 2 of `props'
  /* End recovering plssem properties */
  
	/* Recover the structural option */
  local structural
  local struct "structural"
	local isstruct : list struct in props
  if (!`isstruct') {
    local structural "nostructural"
  }
	local num_ind: word count `allindicators'
	local num_lv_A: word count `modeA'
	/* End of recovering the structural option */

	/* Recover the digits option */
  local digits
  gettoken part rest : options
  while `"`part'"' != "" {
    local lopt = length("`part'")
    if (substr("`part'", 1, min(3, `lopt')) == "dig") {
      tokenize `"`part'"', parse("()")
      local digits = `3'
      continue, break
    }
    gettoken part rest : rest
  }
  if ("`digits'" == "") {
    local digits = 3
  }
	/* End of recovering the digits option */
  
  /* Recover the boot option */
  local bootstrap "bootstrap"
	local isboot : list bootstrap in props
  local boot
  if (`isboot') {
    local boot "bootstrap"
  }
	/* End of recovering the boot option */
	
	/* Recover the rawsum option */
  local rawsum
  gettoken part rest : options
  while `"`part'"' != "" {
    local lopt = length("`part'")
    if (substr("`part'", 1, min(3, `lopt')) == "raw") {
      local rawsum = "rawsum"
      continue, break
    }
    gettoken part rest : rest
  }
	/* End of recovering the rawsum option */
	
	/* Recover the noheader option */
  local header
  gettoken part rest : options
  while `"`part'"' != "" {
    local lopt = length("`part'")
    if (substr("`part'", 3, min(4, `lopt')) == "head") {
      local header = "noheader"
      continue, break
    }
    gettoken part rest : rest
  }
	/* End of recovering the noheader option */

	/* Recover the nomeastable option */
  local meastable
  gettoken part rest : options
  while `"`part'"' != "" {
    local lopt = length("`part'")
    if (substr("`part'", 3, min(4, `lopt')) == "meas") {
      local meastable = "nomeastable"
      continue, break
    }
    gettoken part rest : rest
  }
	/* End of recovering the nomeastable option */

	/* Recover the nostructtable option */
  local structtable
  gettoken part rest : options
  while `"`part'"' != "" {
    local lopt = length("`part'")
    if (substr("`part'", 3, min(6, `lopt')) == "struct") {
      local structtable = "nostructtable"
      continue, break
    }
    gettoken part rest : rest
  }
	/* End of recovering the nostructtable option */

	/* Recover the loadpval option */
  local loadpval
  gettoken part rest : options
  while `"`part'"' != "" {
    local lopt = length("`part'")
    if (substr("`part'", 1, min(5, `lopt')) == "loadp") {
      local loadpval = "loadpval"
      continue, break
    }
    gettoken part rest : rest
  }
	/* End of recovering the loadpval option */

  if ("`midisplay'" == "") {
    /* Multiple imputation analysis */
    mi estimate, cmdok noheader notable `mioptions': plssemc `anything', `options'
    /* End of multiple imputation analysis */
    
    /* Display results */
    Display_c_mi, wscheme(`wscheme') init(`init') `structural' ///
      digits(`digits') boot(`boot') `header' `meastable' `structtable' ///
      `loadpval' binary(`binary') `rawsum'
    /* End of displaying results */
  }
  else {
    /* Multiple imputation analysis */
    mi estimate, cmdok `mioptions': plssemc `anything', `options'
    /* End of multiple imputation analysis */
  }

	/* Clean up */
	if ("`cleanup'" == "") {
		capture mata: cleanup()
    // capture mata: cleanup(st_local("tempnamelist"))
	}
	/* End of cleaning up */
end

program Display_c_mi
	version 15.1
	syntax [, WScheme(string) INit(string) noSTRuctural DIGits(integer 3) ///
    Boot(string) noHEADer noMEAStable noSTRUCTtable LOADPval ///
    BINary(namelist min=1) RAWsum * ]
  
  if (`digits' < 0) {
		display as error "number of digits to display must be a nonnegative integer"
		exit
	}
	
	if ("`header'" == "") {
		if ("`structural'" != "nostructural") {
      mkheader_mi, wscheme(`wscheme') init(`init') digits(5)
		}
		else {
			mkheader_mi, wscheme(`wscheme') init(`init') digits(5) nostructural
		}
	}

	if ("`meastable'" == "") {
		tempname adj_meas b_mi V_mi loadings
		mata: `adj_meas' = st_matrix("e(adj_meas)")
		mata: `b_mi' = st_matrix("e(b_mi)")
		mata: `b_mi' = `b_mi'[1, 1..rows(`adj_meas')]
		mata: st_matrix("`loadings'", editvalue(`adj_meas' :* `b_mi'', 0, .))
    matrix rownames `loadings' = `: rownames e(adj_meas)'
    matrix colnames `loadings' = `: colfullnames e(adj_meas)'
		local num_lv = colsof(`loadings')
		local allformative = e(formative)
		local num_lv_B : word count `allformative'
		local num_lv_A = `num_lv' - `num_lv_B'
		local num_ind = rowsof(`loadings')
		local title_meas "Measurement model - Standardized loadings"
		mktable, matrix(`loadings') digits(`digits') firstcolname("") ///
			title(`title_meas') firstcolwidth(14) colwidth(14) ///
			hlines(`num_ind') novlines
		
		if ("`loadpval'" != "") {
			tempname SE_mi loadings_se nindblock loadings_df loadings_pval
      mata: `SE_mi' = sqrt(diagonal(st_matrix("e(V_mi)")))'
      mata: `SE_mi' = `SE_mi'[1, 1..rows(`adj_meas')]
      mata: st_matrix("`loadings_se'", editvalue(`adj_meas' :* `SE_mi'', 0, .))
      matrix rownames `loadings_se' = `: rownames e(adj_meas)'
      matrix colnames `loadings_se' = `: colfullnames e(adj_meas)'
			mata: st_matrix("`nindblock'", colsum(`adj_meas'))
			mata: st_matrix("`loadings_df'", ///
				`adj_meas'*(st_numscalar("e(N)") - 1))
			local alllatents = e(lvs)
			if (`num_lv_B' > 0) {
				local indi = 1
				local p = 1
				foreach lv in `alllatents' {
					local nbl = `nindblock'[1, `p']
					forvalues k = 1/`nbl' {
						if (`: list lv in allformative') {
							matrix `loadings_df'[`indi', `p'] = `loadings_df'[`indi', `p'] - ///
								`nindblock'[1, `p']
						}
						else {
							matrix `loadings_df'[`indi', `p'] = `loadings_df'[`indi', `p'] - 1
						}
						local ++indi
					}
					local ++p
				}
			}
			mata: st_matrix("`loadings_pval'", 2*ttail(st_matrix("`loadings_df'"), ///
				abs(st_matrix("`loadings'") :/ st_matrix("`loadings_se'"))))
      // mata: _editmissing(st_matrix("`loadings_pval'"), 0)
			matrix rownames `loadings_pval' = `: rowfullnames `loadings''
			matrix colnames `loadings_pval' = `: colfullnames `loadings''
			if ("`boot'" == "") {
				local title_meas "Measurement model - Standardized loading p-values"
			}
			else {
				local title_meas "Measurement model - Standardized loading p-values (Bootstrap)"
			}
			mktable, matrix(`loadings_pval') digits(`digits') firstcolname("") ///
				title(`title_meas') firstcolwidth(14) colwidth(14) ///
				hlines(`num_ind') novlines
		}
	}

	if ("`structural'" != "nostructural") {
		if ("`structtable'" == "") {
      tempname adj_meas adj_struct path_mi pathSE_mi path_pval_mi ///
        pathtab endo exo
      mata: `adj_meas' = st_matrix("e(adj_meas)")
      mata: `adj_struct' = st_matrix("e(adj_struct)")
      mata: `path_mi' = st_matrix("e(b_mi)")
      mata: `path_mi' = ///
        `path_mi'[1, (sum(`adj_meas') + 1)..(sum(`adj_meas') + sum(`adj_struct'))]'
      mata: `path_mi' = plssem_vec_to_matrix(`path_mi', `adj_struct')
      mata: `pathSE_mi' = sqrt(diagonal(st_matrix("e(V_mi)")))
      mata: `pathSE_mi' = ///
        `pathSE_mi'[(sum(`adj_meas') + 1)..(sum(`adj_meas') + sum(`adj_struct')), 1]
      mata: `pathSE_mi' = plssem_vec_to_matrix(`pathSE_mi', `adj_struct')
      local alllatents = e(lvs)
      mata: `path_pval_mi' = ///
        plssem_pval( ///
          `path_mi', ///
          `pathSE_mi', ///
          "`alllatents'", ///
          "`binary'", ///
          st_numscalar("e(N_mi)"), ///
          ("`boot'" != ""))
      mata: st_matrix("`pathtab'", plssem_pathtab(`path_mi', `path_pval_mi'))
      mata: st_matrix("`endo'", colsum(`adj_struct') :> 0)
      mata: st_matrix("`exo'", rowsum(`adj_struct') :> 0)
      mata: st_local("num_lv", strofreal(cols(`adj_struct')))
      local pathtab_rown ""
      local idx = 1
      foreach lv in `: rownames e(adj_struct)' {
        if (`exo'[`idx', 1]) {
          local pathtab_rown "`pathtab_rown' `lv' ."
        }
        local ++idx
      }
      local pathtab_coln ""
      local idx = 1
      foreach lv in `: colnames e(adj_struct)' {
        if (`endo'[1, `idx']) {
          local pathtab_coln "`pathtab_coln' `lv'"
        }
        local ++idx
      }
      matrix rownames `pathtab' = `pathtab_rown'
      matrix colnames `pathtab' = `pathtab_coln'
      
			local hline_path = rowsof(`pathtab')
			if ("`boot'" == "") {
				local title_st "Structural model - Standardized path coefficients"
			}
			else {
				local title_st "Structural model - Standardized path coefficients (Bootstrap)"
			}
			mktable, matrix(`pathtab') digits(`digits') firstcolname("Variable") ///
				title(`title_st') firstcolwidth(14) colwidth(14) ///
				hlines(`hline_path') novlines path binary(`binary')
		}
	}
end
