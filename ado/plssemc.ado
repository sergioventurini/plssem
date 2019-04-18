*!plssemc version 0.3.0
*!Written 18Apr2019
*!Written by Sergio Venturini and Mehmet Mehmetoglu
*!The following code is distributed under GNU General Public License version 3 (GPL-3)

program plssemc, byable(onecall)
	version 15.1
	syntax [anything] [if] [in] [, * ]
	
	if replay() {
		if ("`e(cmd)'" != "plssem") {
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
		if (`hasstruct' == 0) {
			local nostructural "nostructural"
		}
		
		if (`"`options'"' != "") {
			Display, `options' `nostructural'
		}
		else {
			Display, `old_options' `nostructural'
		}
		exit
	}
	
	if _by() {
		local BY `"by `_byvars' `_byrc0':"'
	}
	if (_caller() < 8) {
		local version : display "version " string(_caller()) ", missing :"
	}
	else {
		local version : display "version " string(_caller()) " :"
	}
	
	`version' `BY' Estimate_c `0'		// version is not necessary
end

program Estimate_c, eclass byable(recall)
	version 15.1
	syntax anything(id="Measurement model" name=blocks) [if] [in], ///
		[ STRuctural(string) Wscheme(string) BINary(namelist min=1) ///
		Boot(numlist integer >0 max=1) Seed(numlist max=1) Tol(real 1e-7) ///
		MAXiter(integer 100) MISsing(string) k(numlist integer >0 max=1) ///
		INIT(string) DIGits(integer 3) noHEADer noMEAStable noDISCRIMtable ///
		noSTRUCTtable LOADPval STATs CORRelate(string) RAWsum noSCale ///
		CONVcrit(string) noCLeanup ]
	
	/* Options:
	   --------
		 structural(string)							--> structural model specification
		 wscheme(string)								--> weighting scheme ("centroid", "factorial"
																				or "path")
		 binary(namelist min=1)					--> latent variables to fit using logit
		 boot(numlist integer >0 max=1)	--> bootstrap estimation (# of replications)
		 seed(numlist max=1)						--> bootstrap seed number
		 tol(real 1e-7)									--> tolerance (default 1e-7)
		 maxiter(integer 100)						--> maximum number of iterations (default
																				100)
		 missing(string)								--> imputation method ("mean" or "knn")
		 k(numlist integer >0 max=1)		--> neighbor size; to use with
																				missing("knn") (default 5)
		 init(string)										--> initialization method ("eigen" or "indsum")
		 digits(integer 3)							--> number of digits to display (default 3)
		 noheader												--> do not show the header
		 nomeastable										--> do not show the measurement table
		 nodiscrimtable									--> do not show the discriminant validity
																				table
		 nostructtable									--> do not show the structural table
		 loadpval												--> show the loadings' p-values
		 stats													--> print a table of summary statistics for
																				the indicators
		 correlate(string)							--> report the correlation among the
																				indicators, the LVs as well as the cross
																				loadings; accepts any combination of 'mv',
																				'lv' and 'cross' (i.e. cross-loadings);
																				suboption 'cutoff()' can be used to avoid
																				showing correlations smaller than the
																				specified value
		 rawsum													--> estimate the latent scores as the raw
																				sum of the indicators
		 noscale												--> manifest variables are not standardized
		 convcrit												--> convergence criterion (either 'relative'
																				or 'square')
		 nocleanup											--> Mata temporary objects are not removed
																				(undocumented)
	 */
	
	/* Macros:
	 * 	LV`i'							- latent variable name
	 * 	i`i'							- indicators for latent, by latent index
	 * 	istd`i'						- indicators for latent, by latent index (standardized)
	 * 	i`LV`i''					- indicators for latent, by latent name
	 * 	istd`LV`i''				- indicators for latent, by latent name (standardized)
	 *	allindicators			- all indicators
	 *	allstdindicators	- all standardized indicators
	 *	alllatents				- all latents
	 */

	local cmdline : list clean 0

	/* Parse the blocks */
	local blocks : list clean blocks
	tokenize `"`blocks'"', parse(" ()<>")
	
	tempname inblock startblock
	scalar `inblock' = 0
	scalar `startblock' = 0
	
	local j = 0
	local tok_i = 1
	
	while ("``tok_i''" != "") {
		if ("``tok_i''" == "(") {
			if (`inblock') {
				display as error "unexpected ("
				error 197
			}
			scalar `inblock' = 1
			scalar `startblock' = 1
			local ++j
		}
		else if (`inblock') {
			if ("``tok_i''" == ")") {
				if ("LV`j'" == "" | "i`j'" == "") {
					display as error "incomplete measurement model specification"
					error 197
				}
				else {
					scalar `inblock' = 0
					local i`LV`j'' `i`j''
					local allindicators "`allindicators' `i`j''"
					local alllatents "`alllatents' `LV`j''"
				}
			}
			else if ("``tok_i''" == "<" | "``tok_i''" == ">") {
				scalar `startblock' = 0
				if ("``tok_i''" == ">") {
					local modeA "`modeA' `LV`j''"
				}
				else if ("``tok_i''" == "<") {
					local modeB "`modeB' `LV`j''"
				}
			}
			else if (`startblock') {
				if ("`LV`j''" != "") {
					display as error "missing ="
					error 197
				}
				if (_byindex() == 1) {  // this provides the LVs predictions when by'd
					capture confirm new variable ``tok_i''
					if (_rc == 110) {
						capture quietly drop ``tok_i''
					}
				}
				local LV`j' ``tok_i''
			}
			else {
				foreach var of varlist ``tok_i'' {
					local i`j' "`i`j'' `var'"
				}
			}
		}
		else {
			error 197
		}
		local ++tok_i
	}
	local modeA : list clean modeA
	local modeB : list clean modeB
	local modeA : list uniq modeA
	local modeB : list uniq modeB
	
	if (`inblock') {
		display as error "missing )"
		error 197
	}
	
	local allindicators : list clean allindicators
	local alllatents : list clean alllatents
	/* End of parsing the blocks */
	
	/* Parse the binary() option */
	if ("`binary'" != "") {
		noisily {
			display
			display as error "WARNING: the use of binary latent variables " _continue
			display as error "goes beyond the original scopes of PLS-SEM. "
			display as error "         They are provided here only for exploratory " _continue
			display as error "purposes and we suggest not to "
			display as error "         report the corresponding results " _continue
			display as error "in any published study or report."
		}
		local binary : list clean binary
		foreach var in `binary' {
			local nind_bin : word count `i`var''
			if (`nind_bin' > 1) {
				display as error "latent variables in 'binary()' option must be defined through a single binary indicator"
				exit
			}
		}
	}
	/* End of parsing the binary() option */
	
	/* Set obs to use (1/2) */
	marksample touse
	/* End of setting observations to use (1/2) */
	
	/* Save original data set */
	if ("`missing'" != "") {
		tempname original_data
		mata: `original_data' = st_data(., "`: list uniq allindicators'", "`touse'")
	}
	/* End saving original data set */
	
	/* Imputation of missing values */
	if ("`missing'" != "") {
		local missing : list clean missing
		foreach var in `binary' {
			local indbin "`indbin' `i`var''"
		}
		local indbin : list clean indbin
		
		if ("`missing'" == "mean") {
			mata: meanimp( ///
				st_data(., "`allindicators'", "`touse'"), ///
				"`allindicators'", ///
				"`indbin'", ///
				"`touse'")
		}
		else if ("`missing'" == "knn") {
			tempname ind_mata allmiss
			mata: `ind_mata' = st_data(., "`allindicators'", "`touse'")
			mata: st_numscalar("`allmiss'", anyof(rowmissing(`ind_mata'), ///
				cols(`ind_mata')))
			if (`allmiss') {
				display as error "some observations have all indicators missing; " _continue
				display as error "remove them manually before proceeding"
				exit
			}
			
			if ("`k'" == "") {
				local k = 5
			}
			mata: knnimp( ///
				st_data(., "`allindicators'", "`touse'"), ///
				"`allindicators'", ///
				strtoreal("`k'"), ///
				"`indbin'", ///
				"`touse'", ///
				1)
		}
		else {
			display as error "missing can be either 'mean' or 'knn'"
			exit
		}
	}
	/* End of imputing the missing values */
	
	/* Set obs to use (2/2) */
	tempvar touse_nomiss
	quietly generate `touse_nomiss' = `touse'
	markout `touse' `allindicators'
	
	quietly count if `touse'
	if (r(N) == 0) {
		error 2000
	}
	else if (r(N) < 10) {
		display as error "at least 10 complete observations required"
		error 2001
	}

	quietly misstable patterns `allindicators' if `touse'
	local nobs = r(N_complete)
	/* End of setting observations to use (2/2) */
	
	/* Check that there are no "zero-variance" indicators */
	tempname zerovar
	mata: st_numscalar("`zerovar'", any(selectindex( ///
			sd(st_data(., "`allindicators'", "`touse'")) :== 0)))
	if (`zerovar') {
		display as error "at least one indicator has zero variance in one of the groups"
		exit
	}
	/* End of checking that there are no "zero-variance" indicators */
	
	/* Parse the inner weight scheme */	
	if ("`wscheme'" == "") {
		local wscheme "path"
	}
	if !("`wscheme'" == "centroid" | "`wscheme'" == "factorial" | "`wscheme'" == "path") {
		display as error "scheme can be either 'centroid', 'factorial', or 'path'"
		exit
	}
	/* End of parsing the inner weight scheme */

	/* Check convergence parameters */
	if (`maxiter' <= 0) {
		display as error "maximum number of iterations must be a positive integer"
		exit
	}
	if (`tol' <= 0) {
		display as error "tolerance must be a strictly positive number"
		exit
	}
	if ("`convcrit'" == "") {
		local convcrit "relative"
	}
	if !("`convcrit'" == "relative" | "`convcrit'" == "square") {
		display as error "convergence criterion can be either 'relative' or 'square'"
		exit
	}
	/* End of checking convergence parameters */

	/* Check that digits is nonnegative */
	if (`digits' < 0) {
		display as error "number of digits to display must be a nonnegative integer"
		exit
	}
	/* End of checking that digits is nonnegative */
	
	/* Check the initialization method */
	if ("`init'" == "") {
		local init "indsum"
	}
	if !("`init'" == "eigen" | "`init'" == "indsum") {
		display as error "initialization method can be either 'eigen' or 'indsum'"
		error 198
	}
	if (("`structural'" == "") & ("`init'" != "eigen")) {
		display as error "'init()' option must be set to 'eigen' when no structural model is provided"
		error 198
	}
	/* End of checking the initialization method */

	/* Parse the correlate() option */
	if ("`correlate'" != "") {
		local tok_i = 1
		tokenize "`correlate'", parse(",() ")
		while ("``tok_i''" != "") {
			if ("``tok_i''" == "mv") {
				local corrind = "corrind"
			}
			if ("``tok_i''" == "lv") {
				local corrlv = "corrlv"
			}
			if ("``tok_i''" == "cross") {
				local crossload = "crossload"
			}
			if ("``tok_i''" == "cutoff") {
				local tok_tmp = `tok_i' + 2
				local cutoff = "``tok_tmp''"
				if (`cutoff' < 0 | `cutoff' > 1) {
					display as error "cutoff must be in between 0 and 1"
					exit
				}
			}
			local ++tok_i
		}
		if ("`corrind'" == "" & "`corrlv'" == "" & "`crossload'" == "") {
			display as error "'correlate()' option requires at least one of 'mv', 'lv' or 'cross'"
			exit
		}
	}
	/* End of parsing the correlate() option */
	
	/* Parse inner relationships */
	local num_lv: word count `alllatents'
	
	if ("`structural'" != "") {
		tokenize `"`structural'"', parse(",")
		
		local reg3eqs = "("
		local tok_i = 1
		while ("``tok_i''" != "") {
			// Parse the interactions
			if (strpos("``tok_i''", "*")) {
				local inter "``tok_i''"
				gettoken depvar : inter
				if (strpos("`depvar'", "*")) {
					display as error "`depvar': interactions cannot involve a dependent variable"
					exit
				}
				local new_tok "`depvar'"
				gettoken blk inter : inter
				local inter : list clean inter
				while ("`inter'" != "") {
					gettoken toadd inter : inter
					local inter : list clean inter
					if (strpos("`toadd'", "*")) {
						local inter_ind ""
						gettoken var1 toadd : toadd, parse("*")
						gettoken ast toadd : toadd, parse("*")
						gettoken var2 toadd : toadd, parse("*")
						if (strpos("`toadd'", "*")) {
							display as error "interactions can involve only two variables"
							exit
						}

						local internm = "`var1'`var2'"
						if (_byindex() == 1) {  // this provides the LVs predictions when by'd
							capture confirm new variable `internm'
							if (_rc == 110) {
								capture quietly drop `internm'
							}
						}

						// Add the interactions to the list of LVs
						local alllatents "`alllatents' `internm'"
						local modeA "`modeA' `internm'"
						local num_lv: word count `alllatents'
						local LV`num_lv' `internm'
						
						// Check that the new interaction is reflective
						local int_in_modeB : list internm in modeB
						if (`int_in_modeB') {
							display as error "interactions of latent variables must be reflective"
							exit
						}

						// Generate the product indicators
						foreach ind1 in `i`var1'' {
							foreach ind2 in `i`var2'' {
								local indnm = "`ind1'`ind2'"
								if (_byindex() == 1) {  // this provides the indicators when by'd
									capture confirm new variable `indnm'
									if (_rc == 110) {
										capture quietly drop `indnm'
									}
									quietly generate `indnm' = `ind1'*`ind2' if `touse'
								}
								else {
									quietly replace `indnm' = `ind1'*`ind2' if `touse'
								}
								
								// Add the product indicator to the list of indicators
								local allindicators "`allindicators' `indnm'"
								local inter_ind "`inter_ind' `indnm'"
								
								// Create a macro for the product indicators to use for scaling
								local inter_for_scale "`inter_for_scale' `ind1'%£%]_[%£%`ind2'"
							}
						}
						local i`num_lv' `inter_ind'
						local i`internm' `inter_ind'

						local new_tok "`new_tok' `var1' `var2' `internm'"
					}
					else {
						local new_tok "`new_tok' `toadd'"
					}
				}
				local reg3eqs "`reg3eqs'`new_tok') ("
				local 0 `new_tok'
			}
			else {
				local reg3eqs "`reg3eqs'``tok_i'') ("
				local 0 ``tok_i''
			}
			syntax namelist(min=2)
			
			// Use macros r`var' and c`latentname' to store whether the adjacency is
			// treated as directional in path weighting scheme.
			// The regest`var' macro is needed to store the path coefficients.
			
			local dvar : word 1 of `namelist'
			local ivars : list namelist - dvar
			local r`dvar' `r`dvar'' `ivars' // r`dvar' contains the predecessors of each LV
			local regest`dvar' `regest`dvar'' `ivars'
			
			foreach iv in `ivars' {
				local c`iv' `c`iv'' `dvar' // c`iv' contains the successors of each LV
			}
			
			local tok_i = `tok_i' + 2
		}
		local parenthesis = "("
		local reg3eqs : list reg3eqs - parenthesis
		local reg3eqs : list clean reg3eqs

		foreach var in `alllatents' {
			if ("`r`var''`c`var''" == "") {
				display as error "latent `var' is not adjacent to any other latent"
				exit
			}
			
			local regest`var' : list uniq regest`var'
			
			if ("`wscheme'" == "path") {
				local c`var' : list uniq c`var'
				local r`var' : list uniq r`var'
				local c`var' : list c`var' - r`var'
			}
			else {
				local c`var' `c`var'' `r`var''
				local c`var' : list uniq c`var'
				local r`var' ""
			}
		}
	}
	/* End of parsing the inner relationships */
	
	/* Create other macros and scalars */
	foreach var in `allindicators' {
		local allstdindicators "`allstdindicators' std`var'"
	}
	local allstdindicators : list clean allstdindicators
	forvalues k = 1/`num_lv' {
		foreach var in `i`k'' {
			local istd`k' "`istd`k'' std`var'"
			local istd`LV`k'' `istd`k''
		}
	}

	tempname rawsum_sc struct_sc
	if ("`rawsum'" == "") {
		scalar `rawsum_sc' = 0
	}
	else {
		scalar `rawsum_sc' = 1
	}
	if ("`structural'" == "") {
		scalar `struct_sc' = 0
	}
	else {
		scalar `struct_sc' = 1
	}
	/* End of creating other macros */
	
	/* Create the adjacency matrices */
	tempname modes adj_meas adj_struct
	local num_ind : word count `allindicators'

	matrix `modes' = J(`num_lv', 1, 1)
	local i = 1
	foreach var in `alllatents' {
		if (`: list var in modeA') {
			matrix `modes'[`i', 1] = 0
		}
		local ++i
	}
	
	matrix `adj_meas' = J(`num_ind', `num_lv', 0)	
	local i = 1
	local j = 1
	foreach var in `alllatents' {
		foreach var2 in `i`var'' {
			matrix `adj_meas'[`i', `j'] = 1
			local ++i
		}
		if (`: list var in modeA') {
			local loadcolnames `loadcolnames' "Reflective:`var'"
		}
		else {
			local loadcolnames `loadcolnames' "Formative:`var'"
		}
		local loadrownames `loadrownames' `i`var''
		local ++j
	}
	matrix rownames `adj_meas' = `loadrownames'
	matrix colnames `adj_meas' = `loadcolnames'

	matrix `adj_struct' = J(`num_lv', `num_lv', 0)
	local B_j = 1
	foreach var in `alllatents' {
		if ("`regest`var''" != "") {
			forvalues lv_i = 1/`num_lv' {
				local xn : word `lv_i' of `alllatents'
				local isthere : list xn in regest`var'
				if (`isthere') {
					local B_i : list posof "`xn'" in alllatents
					matrix `adj_struct'[`B_i', `B_j'] = 1
				}
			}
		}
		local ++B_j
	}
	matrix rownames `adj_struct' = `alllatents'
	matrix colnames `adj_struct' = `alllatents'
	/* End of creating the adjacency matrices */
	
	capture noisily {
		/* Standardize the MVs (if requested) */
		mata: ///
			plssem_scale( ///
				st_data(., "`allindicators'", "`touse'"), ///
				"`allstdindicators'", ///
				"`touse'", ///
				"`scale'")
		/* End of standardizing the MVs */
		
		/* Initialize LVs */
		mata: ///
			plssem_init( ///
				st_data(., "`allstdindicators'", "`touse'"), ///
				st_matrix("`adj_meas'"), ///
				"`allindicators'", ///
				"`allstdindicators'", ///
				"`alllatents'", ///
				"`touse'", ///
				st_numscalar("`rawsum_sc'"), ///
				"`init'")
		/* End of initializing the LVs */
		
		/* Run the PLS algorithm */
		tempname res converged outerW iter matreldiff Whistory
		
		mata: `res' = ///
			plssem_base( ///
				st_data(., "`allstdindicators'", "`touse'"), ///
				st_data(., "`alllatents'", "`touse'"), ///
				st_matrix("`adj_meas'"), ///
				st_matrix("`adj_struct'"), ///
				st_matrix("`modes'"), ///
				"`alllatents'", ///
				"`binary'", ///
				strtoreal("`tol'"), ///
				strtoreal("`maxiter'"), ///
				"`touse'", ///
				"`wscheme'", ///
				"`convcrit'", ///
				st_numscalar("`struct_sc'"), ///
				st_numscalar("`rawsum_sc'"), ///
				1)
		
		mata: st_numscalar("`converged'", `res'.converged)
		mata: st_numscalar("`iter'", `res'.niter)
		if (("`structural'" != "") & ("`rawsum'" == "")) {
			mata: st_matrix("`matreldiff'", `res'.diff)
			mata: st_matrix("`Whistory'", `res'.evo)
			foreach var in `alllatents' {
				foreach var2 in `i`var'' {
					local Whistcolnames `Whistcolnames' "`var2':`var'"
				}
			}
			matrix colnames `Whistory' = `Whistcolnames'
			mata: st_matrix("`outerW'", `res'.outer_weights)
			matrix rownames `outerW' = `loadrownames'
			matrix colnames `outerW' = `loadcolnames'
		}
		/* End of the PLS algorithm */
		
		/* Label the LVs */
		local now "`c(current_date)', `c(current_time)'"
		local now : list clean now
		foreach var of varlist `alllatents' {
			label variable `var' "Scores of `var' latent variable [`now']"
		}
		/* End of labeling the LVs */
		
		/* Compute the model parameters' variances */
		tempname xload_v path_v
		
		if ("`boot'" == "") {
			// in PLSc variances are available only using bootstrap
			mata: `xload_v' = J(rows(`res'.loadings_c), cols(`res'.loadings_c), .)
			mata: `path_v' = J(rows(`res'.path_c), cols(`res'.path_c), .)
		}
		else {
			tempname res_bs
			
			mata: `res_bs' = ///
				plssem_boot( ///
					st_data(., "`allstdindicators'", "`touse'"), ///
					st_data(., "`alllatents'", "`touse'"), ///
					st_matrix("`adj_meas'"), ///
					st_matrix("`adj_struct'"), ///
					st_matrix("`modes'"), ///
					"`alllatents'", ///
					"`binary'", ///
					strtoreal("`tol'"), ///
					strtoreal("`maxiter'"), ///
					"`touse'", ///
					"`wscheme'", ///
					"`convcrit'", ///
					st_numscalar("`struct_sc'"), ///
					st_numscalar("`rawsum_sc'"), ///
					1, ///
					strtoreal("`boot'"), ///
					strtoreal("`seed'"), ///
					1)
			mata: `xload_v' = `res_bs'.xloadings_v
			if ("`structural'" != "") {
				mata: `path_v' = `res_bs'.path_v
			}
		}
		/* End of computing the model parameters' variances */
	} // end of -capture-
	local rc = _rc
	if (`rc' == 1) {
		display
		display as error "you pressed the Break key; " _continue
		display as error "calculation interrupted"
	}
	if (`rc' >= 1) {
		/* Clean up */
		foreach var in `allstdindicators' {
			capture quietly drop `var'
		}
		if ("`rawsum'" != "") {
			foreach var in `alllatents' {
				quietly replace `var' = rs_`var'
				capture quietly drop rs_`var'
			}
		}
		if ("`missing'" != "") {
			tempvar __touse__
			quietly generate `__touse__' = e(sample)
			mata: st_store(., tokens("`: list uniq allindicators'"), "`__touse__'", ///
				`original_data')
		}
		if ("`cleanup'" == "") {
			capture mata: cleanup()
		}
		/* End of cleaning up */
		
		error `rc'
	}
	
	/* Compute the table of measurement loadings */
	tempname loadings xloadings annihilate_se loadings_se xloadings_se

	mata: st_matrix("`loadings'", `res'.loadings_c)
	matrix rownames `loadings' = `loadrownames'
	matrix colnames `loadings' = `loadcolnames'

	mata: st_matrix("`xloadings'", `res'.xloadings_c)
	matrix rownames `xloadings' = `loadrownames'
	matrix colnames `xloadings' = `loadcolnames'

	if ("`boot'" != "") {
		tempname loadings_bs xloadings_bs
		
		mata: st_matrix("`loadings_bs'", `res_bs'.loadings_bs)
		matrix rownames `loadings_bs' = `loadrownames'
		matrix colnames `loadings_bs' = `loadcolnames'

		mata: st_matrix("`xloadings_bs'", `res_bs'.xloadings_bs)
		matrix rownames `xloadings_bs' = `loadrownames'
		matrix colnames `xloadings_bs' = `loadcolnames'
	}

	mata: `annihilate_se' = editvalue(st_matrix("`adj_meas'"), 0, .)
	mata: st_matrix("`loadings_se'", sqrt(`xload_v') :* `annihilate_se')
	matrix rownames `loadings_se' = `loadrownames'
	matrix colnames `loadings_se' = `loadcolnames'

	mata: st_matrix("`xloadings_se'", sqrt(`xload_v'))
	matrix rownames `xloadings_se' = `loadrownames'
	matrix colnames `xloadings_se' = `loadcolnames'
	/* End of computing the table of measurement loadings */
	
	/* Compute the table of average variance extracted (AVE) */
	local num_ind: word count `allindicators'
	local num_lv_A: word count `modeA'
	
	if (`num_lv_A' > 0) {
		tempname sqcorr
		
		/*
		mata: st_matrix("`sqcorr'", ///
			correlation(st_data(., "`modeA'", "`touse'")) :^ 2)
		matrix rownames `sqcorr' = `modeA'
		matrix colnames `sqcorr' = `modeA'
		*/
		mata: st_matrix("`sqcorr'", `res'.R_c :^ 2)
		matrix rownames `sqcorr' = `alllatents'
		matrix colnames `sqcorr' = `alllatents'
	}

	tempname ave_num ave_mata whichB ave
	
	mata: st_matrix("`ave_num'", colsum(st_matrix("`adj_meas'")))
	mata: `ave_mata' = ///
		colsum(st_matrix("`loadings'") :^ 2) :/ st_matrix("`ave_num'")
	if ("`modeB'" != "") {
		mata: `whichB' = colsum(strmatch(tokens("`alllatents'"), tokens("`modeB'")'))
		mata: `ave_mata'[., selectindex(`whichB')] = J(1, sum(`whichB'), .)
	}
	mata: st_matrix("`ave'", `ave_mata')
	matrix colnames `ave_num' = `alllatents'
	matrix rownames `ave' = "AVE"
	matrix colnames `ave' = `alllatents'
	/* End of computing the AVE table */
	
	/* Compute the table of summary statistics for indicators */
	// ("if `touse'" is not used here to correctly account for the missing values)
	if ("`stats'" != "") {
		tempname indstats
		local kk = 1
		matrix `indstats' = J(`num_ind', 7, .)
		foreach ind in `allindicators' {
			quietly summarize `ind' if `touse_nomiss', detail
			matrix `indstats'[`kk', 1] = r(mean)
			matrix `indstats'[`kk', 2] = r(sd)
			matrix `indstats'[`kk', 3] = r(p50)
			matrix `indstats'[`kk', 4] = r(min)
			matrix `indstats'[`kk', 5] = r(max)
			mata: st_local("indstats_mata", ///
				strofreal(nonmissing(st_data(., "`ind'", "`touse_nomiss'"))))
			matrix `indstats'[`kk', 6] = real("`indstats_mata'")
			mata: st_local("indstats_mata", ///
				strofreal(missing(st_data(., "`ind'", "`touse_nomiss'"))))
			matrix `indstats'[`kk', 7] = real("`indstats_mata'")
			local ++kk
		}
		matrix rownames `indstats' = `allindicators'
		matrix colnames `indstats' = "mean" "sd" "median" "min" "max" "N" "missing"
	}
	/* End of computing the table of summary statistics for indicators */
	
	/* Compute the table of structural path coefficients */
	if ("`structural'" != "") {
		tempname path path_se rsquared redundancy path_pval pathtab path_toteff
		
		mata: st_matrix("`path'", `res'.path_c)
		matrix rownames `path' = `alllatents'
		matrix colnames `path' = `alllatents'

		if ("`boot'" != "") {
			tempname path_bs
			
			mata: st_matrix("`path_bs'", `res_bs'.path_bs)
			matrix rownames `path_bs' = `alllatents'
			matrix colnames `path_bs' = `alllatents'
		}
		
		mata: st_matrix("`path_se'", sqrt(`path_v'))
		matrix rownames `path_se' = `alllatents'
		matrix colnames `path_se' = `alllatents'

		mata: st_matrix("`rsquared'", `res'.r2_c)
		matrix rownames `rsquared' = "r2"
		matrix colnames `rsquared' = `alllatents'
		
		mata: st_matrix("`redundancy'", `res'.r2_c :* st_matrix("`ave'"))
		matrix rownames `redundancy' = "redundancy"
		matrix colnames `redundancy' = `alllatents'
		
		mata: st_local("regest_all", ///
			invtokens(tokens("`alllatents'")[selectindex(rowsum(`res'.S))]))
		mata: st_local("lv_regest_all", ///
			invtokens(tokens("`alllatents'")[selectindex(colsum(`res'.S))]))
		local num_lv_path : word count `lv_regest_all'
		
		mata: `path_pval' = ///
			plssem_pval( ///
				st_matrix("`path'"), ///
				st_matrix("`path_se'"), ///
				"`alllatents'", ///
				"`binary'", ///
				strtoreal("`nobs'"), ///
				("`boot'" != ""))

		mata: st_matrix("`pathtab'", ///
			plssem_pathtab(st_matrix("`path'"), `path_pval', `res'.r2_a_c))
		foreach var in `regest_all' {
			local pathtab_rownames "`pathtab_rownames' `var' ."
		}
		matrix rownames `pathtab' = `pathtab_rownames' "r2_a"
		matrix colnames `pathtab' = `lv_regest_all'

		mata: st_matrix("`path_toteff'", `res'.total_effects_c)
		matrix rownames `path_toteff' = `alllatents'
		matrix colnames `path_toteff' = `alllatents'
	}
	/* End of computing the table of structural path coefficients */
	
	/* Compute the reliability coefficients */
	tempname relcoef
	
	mata: st_matrix("`relcoef'", ///
		plssem_reliability( ///
			st_data(., "`allindicators'", "`touse'"), ///
			st_matrix("`loadings'"), ///
			st_matrix("`modes'"), ///
			st_matrix("`outerW'")))
	matrix rownames `relcoef' = "Cronbach" "DG" "rho_A"
	matrix colnames `relcoef' = `loadcolnames'
	/* End of computing the reliability coefficients */
	
	/* Setup the tables of structural coefficients */
	if ("`structural'" != "") {
		tempname path_mata strcoef path_se_mata strse path_vif_mata strvif ///
			rowi_path coli_path R_Y
		
		mata: `path_mata' = st_matrix("`path'")
		mata: `rowi_path' = selectindex(rownonmissing(`path_mata'))
		mata: `coli_path' = selectindex(colnonmissing(`path_mata'))
		mata: st_matrix("`strcoef'", `path_mata'[`rowi_path', `coli_path'])
		foreach var in `regest_all' {
			local struct_rownames "`struct_rownames' `var'"
		}
		matrix rownames `strcoef' = `struct_rownames'
		matrix colnames `strcoef' = `lv_regest_all'
		
		mata: `path_se_mata' = st_matrix("`path_se'")
		mata: st_matrix("`strse'", `path_se_mata'[`rowi_path', `coli_path'])
		matrix rownames `strse' = `struct_rownames'
		matrix colnames `strse' = `lv_regest_all'
		
		mata: st_matrix("`R_Y'", `res'.R_c)
		matrix rownames `R_Y' = `alllatents'
		matrix colnames `R_Y' = `alllatents'
		mata: `path_vif_mata' = ///
			plssem_vif(st_matrix("`R_Y'"), st_matrix("`adj_struct'"))
		mata: st_matrix("`strvif'", ///
			`path_vif_mata'[selectindex(rownonmissing(`path_vif_mata')), ///
				selectindex(colnonmissing(`path_vif_mata'))])
		matrix rownames `strvif' = `struct_rownames'
		matrix colnames `strvif' = `lv_regest_all'
	}
	/* End of setting up the tables of structural coefficients */
	
	/* Calculate the model assessment indexes (GOF, etc.) */
	if (("`structural'" != "") & (`num_lv_A' > 0)) {
		tempname gof avrsq avave avred mod_assessment tmp tmp_w
		matrix `tmp' = (`rsquared')'
		mata: st_numscalar("`avrsq'", mean(st_matrix("`tmp'")))
		matrix `tmp' = (`ave')'
		matrix `tmp_w' = (`ave_num')'
		mata: st_numscalar("`avave'", mean(st_matrix("`tmp'"), ///
			st_matrix("`tmp_w'")))
		matrix `tmp' = (`redundancy')'
		mata: st_numscalar("`avred'", mean(st_matrix("`tmp'")))
		scalar `gof' = sqrt(`avrsq' * `avave')
		matrix `mod_assessment' = J(1, 4, .)
		matrix `mod_assessment'[1, 1] = `avrsq'
		matrix `mod_assessment'[1, 2] = `avave'
		matrix `mod_assessment'[1, 3] = `gof'
		matrix `mod_assessment'[1, 4] = `avred'
		matrix colnames `mod_assessment' = "Average R-squared" ///
			"Average communality (AVE)" "GoF" "Average redundancy"
	}
	/* End of calculating the model assessment indexes (GOF, etc.) */
	
	/* Return values */
	local props "`init'"
	if ("`structural'" != "") {
		local props "`props' `wscheme'"
	}
	local props "`props' `missing'"
	if ("`boot'" != "") {
		local props "`props' bootstrap"
	}
	if ("`structural'" != "") {
		local props "`props' structural"
	}
	if ("`rawsum'" != "") {
		local props "`props' rawsum"
	}
	if ("`scale'" == "") {
		local props "`props' scaled"
	}
	else {
		local props "`props' unscaled"
	}
	local props "`props' `convcrit'"
	ereturn post, obs(`nobs') esample(`touse') properties(`props')
	if ("`boot'" != "") {
		ereturn scalar reps = `boot'
	}
	if ("`structural'" != "") {
		if ("`rawsum'" == "") {
			ereturn scalar maxiter = `maxiter'
			ereturn scalar iterations = `iter'
		}
		else {
			ereturn scalar maxiter = .
		}
		ereturn scalar tolerance = `tol'
	}
	ereturn scalar converged = `converged'

	ereturn local struct_eqs `"`reg3eqs'"'
	ereturn local formative `"`modeB'"'
	ereturn local reflective `"`modeA'"'
	if ("`binary'" != "") {
		ereturn local binarylvs `"`binary'"'
	}
	ereturn local lvs `"`alllatents'"'
	ereturn local mvs `"`allindicators'"'
	local varlist `e(mvs)' `e(lvs)'
	signestimationsample `varlist'
	ereturn local title "Consistent partial least squares (PLSc) structural equation modeling"
	ereturn local predict plssem_p
	ereturn local estat_cmd "plssemc_estat"
	ereturn local cmdline "plssemc `cmdline'"
	ereturn local cmd "plssemc"

	if ("`missing'" != "") {
		tempname __touse__ imp_data
		quietly generate `__touse__' = e(sample)
		mata: st_matrix("`imp_data'", ///
			st_data(., tokens("`: list uniq allindicators'"), "`__touse__'"))
		matrix colnames `imp_data' = `: list uniq allindicators'
		ereturn matrix imputed_data = `imp_data'
	}
	if ("`structural'" != "") {
		ereturn matrix R = `R_Y'
		if ("`rawsum'" == "") {
			ereturn matrix reldiff = `matreldiff'
			ereturn matrix outerweights =  `outerW'
			ereturn matrix ow_history =  `Whistory'
		}
		if (`num_lv_A' > 0) {
			ereturn matrix assessment = `mod_assessment'
			ereturn matrix redundancy = `redundancy'
			ereturn matrix rsquared = `rsquared'
		}
		ereturn matrix total_effects = `path_toteff'
		ereturn matrix adj_struct = `adj_struct'
		if ("`boot'" != "") {
			mata: st_matrix("`path_bs'", editmissing(st_matrix("`path_bs'"), 0))
			matrix rownames `path_bs' = `alllatents'
			matrix colnames `path_bs' = `alllatents'
			ereturn matrix pathcoef_bs = `path_bs'
		}
		mata: st_matrix("`path'", editmissing(st_matrix("`path'"), 0))
		matrix rownames `path' = `alllatents'
		matrix colnames `path' = `alllatents'
		ereturn matrix pathcoef = `path'
		ereturn matrix struct_table = `pathtab'
		if ("`boot'" == "") {
			ereturn matrix struct_vif = `strvif'
		}
		ereturn matrix struct_se = `strse'
		ereturn matrix struct_b = `strcoef'
	}
	if (`num_lv_A' > 0) {
		ereturn matrix ave = `ave'
		ereturn matrix sqcorr = `sqcorr'
		ereturn matrix relcoef = `relcoef'
	}
	ereturn matrix adj_meas = `adj_meas'
	ereturn matrix cross_loadings_se = `xloadings_se'
	if ("`boot'" != "") {
		ereturn matrix cross_loadings_bs = `xloadings_bs'
	}
	ereturn matrix cross_loadings = `xloadings'
	ereturn matrix loadings_se = `loadings_se'
	if ("`boot'" != "") {
		ereturn matrix loadings_bs = `loadings_bs'
	}
	ereturn matrix loadings = `loadings'
	if ("`stats'" != "") {
		ereturn matrix indstats = `indstats'
	}
	/* End of returning values */
	
	/* Display results */
	if ("`display'" == "") {
		if ("`structural'" == "") {
			local nostructural "nostructural"
		}
		if (`num_lv_A' == 0) {
			local nodiscrim "nodiscrimtable"
		}
		else {
			local nodiscrim `discrimtable'
		}
		Display_c, `nostructural' digits(`digits') boot(`boot') ///
			`header' `meastable' `nodiscrim' `structtable' `loadpval' `stats' ///
			`corrind' `corrlv' `crossload' cutoff(`cutoff') binary(`binary') `rawsum'
	}
	/* End of displaying results */

	/* Maximum number of iterations reached */
	if (("`structural'" != "") & ("`rawsum'" == "")) {
		if (!`converged') {
			display as error "warning: convergence not achieved --> " _continue
			display as error "maximum number of iterations reached"
			display as error _skip(9) "the solution provided may not be acceptable; " _continue
			display as error "try to increase the 'maxiter' option"
		}
	}
	/* End of 'convergence not attained' message */
	
	/* Clean up */
	foreach var in `allstdindicators' {
		capture quietly drop `var'
	}
	if ("`rawsum'" != "") {
		foreach var in `alllatents' {
			quietly replace `var' = rs_`var'
			capture quietly drop rs_`var'
		}
	}
	if ("`missing'" != "") {
		mata: st_store(., tokens("`: list uniq allindicators'"), "`__touse__'", ///
			`original_data')
	}
	if ("`cleanup'" == "") {
		capture mata: cleanup()
	}
	/* End of cleaning up */
end

program Display_c
	version 15.1
	syntax [, noSTRuctural DIGits(integer 3) noHEADer noMEAStable ///
		noDISCRIMtable noSTRUCTtable LOADPval stats corrind corrlv crossload ///
		CUToff(real 0) BINary(namelist min=1) RAWsum * ]
	
	if (`digits' < 0) {
		display as error "number of digits to display must be a nonnegative integer"
		exit
	}
	
	local props = e(properties)
	local boot_lbl "bootstrap"
	local isboot : list boot_lbl in props
	
	tempvar __touse__
	quietly generate `__touse__' = e(sample)
	
	if ("`header'" == "") {
		if ("`structural'" != "nostructural") {
			if ("`rawsum'" == "") {
				tempname mod_assessment reldiff
				matrix `mod_assessment' = e(assessment)
				matrix `reldiff' = e(reldiff)
				mkheader, matrix1(`mod_assessment') matrix2(`reldiff') digits(5) ///
					nogroup consistent
			}
			else {
				tempname mod_assessment
				matrix `mod_assessment' = e(assessment)
				mkheader, matrix1(`mod_assessment') digits(5) nogroup rawsum consistent
			}
		}
		else {
			mkheader, digits(5) nogroup nostructural consistent
		}
	}

	if ("`stats'" != "") {
		tempname indstats
		matrix `indstats' = e(indstats)
		if (matmissing(`indstats')) {
			local allindicators = e(mvs)
			local num_ind: word count `allindicators'
			local k = 1
			matrix `indstats' = J(`num_ind', 7, .)
			foreach ind in `allindicators' {
				// quietly summarize `ind' if `__touse__', detail
				quietly summarize `ind' , detail
				matrix `indstats'[`k', 1] = r(mean)
				matrix `indstats'[`k', 2] = r(sd)
				matrix `indstats'[`k', 3] = r(p50)
				matrix `indstats'[`k', 4] = r(min)
				matrix `indstats'[`k', 5] = r(max)
				local nonmissN = r(N)
				// quietly misstable summarize `ind' if `__touse__'
				quietly misstable summarize `ind'
				if (missing(r(N_eq_dot))) {
					matrix `indstats'[`k', 6] = `nonmissN'
					matrix `indstats'[`k', 7] = 0
				}
				else {
					matrix `indstats'[`k', 6] = r(N_lt_dot)
					matrix `indstats'[`k', 7] = r(N_eq_dot)
				}
				local ++k
			}
			matrix rownames `indstats' = `allindicators'
			matrix colnames `indstats' = "mean" "sd" "median" "min" "max" "N" "missing"
		}
		local num_ind = rowsof(`indstats')
		mktable, matrix(`indstats') digits(`digits') firstcolname("Indicator") ///
			title("Table of summary statistics for indicator variables") ///
			firstcolwidth(14) colwidth(12) hlines(`num_ind') novlines stats
	}

	if ("`meastable'" == "") {
		tempname loadings loadings_d
		matrix `loadings' = e(loadings)
		local num_lv = colsof(`loadings')
		local allformative = e(formative)
		local num_lv_B : word count `allformative'
		local num_lv_A = `num_lv' - `num_lv_B'
		local num_ind = rowsof(`loadings')
		if (`num_lv_A' == 0) {
			matrix `loadings_d' = `loadings'
		}
		else {
			local num_ind_plus_3 = `num_ind' + 3
			matrix `loadings_d' = J(`num_ind_plus_3', `num_lv', .)
			matrix `loadings_d'[1, 1] = `loadings'
			matrix `loadings_d'[`num_ind' + 1, 1] = e(relcoef)
			matrix rownames `loadings_d' = `:rowfullnames `loadings'' "Cronbach" ///
				"DG" "rho_A"
			matrix colnames `loadings_d' = `:colfullnames `loadings''
		}
		local title_meas "Measurement model - Standardized loadings"
		mktable, matrix(`loadings_d') digits(`digits') firstcolname("") ///
			title(`title_meas') firstcolwidth(14) colwidth(14) ///
			hlines(`num_ind' `num_ind_plus_3') novlines
		
		if ("`loadpval'" != "") {
			if (`isboot') {
				tempname adj_meas loadings_se nindblock loadings_df loadings_pval
				matrix `adj_meas' = e(adj_meas)
				matrix `loadings_se' = e(loadings_se)
				mata: st_matrix("`nindblock'", colsum(st_matrix("`adj_meas'")))
				mata: st_matrix("`loadings_df'", ///
					st_matrix("`adj_meas'")*(st_numscalar("e(N)") - 1))
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
				matrix rownames `loadings_pval' = `:rowfullnames `loadings''
				matrix colnames `loadings_pval' = `:colfullnames `loadings''
				if (!`isboot') {
					local title_meas "Measurement model - Standardized loading p-values"
				}
				else {
					local title_meas "Measurement model - Standardized loading p-values (Bootstrap)"
				}
				mktable, matrix(`loadings_pval') digits(`digits') firstcolname("") ///
					title(`title_meas') firstcolwidth(14) colwidth(14) ///
					hlines(`num_ind') novlines
			}
			else {
				// when using PLSc p-values are shown only with bootstrap
			}
		}
	}

	if ("`discrimtable'" == "") {
		tempname sqcorr sqcorr_d
		matrix `sqcorr' = e(sqcorr)
		local num_lv_A = colsof(`sqcorr')
		local num_lv_A_plus_one = `num_lv_A' + 1
		matrix `sqcorr_d' = J(`num_lv_A_plus_one', `num_lv_A', .)
		matrix `sqcorr_d'[1, 1] = `sqcorr'
		matrix `sqcorr_d'[`num_lv_A_plus_one', 1] = e(ave)
		matrix rownames `sqcorr_d' = `:rowfullnames `sqcorr'' "AVE"
		matrix colnames `sqcorr_d' = `:colfullnames `sqcorr''
		mktable, matrix(`sqcorr_d') digits(`digits') firstcolname("") ///
			title("Discriminant validity - Squared interfactor correlation vs. Average variance extracted (AVE)") ///
			firstcolwidth(14) colwidth(14) hlines(`num_lv_A' `num_lv_A_plus_one') novlines
	}
	
	if ("`structural'" != "nostructural") {
		if ("`structtable'" == "") {
			tempname pathtab
			matrix `pathtab' = e(struct_table)
			local hline_path_2 = rowsof(`pathtab')
			local hline_path_1 = `hline_path_2' - 1
			if (!`isboot') {
				local title_st "Structural model - Standardized path coefficients"
			}
			else {
				local title_st "Structural model - Standardized path coefficients (Bootstrap)"
			}
			mktable, matrix(`pathtab') digits(`digits') firstcolname("Variable") ///
				title(`title_st') firstcolwidth(14) colwidth(14) ///
				hlines(`hline_path_1' `hline_path_2') novlines path binary(`binary') ///
				consistent
		}
	}

	if ("`corrind'" != "") {
		tempname loadings
		matrix `loadings' = e(loadings)
		local num_ind = rowsof(`loadings')
		quietly correlate `e(mvs)' if `__touse__'
		tempname C
		matrix `C' = r(C)
		forvalues i = 1/`num_ind' {
			local ip1 = `i' + 1
			forvalues j = `ip1'/`num_ind' {
				matrix `C'[`i', `j'] = .
			}
		}
		mktable_corr, matrix(`C') title("Correlation of indicators") cutoff(`cutoff')
	}

	if ("`corrlv'" != "") {
		tempname loadings
		matrix `loadings' = e(loadings)
		local num_lv = colsof(`loadings')
		tempname C
		matrix `C' = e(R)
		forvalues i = 1/`num_lv' {
			local ip1 = `i' + 1
			forvalues j = `ip1'/`num_lv' {
				matrix `C'[`i', `j'] = .
			}
		}
		mktable_corr, matrix(`C') title("Correlation of latent variables") ///
			cutoff(`cutoff')
	}

	if ("`crossload'" != "") {
		tempname xload
		matrix `xload' = e(cross_loadings)
		local num_ind = rowsof(`xload')
		matrix rownames `xload' = `e(mvs)'
		matrix colnames `xload' = `e(lvs)'
		matrix coleq `xload' = ""
		mktable, matrix(`xload') digits(4) firstcolname("") ///
			title("Cross loadings") firstcolwidth(12) colwidth(9) ///
			hlines(`num_ind') novlines corr cutoff(`cutoff')
	}
end
