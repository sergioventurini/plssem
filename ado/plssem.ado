*!plssem version 0.3.0
*!Written 04Sep2017
*!Written by Sergio Venturini and Mehmet Mehmetoglu
*!The following code is distributed under GNU General Public License version 3 (GPL-3)

program plssem, byable(onecall)
	version 14.2
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

	local isgroup = strpos(`"`options'"', "gr")
	if (`isgroup') {
		if _by() {
			display as error "the 'group()' option is not allowed with by"
			exit
		}
		else {
			`version' Compare `0'  // version is not necessary
		}
	}
	else {
		`version' `BY' Estimate `0'  // version is not necessary
	}
end

program Estimate, eclass byable(recall)
	version 14.2
	syntax anything(id="Measurement model" name=blocks) [if] [in], ///
		[ STRuctural(string) Wscheme(string) BINary(namelist min=1) ///
		Boot(numlist integer >0 max=1) Seed(numlist max=1) Tol(real 1e-7) ///
		MAXiter(integer 100) MISsing(string) k(numlist integer >0 max=1) ///
		INIT(string) DIGits(integer 3) noHEADer noMEAStable noDISCRIMtable ///
		noSTRUCTtable STATs GRoup(string) CORRelate(string) RAWsum noSCale ///
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
		 stats													--> print a table of summary statistics for
																				the indicators
		 group(string)									--> perform multigroup analysis; accepts
																				the suboptions reps(#), method(), plot,
																				alpha, what() and groupseed(numlist max=1);
																				method() accepts either 'permutation',
																				'bootstrap' or 'normal'; groupseed(numlist
																				max=1) accepts an optional seed for
																				multigroup analysis; what() accepts either
																				'path' or 'loadings'
		 correlate(string)							--> report the correlation among the
																				indicators, the LVs as well as the cross
																				loadings; accepts any combination of 'mv',
																				'lv' and 'cross' (i.e. cross-loadings);
																				suboption 'cutoff()' can be used to avoid
																				showing correlations smaller than the
																				specified value
		 rawsum													--> estimate the latent scores as the raw
																				sum of the indicators
		 noscale												--> manifest variables are not scaled but
																				only centered
		 convcrit												--> convergence criterion (either 'relative'
																				or 'square')
		 nocleanup											--> Mata temporary objects are not removed
																				(undocumented)
	 */
	
	/* Parse the specified blocks.
	 * Macros:
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
			mata: `ind_mata' = st_data(., "`allindicators'")
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
		exit
	}
	if (("`structural'" == "") & ("`init'" != "eigen")) {
		display as error "'init()' option must be set to 'eigen' when no structural model is provided"
		exit
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
			st_numscalar("`rawsum_sc'"))
	
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
		mata: `xload_v' = ///
			plssem_lv( ///
				st_data(., "`allstdindicators'", "`touse'"), ///
				st_data(., "`alllatents'", "`touse'"), ///
				"`alllatents'", ///
				"`binary'")
		if ("`structural'" != "") {
			mata: `path_v' = `res'.path_v
		}
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
				strtoreal("`boot'"), ///
				strtoreal("`seed'"), ///
				1)
		mata: `xload_v' = `res_bs'.xloadings_v
		if ("`structural'" != "") {
			mata: `path_v' = `res_bs'.path_v
		}
	}
	/* End of computing the model parameters' variances */

	/* Compute the table of measurement loadings */
	tempname loadings xloadings annihilate_se loadings_se xloadings_se

	mata: st_matrix("`loadings'", `res'.loadings)
	matrix rownames `loadings' = `loadrownames'
	matrix colnames `loadings' = `loadcolnames'

	mata: st_matrix("`xloadings'", `res'.xloadings)
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
		mata: st_matrix("`sqcorr'", ///
			correlation(st_data(., "`alllatents'", "`touse'")) :^ 2)
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
			quietly summarize `ind' if `touse', detail
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
		
		mata: st_matrix("`path'", `res'.path)
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

		mata: st_matrix("`rsquared'", `res'.r2)
		matrix rownames `rsquared' = "r2"
		matrix colnames `rsquared' = `alllatents'
		
		mata: st_matrix("`redundancy'", `res'.r2 :* st_matrix("`ave'"))
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
			plssem_pathtab(st_matrix("`path'"), `path_pval', `res'.r2_a))
		foreach var in `regest_all' {
			local pathtab_rownames "`pathtab_rownames' `var' ."
		}
		matrix rownames `pathtab' = `pathtab_rownames' "r2_a"
		matrix colnames `pathtab' = `lv_regest_all'

		mata: st_matrix("`path_toteff'", `res'.total_effects)
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
			st_matrix("`ave'"), ///
			st_matrix("`modes'")))
	matrix rownames `relcoef' = "Cronbach" "DG"
	matrix colnames `relcoef' = `loadcolnames'
	/* End of computing the reliability coefficients */
	
	/* Setup the tables of structural coefficients */
	if ("`structural'" != "") {
		tempname path_mata strcoef path_se_mata strse path_vif_mata strvif
		
		mata: `path_mata' = st_matrix("`path'")
		mata: st_matrix("`strcoef'", ///
			`path_mata'[selectindex(rownonmissing(`path_mata')), ///
				selectindex(colnonmissing(`path_mata'))])
		foreach var in `regest_all' {
			local struct_rownames "`struct_rownames' `var'"
		}
		matrix rownames `strcoef' = `struct_rownames'
		matrix colnames `strcoef' = `lv_regest_all'
		
		mata: `path_se_mata' = st_matrix("`path_se'")
		mata: st_matrix("`strse'", ///
			`path_se_mata'[selectindex(rownonmissing(`path_se_mata')), ///
				selectindex(colnonmissing(`path_se_mata'))])
		matrix rownames `strse' = `struct_rownames'
		matrix colnames `strse' = `lv_regest_all'
		
		mata: `path_vif_mata' = ///
			plssem_vif( ///
				st_data(., "`alllatents'", "`touse'"), ///
				st_matrix("`adj_struct'"))
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
	ereturn local title "Partial least squares structural equation modeling"
	ereturn local predict plssem_p
	ereturn local estat_cmd "plssem_estat"
	ereturn local cmdline "plssem `cmdline'"
	ereturn local cmd "plssem"

	if ("`structural'" != "") {
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
	if (("`display'" == "") & ("`group'" == "")) {  // this saves some time when used with 'group' 
		if ("`structural'" == "") {
			local nostructural "nostructural"
		}
		if (`num_lv_A' == 0) {
			local nodiscrim "nodiscrimtable"
		}
		else {
			local nodiscrim `discrimtable'
		}
		Display, `nostructural' digits(`digits') boot(`boot') ///
			`header' `meastable' `nodiscrim' `structtable' `stats' `corrind' ///
			`corrlv' `crossload' cutoff(`cutoff') binary(`binary') `rawsum'
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
	if ("`cleanup'" == "") {
		capture mata: cleanup()
	}
	/* End of cleaning up */
	}
end

program Compare, eclass sortpreserve
	version 14.2
	syntax anything(id="Measurement model" name=blocks) [if] [in], ///
		[ STRuctural(string) Wscheme(string) BINary(namelist min=1) ///
		Boot(numlist integer >0 max=1) Seed(numlist max=1) Tol(real 1e-7) ///
		MAXiter(integer 100) MISsing(string) k(numlist integer >0 max=1) ///
		INIT(string) DIGits(integer 3) noHEADer noMEAStable noDISCRIMtable ///
		noSTRUCTtable STATs GRoup(string) CORRelate(string) RAWsum noSCale ///
		CONVcrit(string) noCLeanup ]

	/* Options:
	   --------
		 structural(string)							--> structural model specification
		 wscheme(string)								--> weighting scheme ("centroid", "factorial"
																				or "path")
		 binary(namelist min=1)					--> indicators to fit using logit
		 boot(numlist integer >0 max=1)	--> bootstrap estimation (# of repetions)
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
		 stats													--> print a table of summary statistics for
																				the indicators
		 group(string)									--> performs multigroup analysis; accepts
																				the suboptions reps(#), method(), plot,
																				alpha, what() and groupseed(numlist max=1);
																				method() accepts either 'permutation',
																				'bootstrap' or 'normal'; groupseed(numlist
																				max=1) accepts an optional seed for
																				multigroup analysis; what() accepts either
																				'path' or 'loadings'
		 correlate(string)							--> reports the correlation among the
																				indicators, the LVs as well as the cross
																				loadings; accepts any combination of 'mv',
																				'lv' and 'cross' (i.e. cross-loadings);
																				suboption 'cutoff()' can be used to avoid
																				showing correlations smaller than the
																				specified value
		 rawsum													--> estimates the latent scores as the raw
																				sum of the indicators
		 noscale												--> manifest variables are not scaled but
																				only centered
		 convcrit												--> convergence criterion (either 'relative'
																				or 'square')
		 nocleanup											--> Mata temporary objects are not removed
																				(undocumented)
	 */
	
	/* Warning */
	if ("`boot'" != "") {
		display as error "using the 'boot()' option together with 'group()' will slow down excessively the calculation"
		display as error "try removing 'boot()'"
		exit
	}
	
	/* Check that digits is nonnegative */
	if (`digits' < 0) {
		display as error "number of digits to display must be a nonnegative integer"
		exit
	}
	
	/* Set macros */
	if ("`wscheme'" == "") {
		local wscheme = "path"
	}
	if ("`reps'" == "") {
		local reps = 100
	}
	if ("`method'" == "") {
		local method = "normal"
	}
	if ("`alpha'" == "") {
		local alpha = 0.05
	}
	if ("`what'" == "") {
		local what = "path"
	}
	if (("`structural'" == "") & ("`what'" == "path")) {
		display as error "choosing 'what(path)' requires a structural part"
		exit
	}
	if ("`init'" == "") {
		local init "indsum"
	}
	if ("`convcrit'" == "") {
		local convcrit "relative"
	}
	local skip1 = 1
	local skip3 = 3

	/* Parse the group() option */
	tokenize `"`group'"', parse(",")
	local ngroupvars : word count `1'
	if (`ngroupvars' > 1) {
		display as error "a single variable allowed in the 'group()' option"
		exit
	}
	local groupvar `"`1'"'
	local rest `"`3'"'
	local tok_i = 1
	tokenize `"`rest'"', parse("() ")
	while ("``tok_i''" != "") {
		if (("``tok_i''" != "(") | ("``tok_i''" != ")")) {
			if ("``tok_i''" == "reps") {
				local tok_tmp = `tok_i' + 2
				local reps = int(``tok_tmp'')
				if (`reps' < 2) {
					display as error "'reps()' suboption requires an integer larger than 1"
					exit
				}
			}
			else if ("``tok_i''" == "method") {
				local tok_tmp = `tok_i' + 2
				local method = "``tok_tmp''"
				if (("`method'" != "normal") & ("`method'" != "permutation") & ("`method'" != "bootstrap")) {
					display as error "'method()' suboption accepts either 'permutation', 'bootstrap' or 'normal'"
					exit
				}
			}
			else if ("``tok_i''" == "plot") {
				local plot = "``tok_i''"
			}
			else if ("``tok_i''" == "alpha") {
				local tok_tmp = `tok_i' + 2
				local alpha = ``tok_tmp''
				if (`alpha' <= 0 | `alpha' >= 1) {
					display as error "'alpha()' suboption requires a real scalar in between 0 and 1"
					exit
				}
			}
			else if ("``tok_i''" == "groupseed") {
				local tok_tmp = `tok_i' + 2
				local groupseed = ``tok_tmp''
				if (`groupseed' < 0 | `groupseed' >  2^31-1) {
					display as error "'groupseed()' suboption requires a value between 0 and 2^31-1"
					exit
				}
			}
			else if ("``tok_i''" == "what") {
				local tok_tmp = `tok_i' + 2
				local what = "``tok_tmp''"
				if (("`what'" != "path") & ("`what'" != "loadings")) {
					display as error "'what()' suboption accepts either 'path' or 'loadings'"
					exit
				}
			}
		}
		local ++tok_i
	}
	/* End of parsing the group() option */
	
	/* Estimation using the whole sample */
	tempname results whole adj_meas adj_struct modes rawsum_sc struct_sc ehold tmp
	quietly Estimate `0'
	local alllatents = e(lvs)
	local allindicators = e(mvs)
	local allreflective = e(reflective)
	local nlv : word count `alllatents'
	if ("`what'" == "path") {
		matrix `tmp' = e(struct_b)
		local loadcolnames : colnames `tmp'
		matrix `whole' = vec(`tmp')
	}
	else if ("`what'" == "loadings") {
		matrix `tmp' = e(loadings)
		local loadrownames : rownames `tmp'
		local loadcolnames : colnames `tmp'
		matrix `whole' = vec(`tmp')
	}
	matrix `adj_meas' = e(adj_meas)
	if ("`structural'" != "") {
		matrix `adj_struct' = e(adj_struct)
	}
	else {
		matrix `adj_struct' = J(`nlv', `nlv', 0)
	}
	matrix `modes' = J(`nlv', 1, 1)
	local i = 1
	foreach var in `alllatents' {
		if (`: list var in allreflective') {
			matrix `modes'[`i', 1] = 0
		}
		local ++i
	}
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
	
	/* Set obs to use */
	tempvar touse
	quietly generate byte `touse' = e(sample)
	quietly count if `touse'
	if (r(N) == 0) {
		error 2000
	}
	else if (r(N) < 10) {
		display as error "at least 10 complete observations required"
		error 2001
	}

	_estimates hold `ehold', restore
	
	/* Compute group values and sizes */
	tempname grpvar_tmp groupvar_m groupvals groupsizes alln
	mata: `grpvar_tmp' = st_data(., "`groupvar'", "`touse'")
	mata: st_matrix("`groupvar_m'", ///
		uniqrows(select(`grpvar_tmp', `grpvar_tmp' :!= .), 1))
	local ngroups = rowsof(`groupvar_m')
	if (`ngroups' == 1) {
		display as error "the group() option requires at least two groups to compare"
		exit
	}
	matrix `groupvals' = `groupvar_m'[1..., 1]
	matrix `groupsizes' = `groupvar_m'[1..., 2]
	mata: st_numscalar("`alln'", sum(st_matrix("`groupsizes'")))
	
	/* Calculate observed test statistic value */
	tempname idx res
	mata: `idx' = selectindex(rowsum(st_matrix("`adj_struct'")))
	mata: st_local("regest_all", invtokens(tokens("`alllatents'")[`idx']))
	mata: `idx' = selectindex(colsum(st_matrix("`adj_struct'")))
	mata: st_local("lv_regest_all", invtokens(tokens("`alllatents'")[`idx']))
	foreach var in `regest_all' {
		local struct_rownames "`struct_rownames' `var'"
	}
	local struct_rownames : list clean struct_rownames
	foreach var in `allindicators' {
		local allstdindicators "`allstdindicators' std`var'"
	}
	local allstdindicators : list clean allstdindicators
	
	tempname zerovar tmp_mata
	forvalues ng = 1/`ngroups' {
		preserve
		
		tempvar touse_gr
		quietly generate `touse_gr' = `touse'
		quietly replace `touse_gr' = 0 if (`groupvar' != `groupvals'[`ng', 1])
		
		/* Check that there are no "zero-variance" indicators in any group */
		tempvar touse_gr
		quietly generate `touse_gr' = `touse'
		quietly replace `touse_gr' = 0 if (`groupvar' != `groupvals'[`ng', 1])
		
		mata: st_numscalar("`zerovar'", ///
			any(selectindex(sd(st_data(., "`allindicators'", "`touse_gr'")) :== 0)))
		if (`zerovar') {
			continue, break
		}
		
		/* Standardize the MVs (if requested) */
		mata: ///
			plssem_scale( ///
				st_data(., "`allindicators'", "`touse_gr'"), ///
				"`allstdindicators'", ///
				"`touse_gr'", ///
				"`scale'")
		
		/* Initialize LVs */
		mata: ///
			plssem_init( ///
				st_data(., "`allstdindicators'", "`touse_gr'"), ///
				st_matrix("`adj_meas'"), ///
				"`allindicators'", ///
				"`allstdindicators'", ///
				"`alllatents'", ///
				"`touse_gr'", ///
				st_numscalar("`rawsum_sc'"), ///
				"`init'")
		
		/* Run the PLS algorithm */
		mata: `res' = ///
			plssem_base( ///
				st_data(., "`allstdindicators'", "`touse_gr'"), ///
				st_data(., "`alllatents'", "`touse_gr'"), ///
				st_matrix("`adj_meas'"), ///
				st_matrix("`adj_struct'"), ///
				st_matrix("`modes'"), ///
				"`alllatents'", ///
				"`binary'", ///
				strtoreal("`tol'"), ///
				strtoreal("`maxiter'"), ///
				"`touse_gr'", ///
				"`wscheme'", ///
				"`convcrit'", ///
				st_numscalar("`struct_sc'"), ///
				st_numscalar("`rawsum_sc'"))
		
		tempname strb_`ng' strvar_`ng' group_`ng'
		if ("`what'" == "path") {
			mata: `tmp_mata' = `res'.path
			mata: `tmp_mata' = `tmp_mata'[selectindex(rownonmissing(`tmp_mata')), ///
				selectindex(colnonmissing(`tmp_mata'))]
			mata: st_matrix("`strb_`ng''", `tmp_mata')
			matrix rownames `strb_`ng'' = `struct_rownames'
			matrix colnames `strb_`ng'' = `lv_regest_all'
			matrix `strb_`ng'' = vec(`strb_`ng'')
			
			mata: `tmp_mata' = sqrt(`res'.path_v)
			mata: `tmp_mata' = `tmp_mata'[selectindex(rownonmissing(`tmp_mata')), ///
				selectindex(colnonmissing(`tmp_mata'))]
			mata: st_matrix("`strvar_`ng''", `tmp_mata')
			matrix rownames `strvar_`ng'' = `struct_rownames'
			matrix colnames `strvar_`ng'' = `lv_regest_all'
			matrix `strvar_`ng'' = vec(`strvar_`ng'')
		}
		else if ("`what'" == "loadings") {
			mata: st_matrix("`strb_`ng''", `res'.loadings)
			matrix rownames `strb_`ng'' = `loadrownames'
			matrix colnames `strb_`ng'' = `loadcolnames'
			matrix `strb_`ng'' = vec(`strb_`ng'')

			mata: `tmp_mata' = ///
				plssem_lv( ///
					st_data(., "`allstdindicators'", "`touse_gr'"), ///
					st_data(., "`alllatents'", "`touse_gr'"), ///
					"`alllatents'", ///
					"`binary'")
			mata: `tmp_mata' = sqrt(`tmp_mata') :* editvalue(`res'.M, 0, .)
			mata: st_matrix("`strvar_`ng''", `tmp_mata')
			matrix rownames `strvar_`ng'' = `loadrownames'
			matrix colnames `strvar_`ng'' = `loadcolnames'
			matrix `strvar_`ng'' = vec(`strvar_`ng'')
		}
		
		restore
	}
	if (`zerovar') {
		restore
		_estimates unhold `ehold'
		display as error "at least one indicator has zero variance in one of the groups"
		exit
	}
	
	local nrows = rowsof(`tmp')
	local ncols = colsof(`tmp')
	local totrows = `nrows'*`ncols'
	tempname nonmiss
	matrix `nonmiss' = J(`totrows', 1, 0)
	local strb_rn : rowfullnames `strb_1'
	forvalues i = 1/`totrows' {
		if (!missing(`strb_1'[`i', 1])) {
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
			if ("`what'" == "path") {
				local strbok_rn "`strbok_rn' `nm_X':`nm_Y'"
			}
			else if ("`what'" == "loadings") {
				local lv_is_ref : list nm_Y in allreflective
				if (`lv_is_ref') {
					local strbok_rn "`strbok_rn' `nm_Y':`nm_X'"
				}
				else {
					local strbok_rn "`strbok_rn' `nm_X':`nm_Y'"
				}
			}
		}
	}
	forvalues ng = 1/`ngroups' {
		tempname strbok_`ng' strvarok_`ng' strb_var_`ng'
		mata: st_matrix("`strbok_`ng''", select(st_matrix("`strb_`ng''"), ///
			st_matrix("`nonmiss'")))
		mata: st_matrix("`strvarok_`ng''", select(st_matrix("`strvar_`ng''"), ///
			st_matrix("`nonmiss'")))
		mata: st_matrix("`strvarok_`ng''", st_matrix("`strvarok_`ng''"):^2)
		matrix `group_`ng'' = `strbok_`ng''
	}
	tempname numeff obstest_tmp obstest
	mata: st_numscalar("`numeff'", colsum(st_matrix("`nonmiss'")))
	local neff = `numeff'
	matrix `obstest' = J(`neff', `ngroups' - 1, .)
	forvalues ng = 2/`ngroups' {
		mata: st_matrix("`obstest_tmp'", abs(st_matrix("`strbok_1'") - ///
			st_matrix("`strbok_`ng''")))
		matrix `obstest'[1, `ng' - 1] = `obstest_tmp'
		local obstest_cn "`obstest_cn' `ng'vs1"
	}
	matrix rownames `obstest' = `strbok_rn'
	matrix colnames `obstest' = `obstest_cn'
	
	/* Standardize the MVs (if requested) */
	mata: ///
		plssem_scale( ///
			st_data(., "`allindicators'", "`touse'"), ///
			"`allstdindicators'", ///
			"`touse'", ///
			"`scale'")
	
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
	
	tempname dtest
	matrix `dtest' = J(`ngroups' - 1, `neff', .)
	if ("`method'" != "normal") {
		tempname res_mga

		capture noisily {
			if ("`method'" == "permutation") {
				local title "Multigroup comparison (`groupvar') - Permutation test"
				
				/* Compute the permutation distribution */
				mata: `res_mga' = ///
					plssem_mga_perm( ///
						st_data(., "`allstdindicators'"), ///		 note: `touse' not used here
						st_data(., "`alllatents'"), ///					 note: `touse' not used here
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
						st_data(., "`groupvar'"), ///						 note: `touse' not used here
						strtoreal("`reps'"), ///
						strtoreal("`groupseed'"), ///
						1)
				
				mata: st_matrix("`dtest'", plssem_mga_perm_diff(`res_mga', ///
					st_matrix("`obstest'"), "`what'"))
			}
			else if ("`method'" == "bootstrap") {
				local title "Multigroup comparison (`groupvar') - Bootstrap t-test"
				
				/* Compute the bootstrap distribution */
				mata: `res_mga' = ///
					plssem_mga_boot( ///
						st_data(., "`allstdindicators'"), ///		 note: `touse' not used here
						st_data(., "`alllatents'"), ///					 note: `touse' not used here
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
						st_data(., "`groupvar'"), ///						 note: `touse' not used here
						strtoreal("`reps'"), ///
						strtoreal("`groupseed'"), ///
						1)
				
				mata: st_matrix("`dtest'", plssem_mga_boot_diff(`res_mga', ///
					st_matrix("`groupsizes'"), strtoreal("`neff'"), "`what'"))
			}
		} // end of -capture-
		local rc = _rc
		_estimates unhold `ehold'
		if (`rc' == 1) {
			display
			display as error "you pressed the Break key;" _continue
			display as error "calculation interrupted"
			error `rc'
		}
		else if (`rc' > 1) {
			error `rc'
		}
	}
	else {
		local title "Multigroup comparison (`groupvar') - Normal-based t-test"
		
		tempname k0 k1 k2 k3
		scalar `k0' = `alln' - 2
		scalar `k1' = ((`groupsizes'[1, 1] - 1)^2)/`k0'
		forvalues ng = 2/`ngroups' {
			scalar `k2' = ((`groupsizes'[`ng', 1] - 1)^2)/`k0'
			scalar `k3' = sqrt(1/`groupsizes'[1, 1] + 1/`groupsizes'[`ng', 1])
			forvalues j = 1/`neff' {
				matrix `dtest'[`ng' - 1, `j'] = `obstest'[`j', `ng' - 1]/(sqrt( ///
					`k1'*`strvarok_1'[`j', 1] + `k2'*`strvarok_`ng''[`j', 1])*`k3')
			}
		}
	}
	matrix rownames `dtest' = `obstest_cn'
	matrix colnames `dtest' = `strbok_rn'
	
	matrix `results' = J(`neff', 1 + `ngroups' + (`ngroups' - 1)*3, .)
	local resnm : subinstr local strbok_rn ":" "->", all
	matrix rownames `results' = `resnm'
	if (`ngroups' > 2) {
		forvalues ng = 1/`ngroups' {
			local grp_cn `grp_cn' "Group_`ng'"
			if (`ng' > 1) {
				local absd_cn `absd_cn' "AD_`ng'vs1"
				local stat_cn `stat_cn' "S_`ng'vs1"
				local pval_cn `pval_cn' "P_`ng'vs1"
			}
		}
		matrix colnames `results' = "Global" `grp_cn' `absd_cn' `stat_cn' `pval_cn'
	}
	else {
		matrix colnames `results' = "Global" "Group_1" "Group_2" "Abs_Diff" "Statistic" "P-value"
	}

	tempname alllen anysig
	matrix `alllen' = J(`neff', 1, .)
	mata: st_matrix("`whole'", select(st_matrix("`whole'"), ///
		st_matrix("`nonmiss'")))
	local anysig = 0
	forvalues j = 1/`neff' {
		matrix `results'[`j', 1] = `whole'[`j', 1]
		forvalues ng = 1/`ngroups' {
			matrix `results'[`j', 1 + `ng'] = `group_`ng''[`j', 1]
			if (`ng' > 1) {
				matrix `results'[`j', 1 + `ngroups' + (`ng' - 1)] = ///
					`obstest'[`j', `ng' - 1]
				matrix `results'[`j', 1 + `ngroups' + (`ngroups' - 1) + (`ng' - 1)] = ///
					`dtest'[`ng' - 1, `j']
				if ("`method'" == "permutation") {
					matrix `results'[`j', 1 + `ngroups' + (`ngroups' - 1)*2 + (`ng' - 1)] = ///
						(`dtest'[`ng' - 1, `j'] + 1)/(`reps' + 1)
				}
				else if (("`method'" == "bootstrap") | ("`method'" == "normal")) {
					matrix `results'[`j', 1 + `ngroups' + (`ngroups' - 1)*2 + (`ng' - 1)] = ///
						2*ttail(`alln' - 2, `dtest'[`ng' - 1, `j'])
				}
			}
		}
		local nm : word `j' of `resnm'
		local nm_len : strlen local nm
		matrix `alllen'[`j', 1] = `nm_len' + 2
		local nm_new = subinstr("`nm'", "->", " -> ", 1)
		forvalues ng = 2/`ngroups' {
			if (`results'[`j', 1 + `ngroups' + (`ngroups' - 1)*2 + (`ng' - 1)] <= `alpha') {
				local nm_new "(*) `nm_new'"
				local anysig = 1
				continue, break
			}
		}
		local lblbar "`lblbar' `j' "`nm_new'""
	}
	
	if ("`rawsum'" == "") {
		if ("`structural'" != "") {
			mkheader, digits(5)
		}
		else {
			mkheader, digits(5) nostructural
		}
	}
	else {
		if ("`structural'" != "") {
			mkheader, digits(5) rawsum
		}
		else {
			mkheader, digits(5) nostructural rawsum
		}
	}
	
	if ("`what'" == "path") {
		local firstcollbl "Structural effect"
	}
	else {
		local firstcollbl "Measurement effect"
	}
	tempname maxlen
	mata: st_numscalar("`maxlen'", max(st_matrix("`alllen'")))
	local `maxlen' = max(strlen("`firstcollbl'"), `maxlen') + 2
	if (`ngroups' == 2) {
		local colw = 11
	}
	else {
		local colw = 9
	}
	mktable, matrix(`results') digits(`digits') firstcolname(`firstcollbl') ///
		title(`title') firstcolwidth(``maxlen'') colwidth(`colw') hlines(`neff') ///
		novlines total
	if ("`method'" != "normal") {
		display as text _skip(`skip1') "number of replications: `reps'"
	}
	if (`ngroups' > 2) {
		display as text _skip(`skip1') "legend:"
		display as text _skip(`skip3') "AD: absolute difference"
		display as text _skip(`skip3') "S: statistic"
		display as text _skip(`skip3') "P: p-value"
	}
	display as text _skip(`skip1') "group labels:"
	local grplbl : value label `groupvar'
	forvalues ng = 1/`ngroups' {
		local grpvarlbl`ng' = `groupvals'[`ng', 1]
		if ("`grplbl'" != "") {
			local grpvarlbl`ng' : label `grplbl' `grpvarlbl`ng''
		}
		display as text _skip(`skip3') "Group `ng': `grpvarlbl`ng''"
	}
	
	if ("`plot'" != "") {
		local reslbl "__result__"
		quietly svmat `results', names(`reslbl')
		generate id = _n
		local grplbl : value label `groupvar'
		forvalues ng = 1/`ngroups' {
			local grpvarlbl`ng' = `groupvals'[`ng', 1]
			if ("`grplbl'" != "") {
				local grpvarlbl`ng' : label `grplbl' `grpvarlbl`ng''
			}
			local ngp1 = `ng' +  1
			local reslbl_toplot "`reslbl_toplot' `reslbl'`ngp1'"
			local res_tolbl "`res_tolbl' label(`ng' `grpvarlbl`ng'')"
		}
		if ("`what'" == "path") {
			local whatplot "Path Coefficients"
		}
		else if ("`what'" == "loadings") {
			local whatplot "Loadings"
		}
		local bartitle "`whatplot' Comparison across Groups (`groupvar')"
		if ("`method'" != "normal") {
			local barnote "Method: `method' - number of replications: `reps'"
		}
		else {
			local barnote "Method: `method'"
		}
		if (`anysig') {
			local alpha_perc = string(`alpha'*100, "%5.0g")
			if (`ngroups' == 2) {
				local barcaption "(*) Difference significant at {&alpha} = `alpha_perc'%"
			}
			else {
				local barcaption "(*) At least one of the differences is significant at {&alpha} = `alpha_perc'%"
			}
		}
		else {
			local barcaption ""
		}

		graph bar (asis) `reslbl_toplot', over(id, relabel(`lblbar') ///
			label(angle(90) labsize(small))) nofill legend(pos(12) ///
			region(style(none)) `res_tolbl' rows(1)) ///
			ylabel(, nogrid) yline(0, lwidth(vvthin)) title(`bartitle') scheme(sj) ///
			note(`barnote') caption(`barcaption', size(small))
		capture quietly drop `reslbl'* id
	}
	
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
	if ("`cleanup'" == "") {
		capture mata: cleanup()
	}

	/* Return */
	// ereturn local method = "`method'"
	// ereturn matrix comparison = `results'
end

program Display
	version 14.2
	syntax [, noSTRuctural DIGits(integer 3) noHEADer noMEAStable ///
		noDISCRIMtable noSTRUCTtable stats corrind corrlv crossload ///
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
				mkheader, matrix1(`mod_assessment') matrix2(`reldiff') digits(5) nogroup
			}
			else {
				tempname mod_assessment
				matrix `mod_assessment' = e(assessment)
				mkheader, matrix1(`mod_assessment') digits(5) nogroup rawsum
			}
		}
		else {
			mkheader, digits(5) nogroup nostructural
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
			local num_ind_plus_2 = `num_ind' + 2
			matrix `loadings_d' = J(`num_ind_plus_2', `num_lv', .)
			matrix `loadings_d'[1, 1] = `loadings'
			matrix `loadings_d'[`num_ind' + 1, 1] = e(relcoef)
			matrix rownames `loadings_d' = `:rowfullnames `loadings'' "Cronbach" "DG"
			matrix colnames `loadings_d' = `:colfullnames `loadings''
		}
		//if (!`isboot') {
			local title_meas "Measurement model - Standardized loadings"
		//}
		//else {
		//	local title_meas "Measurement model - Standardized loadings (Bootstrap)"
		//}
		mktable, matrix(`loadings_d') digits(`digits') firstcolname("") ///
			title(`title_meas') firstcolwidth(14) colwidth(14) ///
			hlines(`num_ind' `num_ind_plus_2') novlines
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
				hlines(`hline_path_1' `hline_path_2') novlines path binary(`binary')
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
		quietly correlate `e(lvs)' if `__touse__'
		tempname C
		matrix `C' = r(C)
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
		tempname loadings
		matrix `loadings' = e(loadings)
		local num_ind = rowsof(`loadings')
		local num_lv = colsof(`loadings')
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
			hlines(`num_ind') novlines corr cutoff(`cutoff')
	}
end
