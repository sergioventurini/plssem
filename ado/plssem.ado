*!plssem version 0.1
*!Written 06Jan2017
*!Written by Sergio Venturini and Mehmet Mehmetoglu
*!The following code is distributed under GNU General Public License version 3 (GPL-3)

program plssem, byable(onecall)
	version 10
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
			error 198
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
	version 10
	syntax anything(id="Indicator blocks" name=blocks) [if] [in], ///
		[ STRuctural(string) Wscheme(string) BINary(namelist min=1) ///
		Boot(numlist integer >0 max=1) Seed(numlist max=1) Tol(real 1e-7) ///
		MAXiter(integer 1000) INIT(string) DIGits(integer 3) noHEADer ///
		noMEAStable noDISCRIMtable noSTRUCTtable STATs GRoup(string) ///
		CORRelate(string) RAWsum ]
	
	/* Options:
	   --------
		 structural(string)							--> structural model specification
		 wscheme(string)								--> weighting scheme ("centroid", "factorial"
																				or "path")
		 binary(namelist min=1)					--> latent variables to fit using logit
		 boot(numlist integer >0 max=1)	--> bootstrap estimation (# of replications)
		 seed(numlist max=1)						--> bootstrap seed number
		 tol(real 1e-7)									--> tolerance (default 1e-7)
		 maxiter(integer 1000)					--> maximum number of iterations (default
																				1000)
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
	 */
	
	/* Parse the specified blocks.
	 * Macros:
	 * 	LV`i'						- latent variable name
	 * 	i`i'						- indicators for latent, by latent index
	 * 	istd`i'					- indicators for latent, by latent index (standardized)
	 * 	i`LV`i''				- indicators for latent, by latent name
	 * 	istd`LV`i''			- indicators for latent, by latent name (standardized)
	 *	allindicators		- all indicators
	 *	alllatents			- all latents
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
						quietly drop ``tok_i''
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

	/* Set obs to use */
	marksample touse
	markout `touse' `allindicators'

	quietly count if `touse'
	if (r(N) == 0) {
		error 2000
	}

	quietly misstable patterns `allindicators' if `touse'
	local nobs = r(N_complete)

	/* Parse the inner weight scheme */	
	if ("`wscheme'" == "") {
		local wscheme "path"
	}
	if !("`wscheme'" == "centroid" | "`wscheme'" == "factorial"  | "`wscheme'" == "path") {
		display as error "scheme can be either 'centroid', 'factorial', or 'path'"
		error 198
	}
	/* End of parsing the inner weight scheme */

	/* Parse the binary indicators */
	if ("`binary'" != "") {
		foreach var in `binary' {
			local nind_bin : word count `i`var''
			if (`nind_bin' > 1) {
				display as error "latent variables in 'binary()' option must be defined through a single binary indicator"
				error 198
			}
		}
	}
	/* End of parsing the binary indicators */
	
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
					error 198
				}
			}
			local ++tok_i
		}
		if ("`corrind'" == "" & "`corrlv'" == "" & "`crossload'" == "") {
			display as error "'correlate()' option requires at least one of 'mv', 'lv' or 'cross'"
			error 198
		}
	}
	/* End of parsing the correlate() option */

	/* Standardize the indicators */
	foreach var in `allindicators' {
		// quietly tabulate `var' if `touse'
		// if (r(r) > 2) {
			capture confirm new variable std`var'
			if (_rc == 0) {
				quietly summarize `var' if `touse'
				quietly generate std`var' = (`var' - r(mean))/r(sd) if `touse'
			}
		// }
		// else {
		// 	capture confirm new variable std`var'
		// 	if (_rc == 0) {
		// 		quietly generate std`var' = `var' if `touse'  // do not standardize binary indicators
		// 	}
		// }
		if (_rc == 0) {
			local allstdindicators "`allstdindicators' std`var'"
		}
	}
	local num_lv: word count `alllatents'
	forvalues k = 1/`num_lv' {
		foreach var in `i`k'' {
			local istd`k' "`istd`k'' std`var'"
			local istd`LV`k'' `istd`k''
		}
	}
	
	/* Initialize latents with equal weights (STEP 1) */
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
	
	if ("`rawsum'" == "") {
		if ("`init'" == "indsum") {
			foreach k of numlist 1/`num_lv' {
				tokenize `"`istd`k''"'
				if (_byindex() > 1) {  // this provides all the LVs predictions when by'd
					quietly replace `LV`k'' = cond(`1' >= ., 0, `1') if `touse' // equation (4)
				}
				else {
					quietly generate `LV`k'' = cond(`1' >= ., 0, `1') if `touse' // equation (4)
				}
				macro shift
				while ("`1'" != "") {
					quietly replace `LV`k'' = `LV`k'' + cond(`1' >= ., 0, `1') if `touse'
					macro shift
				}
				quietly summarize `LV`k'' if `touse'
				local lv : word `k' of `alllatents'
				local isbinary : list lv in binary
				if (!`isbinary') {
					quietly replace `LV`k'' = (`LV`k'' - r(mean))/r(sd) if `touse' // equation (5)
				}
			}
		}
		else if ("`init'" == "eigen") {
			tempvar tmpscores
			foreach k of numlist 1/`num_lv' {
				if (_byindex() > 1) {  // this provides all the LVs predictions when by'd
					if (`: word count `istd`k''' > 1) {
						quietly factor `istd`k'' if `touse', factors(1) pcf
						quietly predict `tmpscores' if `touse' // factor scores
						quietly replace `LV`k'' = `tmpscores' if `touse'
						quietly label variable `LV`k''
						quietly drop `tmpscores'
					}
					else {
						quietly replace `LV`k'' = `: word 1 of `istd`k''' if `touse'
					}
				}
				else {
					if (`: word count `istd`k''' > 1) {
						quietly factor `istd`k'' if `touse', factors(1) pcf
						quietly predict `LV`k'' if `touse' // component scores
						quietly label variable `LV`k''
					}
					else {
						quietly generate `LV`k'' = `: word 1 of `istd`k''' if `touse'
					}
				}
				quietly summarize `LV`k'' if `touse'
				quietly replace `LV`k'' = (`LV`k'' - r(mean))/r(sd) if `touse' // equation (5)
			}
		}
	}
	else {
		foreach k of numlist 1/`num_lv' {
			tokenize `"`i`k''"'
			if (_byindex() > 1) {  // this provides all the LVs predictions when by'd
				quietly replace `LV`k'' = cond(`1' >= ., 0, `1') if `touse' // equation (4)
				// rs_* are the unstandardized versions of the latent variables
				quietly replace rs_`LV`k'' = cond(`1' >= ., 0, `1') if `touse' // equation (4)
			}
			else {
				quietly generate `LV`k'' = cond(`1' >= ., 0, `1') if `touse' // equation (4)
				quietly generate rs_`LV`k'' = cond(`1' >= ., 0, `1') if `touse' // equation (4)
			}
			macro shift
			while ("`1'" != "") {
				quietly replace `LV`k'' = `LV`k'' + cond(`1' >= ., 0, `1') if `touse'
				quietly replace rs_`LV`k'' = rs_`LV`k'' + cond(`1' >= ., 0, `1') if `touse'
				macro shift
			}
			quietly summarize `LV`k'' if `touse'
			local lv : word `k' of `alllatents'
			local isbinary : list lv in binary
			if (!`isbinary') {
				quietly replace `LV`k'' = (`LV`k'' - r(mean))/r(sd) if `touse' // equation (5)
			}
		}
	}

	/* Parse inner relationship */
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
					error 197
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
							error 197
						}

						local internm = "`var1'`var2'"
						if (_byindex() == 1) {  // this provides the LVs predictions when by'd
							capture confirm new variable `internm'
							if (_rc == 110) {
								quietly drop `internm'
							}
						}

						// Add the interaction to the list of latent variables
						local alllatents "`alllatents' `internm'"
						local modeA "`modeA' `internm'"
						local num_lv: word count `alllatents'
						local LV`num_lv' `internm'
						
						// Check that the new interaction is reflective
						local int_in_modeB : list internm in modeB
						if (`int_in_modeB') {
							display as error "interactions of latent variables must be reflective"
							error 197
						}

						// Generate the product indicators
						foreach ind1 in `i`var1'' {
							foreach ind2 in `i`var2'' {
								local indnm = "`ind1'`ind2'"
								if (_byindex() == 1) {  // this provides the indicators when by'd
									capture confirm new variable `indnm'
									if (_rc == 110) {
										quietly drop `indnm'
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
						
						// Standardize the product indicators
						foreach var in `inter_ind' {
							// quietly tabulate `var' if `touse'
							// if (r(r) > 2) {
								quietly summarize `var' if `touse'
								quietly generate std`var' = (`var' - r(mean))/r(sd) if `touse'
							// }
							// else {
							// 	quietly generate std`var' = `var' if `touse'  // do not standardize binary indicators
							// }
							local allstdindicators "`allstdindicators' std`var'"
						}
						foreach var in `i`num_lv'' {
							local istd`num_lv' "`istd`num_lv'' std`var'"
							local istd`LV`num_lv'' `istd`num_lv''
						}
						
						// Initialize the interactions
						if ("`init'" == "indsum") {
							local init_tmp `istd`num_lv''
							gettoken init_var init_tmp : init_tmp
							if (_byindex() > 1) {  // this provides all the LVs predictions when by'd
								quietly replace `LV`num_lv'' = cond(`init_var' >= ., 0, `init_var') if `touse' // equation (4)
							}
							else {
								quietly generate `LV`num_lv'' = cond(`init_var' >= ., 0, `init_var') if `touse' // equation (4)
							}
							gettoken init_var init_tmp : init_tmp
							while ("`init_var'" != "") {
								quietly replace `LV`num_lv'' = `LV`num_lv'' + cond(`init_var' >= ., 0, `init_var') if `touse'
								gettoken init_var init_tmp : init_tmp
							}
							quietly summarize `LV`num_lv'' if `touse'
							quietly replace `LV`num_lv'' = (`LV`num_lv'' - r(mean))/r(sd) if `touse' // equation (5)
						}
						else if ("`init'" == "eigen") {
							tempvar tmpscores
							if (_byindex() > 1) {  // this provides all the LVs predictions when by'd
								if (`: word count `istd`num_lv''' > 1) {
									quietly factor `istd`num_lv'' if `touse', factors(1) pcf
									quietly predict `tmpscores' if `touse' // component scores
									quietly replace `LV`num_lv'' = `tmpscores' if `touse'
									quietly label variable `LV`num_lv''
									quietly drop `tmpscores'
								}
								else {
									quietly replace `LV`num_lv'' = `: word 1 of `istd`num_lv''' if `touse'
								}
							}
							else {
								if (`: word count `istd`num_lv''' > 1) {
									quietly factor `istd`num_lv'' if `touse', factors(1) pcf
									quietly predict `LV`num_lv'' if `touse' // component scores
									quietly label variable `LV`num_lv''
								}
								else {
									quietly generate `LV`num_lv'' = `: word 1 of `istd`num_lv''' if `touse'
								}
							}
							quietly summarize `LV`num_lv'' if `touse'
							quietly replace `LV`num_lv'' = (`LV`num_lv'' - r(mean))/r(sd) if `touse' // equation (5)
						}
						
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
			syntax varlist(min=2)
			
			// Use macros r`var' and c`latentname' to store whether
			// the adjacency is treated as directional in path weighting scheme
			// The regest`var' macro is needed to store the path coefficients
			
			local dvar : word 1 of `varlist'
			local ivars : list varlist - dvar
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
				error 198
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
	
	/* Start of the PLS algorithm */	
	tempname converged iteration reldiff
	
	scalar `converged' = 0
	scalar `iteration' = 0
	
	if ("`structural'" != "") {
		tempname C W Wold mat0 matreldiff
		
		if ("`rawsum'" == "") {
			while (!`converged') {
				// Inner estimation. The three commonly used schemes are:
				// - centroid	==> sign of correlations
				// - factor		==> correlations
				// - path			==> correlations and regresions
				
				// Update the latents as weighted sums of adjacent latents
				// These are stored as separate temporary variables
				
				// (STEP 2)
				foreach var in `alllatents' {
					tempvar t`var'
					
					/* Inner estimation with regression weights (path) */
					if ("`r`var''" != "") {
						quietly regress `var' `r`var'' if `touse'
						quietly predict `t`var'' if `touse'
					}
					else {
						quietly generate `t`var'' = 0 if `touse'
					}
					
					/* Inner estimation with correlational weights (all schemes) */
					foreach var2 in `c`var'' {
						quietly correlate `var' `var2' if `touse'
						matrix `C' = r(C)
						if ("`wscheme'" == "centroid") {
							quietly replace `t`var'' = `t`var'' + `var2' * `C'[1,2]/abs(`C'[1,2]) ///
								if `touse' // Y~ = Y^*E for centroid scheme + equation (13)
						}
						else {
							quietly replace `t`var'' = `t`var'' + `var2' * `C'[1,2] ///
								if `touse' // Y~ = Y^*E for factorial and path schemes + equation (14)
						}
					}
					
					/* Standardize the LVs */
					local isbinary : list var in binary
					if (!`isbinary') {
						quietly summarize `var' if `touse'
						quietly replace `var' = (`var' - r(mean))/r(sd) if `touse' // equation (10)
						quietly summarize `t`var'' if `touse'
						quietly replace `t`var'' = (`t`var'' - r(mean))/r(sd) if `touse' // equation (10)
					}
				}
				
				// (STEP 3&4)
				// Store weights in matrix. These are unscaled and only used for
				// convergence check
				matrix `W' = 1
				
				/* Outer estimation (Mode A) */
				foreach var in `modeA' {
					quietly replace `var' = 0 if `touse'
					foreach var2 in `istd`var'' {
						quietly correlate `t`var'' `var2' if `touse' // equation (7)
						matrix `C' = r(C)
						matrix `W' = (`W', `C'[1, 2])
						quietly replace `var' = `var' + `var2' * `C'[1, 2] ///
							if `touse' // equation (9) for mode A LVs
					}
				}
				
				/* Outer estimation (Mode B) */
				foreach var in `modeB' {
					tempvar tv
					quietly regress `t`var'' `istd`var'' if `touse' // equation (8)
					matrix `W' = (`W', e(b))
					quietly predict `tv' if `touse'
					quietly replace `var' = `tv' if `touse' // equation (9) for mode B LVs
					quietly drop `tv'
				}

				/* Standardize the LVs */
				foreach var of varlist `alllatents' {
					local isbinary : list var in binary
					if (!`isbinary') {
						quietly summarize `var' if `touse'
						quietly replace `var' = (`var' - r(mean))/r(sd) if `touse' // equation (10)
					}
				}
				
				// (STEP 5)
				// Convergence check: compare new weights (W) with those from previous
				// iteration (Wold)

				if (`tol' <= 0) {
					display as error "tolerance must be a strictly positive number"
					error 197
				}
				if (`iteration' == 0) {
					matrix `mat0' = J(1, colsof(`W'), 0)
					scalar `reldiff' = mreldif(`W', `mat0')
					matrix `matreldiff' = `reldiff'
				}
				if (`iteration' > 0) {
					scalar `reldiff' = mreldif(`W', `Wold')
					matrix `matreldiff' = (`matreldiff', `reldiff')
					scalar `converged' = `reldiff' < `tol' // equation (11)
				}
				matrix `Wold' = `W'
			
				scalar `iteration' = `iteration' + 1

				/* No convergence */
				if (`maxiter' <= 0) {
					display as error "maximum number of iterations must be a positive integer"
					error 197
				}
				if (`iteration' > `maxiter') {
					error 430
				}
			}
		}
	}
	else {
		tempname C
		
		/* Outer estimation (Mode A) */
		foreach var in `modeA' {
			tempvar t`var'
			capture confirm new variable `t`var''
			if (_rc == 0) {
				quietly generate `t`var'' = `var' if `touse'
			}
			else {
				quietly replace `t`var'' = `var' if `touse'
			}
			quietly replace `var' = 0 if `touse'
			foreach var2 in `istd`var'' {
				quietly correlate `t`var'' `var2' if `touse' // equation (7)
				matrix `C' = r(C)
				quietly replace `var' = `var' + `var2' * `C'[1, 2] ///
					if `touse' // equation (9) for mode A LVs
			}
		}
		
		/* Outer estimation (Mode B) */
		foreach var in `modeB' {
			tempvar t`var'
			capture confirm new variable `t`var''
			if (_rc == 0) {
				quietly generate `t`var'' = `var' if `touse'
			}
			else {
				quietly replace `t`var'' = `var' if `touse'
			}
			tempvar tv
			quietly regress `t`var'' `istd`var'' if `touse' // equation (8)
			quietly predict `tv' if `touse'
			quietly replace `var' = `tv' if `touse' // equation (9) for mode B LVs
			quietly drop `tv'
		}

		/* Standardize the LVs */
		foreach var of varlist `alllatents' {
			local isbinary : list var in binary
			if (!`isbinary') {
				quietly summarize `var' if `touse'
				quietly replace `var' = (`var' - r(mean))/r(sd) if `touse' // equation (10)
			}
		}
	}
	
	/* Make the binary latent equal to 0 or 1 */
	foreach var in `binary' {
		quietly summarize `var' if `touse'
		local minval r(min)
		local maxval r(max)
		if (`minval' < `maxval') {
			quietly replace `var' = 0 if (`var' == `minval') & `touse'
			quietly replace `var' = 1 if (`var' == `maxval') & `touse'
		}
		else {
			quietly replace `var' = 1 if (`var' == `minval') & `touse'
			quietly replace `var' = 0 if (`var' == `maxval') & `touse'
		}
	}
	
	/* Compute the table of measurement loadings */
	tempname loadings loadings_se adj_meas Cor Cor_var ave ave_num
	if ("`boot'" != "") {
		tempname loadings_bs Cor_bs
	}
	local num_ind: word count `allindicators'
	local num_lv_A: word count `modeA'
	local num_lv_B: word count `modeB'
	matrix `loadings' = J(`num_ind', `num_lv', .)
	if ("`boot'" != "") {
		matrix `loadings_bs' = J(`num_ind', `num_lv', .)
	}
	matrix `loadings_se' = J(`num_ind', `num_lv', .)
	matrix `adj_meas' = J(`num_ind', `num_lv', 0)
	if (`num_lv_A' > 0) {
		matrix `ave' = J(1, `num_lv_A', .)
		matrix `ave_num' = J(1, `num_lv_A', .)
	}
	foreach var in `modeA' {
		local loadcolnames `loadcolnames' "Reflective:`var'"
		local loadrownames `loadrownames' `i`var''
	}
	foreach var in `modeB' {
		local loadcolnames `loadcolnames' "Formative:`var'"
		local loadrownames `loadrownames' `i`var''
	}
	matrix rownames `loadings' = `loadrownames'
	matrix colnames `loadings' = `loadcolnames'
	if ("`boot'" != "") {
		matrix rownames `loadings_bs' = `loadrownames'
		matrix colnames `loadings_bs' = `loadcolnames'
	}
	matrix rownames `adj_meas' = `loadrownames'
	matrix colnames `adj_meas' = `loadcolnames'
	if (`num_lv_A' > 0) {
		matrix rownames `ave' = "AVE"
		matrix colnames `ave' = `modeA'
		matrix colnames `ave_num' = `modeA'
	}
	
	local start = 1
	local i = 1
	tempname ave_tmp
	tempvar var_std
	quietly generate `var_std' = . if `touse'
	foreach var in `modeA' {
		scalar `ave_tmp' = 0
		local to_use: word count `i`var''
		local isbinary : list var in binary
		foreach var2 in `istd`var'' {
			if ("`boot'" == "") {
				if ("`rawsum'" == "") {
					quietly regress `var' `var2' if `touse'
				}
				else {
					quietly summarize `var' if `touse'
					quietly replace `var_std' = (`var' - r(mean))/r(sd) if `touse'
					quietly regress `var_std' `var2' if `touse'
				}
				matrix `Cor' = e(b)
				matrix `Cor_var' = e(V)
			}
			else {
				if ("`rawsum'" == "") {
					quietly bootstrap, reps(`boot') seed(`seed'): regress `var' `var2' if `touse'
				}
				else {
					quietly summarize `var' if `touse'
					quietly replace `var_std' = (`var' - r(mean))/r(sd) if `touse'
					quietly bootstrap, reps(`boot') seed(`seed'): regress `var_std' `var2' if `touse'
				}
				matrix `Cor' = e(b)
				if ("`boot'" != "") {
					matrix `Cor_bs' = e(b_bs)
				}
				matrix `Cor_var' = e(V)
			}
			if (!`isbinary') {
				matrix `loadings'[`start', `i'] = `Cor'[1, 1]
				if ("`boot'" != "") {
					matrix `loadings_bs'[`start', `i'] = `Cor_bs'[1, 1]
				}
				matrix `loadings_se'[`start', `i'] = `Cor_var'[1, 1]
			}
			else {
				matrix `loadings'[`start', `i'] = 1
				if ("`boot'" != "") {
					matrix `loadings_bs'[`start', `i'] = 1
				}
				matrix `loadings_se'[`start', `i'] = 0
			}
			local ++start
*			if (!`isbinary') {
				scalar `ave_tmp' = `ave_tmp' + `Cor'[1, 1]^2/`to_use'
*			}
*			else {
*				scalar `ave_tmp' = 1
*			}
		}
		matrix `ave_num'[1, `i'] = `to_use'
		matrix `ave'[1, `i'] = `ave_tmp'
		local end = `start' + `to_use' - 1
		local ++i
	}
	
	tempname modeB_b modeB_var
	if ("`boot'" != "") {
		tempname modeB_b_bs
	}
	foreach var in `modeB' {
		local isbinary : list var in binary
		if ("`boot'" == "") {
			if ("`rawsum'" == "") {
				quietly regress `var' `istd`var'' if `touse'
			}
			else {
				quietly summarize `var' if `touse'
				quietly replace `var_std' = (`var' - r(mean))/r(sd) if `touse'
				quietly regress `var_std' `istd`var'' if `touse'
			}
			matrix `modeB_b' = e(b)
			matrix `modeB_var' = e(V)
		}
		else {
			if ("`rawsum'" == "") {
				quietly bootstrap, reps(`boot') seed(`seed'): regress `var' `istd`var'' if `touse'
			}
			else {
				quietly summarize `var' if `touse'
				quietly replace `var_std' = (`var' - r(mean))/r(sd) if `touse'
				quietly bootstrap, reps(`boot') seed(`seed'): regress `var_std' `istd`var'' if `touse'
			}
			matrix `modeB_b' = e(b)
			matrix `modeB_b_bs' = e(b_bs)
			matrix `modeB_var' = e(V)
		}
		local to_use: word count `i`var''
		local end = `start' + `to_use' - 1
		matrix `modeB_b' = `modeB_b''
		if ("`boot'" != "") {
			matrix `modeB_b_bs' = `modeB_b_bs''
		}
		matrix `modeB_var' = vecdiag(`modeB_var')
		matrix `modeB_var' = `modeB_var''
		if (!`isbinary') {
			matrix `loadings'[`start', `i'] = `modeB_b'[1..`to_use', 1]
			if ("`boot'" != "") {
				matrix `loadings_bs'[`start', `i'] = `modeB_b_bs'[1..`to_use', 1]
			}
			matrix `loadings_se'[`start', `i'] = `modeB_var'[1..`to_use', 1]
		}
		else {
			matrix `loadings'[`start', `i'] = 1
			if ("`boot'" != "") {
				matrix `loadings_bs'[`start', `i'] = 1
			}
			matrix `loadings_se'[`start', `i'] = 0
		}
		local start = `end' + 1
		local ++i
	}
	mata: st_matrix("`loadings_se'", sqrt(st_matrix("`loadings_se'")))
	matrix rownames `loadings_se' = `loadrownames'
	matrix colnames `loadings_se' = `loadcolnames'

	forvalues i = 1/`num_ind' {
		forvalues j = 1/`num_lv' {
			if (!missing(`loadings'[`i', `j'])) {
				matrix `adj_meas'[`i', `j'] = 1
			}
		}
	}
	
	/* Compute the table of average variance extracted */
	if (`num_lv_A' > 0) {
		tempname sqcorr
		matrix `sqcorr' = J(`num_lv_A', `num_lv_A', .)
		local count = 2
		local j = 1
		foreach var in `modeA' {
			local isbinary : list var in binary
			local resttmp `resttmp' `var'
			local rest : list modeA - resttmp
			local i = `count'
			foreach var2 in `rest' {
				if ("`boot'" == "") {
					quietly regress `var' `var2' if `touse'
					matrix `Cor' = e(b)
				}
				else {
					quietly bootstrap, reps(`boot') seed(`seed'): regress `var' `var2' if `touse'
					matrix `Cor' = e(b)
				}
				matrix `Cor' = `Cor'[1, 1]
	*			if (!`isbinary') {
					matrix `sqcorr'[`i', `j'] = `Cor'[1, 1]^2
	*			}
	*			else {
	*				matrix `sqcorr'[`i', `j'] = 1
	*			}
				local ++i
			}
			local ++j
			local ++count
		}
		forvalues ii = 1/`num_lv_A' {
			matrix `sqcorr'[`ii', `ii'] = 1
		}
		matrix rownames `sqcorr' = `modeA'
		matrix colnames `sqcorr' = `modeA'
	}
	
	/* Compute the table of summary statistics for indicators */
	// ("if `touse'" is not used here to correctly count missing values)
	if ("`stats'" != "") {
		tempname indstats
		local k = 1
		matrix `indstats' = J(`num_ind', 7, .)
		foreach ind in `allindicators' {
			// quietly summarize `ind' if `touse', detail
			quietly summarize `ind', detail
			matrix `indstats'[`k', 1] = r(mean)
			matrix `indstats'[`k', 2] = r(sd)
			matrix `indstats'[`k', 3] = r(p50)
			matrix `indstats'[`k', 4] = r(min)
			matrix `indstats'[`k', 5] = r(max)
			local nonmissN = r(N)
			// quietly misstable summarize `ind' if `touse'
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
	
	/* Compute the table of structural path coefficients */
	if ("`structural'" != "") {
		foreach var in `alllatents' {
			if ("`regest`var''" != "") {
				local lv_regest_all `lv_regest_all' `var'
			}
			local regest_all `regest_all' `regest`var''
		}
		local lv_regest_all : list uniq lv_regest_all
		local regest_all : list uniq regest_all
		local regest_all = strtrim("`regest_all'")
		foreach var in `regest_all' {
			local pathtab_rownames "`pathtab_rownames' `var' ."
		}
		local num_lv_path : word count `lv_regest_all'
		local num_path : word count `regest_all'
		tempname pathtab rsquared redundancy tmp_b tmp_var tmp_se tmp_df
		if ("`boot'" != "") {
			tempname tmp_b_bs
		}
		matrix `pathtab' = J(2*`num_path' + 1, `num_lv_path', .)
		if (`num_lv_A' > 0) {
			matrix `rsquared' = J(1, `num_lv_A', .)
			matrix `redundancy' = J(1, `num_lv_A', .)
		}
		local lv_j = 1
		foreach var in `alllatents' {
			local isbinary : list var in binary
			if ("`regest`var''" != "") {
				local numind_lv : word count `regest`var''
				if ("`boot'" == "") {
					if (!`isbinary') {
						quietly regress `var' `regest`var'' if `touse'
					}
					else {
						quietly logit `i`var'' `regest`var'' if `touse'
					}
					matrix `tmp_b' = e(b)
				}
				else {
					if (!`isbinary') {
						quietly bootstrap, reps(`boot') seed(`seed'): regress `var' `regest`var'' if `touse'
					}
					else {
						quietly bootstrap, reps(`boot') seed(`seed'): logit `i`var'' `regest`var'' if `touse'
					}
					matrix `tmp_b' = e(b)
					matrix `tmp_b_bs' = e(b_bs)
				}
				matrix `tmp_var' = e(V)
				scalar `tmp_df' = e(N) - e(df_m) - 1
				local b_colnames : colnames `tmp_b'
				forvalues lv_i = 1/`numind_lv' {
					local xn : word `lv_i' of `b_colnames'
					local xi : list posof "`xn'" in pathtab_rownames
					matrix `pathtab'[`xi', `lv_j'] = `tmp_b'[1, `lv_i']
					matrix `pathtab'[`xi' + 1, `lv_j'] = 2*ttail(`tmp_df', ///
						abs(`tmp_b'[1, `lv_i']/sqrt(`tmp_var'[`lv_i', `lv_i'])))
				}
				if (!`isbinary') {
					matrix `pathtab'[2*`num_path' + 1, `lv_j'] = e(r2_a)
				}
				else {
					matrix `pathtab'[2*`num_path' + 1, `lv_j'] = e(r2_p)
				}
				local ismodeA : list var in modeA
				local isbinary : list var in binary
				if (`ismodeA' & !`isbinary') {
					if (`num_lv_A' > 0) {
					local vari : list posof "`var'" in modeA
					matrix `rsquared'[1, `vari'] = e(r2)
					matrix `redundancy'[1, `vari'] = `ave'[1, `vari'] * `rsquared'[1, `vari']
					}
				}
				estimates store `var'
				local regest_tbl `regest_tbl' `var'
				local ++lv_j
			}
		}
		matrix rownames `pathtab' = `pathtab_rownames' "r2_a"
		matrix colnames `pathtab' = `lv_regest_all'
		if (`num_lv_A' > 0) {
			matrix rownames `rsquared' = "r2"
			matrix colnames `rsquared' = `modeA'
			matrix rownames `redundancy' = "redundancy"
			matrix colnames `redundancy' = `modeA'
		}
	}
	
	/* Compute reliability coefficients */
	local rr = 1
	tempname dgnum dgden relcoef
	matrix `relcoef' = J(2, `num_lv', .)
	local start = 1
	foreach var in `modeA' {
		local num_ind_A = `: word count `i`var'''
		local sqsumload = 0
		if (`num_ind_A' > 1) {
			quietly alpha `i`var'' if `touse', std
			matrix `relcoef'[1, `rr'] = r(alpha)
		}
		else {
			matrix `relcoef'[1, `rr'] = 1
		}
		forvalues hh = 1/`num_ind_A' {
			local sqsumload = `sqsumload' + `loadings'[`start' + `hh' - 1, `rr']
		}
		scalar `dgnum' = `sqsumload'^2
		scalar `dgden' = `dgnum' + `num_ind_A'*(1 - `ave'[1, `rr'])
		matrix `relcoef'[2, `rr'] = `dgnum'/`dgden'
		local start = `start' + `num_ind_A'
		local ++rr
	}
	matrix rownames `relcoef' = "Cronbach" "DG"
	matrix colnames `relcoef' = `loadcolnames'
	
	/* Setup the table of structural coefficients */
	if ("`structural'" != "") {
		tempname strcoef_tmp strcoef strse strstats strvif B adj_struct
		if ("`boot'" != "") {
			tempname B_bs
		}
		quietly estimates table `regest_tbl', b(%12.`digits'f) p(%12.`digits'f) ///
			drop(_cons)
		matrix `strcoef_tmp' = r(coef)
		matrix `strstats' = r(stats)
		local struct_rownames : rownames `strcoef_tmp'
		local num_reg_coef = rowsof(`strcoef_tmp')
		local num_col_coef = colsof(`strcoef_tmp')
		local g = 1
		local q = 1
		matrix `strcoef' = J(`num_reg_coef', `num_lv_path', .)
		matrix rownames `strcoef' = `struct_rownames'
		matrix colnames `strcoef' = `lv_regest_all'
		matrix `strse' = J(`num_reg_coef', `num_lv_path', .)
		matrix rownames `strse' = `struct_rownames'
		matrix colnames `strse' = `lv_regest_all'
		forvalues j = 1/`num_col_coef' {
			forvalues h = 1/`num_reg_coef' {
				if (mod(`j', 2) != 0) {
					if (!missing(`strcoef_tmp'[`h', `j'])) {
						matrix `strcoef'[`h', `g'] = `strcoef_tmp'[`h', `j']
					}
				}
				else {
					if (!missing(`strcoef_tmp'[`h', `j'])) {
						matrix `strse'[`h', `q'] = sqrt(`strcoef_tmp'[`h', `j'])
					}
				}
			}
			if (mod(`j', 2) != 0) {
				local ++g
			}
			else {
				local ++q
			}
		}
		matrix `strvif' = J(`num_reg_coef', `num_lv_path', .)
		matrix rownames `strvif' = `struct_rownames'
		matrix colnames `strvif' = `lv_regest_all'
		matrix `B' = J(`num_lv', `num_lv', 0)
		matrix rownames `B' = `alllatents'
		matrix colnames `B' = `alllatents'
		if ("`boot'" != "") {
			matrix `B_bs' = J(`num_lv', `num_lv', 0)
			matrix rownames `B_bs' = `alllatents'
			matrix colnames `B_bs' = `alllatents'
		}
		matrix `adj_struct' = J(`num_lv', `num_lv', 0)
		matrix rownames `adj_struct' = `alllatents'
		matrix colnames `adj_struct' = `alllatents'
		local lv_j = 1
		local B_j = 1
		foreach var in `alllatents' {
			local isbinary : list var in binary
			if ("`regest`var''" != "" & !`isbinary') {
				quietly estimates restore `var'
				matrix `tmp_b' = e(b)
				if ("`boot'" != "") {
					matrix `tmp_b_bs' = e(b_bs)
				}
				local b_colnames : colnames `tmp_b'
				if ("`boot'" == "") {
					local numind_lv : word count `regest`var''
					quietly estat vif
					forvalues lv_i = 1/`numind_lv' {
						local xn : word `lv_i' of `b_colnames'
						local xi : list posof "`xn'" in struct_rownames
						matrix `strvif'[`xi', `lv_j'] = ${S_`lv_i'}
					}
					local ++lv_j
				}
				forvalues lv_i = 1/`num_lv' {
					local xn : word `lv_i' of `alllatents'
					local isthere : list xn in b_colnames
					if (`isthere') {
						local B_i : list posof "`xn'" in alllatents
						local xi : list posof "`xn'" in b_colnames
						matrix `B'[`B_i', `B_j'] = `tmp_b'[1, `xi']
						if ("`boot'" != "") {
							matrix `B_bs'[`B_i', `B_j'] = `tmp_b_bs'[1, `xi']
						}
						matrix `adj_struct'[`B_i', `B_j'] = 1
					}
				}
			}
			local ++B_j
		}
	}
	
	/* Calculate model assessment indexes (GOF, etc.) */
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
	
	/* Return values */
	local props "`init'"
	if ("`structural'" != "") {
		local props "`props' `wscheme'"
	}
	if ("`boot'" != "") {
		local props "`props' bootstrap"
	}
	if ("`structural'" != "") {
		local props "`props' structural"
	}
	ereturn post, obs(`nobs') esample(`touse') properties(`props')
	if ("`boot'" != "") {
		ereturn scalar reps = `boot'
	}
	if ("`structural'" != "") {
		if ("`rawsum'" == "") {
			ereturn scalar iterations = `iteration' - 1
		}
		else {
			ereturn scalar iterations = 0
		}
		ereturn scalar tolerance = `tol'
	}
	ereturn local struct_eqs `reg3eqs'
	ereturn local formative `modeB'
	ereturn local reflective `modeA'
	ereturn local lvs `alllatents'
	ereturn local mvs `allindicators'
	ereturn local estat_cmd "plssem_estat"
	ereturn local cmdline `cmdline'
	ereturn local cmd "plssem"
	if ("`structural'" != "") {
		if ("`rawsum'" == "") {
			ereturn matrix reldiff = `matreldiff'
		}
		if (`num_lv_A' > 0) {
			ereturn matrix assessment = `mod_assessment'
			ereturn matrix redundancy = `redundancy'
			ereturn matrix rsquared = `rsquared'
		}
		ereturn matrix adj_struct = `adj_struct'
		if ("`boot'" != "") {
			ereturn matrix pathcoef_bs = `B_bs'
		}
		ereturn matrix pathcoef = `B'
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
	ereturn matrix loadings_se = `loadings_se'
	if ("`boot'" != "") {
		ereturn matrix loadings_bs = `loadings_bs'
	}
	ereturn matrix loadings = `loadings'
	if ("`stats'" != "") {
		ereturn matrix indstats = `indstats'
	}
	
	/* Clean up */
	foreach var in `allstdindicators' {
		quietly drop `var'
	}
	if ("`structural'" != "") {
		estimates drop `regest_tbl'
	}
	if ("`rawsum'" != "") {
		foreach var in `alllatents' {
			quietly replace `var' = rs_`var'
			quietly drop rs_`var'
		}
	}
	
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
end

program Compare, eclass sortpreserve
	version 10
	syntax anything(id="Indicator blocks" name=blocks equalok) [if] [in], ///
		[ STRuctural(string) Wscheme(string) BINary(namelist min=1) ///
		Boot(numlist integer >0 max=1) Seed(numlist max=1) Tol(real 1e-7) ///
		MAXiter(integer 1000) INIT(string) DIGits(integer 3) noHEADer ///
		noMEAStable noDISCRIMtable noSTRUCTtable STATs GRoup(string) ///
		CORRelate(string) RAWsum ]

	/* Options:
	   --------
		 structural(string)							--> structural model specification
		 wscheme(string)								--> weighting scheme ("centroid", "factorial"
																				or "path")
		 binary(namelist min=1)					--> indicators to fit using logit
		 boot(numlist integer >0 max=1)	--> bootstrap estimation (# of repetions)
		 seed(numlist max=1)						--> bootstrap seed number
		 tol(real 1e-7)									--> tolerance (default 1e-7)
		 maxiter(integer 1000)					--> maximum number of iterations (default
																				1000)
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
	 */
	
	/* Warning */
	if ("`boot'" != "") {
		display as error "using the 'boot()' option together with 'group()' will slow down excessively the calculation"
		display as error "try removing 'boot()'"
		exit
	}
	
	/* Set obs to use */
	marksample touse
	
	quietly count if `touse'
	if (r(N) == 0) {
		error 2000
	}
	
	local skip1 = 1
	local skip3 = 3

	/* Parse the group() option */
	tokenize `"`group'"', parse(",")
	local ngroupvars : word count `1'
	if (`ngroupvars' > 1) {
		display as error "a single variable allowed in the 'group()' option"
		error 198
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
					error 198
				}
			}
			else if ("``tok_i''" == "method") {
				local tok_tmp = `tok_i' + 2
				local method = "``tok_tmp''"
				if (("`method'" != "normal") & ("`method'" != "permutation") & ("`method'" != "bootstrap")) {
					display as error "'method()' suboption accepts either 'permutation', 'bootstrap' or 'normal'"
					error 198
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
					error 198
				}
			}
			else if ("``tok_i''" == "groupseed") {
				local tok_tmp = `tok_i' + 2
				local groupseed = ``tok_tmp''
				if (`groupseed' < 0 | `groupseed' >  2^31-1) {
					display as error "'groupseed()' suboption requires a value between 0 and 2^31-1"
					error 198
				}
			}
			else if ("``tok_i''" == "what") {
				local tok_tmp = `tok_i' + 2
				local what = "``tok_tmp''"
				if (("`what'" != "path") & ("`what'" != "loadings")) {
					display as error "'what()' suboption accepts either 'path' or 'loadings'"
					error 198
				}
			}
		}
		local ++tok_i
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
		error 198
	}
	/* End of parsing the group() option */
	
	if ("`groupseed'" != "") {
		set seed `groupseed'
	}
	
	tempname groupvar_m groupvals groupsizes alln
	mata: st_matrix("`groupvar_m'", uniqrows(st_data(., "`groupvar'", "`touse'"), 1))
	local ngroups = rowsof(`groupvar_m')
	if (`ngroups' == 1) {
		display as error "the group() option requires at least two groups to compare"
		error 198
	}
	matrix `groupvals' = `groupvar_m'[1..., 1]
	matrix `groupsizes' = `groupvar_m'[1..., 2]
	mata: st_numscalar("`alln'", sum(st_matrix("`groupsizes'")))
	
	tempname results whole ehold tmp

	/* Using the whole sample */
	quietly Estimate `0'
	local allindicators = e(mvs)
	local alllatents = e(lvs)
	local allreflective = e(reflective)
	local nind : word count `allindicators'
	local nlv : word count `alllatents'
	local nlv_A : word count `allreflective'
	local nlv_B = `nlv' - `nlv_A'
	if ("`what'" == "path") {
		matrix `tmp' = e(struct_b)
		matrix `whole' = vec(`tmp')
	}
	else if ("`what'" == "loadings") {
		matrix `tmp' = e(loadings)
//		matrix `tmp' = `tmp'[1..`nind', 1..`nlv_A']
		matrix `whole' = vec(`tmp')
	}
	_estimates hold `ehold', restore
	
	/* Calculate observed test statistic value */
	forvalues ng = 1/`ngroups' {
		preserve
		quietly keep if (`groupvar' == `groupvals'[`ng', 1]) & `touse'
		quietly Estimate `0'
		tempname strb_`ng' strvar_`ng' group_`ng'
		if ("`what'" == "path") {
			matrix `tmp' = e(struct_b)
			matrix `strb_`ng'' = vec(`tmp')
			matrix `tmp' = e(struct_se)
			matrix `strvar_`ng'' = vec(`tmp')
		}
		else if ("`what'" == "loadings") {
			matrix `tmp' = e(loadings)
//			matrix `tmp' = `tmp'[1..`nind', 1..`nlv_A']
			matrix `strb_`ng'' = vec(`tmp')
			matrix `tmp' = e(loadings_se)
//			matrix `tmp' = `tmp'[1..`nind', 1..`nlv_A']
			matrix `strvar_`ng'' = vec(`tmp')
		}
		restore
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
		tempname strbok_`ng' strvarok_`ng' strbok_boot_`ng' strbok_var_`ng'
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
	
	tempname dtest dtest_tmp
	matrix `dtest' = J(`ngroups' - 1, `neff', .)
	matrix rownames `dtest' = `obstest_cn'
	matrix colnames `dtest' = `strbok_rn'
	matrix `dtest_tmp' = J(`reps', `neff', 0)
	matrix colnames `dtest_tmp' = `strbok_rn'
	
	if (("`method'" == "permutation") | ("`method'" == "bootstrap")) {
		capture {
			noisily display
			if ("`method'" == "permutation") {
				local title "Multigroup comparison (`groupvar') - Permutation test"

				/* Compute permutation distribution */
				tempvar __id__ __group__
				quietly generate `__id__' = _n if `touse'
				quietly generate `__group__' = . if `touse'

				noisily display as text "Permutation replications (" ///
					as result string(`reps') as text ")"
				noisily display as text ///
					"{hline 4}{c +}{hline 3}" " 1 " ///
					"{hline 3}{c +}{hline 3}" " 2 " ///
					"{hline 3}{c +}{hline 3}" " 3 " ///
					"{hline 3}{c +}{hline 3}" " 4 " ///
					"{hline 3}{c +}{hline 3}" " 5"
				
				tempname diff_tmp dtest_c
				forvalues ng = 1/`ngroups' {
					local ngm1 = `ng' - 1
					if (`ng' == 1) {
						local from`ng' = 1
						local to`ng' = `groupsizes'[`ng', 1]
					}
					else {
						local from`ng' = `to`ngm1'' + 1
						local to`ng' = `to`ngm1'' + `groupsizes'[`ng', 1]
					}
					if (`ng' > 1) {
						tempname diff`ng' dtest_tmp`ng'
						matrix `diff`ng'' = J(`reps', `neff', .)
						matrix colnames `diff`ng'' = `strbok_rn'
						matrix `dtest_tmp`ng'' = J(`reps', `neff', .)
						matrix colnames `dtest_tmp`ng'' = `strbok_rn'
					}
				}
				
				forvalues h = 1/`reps' {
					if (mod(`h', 50) == 0) {
						local todisp = strtrim(string(`h', "%9.0f"))
						if (strlen("`todisp'") < strlen(string(`reps'))) {
							local tmp_skip = strlen(string(`reps')) - strlen("`todisp'")
						}
						noisily display as text "." _skip(3) _skip(`tmp_skip') "`todisp'"
						local tmp_skip = 0
					}
					else {
						noisily display as text "." _continue
					}

					mata: V = st_data(., "`__id__'", "`touse'")
					mata: st_store(., "`__id__'", "`touse'", jumble(V))
					sort `__id__'
					forvalues ng = 1/`ngroups' {
						quietly replace `__group__' = `groupvals'[`ng', 1] if `touse' in `from`ng''/`to`ng''
					}
					
					forvalues ng = 1/`ngroups' {
						preserve
						quietly keep if (`__group__' == `groupvals'[`ng', 1]) & `touse'
						quietly Estimate `0'
						if ("`what'" == "path") {
							matrix `tmp' = e(struct_b)
							matrix `strb_`ng'' = vec(`tmp')
						}
						else if ("`what'" == "loadings") {
							matrix `tmp' = e(loadings)
//							matrix `tmp' = `tmp'[1..`nind', 1..`nlv_A']
							matrix `strb_`ng'' = vec(`tmp')
						}
						restore
						mata: st_matrix("`strbok_`ng''", select(st_matrix("`strb_`ng''"), ///
							st_matrix("`nonmiss'")))

						if (`ng' > 1) {
							mata: st_matrix("`diff_tmp'", abs(st_matrix("`strbok_1'") - ///
								st_matrix("`strbok_`ng''"))')
							matrix `diff`ng''[`h', 1] = `diff_tmp'
							forvalues j = 1/`neff' {
								if (`obstest'[`j', `ng' - 1] < `diff`ng''[`h', `j']) {
									matrix `dtest_tmp`ng''[`h', `j'] = 1
								}
							}
						}
					}
				}

				forvalues ng = 2/`ngroups' {
					mata: st_matrix("`dtest_c'", colsum(st_matrix("`dtest_tmp`ng''")))
					matrix `dtest'[`ng' - 1, 1] = `dtest_c'
				}
			}
			else if ("`method'" == "bootstrap") {
				local title "Multigroup comparison (`groupvar') - Bootstrap t-test"
				
				quietly snapshot save, label("pls_whole_model")
				local snapwhole = r(snapshot)
				
				tempname tmp_b_boot tmp_var diff_mean diff_var k0 k1 k2 k3
				scalar `k0' = `alln' - 2
				scalar `k1' = ((`groupsizes'[1, 1] - 1)^2)/`k0'

				/* Compute bootstrap distribution */
				forvalues ng = 1/`ngroups' {
					matrix `strbok_boot_`ng'' = J(`reps', `neff', .)
					matrix colnames `strbok_boot_`ng'' = `strbok_rn'
					matrix `strbok_var_`ng'' = J(1, `neff', .)
				}
				
				noisily display as text "Bootstrap replications (" ///
					as result string(`reps') as text ")"
				noisily display as text ///
					"{hline 4}{c +}{hline 3}" " 1 " ///
					"{hline 3}{c +}{hline 3}" " 2 " ///
					"{hline 3}{c +}{hline 3}" " 3 " ///
					"{hline 3}{c +}{hline 3}" " 4 " ///
					"{hline 3}{c +}{hline 3}" " 5"

				local tmp_skip = 0
				forvalues h = 1/`reps' {
					if (mod(`h', 50) == 0) {
						local todisp = strtrim(string(`h', "%9.0f"))
						if (strlen("`todisp'") < strlen(string(`reps'))) {
							local tmp_skip = strlen(string(`reps')) - strlen("`todisp'")
						}
						noisily display as text "." _skip(3) _skip(`tmp_skip') "`todisp'"
						local tmp_skip = 0
					}
					else {
						noisily display as text "." _continue
					}
					
					bsample if `touse', strata(`groupvar')
					forvalues ng = 1/`ngroups' {
						preserve
						quietly keep if (`groupvar' == `groupvals'[`ng', 1]) & `touse'
						quietly Estimate `0'
						if ("`what'" == "path") {
							matrix `tmp' = e(struct_b)
							matrix `strb_`ng'' = vec(`tmp')
						}
						else if ("`what'" == "loadings") {
							matrix `tmp' = e(loadings)
//							matrix `tmp' = `tmp'[1..`nind', 1..`nlv_A']
							matrix `strb_`ng'' = vec(`tmp')
						}
						restore
						mata: st_matrix("`strbok_`ng''", select(st_matrix("`strb_`ng''"), ///
							st_matrix("`nonmiss'")))
						matrix `strbok_boot_`ng''[`h', 1] = `strbok_`ng'''
					}
					snapshot restore `snapwhole'
				}
				
				forvalues ng = 1/`ngroups' {
					forvalues j = 1/`neff' {
						matrix `tmp_b_boot' = `strbok_boot_`ng''[1..`reps', `j']
						mata: st_numscalar("`tmp_var'", variance(st_matrix("`tmp_b_boot'")))
						matrix `strbok_var_`ng''[1, `j'] = `tmp_var'
					}
					if (`ng' > 1) {
						mata: st_matrix("`diff_mean'", abs(mean(st_matrix("`strbok_boot_1'")) - ///
							mean(st_matrix("`strbok_boot_`ng''"))))
						scalar `k2' = ((`groupsizes'[`ng', 1] - 1)^2)/`k0'
						scalar `k3' = sqrt(1/`groupsizes'[1, 1] + 1/`groupsizes'[`ng', 1])
						forvalues j = 1/`neff' {
							matrix `dtest'[`ng' - 1, `j'] = `diff_mean'[1, `j']/(sqrt( ///
								`k1'*`strbok_var_1'[1, `j'] + `k2'*`strbok_var_`ng''[1, `j'])*`k3')
						}
					}
				}

				snapshot erase `snapwhole'
			}
		} // end of -capture-
		local rc = _rc
		_estimates unhold `ehold'
		if (`rc' != 0) {
			if ("`method'" == "bootstrap") {
				snapshot erase `snapwhole'
			}
			error `rc'
		}
	}
	else if ("`method'" == "normal") {
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
						2*ttail(`k0', `dtest'[`ng' - 1, `j'])
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
	
	if ("`structural'" != "") {
		mkheader, digits(5)
	}
	else {
		mkheader, digits(5) nostructural
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
		display as txt _skip(`skip1') "number of replications: `reps'"
	}
	if (`ngroups' > 2) {
		display as txt _skip(`skip1') "legend:"
		display as txt _skip(`skip3') "AD: absolute difference"
		display as txt _skip(`skip3') "S: statistic"
		display as txt _skip(`skip3') "P: p-value"
	}
	display as txt _skip(`skip1') "group labels:"
	local grplbl : value label `groupvar'
	forvalues ng = 1/`ngroups' {
		local grpvarlbl`ng' = `groupvals'[`ng', 1]
		if ("`grplbl'" != "") {
			local grpvarlbl`ng' : label `grplbl' `grpvarlbl`ng''
		}
		display as txt _skip(`skip3') "Group `ng': `grpvarlbl`ng''"
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
		quietly drop `reslbl'* id
	}
	
	/* Return */
	// ereturn local method = "`method'"
	// ereturn matrix comparison = `results'
end

program Display
	version 10
	syntax [, noSTRuctural DIGits(integer 3) noHEADer noMEAStable ///
		noDISCRIMtable noSTRUCTtable stats corrind corrlv crossload ///
		CUToff(real 0) BINary(namelist min=1) RAWsum * ]
	
	if (`digits' < 0) {
		display as error "number of digits must be a nonnegative integer"
		error 197
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
