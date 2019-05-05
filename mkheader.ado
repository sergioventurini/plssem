*!mkheader version 0.3.0
*!Written 27Apr2018
*!Written by Sergio Venturini and Mehmet Mehmetoglu
*!The following code is distributed under GNU General Public License version 3 (GPL-3)

program mkheader
	version 15.1
	syntax [, matrix1(string) matrix2(string) DIGits(integer 5) noGRoup ///
		noSTRuctural RAWsum rebus_it(integer -999) rebus_gqi(real 0) ///
		 fimix_it(integer -999) fimix_ll(real 0) GAS CONsistent ]

	/* Options:
	   --------
		 matrix1(string)						--> matrix containing the model assessment
																		indexes
		 matrix2(string)						--> matrix containing the outer weights maximum
																		relative difference history
		 digits(integer 5)					--> number of digits to display (default 5)
		 nogroup										--> multigroup indicator
		 nostructural								--> indicator of whether the model includes a
																		structural part
		 rawsum											--> scores are raw sum of indicators
		 rebus_it										--> number of REBUS-PLS iterations
		 rebus_gqi									--> REBUS-PLS group quality index (GQI)
		 fimix_it										--> number of FIMIX-PLS iterations
		 fimix_ll										--> FIMIX-PLS log-likelihood value attained
		 gas												--> PLS-GAS indicator
		 consistent									--> indicator for consistent PLS (PLSc)
	 */

	local props = e(properties)
	local init : word 1 of `props'
	local wscheme : word 2 of `props'
	local tol = e(tolerance)
	local alllatents = e(lvs)
	local allformative = e(formative)
	local num_lv : word count `alllatents'
	local num_lv_B : word count `allformative'
	local num_lv_A = `num_lv' - `num_lv_B'

	display

	if ("`consistent'" == "") {
		local header "Partial least squares path modeling"
	}
	else {
		local header "Consistent PLS path modeling (PLSc)"
	}
	
	if ((`rebus_it' == -999) & (`fimix_it' == -999) & ("`gas'" == "")) {
		if ("`group'" == "nogroup") {
			if ("`structural'" == "nostructural") {
				local nobs: display _skip(0) "`header'" _col(49) ///
					"Number of obs" _col(69) "=" _skip(5) string(e(N))
				display as text "`nobs'"
				
				display

				local init_head: display _skip(0) "Initialization: `init'"
				display as text "`init_head'"
			}
			else {
				if ("`rawsum'" == "") {
					local niter = colsof(`matrix2')
					forvalues i = 1/`niter' {
						local iter_lbl = "Iteration " + string(`i') + ///
							":   outer weights rel. diff. = " + string(`matrix2'[1, `i'], "%9.`digits'e")
						display as text _skip(0) "`iter_lbl'"
					}
					
					display
				}
				
				local nobs: display _skip(0) "`header'" _col(49) ///
					"Number of obs" _col(69) "=" _skip(5) string(e(N))
				display as text "`nobs'"

				if (`num_lv_A' > 0) {
					local avrsq: display _skip(0) _col(49) "Average R-squared" _col(69) "=" ///
						_skip(5) string(`matrix1'[1, 1], "%9.`digits'f")
					display as text "`avrsq'"
				}
				else {
					display
				}

				if ("`rawsum'" == "") {
					local avave: display _skip(0) "Weighting scheme: `wscheme'" _col(49)
				}
				else {
					local avave: display _skip(0) "Weighting scheme: rawsum" _col(49)
				}
				if (`num_lv_A' > 0) {
					local avave: display "`avave'" "Average communality" _col(69) "=" ///
					_skip(5) string(`matrix1'[1, 2], "%9.`digits'f")
				}
				display as text "`avave'"
				
				if ("`rawsum'" == "") {
					local gof: display _skip(0) "Tolerance: " string(`tol', "%9.`digits'e") ///
						_col(49)
				}
				else {
					local gof: display _skip(0) _col(49)
				}
				if (`num_lv_A' > 0) {
					local gof: display "`gof'" "GoF" _col(69) "=" _skip(5) ///
						string(`matrix1'[1, 3], "%9.`digits'f")
				}
				display as text "`gof'"
				
				if ("`rawsum'" == "") {
					local avred: display _skip(0) "Initialization: `init'" _col(49)
				}
				else {
					local avred: display _skip(0) _col(49)
				}
				if (`num_lv_A' > 0) {
					local avred: display "`avred'" "Average redundancy" _col(69) "=" _skip(5) ///
						string(`matrix1'[1, 4], "%9.`digits'f")
				}
				display as text "`avred'"
			}
		}
		else {
			local title: display _skip(0) "`header'"
			display as text "`title'"
			
			display
			
			if ("`structural'" == "nostructural") {
				local initialize: display _skip(0) "Initialization: `init'"
				display as text "`initialize'"
			}
			else {
				if ("`rawsum'" == "") {
					local wgt: display _skip(0) "Weighting scheme: `wscheme'"
					display as text "`wgt'"
					
					local toler: display _skip(0) "Tolerance: " string(`tol', "%9.`digits'e")
					display as text "`toler'"
				}
				else {
					local wgt: display _skip(0) "Weighting scheme: rawsum"
					display as text "`wgt'"
				}
				
				if ("`rawsum'" == "") {
					local initialize: display _skip(0) "Initialization: `init'"
				}
				else {
					local initialize: display _skip(0)
				}
				display as text "`initialize'"
			}
		}
	}
	else if ((`rebus_it' > 0) & (`fimix_it' == -999) & ("`gas'" == "")) {
		local title: display _skip(0) "Response-based unit segmentation partial least squares (REBUS-PLS)"
		display as text "`title'"
		
		display

		local wgt: display _skip(0) "Weighting scheme: `wscheme'"
		display as text "`wgt'"
		
		local toler: display _skip(0) "Tolerance: " string(`tol', "%9.`digits'e")
		display as text "`toler'"

		local initialize: display _skip(0) "Initialization: `init'"
		display as text "`initialize'"
		
		local rebit: display _skip(0) "Number of REBUS-PLS iterations: `rebus_it'"
		display as text "`rebit'"
		
		local rebgqi: display _skip(0) "Group Quality Index (GQI): " ///
			string(`rebus_gqi', "%9.`digits'f")
		display as text "`rebgqi'"
	}
	else if ((`rebus_it' == -999) & (`fimix_it' > 0) & ("`gas'" == "")) {
		local title: display _skip(0) "Finite mixture partial least squares (FIMIX-PLS)"
		display as text "`title'"
		
		display

		local wgt: display _skip(0) "Weighting scheme: `wscheme'"
		display as text "`wgt'"
		
		local toler: display _skip(0) "Tolerance: " string(`tol', "%9.`digits'e")
		display as text "`toler'"

		local initialize: display _skip(0) "Initialization: `init'"
		display as text "`initialize'"
		
		local fimit: display _skip(0) "Number of FIMIX-PLS iterations: `fimix_it'"
		display as text "`fimit'"
		
		local fimll: display _skip(0) "Log-likelihood value: " ///
			string(`fimix_ll', "%99.`digits'f")
		display as text "`fimll'"
	}
	if ("`gas'" != "") {
		local title: display _skip(0) "Partial least squares genetic algorithm segmentation (PLS-GAS)"
		display as text "`title'"
		
		display

		local wgt: display _skip(0) "Weighting scheme: `wscheme'"
		display as text "`wgt'"
		
		local toler: display _skip(0) "Tolerance: " string(`tol', "%9.`digits'e")
		display as text "`toler'"

		local initialize: display _skip(0) "Initialization: `init'"
		display as text "`initialize'"
	}
end
