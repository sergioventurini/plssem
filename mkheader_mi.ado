*!mkheader_mi version 0.1.0
*!Written 22Apr2023
*!Written by Sergio Venturini and Mehmet Mehmetoglu
*!The following code is distributed under GNU General Public License version 3 (GPL-3)

program mkheader_mi
	version 15.1
	syntax [, Wscheme(string) Init(string) DIGits(integer 5) noSTRuctural ///
    RAWsum CONsistent ]

	/* Options:
	   --------
		 digits(integer 5)					--> number of digits to display (default 5)
		 nostructural								--> indicator of whether the model includes a
																		structural part
		 rawsum											--> scores are raw sum of indicators
		 consistent									--> indicator for consistent PLS (PLSc)
	 */

	local tol = e(tolerance)

	display

	display as text "Multiple-imputation estimates"

	// display

  if ("`consistent'" == "") {
		local header "Partial least squares SEM"
	}
	else {
		local header "Consistent partial least squares (PLSc) SEM"
	}
	
  local nobs: display _skip(0) "`header'" _col(49) ///
    "Number of obs" _col(69) "=" _skip(5) //string(e(N))
  display as text "`nobs'" _continue
  display as result e(N)
	if ("`structural'" == "nostructural") {
		local init_head: display _skip(0) "Initialization: "
    display as text "`init_head'" _continue
    display as result "`init'" _continue
    local init_head: display _skip(0) "`init_head'" _col(49) ///
      "Imputations" _col(69) "=" _skip(5) //string(e(M_mi))
		display as text "`init_head'" _continue
    display as result e(M_mi)
	}
	else {
		if ("`rawsum'" == "") {
			local gof: display _skip(0) "Weighting scheme: " 
      display as text "`gof'" _continue
      display as result "`wscheme'" _col(49) _continue
		}
		else {
			local gof: display _skip(0) "Weighting scheme: " _col(49)
      display as text "`gof'" _continue
      display as result "rawsum" _col(49) _continue
		}
    local gof: display _skip(0) _col(49) ///
      "Imputations" _col(69) "=" _skip(5) //string(e(M_mi))
		display as text "`gof'" _continue
    display as result e(M_mi)
		
		if ("`rawsum'" == "") {
			local gof_rel: display _skip(0) "Tolerance: " //string(`tol', "%9.`digits'e")
      display as text "`gof_rel'" _continue
			display as result string(`tol', "%9.`digits'e") _col(49)
		}
		else {
			local gof_rel: display _skip(0) _col(49)
      display as text "`gof_rel'"
		}
		
		if ("`rawsum'" == "") {
			local avred: display _skip(0) "Initialization: "
      display as text "`avred'" _continue
      display as result "`init'"
		}
		else {
			local avred: display _skip(0)
      display as text "`avred'"
		}
	}
end
