*!mktable_mediate version 0.1.0
*!Written 10Oct2018
*!Written by Sergio Venturini and Mehmet Mehmetoglu
*!The following code is distributed under GNU General Public License version 3 (GPL-3)

program mktable_mediate
	version 15.1
	syntax , matrix(string) matrix_bk(string) [ matrix_zlc(string) ///
		matrix_rit(string) matrix_rid(string) FIRSTCOLName(string) ///
		depv(string) medv(string) indepv(string) FIRSTCOLWidth(integer 25) ///
		COLWidth(integer 15) Title(string) HLines(numlist >0 integer sort) ///
		NOVLines Reps(integer 50) DIGits(integer 3) Level(real 0.95) ]

	/* Options:
	   --------
		 matrix(string)							--> matrix containing the numbers to display
		 matrix_bk(string)					--> matrix containing the baron-kenny results
		 matrix_zlc(string)					--> matrix containing the zhao et al. results
		 matrix_rit(string)					--> matrix containing the RIT results
		 matrix_rid(string)					--> matrix containing the RID results
		 depv(string)								--> dependent variable
		 medv(string)								--> mediation variable
		 indepv(string)							--> independent variable
		 firstcolname(string)				--> first column's name
		 firstcolwidth(integer 25)	--> first column's width (default 25)
		 colwidth(integer 15)				--> other columns' width (default 15)
		 title(string)							--> table's main title
		 hlines											--> rows at which to insert an horizontal line
		 novlines										--> no vertical lines between columns
		 reps(integer 50)						--> number of bootstrap replications
		 digits(integer 3)					--> number of digits to display (default 3)
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
		exit
	}

	local firstcolwidth_p1 = `firstcolwidth' + 1
	local ncols = colsof(`matrix')
	local nrows = rowsof(`matrix')
	local tabrows = 5
	local matcolnames : colnames `matrix'
	local matrownames "Indirect effect" "_" "Standard error" "_" ///
		"Z statistic" "_" "P-value" "_" "Confidence interval"
	local ncols_m1 = `ncols' - 1
	local usable = `colwidth' - 2*`skip1'
	local usable2 = floor(`usable'/2 - 1)
	if (`digits' >= `usable') {
		display as error "the number of digits chosen is too large"
		exit
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
	local todisp "`firstcolname'"
	if (strlen("`todisp'") > `firstcolwidth' - 2*`skip1') {
		local todisp = abbrev("`todisp'", `firstcolwidth' - 2*`skip1')
	}
	local tmp_skip = `firstcolwidth' - strlen("`todisp'") - 2*`skip1'
	local title2: display _skip(`skip1') "`firstcolname'" _skip(`tmp_skip') _col(`firstcolwidth_p1') "{c |}"
	local secondline: display "{hline ""`firstcolwidth'""}{c +}"
	local lastline: display "{hline ""`firstcolwidth'""}{c BT}"

	forvalues j = 1/`ncols_m1' {
		local firstline "`firstline'" "`: display "{hline ""`colwidth'""}`tvlines'"'"
		local todisp : word `j' of `matcolnames'
		if (strlen("`todisp'") > `usable') {
			local todisp = abbrev("`todisp'", `usable')
		}
		local tmp_skip = `usable' - strlen("`todisp'")
		local title2 "`title2'" "`: display _skip(`skip1') _skip(`tmp_skip')'" "`todisp'" "`vlines'"
		local secondline "`secondline'" "`: display "{hline ""`colwidth'""}`mvlines'"'"
		local lastline "`lastline'" "`: display "{hline ""`colwidth'""}`bvlines'"'"
	}
	local firstline "`firstline'" "`: display "{hline ""`colwidth'""}"'"
	local todisp : word `ncols' of `matcolnames'
	if (strlen("`todisp'") > `usable') {
		local todisp = abbrev("`todisp'", `usable')
	}
	local tmp_skip = `usable' - strlen("`todisp'")
	local title2 "`title2'" "`: display _skip(`skip1') _skip(`tmp_skip')'" "`todisp'"
	local secondline "`secondline'" "`: display "{hline ""`colwidth'""}"'"
	local lastline "`lastline'" "`: display "{hline ""`colwidth'""}"'"
	
	display
	display as text _skip(0) "`title'"
	display as text "`firstline'"
	display as text "`title2'"
	display as text "`secondline'"
	
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
		display as text "`rownametodisp'" _continue

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
					display as text "`secondline'"
				}
				else if (`num' == `tabrows') {
					display as text "`lastline'"
				}
			}
			continue, break
		}
	}
	
	if (`level' <= 0 | `level' >= 1) {
		display as error "confidence level must be in the range (0, 1)"
		exit
	}
	local conflevel = strtrim(string(`level'*100, "%9.2g"))
	display as text _skip(`skip1') "confidence level: `conflevel'%"
	display as text _skip(`skip1') "bootstrap replications: `reps'"

	local doi `depv':`indepv'
	local moi `medv':`indepv'
	local dom `depv':`medv'

	/* Baron-Kenny mediation testing adjusted by Iacobucci et al. */
	display
	display as text " Baron & Kenny approach to testing mediation"
	if (`matrix_bk'[5, 1] > 0.05 | `matrix_bk'[6, 1] > 0.05) {
		display as result " STEP 1 - `moi' (X -> M) with b = " %-5.3f `matrix_bk'[2, 1] " and p = " %-5.3f `matrix_bk'[5, 1]
		display as result " STEP 2 - `dom' (M -> Y) with b = " %-5.3f `matrix_bk'[3, 1] " and p = " %-5.3f `matrix_bk'[6, 1]
		display as text "          As either STEP 1 or STEP 2 (or both) are not significant,"
		display as text "          there is no mediation"
	}
	else {
		if (`matrix_bk'[5, 1] < 0.05 & `matrix_bk'[6, 1] < 0.05 & `matrix_bk'[7, 1] < 0.05 & `matrix_bk'[4, 1] > 0.05) {
			display as result " STEP 1 - `moi' (X -> M) with b = " %-5.3f `matrix_bk'[2, 1] " and p = " %-5.3f `matrix_bk'[5, 1]
			display as result " STEP 2 - `dom' (M -> Y) with b = " %-5.3f `matrix_bk'[3, 1] " and p = " %-5.3f `matrix_bk'[6, 1]
			display as result " STEP 3 - `doi' (X -> Y) with b = " %-5.3f `matrix_bk'[1, 1] " and p = " %-5.3f `matrix_bk'[4, 1]
			display as text "          As STEP 1, STEP 2 and the Sobel's test above are significant "
			display as text            "          and STEP 3 is not significant the mediation is complete"
		}
		else {
			if (`matrix_bk'[5, 1] < 0.05 & `matrix_bk'[6, 1] < 0.05 & `matrix_bk'[7, 1] < 0.05 & `matrix_bk'[4, 1] < 0.05) {
				display as result " STEP 1 - `moi' (X -> M) with b = " %-5.3f `matrix_bk'[2, 1] " and p = " %-5.3f `matrix_bk'[5, 1]
				display as result " STEP 2 - `dom' (M -> Y) with b = " %-5.3f `matrix_bk'[3, 1] " and p = " %-5.3f `matrix_bk'[6, 1]
				display as result " STEP 3 - `doi' (X -> Y) with b = " %-5.3f `matrix_bk'[1, 1] " and p = " %-5.3f `matrix_bk'[4, 1]
				display as text "          As STEP 1, STEP 2 and STEP 3 as well as the Sobel's test above"
				display as text            "          are significant the mediation is partial"
			}
			else {
				if (`matrix_bk'[5, 1] < 0.05 & `matrix_bk'[6, 1] < 0.05 & `matrix_bk'[7, 1] > 0.05 & `matrix_bk'[4, 1] < 0.05) {
					display as result " STEP 1 - `moi' (X -> M) with b = " %-5.3f `matrix_bk'[2, 1] " and p = " %-5.3f `matrix_bk'[5, 1]
					display as result " STEP 2 - `dom' (M -> Y) with b = " %-5.3f `matrix_bk'[3, 1] " and p = " %-5.3f `matrix_bk'[6, 1]
					display as result " STEP 3 - `doi' (X -> Y) with b = " %-5.3f `matrix_bk'[1, 1] " and p = " %-5.3f `matrix_bk'[4, 1]
					display as text "          As STEP 1, STEP 2 and STEP 3 are all significant and the"
					display as text            "          Sobel's test above is not significant the mediation is partial"
				}
				else {
					if (`matrix_bk'[5, 1] < 0.05 & `matrix_bk'[6, 1] < 0.05 & `matrix_bk'[7, 1] > 0.05 & `matrix_bk'[4, 1] > 0.05) {
						display as result " STEP 1 - `moi' (X -> M) with b = " %-5.3f `matrix_bk'[2, 1] " and p = " %-5.3f `matrix_bk'[5, 1]
						display as result " STEP 2 - `dom' (M -> Y) with b = " %-5.3f `matrix_bk'[3, 1] " and p = " %-5.3f `matrix_bk'[6, 1]
						display as result " STEP 3 - `doi' (X -> Y) with b = " %-5.3f `matrix_bk'[1, 1] " and p = " %-5.3f `matrix_bk'[4, 1]
						display as text "          As STEP 1 and STEP 2 are significant and neither STEP 3 nor"
						display as text            "          the Sobel's test above is significant the mediation is partial"
					}
				}
			}
		}
	}

	/* Zhao et al. mediation testing */
	if (!matmissing("`matrix_zlc'")) {
		display
		display as text " Zhao, Lynch & Chen's approach to testing mediation"
		if (`matrix_zlc'[1, 1] < 0.05 & `matrix_zlc'[2, 1] > 0.05) {
			display as result " STEP 1 - `doi' (X -> Y) with b = " %-5.3f `matrix_zlc'[3, 1] " and p = " %-5.3f `matrix_zlc'[2, 1]
			display as text "          As the bootstrap test above is significant and STEP 1 is not" 
			display as text "          significant you have indirect-only mediation (full mediation)"
		}
		else {
			if (`matrix_zlc'[1, 1] > 0.05 & `matrix_zlc'[2, 1] < 0.05) {
				display as result " STEP 1 - `doi' (X -> Y) with b = " %-5.3f `matrix_zlc'[3, 1] " and p = " %-5.3f `matrix_zlc'[2, 1] 
				display as text "          As the bootstrap test above is not significant and STEP 1 is" 
				display as text "          significant you have direct-only nonmediation (no mediation)"
			}
			else {
				if (`matrix_zlc'[1, 1] > 0.05 & `matrix_zlc'[2, 1] > 0.05) {
					display as result " STEP 1 - `doi' (X -> Y) with b = " %-5.3f `matrix_zlc'[3, 1] " and p = " %-5.3f `matrix_zlc'[2, 1] 
					display as text "          As the bootstrap test above is not significant and STEP 1 is" 
					display as text "          not significant you have no effect nonmediation (no mediation)"
				}
				else {
					if (`matrix_zlc'[1, 1] < 0.05 & `matrix_zlc'[2, 1] < 0.05 & `matrix_zlc'[4, 1] > 0) {
						display as result " STEP 1 - `doi' (X -> Y) with b = " %-5.3f `matrix_zlc'[3, 1] " and p = " %-5.3f `matrix_zlc'[2, 1] 
						display as text "          As the bootstrap test above is significant, STEP 1 is" 
						display as text "          significant and their coefficients point in same direction,"
						display as text "          you have complementary mediation (partial mediation)" 
					}
					else {
						if (`matrix_zlc'[1, 1] < 0.05 & `matrix_zlc'[2, 1] < 0.05 & `matrix_zlc'[4, 1] < 0) {
							display as result " STEP 1 - `doi' (X -> Y) with b = " %-5.3f `matrix_zlc'[3, 1] " and p = " %-5.3f `matrix_zlc'[2, 1] 
							display as text "          As the bootstrap test above is significant, STEP 1 is" 
							display as text "          significant and their coefficients point in opposite"
							display as text "          direction, you have competitive mediation (partial mediation)" 
						}
					}
				}
			}
		}
	}
	
	if (!matmissing("`matrix_rit'")) {
		display
		display as text " RIT  =   (Indirect effect / Total effect)"
		display as result "          (" %-5.3f `matrix_rit'[1, 1] " / " %-5.3f `matrix_rit'[2, 1] ") = " %-5.3f `matrix_rit'[3, 1]
		if (`matrix_rit'[3, 1] < 1) {
			local fmttodisp = "%2.1f"
		}
		else {
			local fmttodisp = "%3.1f"
		}
		display as text "          Meaning that about " `fmttodisp' `matrix_rit'[3, 1]*100 "% of the effect of `indepv'"
		display as text "          on `depv' is mediated by `medv'"
	}
	
	if (!matmissing("`matrix_rid'")) {
		display
		display as text " RID  =   (Indirect effect / Direct effect)"
		display as result "          (" %-5.3f `matrix_rid'[1, 1] " / " %-5.3f `matrix_rid'[2, 1] ") = " %-5.3f `matrix_rid'[3, 1]
		local fmttodisp = "%5.3f"
		display as text "          That is, the mediated effect is about " `fmttodisp' `matrix_rid'[3, 1] " times as"
		display as text "          large as the direct effect of `indepv' on `depv'" 
	}
	
end
