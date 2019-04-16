*!mktable_indirect version 0.1.0
*!Written 01Oct2017
*!Written by Sergio Venturini and Mehmet Mehmetoglu
*!The following code is distributed under GNU General Public License version 3 (GPL-3)

program mktable_indirect
	version 15.1
	syntax , matrix(string) [ FIRSTCOLName(string) FIRSTCOLWidth(integer 25) ///
		COLWidth(integer 15) Title(string) HLines(numlist >0 integer sort) ///
		NOVLines DIGits(integer 3) Boot(numlist min=1 max=1) Level(real 0.95) ]

	/* Options:
	   --------
		 matrix(string)							--> matrix containing the numbers to display
		 firstcolname(string)				--> first column's name
		 firstcolwidth(integer 25)	--> first column's width (default 25)
		 colwidth(integer 15)				--> other columns' width (default 15)
		 title(string)							--> table's main title
		 hlines											--> rows at which to insert an horizontal line
		 novlines										--> no vertical lines between columns
		 digits(integer 3)					--> number of digits to display (default 3)
		 boot												--> flag for bootstrap
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
	if ("`boot'" == "") {
		local matrownames "Indirect effect" "_" "Standard error" "_" ///
			"Z statistic" "_" "P-value" "_" "Conf. interval"
	}
	else {
		local matrownames "Indirect effect" "_" "Standard error" "_" ///
			"Z statistic" "_" "P-value" "_" "Conf. interval (N)" "_" ///
			"Conf. interval (P)" "_" "Conf. interval (BC)"
	}
	if ("`boot'" != "") {
		local tabrows = `tabrows' + 2
	}
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
	local title1: display _col(`firstcolwidth_p1') "{c |}"
	local todisp "`firstcolname'"
	if (strlen("`todisp'") > `firstcolwidth' - 2*`skip1') {
		local todisp = abbrev("`todisp'", `firstcolwidth' - 2*`skip1')
	}
	local tmp_skip = `firstcolwidth' - strlen("`todisp'") - 2*`skip1'
	local title2: display _skip(`skip1') "`firstcolname'" _skip(`tmp_skip') _col(`firstcolwidth_p1') "{c |}"
	local title3: display _col(`firstcolwidth_p1') "{c |}"
	local secondline: display "{hline ""`firstcolwidth'""}{c +}"
	local lastline: display "{hline ""`firstcolwidth'""}{c BT}"

	forvalues j = 1/`ncols_m1' {
		local firstline "`firstline'" "`: display "{hline ""`colwidth'""}`tvlines'"'"
		local colname_j : word `j' of `matcolnames'
		tokenize `"`colname_j'"', parse("_")
		local dep `1'
		macro shift
		local mod `2'
		macro shift
		local indep `3'
		local todisp "`dep'"
		if (strlen("`todisp'") > `usable' - 3) {
			local todisp = abbrev("`todisp'", `usable' - 3)
		}
		local tmp_skip = `usable' - strlen("`todisp'") - 3
		local title1 "`title1'" "`: display _skip(`skip1') _skip(`tmp_skip')'" "`todisp' <-" "`vlines'"
		local todisp "`mod'"
		if (strlen("`todisp'") > `usable' - 3) {
			local todisp = abbrev("`todisp'", `usable' - 3)
		}
		local tmp_skip = `usable' - strlen("`todisp'") - 3
		local title2 "`title2'" "`: display _skip(`skip1') _skip(`tmp_skip')'" "`todisp' <-" "`vlines'"
		local todisp "`indep'"
		if (strlen("`todisp'") > `usable') {
			local todisp = abbrev("`todisp'", `usable')
		}
		local tmp_skip = `usable' - strlen("`todisp'")
		local title3 "`title3'" "`: display _skip(`skip1') _skip(`tmp_skip')'" "`todisp'" "`vlines'"
		local secondline "`secondline'" "`: display "{hline ""`colwidth'""}`mvlines'"'"
		local lastline "`lastline'" "`: display "{hline ""`colwidth'""}`bvlines'"'"
	}
	local firstline "`firstline'" "`: display "{hline ""`colwidth'""}"'"
	local colname_j : word `ncols' of `matcolnames'
	tokenize `"`colname_j'"', parse("_")
	local dep `1'
	macro shift
	local mod `2'
	macro shift
	local indep `3'
	local todisp "`dep'"
	if (strlen("`todisp'") > `usable' - 3) {
		local todisp = abbrev("`todisp'", `usable' - 3)
	}
	local tmp_skip = `usable' - strlen("`todisp'") - 3
	local title1 "`title1'" "`: display _skip(`skip1') _skip(`tmp_skip')'" "`todisp' <-"
	local todisp "`mod'"
	if (strlen("`todisp'") > `usable' - 3) {
		local todisp = abbrev("`todisp'", `usable' - 3)
	}
	local tmp_skip = `usable' - strlen("`todisp'") - 3
	local title2 "`title2'" "`: display _skip(`skip1') _skip(`tmp_skip')'" "`todisp' <-"
	local todisp "`indep'"
	if (strlen("`todisp'") > `usable') {
		local todisp = abbrev("`todisp'", `usable')
	}
	local tmp_skip = `usable' - strlen("`todisp'")
	local title3 "`title3'" "`: display _skip(`skip1') _skip(`tmp_skip')'" "`todisp'"
	local secondline "`secondline'" "`: display "{hline ""`colwidth'""}"'"
	local lastline "`lastline'" "`: display "{hline ""`colwidth'""}"'"
	
	display
	if ("`boot'" == "") {
		display as text _skip(0) "`title'"
	}
	else {
		display as text _skip(0) "`title'" " (Bootstrap)"
	}
	display as text "`firstline'"
	display as text "`title1'"
	display as text "`title2'"
	display as text "`title3'"
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
	if ("`boot'" != "") {
		display as text _skip(`skip1') "(N)    normal confidence interval"
		display as text _skip(`skip1') "(P)    percentile confidence interval"
		display as text _skip(`skip1') "(BC)   bias-corrected confidence interval"
	}
end
