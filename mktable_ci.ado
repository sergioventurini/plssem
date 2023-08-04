*!mktable_ci version 0.1.0
*!Written 28Jul2023
*!Written by Sergio Venturini and Mehmet Mehmetoglu
*!The following code is distributed under GNU General Public License version 3 (GPL-3)

program mktable_ci
	version 15.1
	syntax , matrix(string) [ FIRSTCOLName(string) FIRSTCOLWidth(integer 25) ///
		COLWidth(numlist >0 integer) Title(string) HLines(numlist >0 integer sort) ///
		NOVLines DIGits(integer 3) Path ]

	/* Options:
	   --------
		 matrix(string)							--> matrix containing the numbers to display
		 firstcolname(string)				--> first column's name
		 firstcolwidth(integer 25)	--> first column's width (default 25)
		 colwidth            				--> other columns' widths (potentially different);
                                    this is different w.r.t. mktable, which uses
                                    the same width for all columns
		 title(string)							--> table's main title
		 hlines											--> rows at which to insert an horizontal line
		 novlines										--> no vertical lines between columns
		 digits(integer 3)					--> number of digits to display (default 3)
		 path												--> indicate that the table is for the path
																		coefficients
	 */
	
	local skip0 = 0
	local skip1 = 1
	local skip2 = 2
	local skip3 = 3
	local skip4 = 4

	if (strlen("`firstcolname'") > `firstcolwidth' - 2*`skip1') {
		local firstcolname = abbrev("`firstcolname'", `firstcolwidth' - 2*`skip1')
	}

	local props = e(properties)
	local boot_lbl "bootstrap"
	local isboot : list boot_lbl in props

	local firstcolwidth_p1 = `firstcolwidth' + 1
	local ncols = colsof(`matrix')
	local nrows = rowsof(`matrix')
	local ncols_m1 = `ncols' - 1
	local matrownames : rownames `matrix'
	local matcolnames : colnames `matrix'
	local mateqnames : coleq `matrix'
	local nulleqnames = 0
	forvalues s = 1/`ncols' {
		if ("`: word `s' of `mateqnames''" == "_") {
			local ++nulleqnames
		}
	}
	if ("`colwidth'" != "") {
		numlist "`colwidth'"
		local colwidthtodisp = r(numlist)
	}
  local usable
  local usable2
  foreach num in `colwidthtodisp' {
    local tmpval = `num' - 2*`skip1'
    local usable "`usable' `tmpval'"
    local tmpval = `num' - 2*`skip2'
    local usable2 "`usable2' `tmpval'"
  }
	tempname maxusable
  mata: st_numscalar("`maxusable'", max(strtoreal(tokens(st_local("usable")))))
  if (`digits' >= `maxusable') {
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
	if (`nulleqnames' < `ncols') {
		local title1_eq: display _col(`firstcolwidth_p1') "{c |}"
	}
	local todisp "`firstcolname'"
	if (strlen("`todisp'") > `firstcolwidth' - 2*`skip1') {
		local todisp = abbrev("`todisp'", `firstcolwidth' - 2*`skip1')
	}
	local tmp_skip = `firstcolwidth' - strlen("`todisp'") - 2*`skip1'
	local title1: display _skip(`skip1') _skip(`tmp_skip') "`firstcolname'" ///
		_col(`firstcolwidth_p1') "{c |}"
	local secondline: display "{hline ""`firstcolwidth'""}{c +}"
	local lastline: display "{hline ""`firstcolwidth'""}{c BT}"

	forvalues j = 1/`ncols_m1' {
    local cw_tmp : word `j' of `colwidthtodisp'
		local firstline "`firstline'" "`: display "{hline ""`cw_tmp'""}`tvlines'"'"
    local usabl_tmp : word `j' of `usable'
		if (`nulleqnames' < `ncols') {
			if ("`: word `j' of `mateqnames''" != "_") {
				local todisp : word `j' of `mateqnames'
				if (strlen("`todisp'") >= `') {
					local todisp = abbrev("`todisp'", `usabl_tmp' - 1)
				}
				local tmp_skip = `usabl_tmp' - strlen("`todisp'") - 1
				local title1_eq "`title1_eq'" "`: display _skip(`skip1') _skip(`tmp_skip')'" "`todisp'" ":`vlines'"
			}
		}
		local todisp : word `j' of `matcolnames'
		if (strpos("`todisp'", "_")) {
			local todisp = subinstr("`todisp'", "_", " ", .)
		}
		if (strlen("`todisp'") > `usabl_tmp') {
			local todisp = abbrev("`todisp'", `usabl_tmp')
		}
		local tmp_skip = `usabl_tmp' - strlen("`todisp'")
		local title1 "`title1'" "`: display _skip(`skip1') _skip(`tmp_skip')'" abbrev("`todisp'", `cw_tmp') "`vlines'"
		local secondline "`secondline'" "`: display "{hline ""`cw_tmp'""}`mvlines'"'"
		local lastline "`lastline'" "`: display "{hline ""`cw_tmp'""}`bvlines'"'"
	}
  local last_cw : word count `colwidthtodisp'
  local cw_tmp : word `last_cw' of `colwidthtodisp'
  local usabl_tmp : word `last_cw' of `usable'
	local firstline "`firstline'" "`: display "{hline ""`cw_tmp'""}"'"
	if (`nulleqnames' < `ncols') {
		if ("`: word `ncols' of `mateqnames''" != "_") {
			local todisp : word `ncols' of `mateqnames'
			if (strlen("`todisp'") >= `usabl_tmp') {
				local todisp = abbrev("`todisp'", `usabl_tmp' - 1)
			}
			local tmp_skip = `usabl_tmp' - strlen("`todisp'") - 1
			local title1_eq "`title1_eq'" "`: display _skip(`skip1') _skip(`tmp_skip')'" "`todisp'" ":"
		}
	}
	local todisp : word `ncols' of `matcolnames'
	if (strpos("`todisp'", "_")) {
		local todisp = subinstr("`todisp'", "_", " ", .)
	}
	if (strlen("`todisp'") > `usabl_tmp') {
		local todisp = abbrev("`todisp'", `usabl_tmp')
	}
	local tmp_skip = `usabl_tmp' - strlen("`todisp'")
	local title1 "`title1'" "`: display _skip(`skip1') _skip(`tmp_skip')'" abbrev("`todisp'", `cw_tmp') ""
	local secondline "`secondline'" "`: display "{hline ""`cw_tmp'""}"'"
	local lastline "`lastline'" "`: display "{hline ""`cw_tmp'""}"'"

  /* Start printing in Results window */
	display
	display as text _skip(`skip0') "`title'"
	display as text "`firstline'"
	if (`nulleqnames' < `ncols') {
		display as text "`title1_eq'"
	}
	display as text "`title1'"
	display as text "`secondline'"
	
	forvalues i = 1/`nrows' {
		local rownametodisp : word `i' of `matrownames'
		if ("`rownametodisp'" != ".") {
			local rownametodisp = subinstr("`rownametodisp'", "->", " -> ", 1)
      local rownametodisp = subinstr("`rownametodisp'", "<-", " <- ", 1)
			if (strlen("`rownametodisp'") > (`firstcolwidth' - 2*`skip1')) {
				local rownametodisp = abbrev("`rownametodisp'", `firstcolwidth' - 2*`skip1')
			}
			local tmp_skip = `firstcolwidth' - strlen("`rownametodisp'") - 2*`skip1'
			local rownametodisp "`: display _skip(`skip1') _skip(`tmp_skip')'" "`rownametodisp'" "`: display _skip(`skip1') "{c |}"'"
		}
		else {
			local rownametodisp "`: display _col(`firstcolwidth_p1') "{c |}"'"
		}
		display as text "`rownametodisp'" _continue

		local rownametodisp : word `i' of `matrownames'
		local rowtodisp ""
		forvalues j = 1/`ncols_m1' {
			if (!missing(`matrix'[`i', `j'])) {
				if ("`rownametodisp'" != ".") {
          local usabl_tmp : word `j' of `usable'
					local todisp = strtrim(string(`matrix'[`i', `j'], "%`usabl_tmp'.`digits'f"))
					if (strlen("`todisp'") > `usabl_tmp') {
						local todisp = strtrim(string(`matrix'[`i', `j'], "%`usabl_tmp'.2e"))
					}
					local tmp_skip = `usabl_tmp' - strlen("`todisp'")
					local rowtodisp "`rowtodisp'" "`: display _skip(`skip1') _skip(`tmp_skip') "`todisp'`vlines'"'"
				}
				else {
          local usabl_tmp : word `j' of `usable2'
					local todisp = strtrim(string(`matrix'[`i', `j'], "%`usabl_tmp'.`digits'f"))
					if (strlen("`todisp'") > `usabl_tmp') {
						local todisp = strtrim(string(`matrix'[`i', `j'], "%`usabl_tmp'.2e"))
					}
					local todisp = "(`todisp')"
					local tmp_skip = `usabl_tmp' - strlen("`todisp'")
					local rowtodisp "`rowtodisp'" "`: display _skip(`skip1') _skip(`tmp_skip') "`todisp'`vlines'"'"
				}
			}
			else {
        local usabl_tmp : word `j' of `usable'
				local rowtodisp "`rowtodisp'" "`: display _skip(`skip1') _skip(`usabl_tmp') "`vlines'"'"
			}
		}
		if (!missing(`matrix'[`i', `ncols'])) {
			if ("`rownametodisp'" != ".") {
        local usabl_tmp : word `last_cw' of `usable'
				local todisp = strtrim(string(`matrix'[`i', `ncols'], "%`usabl_tmp'.`digits'f"))
				if (strlen("`todisp'") > `usabl_tmp') {
					local todisp = strtrim(string(`matrix'[`i', `ncols'], "%`usabl_tmp'.2e"))
				}
				local tmp_skip = `usabl_tmp' - strlen("`todisp'")
				local rowtodisp "`rowtodisp'" "`: display _skip(`skip1') _skip(`tmp_skip') "`todisp'"'"
			}
			else {
				local usabl_tmp : word `last_cw' of `usable2'
        local todisp = strtrim(string(`matrix'[`i', `ncols'], "%`usabl_tmp'.`digits'f"))
				if (strlen("`todisp'") > `usabl_tmp') {
					local todisp = strtrim(string(`matrix'[`i', `ncols'], "%`usabl_tmp'.2e"))
				}
				local todisp = "(`todisp')"
				local tmp_skip = `usabl_tmp' - strlen("`todisp'")
				local rowtodisp "`rowtodisp'" "`: display _skip(`skip1') _skip(`tmp_skip') "`todisp'"'"
			}
		}
		display as result "`rowtodisp'"
		
		foreach num in `hlinestodisp' {
			if (`num' == `i') {
				if (`num' < `nrows') {
					display as text "`secondline'"
				}
				else if (`num' == `nrows') {
					display as text "`lastline'"
				}
				break
			}
		}
	}
end
