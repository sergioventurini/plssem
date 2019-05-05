*!mktable version 0.2.0
*!Written 05May2018
*!Written by Sergio Venturini and Mehmet Mehmetoglu
*!The following code is distributed under GNU General Public License version 3 (GPL-3)

program mktable
	version 15.1
	syntax , matrix(string) [ FIRSTCOLName(string) FIRSTCOLWidth(integer 25) ///
		COLWidth(integer 15) Title(string) HLines(numlist >0 integer sort) ///
		NOVLines DIGits(integer 3) Path STats TOTal CORr CUToff(real 0) ///
		BINary(namelist min=1) REBus FIMix CONsistent ]

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
		 path												--> indicate that the table is for the path
																		coefficients
		 stats											--> table of summary statistics for indicators
		 total											--> table of total effects
		 corr												--> correlation table
		 cutoff(real 0)							--> do not show correlation smaller than cutoff
		 binary(namelist min=1)			--> binary indicators/latent variables
		 rebus											--> indicator that the table refers to REBUS-PLS
																		or FIMIX-PLS results
		 fimix											--> indicator that the table refers to FIMIX-PLS
																		results
		 consistent									--> indicator for consistent PLS (PLSc)
	 */
	
	local skip0 = 0
	local skip1 = 1
	local skip2 = 2
	local skip3 = 3
	local skip4 = 4

	if (strlen("`firstcolname'") > `firstcolwidth' - 2*`skip1') {
		local firstcolname = abbrev("`firstcolname'", `firstcolwidth' - 2*`skip1')
	}
	if (`colwidth' < 9) {
		display as error "'colwidth' option must be larger than 8 to properly show the table"
		exit
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
	local usable = `colwidth' - 2*`skip1'
	local usable2 = `colwidth' - 2*`skip2'
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
	
	if ("`corr'" != "") {
		local nr = rowsof(`matrix')
		local nc = colsof(`matrix')
		forvalues i = 1/`nr' {
			forvalues j = 1/`nc' {
				if (!missing(`matrix'[`i', `j']) & abs(`matrix'[`i', `j']) < `cutoff') {
					matrix `matrix'[`i', `j'] = .
				}
			}
		}
	}

	forvalues j = 1/`ncols_m1' {
		local firstline "`firstline'" "`: display "{hline ""`colwidth'""}`tvlines'"'"
		if (`nulleqnames' < `ncols') {
			if ("`: word `j' of `mateqnames''" != "_") {
				local todisp : word `j' of `mateqnames'
				if (strlen("`todisp'") >= `usable') {
					local todisp = abbrev("`todisp'", `usable' - 1)
				}
				local tmp_skip = `usable' - strlen("`todisp'") - 1
				local title1_eq "`title1_eq'" "`: display _skip(`skip1') _skip(`tmp_skip')'" "`todisp'" ":`vlines'"
			}
		}
		local todisp : word `j' of `matcolnames'
		if (strpos("`todisp'", "_")) {
			//local todisp = subinstr("`todisp'", "_", " ", .)
		}
		if (strlen("`todisp'") > `usable') {
			local todisp = abbrev("`todisp'", `usable')
		}
		local tmp_skip = `usable' - strlen("`todisp'")
		local title1 "`title1'" "`: display _skip(`skip1') _skip(`tmp_skip')'" abbrev("`todisp'", `colwidth') "`vlines'"
		local secondline "`secondline'" "`: display "{hline ""`colwidth'""}`mvlines'"'"
		local lastline "`lastline'" "`: display "{hline ""`colwidth'""}`bvlines'"'"
	}
	local firstline "`firstline'" "`: display "{hline ""`colwidth'""}"'"
	if (`nulleqnames' < `ncols') {
		if ("`: word `ncols' of `mateqnames''" != "_") {
			local todisp : word `ncols' of `mateqnames'
			if (strlen("`todisp'") >= `usable') {
				local todisp = abbrev("`todisp'", `usable' - 1)
			}
			local tmp_skip = `usable' - strlen("`todisp'") - 1
			local title1_eq "`title1_eq'" "`: display _skip(`skip1') _skip(`tmp_skip')'" "`todisp'" ":"
		}
	}
	local todisp : word `ncols' of `matcolnames'
	if (strpos("`todisp'", "_")) {
		//local todisp = subinstr("`todisp'", "_", " ", .)
	}
	if (strlen("`todisp'") > `usable') {
		local todisp = abbrev("`todisp'", `usable')
	}
	local tmp_skip = `usable' - strlen("`todisp'")
	local title1 "`title1'" "`: display _skip(`skip1') _skip(`tmp_skip')'" abbrev("`todisp'", `colwidth') ""
	local secondline "`secondline'" "`: display "{hline ""`colwidth'""}"'"
	local lastline "`lastline'" "`: display "{hline ""`colwidth'""}"'"
	
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
			if ("`fimix'" != "") {
				local rownametodisp = subinstr("`rownametodisp'", "%%", " ", 1)
			}
			if ("`total'" != "") {
				local rownametodisp = subinstr("`rownametodisp'", "->", " -> ", 1)
			}
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
					if (`j' == `ncols_m1' & "`stats'" != "") {
						local todisp = strtrim(string(`matrix'[`i', `j'], "%`usable'.0f"))
						if (strlen("`todisp'") > `usable') {
							local todisp = strtrim(string(`matrix'[`i', `j'], "%`usable'.2e"))
						}
					}
					else {
						if ("`rebus'" != "") & (`i' == 1) {
							local todisp = strtrim(string(`matrix'[`i', `j'], "%`usable'.0f"))
						}
						else {
							local todisp = strtrim(string(`matrix'[`i', `j'], "%`usable'.`digits'f"))
						}
						if (strlen("`todisp'") > `usable') {
							local todisp = strtrim(string(`matrix'[`i', `j'], "%`usable'.2e"))
						}
					}
					local tmp_skip = `usable' - strlen("`todisp'")
					local rowtodisp "`rowtodisp'" "`: display _skip(`skip1') _skip(`tmp_skip') "`todisp'`vlines'"'"
				}
				else {
					if ("`rebus'" != "") & (`i' == 1) {
						local todisp = strtrim(string(`matrix'[`i', `j'], "%`usable2'.0f"))
					}
					else {
						local todisp = strtrim(string(`matrix'[`i', `j'], "%`usable2'.`digits'f"))
					}
					if (strlen("`todisp'") > `usable2') {
						local todisp = strtrim(string(`matrix'[`i', `j'], "%`usable2'.2e"))
					}
					local todisp = "(`todisp')"
					local tmp_skip = `usable' - strlen("`todisp'")
					local rowtodisp "`rowtodisp'" "`: display _skip(`skip1') _skip(`tmp_skip') "`todisp'`vlines'"'"
				}
			}
			else {
				local rowtodisp "`rowtodisp'" "`: display _skip(`skip1') _skip(`usable') "`vlines'"'"
			}
		}
		if (!missing(`matrix'[`i', `ncols'])) {
			if ("`rownametodisp'" != ".") {
				if ("`stats'" != "") {
					local todisp = strtrim(string(`matrix'[`i', `ncols'], "%`usable'.0f"))
					if (strlen("`todisp'") > `usable') {
						local todisp = strtrim(string(`matrix'[`i', `ncols'], "%`usable'.2e"))
					}
				}
				else {
					if ("`rebus'" != "") & (`i' == 1) {
						local todisp = strtrim(string(`matrix'[`i', `ncols'], "%`usable'.0f"))
					}
					else {
						local todisp = strtrim(string(`matrix'[`i', `ncols'], "%`usable'.`digits'f"))
					}
					if (strlen("`todisp'") > `usable') {
						local todisp = strtrim(string(`matrix'[`i', `ncols'], "%`usable'.2e"))
					}
				}
				local tmp_skip = `usable' - strlen("`todisp'")
				local rowtodisp "`rowtodisp'" "`: display _skip(`skip1') _skip(`tmp_skip') "`todisp'"'"
			}
			else {
				if ("`rebus'" != "") & (`i' == 1) {
					local todisp = strtrim(string(`matrix'[`i', `ncols'], "%`usable2'.0f"))
				}
				else {
					local todisp = strtrim(string(`matrix'[`i', `ncols'], "%`usable2'.`digits'f"))
				}
				if (strlen("`todisp'") > `usable2') {
					local todisp = strtrim(string(`matrix'[`i', `ncols'], "%`usable2'.2e"))
				}
				local todisp = "(`todisp')"
				local tmp_skip = `usable' - strlen("`todisp'")
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
			}
			break
		}
	}
	
	if ("`path'" != "") {
		if (("`consistent'" != "") & (!`isboot')) {
			display as text _skip(`skip1') "p-values not shown (use the 'boot' option)"
		}
		else {
			display as text _skip(`skip1') "p-values in parentheses"
		}
	}

	if ("`binary'" != "") {
		local binary : list clean binary
		local numbin : word count `binary'
		local numbinm1 = `numbin' - 1
		local new_binary `: word 1 of `binary''
		if (`numbin' > 1) {
			forvalues i = 2/`numbinm1' {
				local tmpbinary `: word `i' of `binary''
				local new_binary "`new_binary', `tmpbinary'"
			}
			local tmpbinary `: word `numbin' of `binary''
			local new_binary "`new_binary' and `tmpbinary'"
		}
		display as text _skip(`skip1') "results for `new_binary' use a logit model"
	}
end
