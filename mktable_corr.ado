*!mktable_corr version 0.2.0
*!Written 06Feb2024
*!Written by Sergio Venturini and Mehmet Mehmetoglu
*!The following code is distributed under GNU General Public License version 3 (GPL-3)

program mktable_corr
  version 15.1
  syntax , matrix(string) [ Title(string) CUToff(real 0) DIGits(integer 4) ]

  /* Options:
     --------
     matrix(string)             --> matrix containing the numbers to display
     title(string)              --> table's main title
     cutoff(real 0)             --> do not show correlation smaller than cutoff
     digits(integer 4)          --> number of digits to display (default 4)
    */
  
  local skip0 = 0
  local len_cell = 5 + `digits'
  local len_num = 3 + `digits'

  display
  display as text _skip(`skip0') "`title'"

  local allvars : rownames `matrix'
  tokenize `allvars'
  local nvar : word count `allvars'
  local j0 = 1
  while (`j0' <= `nvar') {
    local j1 = min(`j0' + `len_cell', `nvar')
    local j = `j0'
    local l = `len_cell'*(`j1' - `j0' + 1)
    display as text "{hline 13}{c TT}{hline `l'}"
    display as text _skip(13) "{c |}" _continue
    while (`j' <= `j1') {
      display as text %`len_cell's abbrev("``j''", `len_cell' - 1) _continue
      local j = `j' + 1
    }
    display as text _newline "{hline 13}{c +}{hline `l'}"

    local i = `j0'
    while (`i' <= `nvar') {
      display as text %12s abbrev("``i''", 12) " {c |} " _continue
      local j = `j0'
      while (`j' <= min(`j1', `i')) {
        if (!missing(`matrix'[`i', `j']) & abs(`matrix'[`i', `j']) < `cutoff') {
          local c`j' = .
        }
        else {
          local c`j' = `matrix'[`i', `j']
        }
        local j = `j' + 1
      }
      local j = `j0'
      while (`j' <= min(`j1', `i')) {
        if (`c`j'' < .) {
          display as result " " %`len_num'.`digits'f `c`j'' " " _continue
        }
        else {
          display as result _skip(`len_cell') _continue
        }
        local j = `j' + 1
      }
      display
      local i = `i' + 1
    }
    local j0 = `j0' + 10
    display as text "{hline 13}{c BT}{hline `l'}"
  }
end
