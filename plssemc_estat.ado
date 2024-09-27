*!plssemc_estat version 0.6.0
*!Written 29Apr2024
*!Written by Sergio Venturini and Mehmet Mehmetoglu
*!The following code is distributed under GNU General Public License version 3 (GPL-3)

program plssemc_estat, rclass
  version 15.1
  gettoken subcmd rest : 0 , parse(", ")
  local lsubcmd = length("`subcmd'")
  local ismi = (e(cmd) == "mi estimate")

  if (`ismi') {
    display as error "plssemc postestimation is not available after multiple imputation"
    exit
  } 
  if ("`subcmd'" == substr("indirect", 1, max(2, `lsubcmd'))) {
    indirect_c `rest'
  }
  else if ("`subcmd'" == substr("total", 1, max(2, `lsubcmd'))) {
    total_c `rest'
  }
  else if ("`subcmd'" == substr("vif", 1, max(2, `lsubcmd'))) {
    plssem_vif_c `rest'
  }
  else if ("`subcmd'" == substr("mediate", 1, max(2, `lsubcmd'))) {
    mediate_c `rest'
  }
  else if ("`subcmd'" == substr("htmt", 1, max(2, `lsubcmd'))) {
    htmt_c `rest'
  }
  else if ("`subcmd'" == substr("f2", 1, max(2, `lsubcmd'))) {
    f2_c `rest'
  }
  else if ("`subcmd'" == substr("ic", 1, max(2, `lsubcmd'))) {
    ic_c `rest'
  }
  else if ("`subcmd'" == substr("dist", 1, max(2, `lsubcmd'))) {
    dist_c `rest'
  }
  else if ("`subcmd'" == substr("ci", 1, max(2, `lsubcmd'))) {
    ci_c `rest'
  }
  else if ("`subcmd'" == substr("blindfolding", 1, max(2, `lsubcmd'))) {
    blindf_c `rest'
  }
  else {
    // estat_default `0'
    display as error "the `subcmd' postestimation command is not implemented for plssemc"
    exit
  }

  return add
end

program indirect_c, rclass
  version 15.1
  syntax , Effects(string) [ Boot(numlist integer >0 min=1 max=1) ///
    Seed(numlist max=1) Level(real 0.95) DIGits(integer 3) ]
  
  /* Options:
     --------
     effects(string)            --> list of indirect effects
     boot(numlist integer >0 min=1 max=1)
                                --> bootstrap estimation (# of repetions;
                                    default 50)
     seed(numlist max=1)        --> bootstrap seed number
     level(real 0.95)           --> confidence level (default 0.95)
     digits(integer 3)          --> number of digits to display (default 3)
   */
   
   /* Description:
      ------------
      This postestimation command provides the estimates for the indirect
      effects mediated by only one LV.
   */
  
  display as error "estat indirect after plssemc will be available in a future version!"
  exit
  
  if ("`effects'" == "") {
    display as error "effects() option must be provided"
    exit
  }
  if (`digits' < 0) {
    display as error "number of digits must be a nonnegative integer"
    exit
  }
  local struct "structural"
  local props = e(properties)
  local isstruct : list struct in props
  if (!`isstruct') {
    display as error "the fitted plssem model includes only the measurement part"
    exit
  }
  
  tempvar __touse__
  quietly generate `__touse__' = e(sample)
  local reg3eqs = e(struct_eqs)
  if (`level' <= 0 | `level' >= 1) {
    display as error "confidence level must be in the range (0, 1)"
    exit
  }
  tempname normperc alpha_cl
  scalar `alpha_cl' = 1 - ((1 - `level')/2)
  scalar `normperc' = invnormal(`alpha_cl')
  
  /* Parse indirect statements (dep1 med1 indep1, dep2 med2 indep2, etc.) */
  local effects : list clean effects
  tokenize `"`effects'"', parse(",")
  
  local num_effects = 0
  local tok_i = 1
  while ("``tok_i''" != "") {
    local 0 ``tok_i''
    syntax varlist(min=3 max=3)
    
    local ind_dvar `ind_dvar' `: word 1 of `varlist''
    local ind_mvar `ind_mvar' `: word 2 of `varlist''
    local ind_ivar `ind_ivar' `: word 3 of `varlist''

    local ++num_effects

    local tok_i = `tok_i' + 2
  }
  if (`num_effects' > 5) {
    display as error "a maximum of 5 indirect effects are allowed"
    error 103
  }
  /* End parsing indirect statements */
  
  tempname sobel_se sobel_z sobel_pv indirect sobel_lci sobel_uci
  tempname indmeas reg3coef reg3coef_bs reg3var
  tempname dommat dommat2 moimat moimat2 //doimat doimat2
  if ("`boot'" == "") {
    matrix `indmeas' = J(6, `num_effects', .)
  }
  else {
    matrix `indmeas' = J(10, `num_effects', .)
  }
  
  tempname ehold
  _estimates hold `ehold'

  if ("`boot'" != "") {
    display
    display as text "Bootstrapping indirect effects..."
  }
  
  capture {
    forvalues ll = 1/`num_effects' {
      local dep `: word `ll' of `ind_dvar''
      local med `: word `ll' of `ind_mvar''
      local indep `: word `ll' of `ind_ivar''
      local indmeas_colnames "`indmeas_colnames' `dep'_`med'_`indep'"
      if ("`dep'" != "" & "`med'" != "" & "`indep'" != "") {
        local doi `dep':`indep'
        local moi `med':`indep'
        local dom `dep':`med'

        if ("`boot'" == "") {
          quietly reg3 `reg3eqs' if `__touse__', mvreg corr(independent)
          matrix `reg3coef' = e(b)
          matrix `reg3var' = e(V)
          local reg3eq_nm : colfullnames `reg3coef'

          tempname moicoef moivar domcoef domvar
          matrix `moimat' = `reg3coef'[1, "`moi'"]
          scalar `moicoef' = `moimat'[1, 1]
          matrix `moimat2' = `reg3var'["`moi'","`moi'"]
          scalar `moivar' = `moimat2'[1, 1]
          matrix `dommat' = `reg3coef'[1, "`dom'"]
          scalar `domcoef' = `dommat'[1, 1]
          matrix `dommat2' = `reg3var'["`dom'", "`dom'"]
          scalar `domvar' = `dommat2'[1, 1]

          scalar `indirect' = `moicoef'*`domcoef'
          scalar `sobel_se' = sqrt(((`domcoef')^2)*`moivar' + ((`moicoef')^2)*`domvar')
          scalar `sobel_z' = `indirect'/`sobel_se'
          scalar `sobel_pv' =  2*(1 - normal(abs(`sobel_z'))) 
          scalar `sobel_lci' = `indirect' - `normperc'*`sobel_se'
          scalar `sobel_uci' = `indirect' + `normperc'*`sobel_se'
          
          matrix `indmeas'[1, `ll'] = `indirect'
          matrix `indmeas'[2, `ll'] = `sobel_se'
          matrix `indmeas'[3, `ll'] = `sobel_z'
          matrix `indmeas'[4, `ll'] = `sobel_pv'
          matrix `indmeas'[5, `ll'] = `sobel_lci'
          matrix `indmeas'[6, `ll'] = `sobel_uci'
        }
        else {
          tempname reg3ciP reg3ciBC
          quietly bootstrap indeff=(_b[`dom']*_b[`moi']), ties reps(`boot') ///
            seed(`seed'): reg3 `reg3eqs' if `__touse__', mvreg corr(independent)
          quietly estat bootstrap, all
          matrix `reg3coef' = e(b)
          matrix `reg3coef_bs' = e(b_bs) // this is not used for now
          matrix `reg3var' = e(V)
          matrix `reg3ciP' = e(ci_percentile)
          matrix `reg3ciBC' = e(ci_bc)

          tempname iecoef ievar perc_lci perc_uci bc_lci bc_uci
          matrix `moimat' = `reg3coef'[1, 1]
          scalar `iecoef' = `moimat'[1, 1]
          matrix `moimat2' = `reg3var'[1, 1]
          scalar `ievar' = `moimat2'[1, 1]

          scalar `indirect' = `iecoef'
          scalar `sobel_se' = sqrt(`ievar')
          scalar `sobel_z' = `indirect'/`sobel_se'
          scalar `sobel_pv' =  2*(1 - normal(abs(`sobel_z'))) 
          scalar `sobel_lci' = `indirect' - `normperc'*`sobel_se'
          scalar `sobel_uci' = `indirect' + `normperc'*`sobel_se'
          scalar `perc_lci' = `reg3ciP'[1, 1]
          scalar `perc_uci' = `reg3ciP'[2, 1]
          scalar `bc_lci' = `reg3ciBC'[1, 1]
          scalar `bc_uci' = `reg3ciBC'[2, 1]
          
          matrix `indmeas'[1, `ll'] = `indirect'
          matrix `indmeas'[2, `ll'] = `sobel_se'
          matrix `indmeas'[3, `ll'] = `sobel_z'
          matrix `indmeas'[4, `ll'] = `sobel_pv'
          matrix `indmeas'[5, `ll'] = `sobel_lci'
          matrix `indmeas'[6, `ll'] = `sobel_uci'
          matrix `indmeas'[7, `ll'] = `perc_lci'
          matrix `indmeas'[8, `ll'] = `perc_uci'
          matrix `indmeas'[9, `ll'] = `bc_lci'
          matrix `indmeas'[10, `ll'] = `bc_uci'
        }
      }
      else {
        display as error "provide a single dependent, independent and mediator variable for each indirect effect"
        exit
      }
    }
    if ("`boot'" == "") {
      matrix rownames `indmeas' = "Indirect effect" "Standard error" ///
        "Z statistic" "P-value" "Lower CI" "Upper CI"
    }
    else {
      matrix rownames `indmeas' = "Indirect effect" "Standard error" ///
        "Z statistic" "P-value" "Lower CI (N)" "Upper CI (N)" ///
        "Lower CI (P)" "Upper CI (P)" "Lower CI (BC)" "Upper CI (BC)"
    }
    matrix colnames `indmeas' = `indmeas_colnames'
  } // end of -capture-
  local rc = _rc
  _estimates unhold `ehold'
  error `rc'
  
  /* Display results */
  local ind_line = 5
  if ("`boot'" != "") {
    local ind_line = `ind_line' + 2
  }
  mktable_indirect, matrix(`indmeas') digits(`digits') boot(`boot') ///
    title("Significance testing of (standardized) indirect effects") ///
    firstcolname("Statistics") firstcolwidth(24) colwidth(18) ///
    hlines(`ind_line') level(`level') //novlines

  /* Return values */
  return matrix indirect `indmeas'
end

program total_c, rclass
  version 15.1
  syntax [ , DIGits(integer 3) Plot ]
  
  /* Options:
     --------
     digits(integer 3)    --> number of digits to display (default 3)
     plot                 --> bar plot of the effects
   */
  
  local struct "structural"
  local props = e(properties)
  local isstruct : list struct in props
  if (!`isstruct') {
    display as error "the fitted plssem model includes only the measurement part"
    exit
  }
  if (`digits' < 0) {
    display as error "number of digits must be a nonnegative integer"
    exit
  }
    
  tempname B path direct indirect total
  matrix `B' = e(pathcoef)
  local nlv = colsof(`B')
  local nlv_m1 = `nlv' - 1
  matrix `indirect' = J(`nlv', `nlv', 0)
  if (`nlv' <= 2) {
    matrix `total' = `B'
  }
  else {
    matrix `path' = `B'
    forvalues k = 2/`nlv_m1' {
      matrix `path' = `path' * `B'
      matrix `indirect' = `indirect' + `path'
    }
  }
  matrix `total' = `B' + `indirect'
  
  tempname alleffects alllen maxlen alleffects_new
  local nall = `nlv'*(`nlv' - 1)
  matrix `alleffects' = J(`nall', 3, .)
  matrix `alllen' = J(`nall', 1, .)
  local lvnames : colnames `B'
  local i = 1
  foreach var_f in `lvnames' {
    local lvnames_rest : list lvnames - var_f
    foreach var_t in `lvnames_rest' {
      local allnames "`allnames' `var_f'->`var_t'"
      local rowi : list posof "`var_f'" in lvnames
      local coli : list posof "`var_t'" in lvnames
      if (`B'[`rowi', `coli'] == 0) {
        matrix `alleffects'[`i', 1] = .
      }
      else {
        matrix `alleffects'[`i', 1] = `B'[`rowi', `coli']
      }
      if (`indirect'[`rowi', `coli'] == 0) {
        matrix `alleffects'[`i', 2] = .
      }
      else {
        matrix `alleffects'[`i', 2] = `indirect'[`rowi', `coli']
      }
      if (`total'[`rowi', `coli'] == 0) {
        matrix `alleffects'[`i', 3] = .
      }
      else {
        matrix `alleffects'[`i', 3] = `total'[`rowi', `coli']
      }
      local ++i
    }
  }
  local firstnonmiss = 0
  forvalues i = 1/`nall' {
    if (!missing(`alleffects'[`i', 3])) {
      if (`firstnonmiss' == 0) {
        matrix `alleffects_new' = `alleffects'[`i', 1..3]
        local ++firstnonmiss
      }
      else {
        matrix `alleffects_new' = (`alleffects_new' \ `alleffects'[`i', 1..3])
      }
      local efftoleave : word `i' of `allnames'
      local allnames_new "`allnames_new' `efftoleave'"
    }
  }
  matrix rownames `alleffects_new' = `allnames_new'
  matrix colnames `alleffects_new' = "Direct" "Indirect" "Total"

  local nall = rowsof(`alleffects_new')
  forvalues i = 1/`nall' {
    local nm : word `i' of `allnames_new'
    local nm_len : strlen local nm
    matrix `alllen'[`i', 1] = `nm_len' + 2
    local nm_new = subinstr("`nm'", "->", " -> ", 1)
    local lblbar "`lblbar' `i' "`nm_new'""
  }
  mata: st_numscalar("`maxlen'", max(st_matrix("`alllen'")))
  local `maxlen' = `maxlen' + 2
  mktable, matrix(`alleffects_new') digits(`digits') firstcolname("Effect") ///
    title("Direct, Indirect (overall) and Total Effects") ///
    firstcolwidth(``maxlen'') colwidth(12) hlines(`nall') novlines total
  
  if ("`plot'" != "") {
    tempvar __id__ // __touse__
    // quietly generate `__touse__' = e(sample)
    quietly generate `__id__' = _n // if `__touse__'
    quietly svmat `alleffects_new', names(col)
    graph bar (asis) Direct Indirect, over(`__id__', relabel(`lblbar') ///
      label(angle(90) labsize(small))) stack nofill legend(position(12) ///
      region(style(none))) ylabel(, nogrid) yline(0, lwidth(thin)) ///
      title("Direct and Indirect (overall) Effects") scheme(sj)
    quietly drop Direct Indirect Total
  }

  /* Return values */
  return matrix total `alleffects_new'
end

program plssem_vif_c, rclass
  version 15.1
  syntax [ , DIGits(integer 3) ]
  
  /* Options:
     --------
     digits(integer 3)          --> number of digits to display (default 3)
   */
  
  if (`digits' < 0) {
    display as error "number of digits must be a nonnegative integer"
    exit
  }

  local struct "structural"
  local boot "bootstrap"
  local props = e(properties)
  local isstruct : list struct in props
  local isboot : list boot in props
  if (!`isstruct') {
    display as error "the fitted plssem model includes only the measurement part"
    exit
  }
  else if (`isboot') {
    display as error "estat vif is not allowed after bootstrap"
    exit
  }
  else {
    tempname strvif
    matrix `strvif' = e(struct_vif)
    local hline_path = rowsof(`strvif')
    mktable, matrix(`strvif') digits(`digits') firstcolname("Variable") ///
      title("Structural model - Multicollinearity check (VIFs)") ///
      firstcolwidth(14) colwidth(12) hlines(`hline_path') novlines
  }
  
  /* Return values */
  return matrix strvif `strvif'
end

program mediate_c, rclass
  version 15.1
  syntax , INDep(string) MED(string) DEP(string) ///
    [ BReps(numlist integer >0 min=1 max=1) Seed(numlist max=1) ZLC RIT RID ///
    BCa Level(real 0.95) DIGits(integer 3) ]
  
  /* Options:
     --------
     indep(string)              --> independent variable
     med(string)                --> mediator variable
     dep(string)                --> dependent variable
     breps(numlist integer >0 min=1 max=1)
                                --> number of bootstrap replications (default
                                    is 50)
     seed(numlist max=1)        --> bootstrap seed number
     zlc                        --> mediation procedures described by Zhao
                                    et al. (2010)
     rit                        --> ratio of the indirect effect to the total
                                    effect
     rid                        --> ratio of the indirect effect to the direct
                                    effect
     bca                        --> returns bias-corrected accelerated bootstrap
                                    confidence intervals instead of percentile
                                    confidence intervals, which is the default
     level(real 0.95)           --> confidence level (default 0.95)
     digits(integer 3)          --> number of digits to display (default 3)
   */
   
   /* Description:
      ------------
      This postestimation command provides the details for a single indirect
      effect mediated by only one LV.
   */

  tempname coef coef_el
  matrix `coef' = e(struct_b)
  
  local exo_nm : rownames `coef'
  local endo_nm : colnames `coef'
  
  if (!`: list dep in endo_nm') {
    display as error "`dep' is not an endogenous variable"
    exit
  }
  if (!`: list med in endo_nm') {
    display as error "`med' is not an endogenous variable"
    exit
  }
  if (!`: list indep in exo_nm') {
    display as error "`indep' is not an exogenous variable"
    exit
  }
  if (`digits' < 0) {
    display as error "number of digits must be a nonnegative integer"
    exit
  }

  local coef_el = `coef'[rownumb(`coef', "`indep'"), colnumb(`coef', "`dep'")]
  if (missing(`coef_el')) {
    display as error "the `indep' variable has no direct effect on `dep'"
    exit
  }
  local coef_el = `coef'[rownumb(`coef', "`med'"), colnumb(`coef', "`dep'")]
  if (missing(`coef_el')) {
    display as error "the `med' variable has no direct effect on `dep'"
    exit
  }
  local coef_el = `coef'[rownumb(`coef', "`indep'"), colnumb(`coef', "`med'")]
  if (missing(`coef_el')) {
    display as error "the `indep' variable has no direct effect on `med'"
    exit
  }

  local struct "structural"
  local props = e(properties)
  local isstruct : list struct in props
  if (!`isstruct') {
    display as error "the fitted plssem model includes only the measurement part"
    exit
  }

  tempvar __touse__
  quietly generate `__touse__' = e(sample)
  local reg3eqs = e(struct_eqs)
  if (`level' <= 0 | `level' >= 1) {
    display as error "confidence level must be in the range (0, 1)"
    exit
  }

  if ("`breps'" == "") {
    local breps = 50
  }

  tempname normperc alpha_cl
  scalar `alpha_cl' = 1 - ((1 - `level')/2)
  scalar `normperc' = invnormal(`alpha_cl')

  tempname ehold
  _estimates hold `ehold'

  display as text "Bootstrapping mediation effect..."

  capture {
    if ("`dep'" != "" & "`med'" != "" & "`indep'" != "") {
      local doi `dep':`indep'
      local moi `med':`indep'
      local dom `dep':`med'

      tempname reg3coef_b reg3var_b reg3ciP reg3ciBC reg3coef reg3var
      if ("`bca'" == "") {
        quietly bootstrap indeff=(_b[`dom']*_b[`moi']), ties reps(`breps') ///
          seed(`seed'): reg3 `reg3eqs' if `__touse__', mvreg corr(independent)
      }
      else {
        quietly bootstrap indeff=(_b[`dom']*_b[`moi']), bca ties reps(`breps') ///
          seed(`seed'): reg3 `reg3eqs' if `__touse__', mvreg corr(independent)
      }
      quietly estat bootstrap, all
      matrix `reg3coef_b' = e(b)
      matrix `reg3var_b' = e(V)
      matrix `reg3ciP' = e(ci_percentile)
      if ("`bca'" == "") {
        matrix `reg3ciBC' = e(ci_bc)
      }
      else {
        matrix `reg3ciBC' = e(ci_bca)
      }

      tempname perc_lci perc_uci bc_lci bc_uci
      scalar `perc_lci' = `reg3ciP'[1, 1]
      scalar `perc_uci' = `reg3ciP'[2, 1]
      scalar `bc_lci' = `reg3ciBC'[1, 1]
      scalar `bc_uci' = `reg3ciBC'[2, 1]

      quietly reg3 `reg3eqs' if `__touse__', mvreg corr(independent)
      matrix `reg3coef' = e(b)
      matrix `reg3var' = e(V)

      tempname moimat moicoef moimat2 moivar dommat domcoef dommat2 domvar
      matrix `moimat' = `reg3coef'[1, "`moi'"]
      scalar `moicoef' = `moimat'[1, 1]
      matrix `moimat2' = `reg3var'["`moi'", "`moi'"]
      scalar `moivar' = `moimat2'[1, 1]
      matrix `dommat' = `reg3coef'[1, "`dom'"]
      scalar `domcoef' = `dommat'[1, 1]
      matrix `dommat2' = `reg3var'["`dom'", "`dom'"]
      scalar `domvar' = `dommat2'[1, 1]
    }
    else {
      display as error "provide a single dependent, independent and mediator variable for each indirect effect"
      exit
    }
  } // end of -capture-
  local rc = _rc
  error `rc'
  
  /* Description of mediating effect */
  tempname indmeas mat_bk mat_zlc mat_rit mat_rid
  matrix `indmeas' = J(6, 3, .)
  matrix `mat_bk' = J(7, 1, .)
  matrix `mat_zlc' = J(4, 1, .)
  matrix `mat_rit' = J(3, 1, .)
  matrix `mat_rid' = J(3, 1, .)
  
  tempname coef_moi se_moi var_moi coef_dom se_dom var_dom prodterm
  scalar `coef_moi' = `moicoef'
  scalar `var_moi' = `moivar'
  scalar `se_moi' = sqrt(`var_moi')
  scalar `coef_dom' = `domcoef'
  scalar `var_dom' = `domvar'
  scalar `se_dom' = sqrt(`var_dom')
  scalar `prodterm' = `coef_moi'*`coef_dom'

  // Sobel's method
  tempname sobel_se sobel_z sobel_pv sobel_lci sobel_uci
  scalar `sobel_se' = sqrt(((`coef_dom')^2)*`var_moi' + ((`coef_moi')^2)*`var_dom')
  scalar `sobel_z' = `prodterm'/`sobel_se'
  scalar `sobel_pv' =  2*(1 - normal(abs(`sobel_z')))
  scalar `sobel_lci' = `prodterm' - `normperc'*`sobel_se'
  scalar `sobel_uci' = `prodterm' + `normperc'*`sobel_se'

  matrix `indmeas'[1, 1] = `prodterm'
  matrix `indmeas'[2, 1] = `sobel_se'
  matrix `indmeas'[3, 1] = `sobel_z'
  matrix `indmeas'[4, 1] = `sobel_pv'
  matrix `indmeas'[5, 1] = `sobel_lci'
  matrix `indmeas'[6, 1] = `sobel_uci'
  
  // Delta method
  tempname v delta_var delta_se delta_z delta_pv delta_lci delta_uci
  quietly nlcom _b[`moi']*_b[`dom']
  matrix `v' = r(V)
  scalar `delta_var' = `v'[1, 1]
  scalar `delta_se' = sqrt(`delta_var')
  scalar `delta_z' = `prodterm'/`delta_se'
  scalar `delta_pv' =  2*(1 - normal(abs(`delta_z')))
  scalar `delta_lci' = `prodterm' - `normperc'*`delta_se'
  scalar `delta_uci' = `prodterm' + `normperc'*`delta_se'

  matrix `indmeas'[1, 2] = `prodterm'
  matrix `indmeas'[2, 2] = `delta_se'
  matrix `indmeas'[3, 2] = `delta_z'
  matrix `indmeas'[4, 2] = `delta_pv'
  matrix `indmeas'[5, 2] = `delta_lci'
  matrix `indmeas'[6, 2] = `delta_uci'
  
  // Bootstrap method
  tempname boot_prod boot_se boot_z boot_pv
  scalar `boot_prod' = `reg3coef_b'[1, 1]
  scalar `boot_se' = sqrt(`reg3var_b'[1, 1])
  scalar `boot_z' = `boot_prod'/`boot_se'
  scalar `boot_pv' =  2*(1 - normal(abs(`boot_z')))             // CHECK THIS!!!
  
  matrix `indmeas'[1, 3] = `boot_prod'
  matrix `indmeas'[2, 3] = `boot_se'
  matrix `indmeas'[3, 3] = `boot_z'
  matrix `indmeas'[4, 3] = `boot_pv'
  if ("`bca'" == "") {
    matrix `indmeas'[5, 3] = `perc_lci'
    matrix `indmeas'[6, 3] = `perc_uci'
  }
  else {
    matrix `indmeas'[5, 3] = `bc_lci'
    matrix `indmeas'[6, 3] = `bc_uci'
  }
  
  matrix rownames `indmeas' = "Indirect effect" "Standard error" ///
    "Z statistic" "P-value" "Lower CI" "Upper CI"
  matrix colnames `indmeas' = "Sobel" "Delta" "Bootstrap"

  /* Clean up */
  _estimates unhold `ehold'

  /* Baron-Kenny mediation testing adjusted by Iacobucci et al. */
  tempname coef st_tab pval stderr
  matrix `coef' = e(struct_b)
  local struct_b_rn : rownames `coef'
  local struct_b_cn : colnames `coef'
  matrix `st_tab' = e(struct_table)
  local st_nrow : rowsof(`st_tab')
  local st_nrow = (`st_nrow' - 1)/2
  local st_ncol : colsof(`st_tab')
  matrix `pval' = J(`st_nrow', `st_ncol', .)
  forvalues i = 1/`st_nrow' {
    matrix `pval'[`i', 1] = `st_tab'[`i'*2, 1...]
  }
  matrix rownames `pval' = `struct_b_rn'
  matrix colnames `pval' = `struct_b_cn'
  mat `stderr' = e(struct_se)

  // X -> M
  tempname coef_moi_m coef_moi moi_pval_m moi_pval
  matrix `coef_moi_m' = `coef'["`indep'", "`med'"]
  scalar `coef_moi' = `coef_moi_m'[1, 1]
  matrix `moi_pval_m' = `pval'["`indep'", "`med'"]
  scalar `moi_pval' = `moi_pval_m'[1, 1]
  
  // M -> Y
  tempname coef_dom_m coef_dom dom_pval_m dom_pval
  matrix `coef_dom_m' = `coef'["`med'", "`dep'"]
  scalar `coef_dom' = `coef_dom_m'[1, 1]
  matrix `dom_pval_m' = `pval'["`med'", "`dep'"]
  scalar `dom_pval' = `dom_pval_m'[1, 1]
  
  // X -> Y
  tempname coef_doi_m coef_doi doi_pval_m doi_pval
  matrix `coef_doi_m' = `coef'["`indep'", "`dep'"]
  scalar `coef_doi' = `coef_doi_m'[1, 1]
  matrix `doi_pval_m' = `pval'["`indep'", "`dep'"]
  scalar `doi_pval' = `doi_pval_m'[1, 1]
  
  matrix `mat_bk'[1, 1] = `coef_doi'
  matrix `mat_bk'[2, 1] = `coef_moi'
  matrix `mat_bk'[3, 1] = `coef_dom'
  matrix `mat_bk'[4, 1] = `doi_pval'
  matrix `mat_bk'[5, 1] = `moi_pval'
  matrix `mat_bk'[6, 1] = `dom_pval'
  matrix `mat_bk'[7, 1] = `sobel_pv'

  /* Zhao et al. mediation testing */
  if ("`zlc'" != "") {
    tempname axbxc
    scalar `axbxc' = `coef_moi'*`coef_dom'*`coef_doi'
    
    matrix `mat_zlc'[1, 1] = `boot_pv'
    matrix `mat_zlc'[2, 1] = `doi_pval'
    matrix `mat_zlc'[3, 1] = `coef_doi'
    matrix `mat_zlc'[4, 1] = `axbxc'
  }
  
  if ("`rit'" != "") {
    tempname prod toteff
    scalar `prod' = abs(`prodterm')
    scalar `toteff' = abs(`prodterm' + `coef_doi')
    
    matrix `mat_rit'[1, 1] = `prod'
    matrix `mat_rit'[2, 1] = `toteff'
    matrix `mat_rit'[3, 1] = `prod'/`toteff'
  }
  
  if ("`rid'" != "") {
    tempname prod absprod
    scalar `prod' = abs(`prodterm')
    scalar `absprod' = abs(`coef_doi')
    
    matrix `mat_rid'[1, 1] = `prod'
    matrix `mat_rid'[2, 1] = `absprod'
    matrix `mat_rid'[3, 1] = `prod'/`absprod'
  }
  
  /* Display results */
  local med_line = 5
  mktable_mediate, matrix(`indmeas') matrix_bk(`mat_bk') ///
    matrix_zlc(`mat_zlc') matrix_rit(`mat_rit') matrix_rid(`mat_rid') ///
    depv("`dep'") medv("`med'") indepv("`indep'") ///
    title("Significance testing of (standardized) indirect effect") ///
    firstcolname("Statistics") firstcolwidth(24) colwidth(20) ///
    hlines(`med_line') reps(`breps') digits(`digits') level(`level') novlines
  
  /* Return values */
  return matrix mediate `indmeas'
end

program htmt_c, rclass
  version 15.1
  syntax , [ All CUToff(real 0) DIGits(integer 3) ]
  
  /* Options:
     --------
     all                        --> include all latent variables
     cutoff(real 0)             --> do not show correlation smaller than cutoff
     digits(integer 3)          --> number of digits to display (default 3)
  */
   
  /* Description:
     ------------
     This postestimation command provides the heterotrait-monotrait ratio of
     correlations for assessing discriminant validity.
  */

  if (`digits' < 0) {
    display as error "number of digits must be a nonnegative integer"
    exit
  }
  
  local is_all = 0
  if ("`all'" != "") {
    local is_all = 1
  }

  /* Set temporary variables */
  tempvar __touse__
  quietly generate `__touse__' = e(sample)
  local allindicators = e(mvs)
  local alllatents = e(lvs)
  foreach var in `allindicators' {
    local allstdindicators "`allstdindicators' std`var'"
  }
  local allstdindicators : list clean allstdindicators
  local tempnamelist

  /* Compute the heterotrait-monotrait ratios (HTMT) */
  tempname res_htmt mata_htmt res_htmt2 mata_htmt2
  local tempnamelist "`tempnamelist' `mata_htmt' `mata_htmt2'"
  capture noisily {
    mata: `mata_htmt' = ///
      plssem_htmt( ///
        st_matrix("e(ind_vcv)"), ///
        st_matrix("e(adj_meas)"), ///
        "`allindicators'", ///
        "`allstdindicators'", ///
        "`alllatents'", ///
        "`e(binarylvs)'", ///
        "`e(reflective)'", ///
        `is_all')
    mata: `mata_htmt2' = ///
      plssem_htmt2( ///
        st_matrix("e(ind_vcv)"), ///
        st_matrix("e(adj_meas)"), ///
        "`allindicators'", ///
        "`allstdindicators'", ///
        "`alllatents'", ///
        "`e(binarylvs)'", ///
        "`e(reflective)'", ///
        `is_all')
  }

  /* Display results */
  mata: st_matrix("`res_htmt'", `mata_htmt')
  if (`is_all') {
    matrix rownames `res_htmt' = `alllatents'
    matrix colnames `res_htmt' = `alllatents'
  }
  else {
    matrix rownames `res_htmt' = `e(reflective)'
    matrix colnames `res_htmt' = `e(reflective)'
  }
  mata: st_matrix("`res_htmt2'", `mata_htmt2')
  if (`is_all') {
    matrix rownames `res_htmt2' = `alllatents'
    matrix colnames `res_htmt2' = `alllatents'
  }
  else {
    matrix rownames `res_htmt2' = `e(reflective)'
    matrix colnames `res_htmt2' = `e(reflective)'
  }

  tempname res_htmt2_null
  mata: st_numscalar("`res_htmt2_null'", allof(`mata_htmt2', .))

  mktable_corr, matrix(`res_htmt') ///
    title("Discriminant validity - Heterotrait-monotrait ratio of correlations (HTMT)") ///
    cutoff(`cutoff') digits(`digits')
  if (!`res_htmt2_null') {
    mktable_corr, matrix(`res_htmt2') ///
      title("Discriminant validity - Advanced heterotrait-monotrait ratio of correlations (HTMT2)") ///
      cutoff(`cutoff') digits(`digits')
  }
  else {
    display as text ""
    display as text "Note: advanced heterotrait-monotrait ratios (HTMT2) are not available"
  }
  
  /* Return values */
  return matrix htmt `res_htmt'
  return matrix htmt2 `res_htmt2'

  /* Clean up */
  capture mata: cleanup(st_local("tempnamelist"))
end

program ci_c, rclass
  version 15.1
  syntax , Type(string) [ Level(real 0.95) DIGits(integer 3) ]
  
  /* Options:
     --------
     type(string)               --> type of confidence interval; one of
                                    - "standard_z" or "standard_t" ==> standard CI with
                                      bootstrap SE's or OLS SE's; critical quantiles can
                                      be based on either the z or or the t distribution; the
                                      latter case may perform better in small samples but
                                      there is no general consensus on the degrees of freedom
                                      to use
                                    - "percentile" ==> the confidence interval's bounds are
                                      estimated as the alpha/2 and 1-alpha/2 quantiles of the
                                      distribution of the resample estimates; works only when
                                      bootstrap has been used
                                    - "bc" ==> bias corrected (Bc) confidence interval
     level(real 0.95)           --> confidence level (default 0.95)
     digits(integer 3)          --> number of digits to display (default 3)
  */
   
  /* Description:
     ------------
     This postestimation command provides the confidence intervals for all
     model parameters.
  */
  
  /* Set temporary variables */
  tempvar __touse__
  quietly generate `__touse__' = e(sample)
  local tempnamelist

  local isboot = (e(vce) == "bootstrap")
  local N = e(N)
  
  if (`level' <= 0 | `level' >= 1) {
    display as error "confidence level must be in the range (0, 1)"
    exit
  }
  if (`digits' < 0) {
    display as error "number of digits must be a nonnegative integer"
    exit
  }
  /* End of setting temporary variables */

  /* Computing confidence intervals */
  tempname alpha_cl mata_b mata_se perc mata_load mata_path ///
    mata_load_se mata_path_se mata_load_lower mata_load_upper ///
    ci_load mata_path_lower mata_path_upper ci_path
  local tempnamelist "`alpha_cl' `tempnamelist' `mata_b' `mata_se'" ///
    "`perc' `mata_load' `mata_path' `mata_load_se'" ///
    "`mata_path_se' `mata_load_lower' `mata_load_upper' `ci_load'" ///
    "`mata_path_lower' `mata_path_upper' `ci_path'"

  mata: `alpha_cl' = 1 - ((1 - `level')/2)
  mata: st_local("num_indep", strofreal(sum(st_matrix("e(adj_meas)"))))
  mata: `mata_b' = st_matrix("e(b)")'
  mata: `mata_se' = sqrt(diagonal(st_matrix("e(V)")))
  mata: `mata_load' = `mata_b'[1::strtoreal(st_local("num_indep"))]
  mata: `mata_load_se' = `mata_se'[1::strtoreal(st_local("num_indep"))]
  mata: `mata_path' = `mata_b'[(strtoreal(st_local("num_indep")) + 1)::length(`mata_b')]
  mata: `mata_path_se' = `mata_se'[(strtoreal(st_local("num_indep")) + 1)::length(`mata_b')]
  local b_eqnm : coleq e(b)
  local b_nm : colnames e(b)
  mata: st_local("load_eqnm", ///
    invtokens(tokens(st_local("b_eqnm"))[1::strtoreal(st_local("num_indep"))]))
  mata: st_local("load_nm", ///
    invtokens(tokens(st_local("b_nm"))[1::strtoreal(st_local("num_indep"))]))
  mata: st_local("path_eqnm", ///
    invtokens(tokens(st_local("b_eqnm"))[(strtoreal(st_local("num_indep")) + 1)::length(`mata_b')]))
  mata: st_local("path_nm", ///
    invtokens(tokens(st_local("b_nm"))[(strtoreal(st_local("num_indep")) + 1)::length(`mata_b')]))

  local i = 1
  local allreflective = e(reflective)
  local allnames_load
  foreach nm_eq in `load_eqnm' {
    local nm_tmp : word `i' of `load_nm'
    if (`: list nm_eq in allreflective') {
      local allnames_load "`allnames_load' `nm_eq'->`nm_tmp'"
    }
    else {
      local allnames_load "`allnames_load' `nm_eq'<-`nm_tmp'"
    }
    local ++i
  }
  local i = 1
  local allnames_path
  foreach nm_eq in `path_eqnm' {
    local nm_tmp : word `i' of `path_nm'
    local allnames_path "`allnames_path' `nm_tmp'->`nm_eq'"
    local ++i
  }

  if (`isboot') {
    if ("`type'" == "standard_z" | "`type'" == "standard_t") {
      if ("`type'" == "standard_z") {
        mata: `perc' = invnormal(`alpha_cl')
      }
      else if ("`type'" == "standard_t") {
        mata: `perc' = invt(strtoreal(st_local("N")) - 1, `alpha_cl')
      }
      mata: `mata_load_lower' = `mata_load' - `perc':*`mata_load_se'
      mata: `mata_load_upper' = `mata_load' + `perc':*`mata_load_se'
      mata: st_matrix("`ci_load'", ///
        (`mata_load', `mata_load_se', `mata_load_lower', `mata_load_upper'))
      matrix rownames `ci_load' = `allnames_load'
      matrix colnames `ci_load' = "Coefficient" "Std._err." "Lower_CI" "Upper_CI"
      mata: `mata_path_lower' = `mata_path' - `perc':*`mata_path_se'
      mata: `mata_path_upper' = `mata_path' + `perc':*`mata_path_se'
      mata: st_matrix("`ci_path'", ///
        (`mata_path', `mata_path_se', `mata_path_lower', `mata_path_upper'))
      matrix rownames `ci_path' = `allnames_path'
      matrix colnames `ci_path' = "Coefficient" "Std._err." "Lower_CI" "Upper_CI"
    }
    else if ("`type'" == "percentile") {
      tempname load_boot path_boot load_q path_q
      local tempnamelist "`tempnamelist' `load_boot' `path_boot'" ///
        "`load_q' `path_q'"
      mata: `load_boot' = st_matrix("e(loadings_breps)")
      mata: `path_boot' = st_matrix("e(pathcoef_breps)")
      mata: `path_boot' = `path_boot'[., selectindex(colnonmissing(`path_boot'))]
      mata: `load_q' = plssem_quantile(`load_boot', (1 - `alpha_cl' \ `alpha_cl'), 1)
      mata: `path_q' = plssem_quantile(`path_boot', (1 - `alpha_cl' \ `alpha_cl'), 1)
      mata: `mata_load_lower' = `load_q'[1, .]'
      mata: `mata_load_upper' = `load_q'[2, .]'
      mata: st_matrix("`ci_load'", ///
        (`mata_load', `mata_load_se', `mata_load_lower', `mata_load_upper'))
      matrix rownames `ci_load' = `allnames_load'
      matrix colnames `ci_load' = "Coefficient" "Std._err." "Lower_CI" "Upper_CI"
      mata: `mata_path_lower' = `path_q'[1, .]'
      mata: `mata_path_upper' = `path_q'[2, .]'
      mata: st_matrix("`ci_path'", ///
        (`mata_path', `mata_path_se', `mata_path_lower', `mata_path_upper'))
      matrix rownames `ci_path' = `allnames_path'
      matrix colnames `ci_path' = "Coefficient" "Std._err." "Lower_CI" "Upper_CI"
    }
    else if ("`type'" == "bc") {
      tempname load_boot path_boot load_q path_q ci_bc
      local tempnamelist "`tempnamelist' `load_boot' `path_boot' `ci_bc'"
      mata: `load_boot' = st_matrix("e(loadings_breps)")
      mata: `path_boot' = st_matrix("e(pathcoef_breps)")
      mata: `path_boot' = `path_boot'[., selectindex(colnonmissing(`path_boot'))]

      mata: `ci_bc' = plssem_bcCIresample(`load_boot', `mata_load', (1 - `alpha_cl' \ `alpha_cl'))
      mata: `mata_load_lower' = `ci_bc'[1, .]'
      mata: `mata_load_upper' = `ci_bc'[2, .]'
      mata: st_matrix("`ci_load'", ///
        (`mata_load', `mata_load_se', `mata_load_lower', `mata_load_upper'))
      matrix rownames `ci_load' = `allnames_load'
      matrix colnames `ci_load' = "Coefficient" "Std._err." "Lower_CI" "Upper_CI"
      mata: `ci_bc' = plssem_bcCIresample(`path_boot', `mata_path', (1 - `alpha_cl' \ `alpha_cl'))
      mata: `mata_path_lower' = `ci_bc'[1, .]'
      mata: `mata_path_upper' = `ci_bc'[2, .]'
      mata: st_matrix("`ci_path'", ///
        (`mata_path', `mata_path_se', `mata_path_lower', `mata_path_upper'))
      matrix rownames `ci_path' = `allnames_path'
      matrix colnames `ci_path' = "Coefficient" "Std._err." "Lower_CI" "Upper_CI"
    }
    else {
      display as error "the specified confidence interval type is not available"
      exit
    }
  }
  else {
    display as error "the calculation of confidence intervals using PLSc is allowed " _continue
    display as error "only with bootstrap"
    exit
  }
  /* End of computing confidence intervals */
  
  /* Display results */
  display
  display as text "Consistent partial least squares (PLSc) SEM"
  display
  local subtitle "Confidence intervals for model's coefficients"
  display as text "`subtitle'"
  display as text "Standard error estimation method: " _continue
  display as result "`e(vce)'"
  display as text "Confidence interval method: " _continue
  display as result "`type'"
  display as text "Number of obs: " _continue
  display as result e(N)

  local skip1 = 1
  local clevel = `level'*100
  local title "Measurement model - Standardized loadings"
  mata: st_local("firstcolwidth", ///
    strofreal(max(strlen(tokens(st_local("allnames_load")))) + 4))
  local neff = rowsof(`ci_load')
  mktable_ci, matrix(`ci_load') digits(`digits') ///
    firstcolname("Effect") title(`title') ///
    firstcolwidth(`firstcolwidth') colwidth(13 11 11 11) ///
    hlines(`neff') novlines
  display as text _skip(`skip1') "confidence level: " _continue
  display as result "`clevel'" _continue
  display as text "%"

  local title "Structural model - Standardized path coefficients"
  mata: st_local("firstcolwidth", ///
    strofreal(max(strlen(tokens(st_local("allnames_path")))) + 4))
  local neff = rowsof(`ci_path')
  mktable_ci, matrix(`ci_path') digits(`digits') ///
    firstcolname("Effect") title(`title') ///
    firstcolwidth(`firstcolwidth') colwidth(13 11 11 11) ///
    hlines(`neff') novlines
  display as text _skip(`skip1') "confidence level: " _continue
  display as result "`clevel'" _continue
  display as text "%"
  /* End of displaying results */
  
  /* Return values */
  return local ci_type "`type'"

  return scalar level = `level'

  return matrix load_ci `ci_load'
  return matrix path_ci `ci_path'

  /* Clean up */
  capture mata: cleanup(st_local("tempnamelist"))
end

program f2_c, rclass
  version 15.1
  syntax , [ DIGits(integer 3) ]
  
  /* Options:
     --------
     digits(integer 3)          --> number of digits to display (default 3)
  */
   
  /* Description:
     ------------
     This postestimation command provides the Cohen's f2  effect sizes.
  */
  
  local struct "structural"
  local props = e(properties)
  local isstruct : list struct in props
  if (!`isstruct') {
    display as error "the fitted plssem model includes only the measurement part"
    exit
  }

  if (`digits' < 0) {
    display as error "number of digits must be a nonnegative integer"
    exit
  }

  /* Set temporary variables */
  tempvar __touse__
  quietly generate `__touse__' = e(sample)
  local allindicators = e(mvs)
  local alllatents = e(lvs)
  local robust = e(robust)
  local ordinal = e(ordinal)
  local tempnamelist

  tempname modes
  local num_lv: word count `alllatents'
  local modeA = e(reflective)
  matrix `modes' = J(`num_lv', 1, 1)
  local i = 1
  foreach var in `alllatents' {
    if (`: list var in modeA') {
      matrix `modes'[`i', 1] = 0
    }
    local ++i
  }
  
  /* Compute the Cohen's f^2 effect sizes */
  tempname res_f2 mata_f2
  local tempnamelist "`tempnamelist' `mata_f2'"
  capture noisily {
    mata: `mata_f2' = ///
      plssem_f2( ///
        st_data(., "`allindicators'"), ///       note: `__touse__' not used here
        st_data(., "`alllatents'"), ///
        st_matrix("e(adj_meas)"), ///
        st_matrix("e(adj_struct)"), ///
        st_matrix("e(outerweights)"), ///
        st_matrix("`modes'"), ///
        st_matrix("e(rsquared)"), ///
        "`alllatents'", ///
        "`e(binarylvs)'", ///
        "`__touse__'", ///
        1, ///
        "`robust'", ///
        "`allindicators'", ///
        "`ordinal'")
  }

  /* Display results */
  mata: st_matrix("`res_f2'", `mata_f2')
  matrix rownames `res_f2' = `alllatents'
  mata: st_local("endo_nm", ///
    invtokens(tokens("`e(lvs)'")[selectindex(colsum(st_matrix("e(adj_struct)")))]))
  matrix colnames `res_f2' = `endo_nm'

  mktable, matrix(`res_f2') digits(`digits') ///
    title("Effect sizes  - Cohen's f^2") ///
    firstcolwidth(13) colwidth(12) novlines hlines(`: word count `alllatents'')
  
  /* Return values */
  return matrix f2 `res_f2'

  /* Clean up */
  capture mata: cleanup(st_local("tempnamelist"))
end

program dist_c, rclass
  version 15.1
  syntax , [ DIGits(integer 3) ]
  
  /* Options:
     --------
     digits(integer 3)          --> number of digits to display (default 3)
  */
   
  /* Description:
     ------------
     This postestimation command provides some model distance measures.
  */
  
  if (`digits' < 0) {
    display as error "number of digits must be a nonnegative integer"
    exit
  }

  /* Set temporary variables */
  tempvar __touse__
  quietly generate `__touse__' = e(sample)
  local allindicators = e(mvs)
  local alllatents = e(lvs)
  local tempnamelist

  tempname modes
  local num_lv: word count `alllatents'
  local modeA = e(reflective)
  matrix `modes' = J(`num_lv', 1, 1)
  local i = 1
  foreach var in `alllatents' {
    if (`: list var in modeA') {
      matrix `modes'[`i', 1] = 0
    }
    local ++i
  }
  
 /* Compute the geodesic distance */
 tempname mata_DG
  local tempnamelist "`tempnamelist' `mata_DG'"
  capture noisily {
    mata: `mata_DG' = ///
      plssem_DG( ///
        st_matrix("e(ind_vcv)"), ///
        st_matrix("e(adj_meas)"), ///
        st_matrix("e(adj_struct)"), ///
        st_matrix("e(outerweights)"), ///
        st_matrix("`modes'"), ///
        editmissing(st_matrix("e(loadings)"), 0), ///
        st_matrix("e(pathcoef)"), ///
        st_matrix("e(construct_vcv)"), ///
        st_matrix("e(proxy_vcv)"))
  }

 /* Compute the squared Euclidean distance */
 tempname mata_DL
  local tempnamelist "`tempnamelist' `mata_DL'"
  capture noisily {
    mata: `mata_DL' = ///
      plssem_DL( ///
        st_matrix("e(ind_vcv)"), ///
        st_matrix("e(adj_meas)"), ///
        st_matrix("e(adj_struct)"), ///
        st_matrix("e(outerweights)"), ///
        st_matrix("`modes'"), ///
        editmissing(st_matrix("e(loadings)"), 0), ///
        st_matrix("e(pathcoef)"), ///
        st_matrix("e(construct_vcv)"), ///
        st_matrix("e(proxy_vcv)"))
  }

 /* Compute the ML distance */
 tempname mata_DML
  local tempnamelist "`tempnamelist' `mata_DML'"
  capture noisily {
    mata: `mata_DML' = ///
      plssem_DML( ///
        st_matrix("e(ind_vcv)"), ///
        st_matrix("e(adj_meas)"), ///
        st_matrix("e(adj_struct)"), ///
        st_matrix("e(outerweights)"), ///
        st_matrix("`modes'"), ///
        editmissing(st_matrix("e(loadings)"), 0), ///
        st_matrix("e(pathcoef)"), ///
        st_matrix("e(construct_vcv)"), ///
        st_matrix("e(proxy_vcv)"))
  }

  /* Display results */
  tempname res_dist
  mata: st_matrix("`res_dist'", (`mata_DG' \ `mata_DL' \ `mata_DML'))
  matrix rownames `res_dist' = Geodesic Squared_Euclidean ML
  matrix colnames `res_dist' = "Value"

  mktable, matrix(`res_dist') digits(`digits') ///
    firstcolname("Measure") title("Distance measures") ///
    firstcolwidth(20) colwidth(14) novlines hlines(3)
  
  /* Return values */
  return matrix dist `res_dist'

  /* Clean up */
  capture mata: cleanup(st_local("tempnamelist"))
end

program ic_c, rclass
  version 15.1
  syntax , [ DIGits(integer 3) ]
  
  /* Options:
     --------
     digits(integer 3)          --> number of digits to display (default 3)
  */
   
  /* Description:
     ------------
     This postestimation command provides some model selection criteria.
  */
  
  local struct "structural"
  local props = e(properties)
  local isstruct : list struct in props
  if (!`isstruct') {
    display as error "the fitted plssem model includes only the measurement part"
    exit
  }

  if (`digits' < 0) {
    display as error "number of digits must be a nonnegative integer"
    exit
  }

  /* Set temporary variables */
  local allendogenous "`: colnames e(struct_b)'"
  local tempnamelist

  /* Compute the model selection criteria */
  tempname res_crit mata_crit
  local tempnamelist "`tempnamelist' `mata_crit'"
  capture noisily {
    mata: `mata_crit' = ///
      plssem_selcrit( ///
        st_numscalar("e(N)"), ///
        st_matrix("e(adj_struct)"), ///
        st_matrix("e(rsquared)"))
  }

  /* Display results */
  mata: st_matrix("`res_crit'", `mata_crit')
  matrix rownames `res_crit' = `allendogenous'
  matrix colnames `res_crit' = "AIC AICc AICu BIC FPE HQ HQc"

  mktable, matrix(`res_crit') digits(`digits') ///
    title("Model information and selection criteria") ///
    firstcolwidth(13) colwidth(15) novlines hlines(`: word count `allendogenous'')
 
  /* Return values */
  return matrix ic `res_crit'

  /* Clean up */
  capture mata: cleanup(st_local("tempnamelist"))
end

program blindf_c, rclass
  version 15.1
  syntax , DIStance(numlist integer >1 min=1 max=1) [ DIGits(integer 3) ]
  
  /* Options:
     --------
     distance(numlist integer >1 min=1 max=1)
                                --> omission distance
     digits(integer 3)          --> number of digits to display (default 3)
  */
   
  /* Description:
     ------------
     This postestimation command applies the blindfolding approach to a fitted model.
  */
  
  tempvar __touse__
  quietly generate `__touse__' = e(sample)

  local struct "structural"
  local props = e(properties)
  local isstruct : list struct in props
  if (!`isstruct') {
    display as error "the fitted plssem model includes only the measurement part"
    exit
  }

  mata: st_local("endo_nm", ///
    invtokens(tokens("`e(lvs)'")[selectindex(colsum(st_matrix("e(adj_struct)")))]))
  local refl_nm `e(reflective)'
  local endo_refl : list endo_nm & refl_nm
  if ("`endo_refl'" == "") {
    display as error "the blindfolding procedure requires at least one " _continue
    display as error "endogenous reflective latent variable"
    exit
  }

  if (`digits' < 0) {
    display as error "number of digits must be a nonnegative integer"
    exit
  }

  local n e(N)
  if (mod(`n', `distance') == 0) {
    display as error "omission distance must not be a divisor of sample size"
    exit
  }

  /* Compute the q2 effect sizes through blindfolding */
  tempname res_q2 mata_q2 X isbinary i nvals modes Xstd Yinit res_blind
  tempname mata_blind
  local tempnamelist "`tempnamelist' `mata_q2' `X' `isbinary' `i' `nvals'"
  local tempnamelist "`tempnamelist' `Xstd' `Yinit' `mata_blind'"

  mata: `X' = st_data(., "`e(mvs)'", "`__touse__'")   // note: ___touse___
                                                      // already applied
  local mean "mean"
  local ismean : list mean in props
  if (`ismean') {
    mata: `isbinary' = plssem_isbinary(`X')
    mata: `X' = meanimp_mat(`X', `isbinary')
  }
  local alllatents = e(lvs)
  local allindicators = e(mvs)
  local robust = e(robust)
  local ordinal = e(ordinal)
  foreach var in `allindicators' {
    local allstdindicators "`allstdindicators' std`var'"
  }
  local allstdindicators : list clean allstdindicators
  local initmeth "eigen indsum"
  local init: list props & initmeth
  local schemeval "centroid factorial path"
  local scheme: list props & schemeval
  local relative "relative square"
  local convcrit: list props & relative
  local rawsum "rawsum"
  local rawsum_sc : list rawsum in props
  local num_lv: word count `alllatents'
  local modeA = e(reflective)
  matrix `modes' = J(`num_lv', 1, 1)
  local i = 1
  foreach var in `alllatents' {
    if (`: list var in modeA') {
      matrix `modes'[`i', 1] = 0
    }
    local ++i
  }
  
  capture noisily {
    mata: `Xstd' = ///
      plssem_scale_mat( ///
        `X', ///
        ., ///
        "")

    mata: `Yinit' = ///
      plssem_init_mat( ///
        `Xstd', ///
        st_matrix("e(adj_meas)"), ///
        "`allindicators'", ///
        "`allstdindicators'", ///
        "`alllatents'", ///
        "`__touse__'", ///
        strtoreal("`rawsum_sc'"), ///
        "`init'")

    mata: `mata_blind' = ///
      plssem_blindfolding( ///
        `X', ///
        `Yinit', ///
        st_matrix("e(adj_meas)"), ///
        st_matrix("e(adj_struct)"), ///
        st_matrix("`modes'"), ///
        "`alllatents'", ///
        "`e(binary)'", ///
        st_numscalar("e(tolerance)"), ///
        st_numscalar("e(maxiter)"), ///
        "`__touse__'", ///
        "`scheme'", ///
        "`convcrit'", ///
        strtoreal("`isstruct'"), ///
        strtoreal("`rawsum_sc'"), ///
        1, ///
        `distance', ///
        1, ///
        0, ///
        "`robust'", ///
        "`allindicators'", ///
        "`ordinal'")

    mata: `mata_q2' = ///
      plssem_q2( ///
        `X', ///
        `Yinit', ///
        st_matrix("e(adj_meas)"), ///
        st_matrix("e(adj_struct)"), ///
        st_matrix("`modes'"), ///
        "`alllatents'", ///
        "`e(binary)'", ///
        st_numscalar("e(tolerance)"), ///
        st_numscalar("e(maxiter)"), ///
        "`__touse__'", ///
        "`scheme'", ///
        "`convcrit'", ///
        strtoreal("`isstruct'"), ///
        strtoreal("`rawsum_sc'"), ///
        1, ///
        `distance', ///
        0, ///
        "`robust'", ///
        "`allindicators'", ///
        "`ordinal'")
  }

  /* Display results */
  mata: st_matrix("`res_blind'", `mata_blind')
  matrix rownames `res_blind' = `endo_refl'
  matrix colnames `res_blind' = "SSO SSE Q2"

  mktable, matrix(`res_blind') digits(`digits') firstcolname("Variable") ///
    title("Blindfolding - Construct cross-validated redundancy - Q2 measures") ///
    firstcolwidth(14) colwidth(14) novlines hlines(`: word count `endo_refl'')

  mata: st_local("q2isempty", ///
    strofreal(`mata_q2' == J(rows(`mata_q2'), cols(`mata_q2'), .)))
  if (!`q2isempty') {
    mata: st_local("q2_nonmiss", invtokens( ///
      tokens("`alllatents'")[selectindex(colnonmissing(`mata_q2') :!= 0)]))
    local q2_nm : list q2_nonmiss & endo_refl
    mata: st_matrix("`res_q2'", plssem_remove_miss_cols(`mata_q2'))
    matrix rownames `res_q2' = `alllatents'
    matrix colnames `res_q2' = `q2_nm'

    mktable, matrix(`res_q2') digits(`digits') ///
      title("Effect sizes  - q2 measures") ///
      firstcolwidth(14) colwidth(14) novlines hlines(`: word count `alllatents'')
  }
  else {
    display as text "note: the q2 measures are not available for this model"
  }
 
  /* Return values */
  return matrix Q2 `res_blind'
  if (!`q2isempty') {
    return matrix q2 `res_q2'
  }

  /* Clean up */
  capture mata: cleanup(st_local("tempnamelist"))
end
