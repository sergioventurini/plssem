{smcl}
{* *! version 0.6.3  10Oct2024}{...}
{vieweralsosee "plssemc" "help plssemc"}{...}
{vieweralsosee "plssemplot" "help plssemplot"}{...}
{viewerjumpto "Postestimation commands" "plssemc postestimation##description"}{...}
{viewerjumpto "estat" "plssemc postestimation##syntax_estat"}{...}
{viewerjumpto "estat options" "plssemc postestimation##options_estat"}{...}
{viewerjumpto "estat indirect stored results" "plssemc postestimation##results_indirect"}{...}
{viewerjumpto "estat f2 stored results" "plssemc postestimation##results_f2"}{...}
{viewerjumpto "estat ic stored results" "plssemc postestimation##results_ic"}{...}
{viewerjumpto "estat fit stored results" "plssemc postestimation##results_fit"}{...}
{viewerjumpto "estat blindfolding stored results" "plssemc postestimation##results_blind"}{...}
{viewerjumpto "predict" "plssemc postestimation##syntax_predict"}{...}
{viewerjumpto "predict options" "plssemc postestimation##options_predict"}{...}
{viewerjumpto "predict stored results" "plssemc postestimation##results_predict"}{...}
{viewerjumpto "Examples" "plssemc postestimation##examples"}{...}
{viewerjumpto "Authors" "plssemc postestimation##authors"}{...}
{viewerjumpto "References" "plssemc postestimation##references"}{...}
{title:Title}

{p 4 18 2}
{hi:plssemc postestimation} {hline 2} Postestimation tools for {helpb plssemc}


{marker description}{...}
{title:Postestimation commands}

{pstd}
The following postestimation commands are of special interest after
{cmd:plssemc}:

{synoptset 22 tabbed}{...}
{p2coldent:Command}Description{p_end}
{synoptline}
{synopt:{helpb plssemc postestimation##indirect:estat indirect}}estimation and
	inference for indirect effects{p_end}
{synopt:{helpb plssemc postestimation##total:estat total}}decomposition of
	total effects{p_end}
{p2coldent:* {helpb plssemc postestimation##vif:estat vif}}variance inflation
	factors for the structural model equations sample{p_end}
{synopt:{helpb plssemc postestimation##htmt:estat htmt}}heterotrait-monotrait
  ratio of correlations for assessing discriminant validity{p_end}
{synopt:{helpb plssemc postestimation##ci:estat ci}}confidence intervals
  for all model's coefficients{p_end}
{synopt:{helpb plssemc postestimation##f2:estat f2}}Cohen's f^2 effect
  sizes{p_end}
{synopt:{helpb plssemc postestimation##ic:estat ic}}Model's information and
  selection criteria{p_end}
{synopt:{helpb plssemc postestimation##fit:estat fit}}model's distance
  and fit measures{p_end}
{synopt:{helpb plssemc postestimation##blind:estat blindfolding}}blindfolding
  procedure{p_end}
{synoptline}
{p2colreset}{...}
{p 4 6 2}
* {cmd:estat vif} is not available for models fitted using bootstrap.
{p_end}

{pstd}
The following standard postestimation commands are also available:

{synoptset 20 tabbed}{...}
{p2coldent :Command}Description{p_end}
{synoptline}
{synopt :{helpb plssemc postestimation##predict:predict}}fitted values and residuals{p_end}
{synoptline}
{p2colreset}{...}


{marker syntax_estat}{...}
{title:Syntax for estat}

{marker indirect}{...}
{pstd}
Display the estimation results for up to 5 indirect effects

{p 8 14 2}
{cmd:estat} {cmdab:in:direct}{cmd:,}
{cmdab:e:ffects(}{it:efflist}{cmd:)} [{opt b:oot(#)} {opt s:eed(#)} {opt l:evel(#)} {opt dig:its(#)}]


{marker total}{...}
{pstd}
Display the decomposition of the total effects in the corresponding direct and indirect effects

{p 8 14 2}
{cmd:estat} {cmdab:to:tal}
[{cmd:,} {opt dig:its(#)} {opt p:lot}]


{marker vif}{...}
{pstd}
Display the variance inflation factors for the structural model equations

{p 8 14 2}
{cmd:estat} {cmdab:vi:f}
[{cmd:,} {opt dig:its(#)}]


{marker htmt}{...}
{pstd}
Display the assessment of discriminant validity using heterotrait-monotrait 
ratios of correlations

{p 8 14 2}
{cmd:estat} {cmdab:ht:mt},
[{opt cut:off(#)} {opt dig:its(#)}]


{marker ci}{...}
{pstd}
Display the confidence intervals for all model's coefficients

{p 8 14 2}
{cmd:estat} {cmdab:ci},
{cmdab:t:ype(}{it:citype}{cmd:)} [{opt l:evel(#)} {opt dig:its(#)}]


{marker f2}{...}
{pstd}
Display the Cohen's f^2 effect sizes

{p 8 14 2}
{cmd:estat} {cmdab:f2},
[{opt dig:its(#)}]


{marker ic}{...}
{pstd}
Display the model's information and selection criteria

{p 8 14 2}
{cmd:estat} {cmdab:ic},
[{opt dig:its(#)}]


{marker fit}{...}
{pstd}
Display the model's distance and fit measures

{p 8 14 2}
{cmd:estat} {cmdab:fi:t},
[{opt dig:its(#)}]


{marker blind}{...}
{pstd}
Blindfolding procedure

{p 8 14 2}
{cmd:estat} {cmdab:bl:indfolding},
{opt dis:tance(#)} [{opt dig:its(#)}]


{marker desc_estat}{...}
{title:Description for estat}

{pstd}
{cmd:estat indirect}
estimates the (standardized) indirect effects and the corresponding tests of significance
using the Sobel's {it:z} statistic (default) as well as the bootstrap
approach ({help plssemc_postestimation##Sobel1982:Sobel 1982},
{help plssemc_postestimation##BaronKenny1986:Baron and Kenny 1986}, 
{help plssemc_postestimation##VanderWeele2015:VanderWeele 2015}). The command can
estimate up to five different indirect effects at a time. Each of these should
specified by sequentially typing the dependent, mediator and independent variable
from any PLS-SEM model. By adding the sub-option {cmd:boot(#)}, you can obtain
the results based on the bootstrap approach. To facilitate the reproducibility of
bootstrap results, the sub-option {cmd:seed(#)} can further be added. Confidence intervals
({cmd:0.95} is the default) for the estimated indirect effects are also provided. To change
the level of confidence interval, the sub-option {cmd:level(#)} can be added. To change the
number of decimals used to display the model estimates, you can change the default ({cmd:3})
to any other value by adding the sub-option {cmd:digits(#)}.
 
{pstd}
{cmd:estat total}
produces the decomposition of the total effects into standardized direct and indirect effects.
Adding the sub-option {cmd:plot} generates a bar plot of the effects. You can change
the number decimals digits reported by setting the sub-option {cmd:digits(#)}.

{pstd} {cmd:estat vif}
computes the variance inflation factors (VIFs) for the independent variables
of the equations in the structural part of a PLS-SEM model. With the
{cmd:digit(#)} sub-option you change the number decimals digits displayed.

{pstd} {cmd:estat htmt}
assesses discriminant validity using heterotrait-monotrait ratios of correlations; both
the HTMT (arithmetic average of correlations) and HTMT2 (geometric average of correlations)
approaches are available.

{pstd} {cmd:estat ci}
computes the confidence intervals for both the measurement and structural models'
coefficients; different types of confidence intervals are available.

{pstd} {cmd:estat f2}
computes the Cohen's f^2 effect sizes ({help plssemc_postestimation##Cohen1988:Cohen 1988}).

{pstd} {cmd:estat ic}
computes some information and selection criteria for the fitted model for each
structural equation of the model separately, as suggested by
{help plssemc_postestimation##Sharma2019:Sharma et al. 2019}. More specifically,
the command provides the following indexes: Akaike information criterion (AIC),
Corrected AIC (AICc), Unbiased AIC (AICu), Bayesian information criterion (BIC),
final prediction error (FPE), Hannan-Quinn criterion (HQ) and Corrected HQ criterion
(HQc). The definitions used by the commend for these indexes can be found in table B1
of {help plssemc_postestimation##Sharma2019:Sharma et al. 2019}.

{pstd} {cmd:estat fit}
computes some measures of the distance between the empirical and the model-implied
indicator correlation matrices. Currently, the geodesic distance, the squared Euclidean
distance and the the maximum likelihood-based distance function are implemented.
In addition, the command reports many model's fit measures like the Chi square,
Chi square degrees freedom, CFI, CN, GFI, IFI, NFI, NNFI, RMSEA, RMS theta and
SRMR indexes (more details on these indexes can be found in
{help plssemc_postestimation##Henseler2021:Henseler 2021}).

{pstd} {cmd:estat blindfolding}
computes the Q2 and q2 values ({help plssemc_postestimation##Geisser1974:Geisser 1974};
{help plssemc_postestimation##Stone1974:Stone 1974}) to assess the out-of-sample
predictive relevance. These measures are obtained by using the blindfolding procedure
for a specified omission distance.

 
{marker options_estat}{...}
{title:Options for estat}

{phang}
{opt boot(#)},
an option used with {cmd:estat indirect}, allows to estimate the indirect
effects using bootstrap; the number of replications is specified via {#}.

{phang}
{opt seed(#)},
an option used with {cmd:estat indirect}, allows to set the bootstrap seed
number.

{phang}
{opt level(#)},
an option used with {cmd:estat indirect} and {cmd:estat ci}, allows to set
the confidence level to use for indirect effects confidence intervals;
default is {cmd:0.95}.

{phang}
{opt digits(#)},
specifies the number of decimal digits to display in the output.

{phang}
{opt plot},
an option used with {cmd:estat total}, provides a graphical representation of
the total effects decomposition.

{phang}
{opt seed(#)},
allows to set the seed for reproducing results.

{phang}
{opt cutoff(#)},
an option used with {cmd:estat htmt}, specifies the cutoff value to use for
showing the heterotrait-monotrait ratios.

{phang}
{opt type(citype)},
an option used with {cmd:estat ci}, allows to choose the type of confidence
intervals to compute; available methods are {cmd:standard_z} (assumes
a standard normal distribution), {cmd:standard_t} (assumes
a t distribution with {it:n} - 1 degrees of freedom), {cmd:percentile} (utilizes
the quantiles of the distribution of the bootstrap resample estimates), and
{cmd:bc} (bias-corrected confidence interval).

{phang}
{opt distance(#)},
an option used with {cmd:estat blindfolding}, allows to set the omission distance
to use in the blindfolding procedure; it must be an integer larger than 1
and it has to be chosen so that the number of observations used in model
estimation divided by the omission distance is not an integer.


{marker syntax_predict}{...}
{marker predict}{...}
{title:Syntax for predict}

{p 8 16 2}
{cmd:predict} [{cmd:,} {it:statistic} {opt noout:er} {opt noin:ner}]

{synoptset 20 tabbed}{...}
{synopthdr :statistic}
{synoptline}
{syntab :Main}
{synopt :{cmd:xb}}linear predictions{p_end}
{synopt :{opt res:iduals}}residuals{p_end}
{synoptline}


{marker des_predict}{...}
{title:Description for predict}

{pstd}
{cmd:predict} creates new variables containing linear predictions and residuals. These
quantities are provided only for reflective blocks of manifest variables in the
measurement/outer model and for endogenous latent variables in the structural/inner model.

{pstd}
The newly computed predictions will replace those already present in the data set.


{marker options_predict}{...}
{title:Options for predict}

{dlgtab:Main}

{phang}
{opt xb} calculates the linear predictions (fitted values).

{phang}
{opt residuals} calculates the residuals.

{dlgtab:Options}

{phang}
{opt nooouter} fitted values and residuals for the measurement/outer model
are not saved in the data set.

{phang}
{opt nooinner} fitted values and residuals for the structural/inner model
are not saved in the data set.


{marker results_ci}{...}
{title:Stored results for estat ci}

{pstd}
{cmd:estat ci} stores the following in {cmd:r()}:

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Scalars}{p_end}
{synopt:{cmd:e(level)}}confidence level{p_end}

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Macros}{p_end}
{synopt:{cmd:e(ci_type)}}type of confidence intervals{p_end}

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Matrices}{p_end}
{synopt:{cmd:r(path_ci)}}matrix of path coefficients confidence intervals{p_end}
{synopt:{cmd:r(load_ci)}}matrix of loadings confidence intervals{p_end}
{p2colreset}{...}


{marker results_htmt}{...}
{title:Stored results for estat htmt}

{pstd}
{cmd:estat htmt} stores the following in {cmd:r()}:

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Matrices}{p_end}
{synopt:{cmd:r(htmt)}}matrix of heterotrait-monotrait ratios (HTMT){p_end}
{synopt:{cmd:r(htmt2)}}matrix of advanced heterotrait-monotrait ratios (HTMT2){p_end}
{p2colreset}{...}


{marker results_f2}{...}
{title:Stored results for estat f2}

{pstd}
{cmd:estat f2} stores the following in {cmd:r()}:

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Matrices}{p_end}
{synopt:{cmd:r(f2)}}matrix of the Cohen's f^2 effect sizes{p_end}
{p2colreset}{...}


{marker results_ic}{...}
{title:Stored results for estat ic}

{pstd}
{cmd:estat ic} stores the following in {cmd:r()}:

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Matrices}{p_end}
{synopt:{cmd:r(ic)}}matrix of the model's information and selection criteria{p_end}
{p2colreset}{...}


{marker results_fit}{...}
{title:Stored results for estat fit}

{pstd}
{cmd:estat fit} stores the following in {cmd:r()}:

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Matrices}{p_end}
{synopt:{cmd:r(dist)}}matrix of the model's distance measures{p_end}
{synopt:{cmd:r(fit)}}matrix of the model's fit measures{p_end}
{p2colreset}{...}


{marker results_blind}{...}
{title:Stored results for estat blindfolding}

{pstd}
{cmd:estat blindfolding} stores the following in {cmd:r()}:

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Matrices}{p_end}
{synopt:{cmd:r(Q2)}}matrix with the Q2 measures{p_end}
{synopt:{cmd:r(q2)}}matrix with the q2 measures{p_end}
{p2colreset}{...}


{marker results_predict}{...}
{title:Stored results for predict}

{pstd}
{cmd:predict} stores the following in {cmd:r()}:

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Matrices}{p_end}
{synopt:{cmd:r(fitted)}}matrix of fitted values for the outer and inner models{p_end}
{synopt:{cmd:r(residuals)}}matrix of residuals for the outer and inner models{p_end}
{p2colreset}{...}


{marker examples}{...}
{title:Examples}

    {hline}
{pstd}Setup{p_end}
{phang2}{cmd:. sysuse workout2, clear}{p_end}

{pstd}Model estimation{p_end}
{phang2}{cmd:. plssemc (Attractive > face sexy) (Appearance > body appear attract) (Muscle > muscle strength endur) (Weight > lweight calories cweight), structural(Appearance Attractive, Muscle Appearance, Weight Appearance)}{p_end}

{pstd}Discriminant validity assessment (HTMT){p_end}
{phang2}{cmd:. estat htmt}{p_end}

{pstd}Structural model assessment{p_end}
{phang2}{cmd:. estat vif}{p_end}
{phang2}{cmd:. estat f2}{p_end}
{phang2}{cmd:. estat ic}{p_end}

{pstd}Indirect effects{p_end}
{phang2}{cmd:. estat indirect, effects(Muscle Appearance Attractive, Weight Appearance Attractive)}{p_end}

{pstd}Confidence intervals{p_end}
{phang2}{cmd:. estat ci, type(bc) digits(5)}{p_end}

{pstd}Predictions{p_end}
{phang2}{cmd:. predict, xb residuals}{p_end}
{phang2}{cmd:. describe *_hat *_res}{p_end}

    {hline}


{marker authors}{...}
{title:Authors}

{pstd} Sergio Venturini{break}
Department of Economics and Social Sciences{break}
Università Cattolica del Sacro Cuore, Italy{break}
{browse "mailto:sergio.venturini@unicatt.it":sergio.venturini@unicatt.it}{break}

{pstd} Mehmet Mehmetoglu{break}
Department of Psychology{break}
Norwegian University of Science and Technology{break}
{browse "mailto:mehmetm@svt.ntnu.no":mehmetm@svt.ntnu.no}{break}
{p_end}


{marker references}{...}
{title:References} 

{marker BaronKenny1986}{...}
{phang}
Baron, R. M., and Kenny, D. A. 1986. The Moderator-Mediator Variable Distinction in Social Psychological
Research: Conceptual, Strategic, and Statistical Considerations. Journal of
Personality and Social Psychology, 51, 1173-1182.

{marker Cohen1988}{...}
{phang}
Cohen, J. 1988. {it:Statistical Power Analysis for the Behavioral Sciences}. Mahwah, NJ: Erlbaum.

{marker Geisser1974}{...}
{phang}
Geisser, S. 1974. A predictive approach to the random effects model. Biometrika, 61, 101–107.

{marker Hahnetal2002}{...}
{phang}
Hahn, C., Johnson, M. D., Herrmann, A., and Huber, F. 2002. Capturing Customer Heterogeneity Using a
Finite Mixture PLS Approach. Schmalenbach Business Review, 54, 243-269.

{marker Hairetal2022}{...}
{phang}
Hair, J. F., Hult, G. T. M., Ringle, C. M., and Sarstedt, M. 2022. {it:A Primer on Partial Least Squares Structural Equation Modeling (PLS-SEM)}. Third edition. Sage.

{marker Hairetal2018}{...}
{phang}
Hair, J. F., Sarstedt, M., Ringle, C. M., and Gudergan, S. P. 2018. {it:Advanced Issues in Partial Least Squares Structural Equation Modeling}. Sage.

{marker Henseler2021}{...}
{phang}
Henseler, J. 2021. {it:Composite-Based Structural Equation Modeling}. The Guilford Press.

{marker Iacobuccietal2007}{...}
{phang}
Iacobucci, D., Saldanha, N., & Deng, X. 2007. A meditation on mediation: evidence that structural equation models perform better than regressions. Journal of
Consumer Psychology, 17, 140-154.

{marker MehmetogluVenturini2021}{...}
{phang}
Mehmetoglu, M., and Venturini, S. 2021. {it:Structural Equation Modelling with Partial Least Squares Using Stata and R}. CRC Press.

{marker Ringleetal2014}{...}
{phang}
Ringle, C. M., Sarstedt, M., and Schlittgen, R. 2014. Genetic algorithm segmentation
in partial least squares structural equation modeling. OR Spectrum, 36, 251–276.

{marker Sharma2019}{...}
{phang}
Sharma, P., Sarstedt, M., Shmueli, G., Kim, K. H., Thiele, K. O. 2019. PLS-Based Model Selection: The Role of Alternative Explanations in Information Systems Research. Journal of the Association for Information Systems, 20(4), 346-397.

{marker Sobel1982}{...}
{phang}
Sobel, M. N. 1982. Asymptotic Confidence Intervals for Indirect Effects in Structural Equations
Models. In Leinhart, S. (ed.), {it:Sociological Methodology}, pp. 290-312. Jossey-Bass.

{marker Stone1974}{...}
{phang}
Stone, S. 1974. Cross-validatory choice and assessment of statistical pre- dictions.
Journal of the Royal Statistical Society, 36, 111–147.

{marker Trinchera2007}{...}
{phang}
Trinchera, L. 2007. {it:Unobserved Heterogeneity in Structural Equation Models: a new approach to latent class detection in PLS Path Modeling}. Ph.D. Thesis.

{marker VanderWeele2015}{...}
{phang}
VanderWeele, T. J. 2015. {it:Explanation in Causal Inference}. Oxford University Press.

{marker Zhaoetal2010}{...}
{phang}
Zhao, X., Lynch, J. G. Jr., & Chen, Q. 2010. Reconsidering Baron and Kenny: myths and
truths about mediation analysis. Journal of Consumer Research, 37, 197-206.
{p_end}
