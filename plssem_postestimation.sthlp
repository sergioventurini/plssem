{smcl}
{* *! version 0.0.1  31May2019}{...}
{vieweralsosee "plssem" "help plssem"}{...}
{vieweralsosee "plssemplot" "help plssemplot"}{...}
{viewerjumpto "Postestimation commands" "plssem postestimation##description"}{...}
{viewerjumpto "estat" "plssem postestimation##syntax_estat"}{...}
{viewerjumpto "estat options" "plssem postestimation##options_estat"}{...}
{viewerjumpto "estat indirect stored results" "plssem postestimation##results_indirect"}{...}
{viewerjumpto "estat mediate stored results" "plssem postestimation##results_meidate"}{...}
{viewerjumpto "predict" "plssem postestimation##syntax_predict"}{...}
{viewerjumpto "predict options" "plssem postestimation##options_predict"}{...}
{viewerjumpto "predict stored results" "plssem postestimation##results_predict"}{...}
{viewerjumpto "Examples" "plssem postestimation##examples"}{...}
{viewerjumpto "Authors" "plssem postestimation##authors"}{...}
{viewerjumpto "References" "plssem postestimation##references"}{...}
{title:Title}

{p 4 18 2}
{hi:plssem postestimation} {hline 2} Postestimation tools for {helpb plssem}


{marker description}{...}
{title:Postestimation commands}

{pstd}
The following postestimation commands are of special interest after
{cmd:plssem}:

{synoptset 22 tabbed}{...}
{p2coldent:Command}Description{p_end}
{synoptline}
{synopt:{helpb plssem postestimation##indirect:estat indirect}}estimation and
	inference for indirect effects{p_end}
{synopt:{helpb plssem postestimation##total:estat total}}decomposition of
	total effects{p_end}
{synopt:{helpb plssem postestimation##mediate:estat mediate}}testing of
	a mediation effect{p_end}
{p2coldent:* {helpb plssem postestimation##vif:estat vif}}variance inflation
	factors for the structural model equations sample{p_end}
{p2coldent:* {helpb plssem postestimation##unobshet:estat unobshet}}unobserved
	heterogeneity assessment{p_end}
{synoptline}
{p2colreset}{...}
{p 4 6 2}
* {cmd:estat vif} and {cmd:estat unobshet} are not available for models fitted using bootstrap.
{p_end}

{pstd}
The following standard postestimation commands are also available:

{synoptset 20 tabbed}{...}
{p2coldent :Command}Description{p_end}
{synoptline}
{synopt :{helpb plssem postestimation##predict:predict}}fitted values and residuals{p_end}
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


{marker mediate}{...}
{pstd}
Display the mediation analysis using partial least squares structural equation modelling

{p 8 14 2}
{cmd:estat} {cmdab:me:diate}
{cmdab:ind:ep(}{it:varname}{cmd:)} {cmdab:med:(}{it:varname}{cmd:)}
{cmdab:dep:(}{it:varname}{cmd:)} [{opt br:eps(#)} {opt s:eed(#)} {opt zlc:}
{opt rit:} {opt rid:} {opt bc:a} {opt l:evel(#)} {opt dig:its(#)}]


{marker vif}{...}
{pstd}
Display the variance inflation factors for the structural model equations

{p 8 14 2}
{cmd:estat} {cmdab:vi:f}
[{cmd:,} {opt dig:its(#)}]


{marker unobshet}{...}
{pstd}
Display the assessment for the presence of unobserved heterogeneity

{p 8 14 2}
{cmd:estat} {cmdab:un:obshet},
{cmdab:m:ethod(}{it:methodname}{cmd:)} [{opt n:umclass(#)} {opt maxcl:ass(#)}
{opt d:endrogram} {opt maxit:er(#)} {opt s:top(#)}  {opt t:est} {opt r:eps(#)}
{opt res:tart(#)} {opth gr:oups(numlist)} {opt pop:size(#)} {opt numg:en(#)}
{opt pm:ut(#)} {opt pt:ransf(#)} {opt maxitg:as(#)} {opt se:ed(#)}
{opt p:lot} {cmdab:name(}{it:varname}{cmd:)} {opt dig:its(#)}]


{marker desc_estat}{...}
{title:Description for estat}

{pstd}
{cmd:estat indirect}
estimates the (standardized) indirect effects and the corresponding tests of significance
using the Sobel's {it:z} statistic (default) as well as the bootstrap
approach ({help plssem_postestimation##Sobel1982:Sobel 1982},
{help plssem_postestimation##BaronKenny1986:Baron and Kenny 1986}, 
{help plssem_postestimation##VanderWeele2015:VanderWeele 2015}). The command can
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

{pstd}
{cmd:estat mediate}
conducts a mediation analysis based on a fitted PLS-SEM model estimated using the
{helpb plssem} command. Two methods are implemented, the 
{help plssem_postestimation##BaronKenny1986:Baron and Kenny (1986)} approach adjusted by
{help plssem_postestimation##Iacobuccietal2007:Iacobucci et al. 2007},
and that by {help plssem_postestimation##Zhaoetal2010:Zhao et al. 2010}.

{pstd} {cmd:estat vif}
computes the variance inflation factors (VIFs) for the independent variables
of the equations in the structural part of a PLS-SEM model. With the
{cmd:digit(#)} sub-option you change the number decimals digits displayed.

{pstd} {cmd:estat unobshet}
assesses the presence of unobserved heterogeneity in the fitted PLS-SEM model
using the {it:methodname} approach. Currently, the REBUS-PLS
({help plssem_postestimation##Trinchera2007:Trinchera 2007}), FIMIX-PLS
({help plssem_postestimation##Hahnetal2002:Hahn et al. 2007}) and PLS-GAS
({help plssem_postestimation##Ringleetal2014:Ringle et al. 2014}) approaches are implemented.
 
 
{marker options_estat}{...}
{title:Options for estat}

{phang}
{opt boot(#)},
an option used with {cmd:estat indirect}, allows to estimate the indirect
effects using bootstrap; the number of replications is specified via {#}.

{phang}
{opt seed(#)},
an option used with {cmd:estat indirect} and {cmd:estat mediate}, allows to set
the bootstrap seed number.

{phang}
{opt level(#)},
an option used with {cmd:estat indirect} and {cmd:estat mediate}, allows to set
the confidence level to use for indirect effects confidence intervals; default
is {cmd:0.95}.

{phang}
{opt digits(#)},
specifies the number of decimal digits to display in the output; default is {cmd:3}.

{phang}
{opt plot},
an option used with {cmd:estat total}, provides a graphical representation of
the total effects decomposition.

{phang}
{opt method(methodname)},
an option used with {cmd:estat unobshet}, allows choosing the method to use for
assessing the presence of unobserved heterogeneity; available methods are {cmd:rebus},
{cmd:fimix} and {cmd:gas}; default is {cmd:rebus}.

{phang}
{opt numclass(#)},
an option used with {cmd:estat unobshet}, allows to manually set the number of
classes to use in the REBUS-PLS analysis; minimum is 2. If not specified, the
number of classes is automatically chosen based on the Calinski-Harabasz
pseudo-F index stopping rule as implemented in {helpb cluster stop:cluster stop}. In
this case, a Ward hierarchical clustering algorithm is used.

{phang}
{opt maxclass(#)},
an option used with {cmd:estat unobshet}, allows to set the maximum number of
classes for which the automatic stopping rule is to be computed when
{cmd: numclass(#)} is not specified; default is {cmd:20}.

{phang}
{opt dendrogram},
an option used with {cmd:estat unobshet}, allows to visualize the dendrogram
for a Ward hierarchical clustering algorithm on the base of the residuals of
the global model. The dendrogram allows to assess the quality of choice for the
number of classes. This option is available only when {cmd: method(rebus)} is
chosen.

{phang}
{opt maxiter(#)},
an option used with {cmd:estat unobshet}, allows to set the maximum number
of iterations the REBUS-PLS algorithm runs; default is {cmd:50}.

{phang}
{opt stop(#)},
an option used with {cmd:estat unobshet}, allows to set the stopping rule for
the chosen algorithm; for REBUS-PLS this refers to the stability in class
composition from one iteration to the other (i.e. percentage of units changing
class at each iteration); for FIMIX-PLS it refers to the (absolute) change in the
complete loglikelihood value in the EM alogrithm. Default is {cmd:0.005}
(i.e. 0.5%) when REBUS-PLS is chosen and {cmd:1e-5} when FIMIX-PLS is used.

{phang}
{opt test},
an option used with {cmd:estat unobshet}, allows to specify whether a
permutation test for the Global Quality Index (GQI) of a REBUS-PLS solution must
be performed. This option is available only when {cmd: method(rebus)} is
chosen.

{phang}
{opt reps(#)},
an option used with {cmd:estat unobshet}, allows to set the number of replications
of the permutation test on the GQI; default is {cmd:50}. This option is available
only when {cmd: method(rebus)} is chosen.

{phang}
{opt restart(#)},
an option used with {cmd:estat unobshet} when the FIMIX-PLS method is chosen,
allows to set the number of EM algortihm runs in the estimation of the
mixture model's parameters; default is {cmd:10}.

{phang}
{opth groups(numlist)},
an option used with {cmd:estat unobshet} when the FIMIX-PLS method is chosen,
allows to compute the solution for a range of selected group number values and
compare them using a set of fit indices, namely Akaike's information criterion
(AIC), modified AIC with factor 3 (AIC3), modified AIC with factor 4 (AIC4),
Bayesian information criterion (BIC), consistent AIC (CAIC),
Hannan-Quinn criterion (HQ), minimum description length with factor 5 (MDL5),
log-likelihood (LnL), entropy statistic (EN), non-fuzzy index (NFI),
normalized entropy criterion (NEC).

{phang}
{opt popsize(#)},
an option used with {cmd:estat unobshet} when the PLS-GAS method is chosen,
allows to set the size of each generation in the genetic algorithm; deafult
is 100.

{phang}
{opt numgen(#)},
an option used with {cmd:estat unobshet} when the PLS-GAS method is chosen,
allows to set the number of generations to create during the genetic algorithm;
deafult is 1000.

{phang}
{opt pmut(#)},
an option used with {cmd:estat unobshet} when the PLS-GAS method is chosen,
allows to set the probability to mutate for a chromosome (i.e. an individual in a
a generation) in the genetic algorithm; deafult is 0.3.

{phang}
{opt ptransf(#)},
an option used with {cmd:estat unobshet} when the PLS-GAS method is chosen,
allows to set the probability to mutate for a single gene in a chromosome during
the genetic algorithm calculation; deafult is 0.1.

{phang}
{opt maxitgas(#)},
an option used with {cmd:estat unobshet}, allows to set the maximum number
of iterations the PLS-GAS algorithm; default is {cmd:30}.

{phang}
{opt seed(#)},
allows to set the seed for reproducing results.

{phang}
{opt plot},
an option used with {cmd:estat unobshet} when the REBUS-PLS method is chosen,
allows to visualize the empirical distribution (i.e. the histogram)
corresponding to the replications of the permutation test on the GQI.

{phang}
{cmdab:name(}{it:varname}{cmd:)},
an option used with {cmd:estat unobshet}, allows to set the name of the variable
that will contain the final classification obtained.

{phang}
{cmdab:dep(}{it:varname}{cmd:)},
an option used with {cmd:estat mediate}, specifies the name of the dependent
variable to use in the mediation analysis.

{phang}
{cmdab:med(}{it:varname}{cmd:)},
an option used with {cmd:estat mediate}, specifies the name of the mediator
variable to use in the mediation analysis.

{phang}
{cmdab:indep(}{it:varname}{cmd:)},
an option used with {cmd:estat mediate}, specifies the name of the independent
variable to use in the mediation analysis.

{phang}
{opt breps(#)},
an option used with {cmd:estat mediate}, specifies the number of bootstrap
replications to be performed. The default is {cmd:50}.

{phang}
{opt zlc},
an option used with {cmd:estat mediate} providing the approach by
{help plssem_postestimation##Zhaoetal2010:Zhao et al. 2010} together with
the {help plssem_postestimation##BaronKenny1986:Baron and Kenny (1986)} approach
adjusted by {help plssem_postestimation##Iacobuccietal2007:Iacobucci et al. 2007}
(default).

{phang}
{opt rit},
an option used with {cmd:estat mediate} providing the ratio of the indirect
effect to the total effect.

{phang}
{opt rid},
an option used with {cmd:estat mediate} providing the ratio of the indirect
effect to the direct effect.

{phang}
{opt bca},
an option used with {cmd:estat mediate} providing the bias-corrected accelerated
(BCa) bootstrap confidence intervals instead of the percentile confidence
intervals (default).


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


{marker results_indirect}{...}
{title:Stored results for estat indirect}

{pstd}
{cmd:estat indirect} stores the following in {cmd:r()}:

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Matrices}{p_end}
{synopt:{cmd:r(indirect)}}matrix of indirect effect testing results{p_end}
{p2colreset}{...}


{marker results_mediate}{...}
{title:Stored results for estat mediate}

{pstd}
{cmd:estat mediate} stores the following in {cmd:r()}:

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Matrices}{p_end}
{synopt:{cmd:r(mediate)}}matrix of indirect effect testing results{p_end}
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
{phang2}{cmd:. plssem (Attractive > face sexy) (Appearance > body appear attract) (Muscle > muscle strength endur) (Weight > lweight calories cweight), structural(Appearance Attractive, Muscle Appearance, Weight Appearance)}{p_end}

{pstd}Multicollinearity assessment{p_end}
{phang2}{cmd:. estat vif}{p_end}

{pstd}Indirect effects{p_end}
{phang2}{cmd:. estat indirect, effects(Muscle Appearance Attractive, Weight Appearance Attractive)}{p_end}

{pstd}Predictions{p_end}
{phang2}{cmd:. predict, xb residuals}{p_end}
{phang2}{cmd:. describe *_hat *_res}{p_end}

{pstd}Assessment of unobserved heterogeneity using REBUS-PLS{p_end}
{phang2}{cmd:. estat unobshet, test reps(200) plot}{p_end}

{pstd}Assessment of unobserved heterogeneity using FIMIX-PLS{p_end}
{phang2}{cmd:. estat unobshet, method(fimix) groups(1/5) stop(1e-5) restart(3)}{p_end}

    {hline}


{marker authors}{...}
{title:Authors}

{pstd} Sergio Venturini{break}
Department of Decision Sciences{break}
Università Bocconi, Italy{break}
{browse "mailto:sergio.venturini@unibocconi.it":sergio.venturini@unibocconi.it}{break}

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

{marker Hahnetal2002}{...}
{phang}
Hahn, C., Johnson, M. D., Herrmann, A., and Huber, F. 2002. Capturing Customer Heterogeneity Using a
Finite Mixture PLS Approach. Schmalenbach Business Review, 54, 243-269.

{marker Hairetal2017}{...}
{phang}
Hair, J. F., Hult, G. T. M., Ringle, C. M., and Sarstedt, M. 2017. {it:A Primer on Partial Least Squares Structural Equation Modeling (PLS-SEM)}. Second edition. Sage.

{marker Hairetal2018}{...}
{phang}
Hair, J. F., Sarstedt, M., Ringle, C. M., and Gudergan, S. P. 2018. {it:Advanced Issues in Partial Least Squares Structural Equation Modeling}. Sage.

{marker Iacobuccietal2007}{...}
{phang}
Iacobucci, D., Saldanha, N., & Deng, X. 2007. A meditation on mediation: evidence
that structural equation models perform better than regressions. Journal of
Consumer Psychology, 17, 140-154.

{marker Ringleetal2014}{...}
{phang}
Ringle, C. M., Sarstedt, M., and Schlittgen, R. 2014. Genetic algorithm segmentation
in partial least squares structural equation modeling. OR Spectrum, 36, 251–276.

{marker Sobel1982}{...}
{phang}
Sobel, M. N. 1982. Asymptotic Confidence Intervals for Indirect Effects in Structural Equations
Models. In Leinhart, S. (ed.), {it:Sociological Methodology}, pp. 290-312. Jossey-Bass.

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
