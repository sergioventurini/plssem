{smcl}
{* *! version 0.4.0  26Apr2023}{...}
{vieweralsosee "plssemc" "help plssemc"}{...}
{vieweralsosee "plssemplot" "help plssemplot"}{...}
{viewerjumpto "Postestimation commands" "plssemc postestimation##description"}{...}
{viewerjumpto "estat" "plssemc postestimation##syntax_estat"}{...}
{viewerjumpto "estat options" "plssemc postestimation##options_estat"}{...}
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
[{opt cut:off(#)}]


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

{pstd} {cmd:estat vif}
computes the variance inflation factors (VIFs) for the independent variables
of the equations in the structural part of a PLS-SEM model. With the
{cmd:digit(#)} sub-option you change the number decimals digits displayed.

{pstd} {cmd:estat htmt}
assesses discriminant validity using heterotrait-monotrait ratios of correlations; both
the HTMT (arithmetic average of correlations) and HTMT2 (geometric average of correlations)
approaches are available.

 
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
an option used with {cmd:estat indirect}, allows to set the confidence level
to use for indirect effects confidence intervals; default is {cmd:0.95}.

{phang}
{opt digits(#)},
specifies the number of decimal digits to display in the output; default is {cmd:3}.

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

{pstd}Multicollinearity assessment{p_end}
{phang2}{cmd:. estat vif}{p_end}

{pstd}Indirect effects{p_end}
{phang2}{cmd:. estat indirect, effects(Muscle Appearance Attractive, Weight Appearance Attractive)}{p_end}

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
{p_end}
