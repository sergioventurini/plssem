{smcl}
{* *! version 0.0.1  30Dec2016}{...}
{vieweralsosee "plssem" "help plssem"}{...}
{vieweralsosee "plssemplot" "help plssemplot"}{...}
{viewerjumpto "Postestimation commands" "plssem postestimation##description"}{...}
{viewerjumpto "estat" "plssem postestimation##syntax_estat"}{...}
{viewerjumpto "estat options" "plssem postestimation##options_estat"}{...}
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
{p2coldent:* {helpb plssem postestimation##vif:estat vif}}variance inflation
	factors for the structural model equations sample{p_end}
{synoptline}
{p2colreset}{...}
{p 4 6 2}
* {cmd:estat vif} is not available after models fitted using bootstrap.
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


{marker vif}{...}
{pstd}
Display the variance inflation factors for the structural model equations

{p 8 14 2}
{cmd:estat} {cmdab:vi:f}
[{cmd:,} {opt dig:its(#)}]


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
{cmd:predict} creates new variables containing linear predictions and residuals.
These quantities are provided only for reflective blocks of manifest variables
in the measurement/outer model and for endogenous latent variables in the
structural/inner model.

{pstd}
The computed variables will replace those already present in the data set.


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

{phang}
{stata sysuse workout2, clear}{p_end}

{phang}
{stata plssem (Attractive > face sexy) (Appearance > body appear attract) (Muscle > muscle strength endur) (Weight > lweight calories cweight), structural(Appearance Attractive, Muscle Appearance, Weight Appearance)}{p_end}

{phang}
{stata estat indirect, effects(Muscle Appearance Attractive, Weight Appearance Attractive)}{p_end}

{phang}
{stata estat total, plot}{p_end}

{phang}
{stata predict, xb residuals}{p_end}
{phang}
{stata describe *_hat *_res}{p_end}


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

{marker Hairetal2017}{...}
{phang}
Hair, J. F., Hult, G. T. M., Ringle, C. M., and Sarstedt, M. 2017. {it:A Primer on Partial Least Squares Structural Equation Modeling (PLS-SEM)}. 2nd edition. Sage.

{marker Sobel1982}{...}
{phang}
Sobel, M. N. 1982. Asymptotic Confidence Intervals for Indirect Effects in Structural Equations
Models. In Leinhart, S. (ed.), {it:Sociological Methodology}, pp. 290-312. Jossey-Bass.

{marker VanderWeele2015}{...}
{phang}
VanderWeele, T. J. 2015. {it:Explanation in Causal Inference}. Oxford University Press.
{p_end}
