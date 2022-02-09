{smcl}
{* *! version 0.0.1  06Oct2017}{...}
{vieweralsosee "plssem" "help plssem"}{...}
{vieweralsosee "plssem postestimation" "help plssem postestimation"}{...}
{viewerjumpto "Syntax" "plssemplot##syntax"}{...}
{viewerjumpto "Description" "plssemplot##description"}{...}
{viewerjumpto "Options" "plssemplot##options"}{...}
{viewerjumpto "Examples" "plssemplot##examples"}{...}
{viewerjumpto "Authors" "plssemplot##authors"}{...}
{title:Title}

{p 4 18 2}
{hi:plssemplot} {hline 2} Graph results from {helpb plssem} (loadings plots, etc.){p_end}


{marker syntax}{...}
{title:Syntax}

{p 8 21 2}{cmd:plssemplot} [{cmd:,} {help plssemplot##plssemplotopts:options}]

{synoptset 20 tabbed}{...}
{marker plssemplotopts}{...}
{synopthdr}
{synoptline}
{synopt:{cmdab:in:nermodel}}plot the structural (inner) model relationships{p_end}
{synopt:{cmdab:out:ermodel}}plot the measurement (outer) model relationships (not yet implemented){p_end}
{synopt:{opth st:ats(syntax##description_of_varlist:varname)}}produce the scatterplot matrix for the indicators of the selected
	latent variables{p_end}
{synopt:{cmdab:sc:ores}}produce the scatterplot matrix of the latent variables scores{p_end}
{synopt:{cmdab:c:rossloadings}}plot the model's cross loadings{p_end}
{synopt:{cmdab:l:oadings}}plot the model's outer loadings{p_end}
{synopt:{cmdab:outerw:eights}}plot the model's outer weights convergence path{p_end}
{synoptline}

{marker description}{...}
{title:Description}

{pstd}
{cmd:plssemplot} graphs the results of the immediately preceding {helpb plssem} command.


{marker options}{...}
{title:Options}

{phang}
{opt plssemplot, crossloadings} provides bar plots of the loadings of indicators for not only
their respective but all the other latent variables (i.e. cross loadings).

{phang}
{opt plssemplot, innermodel} produces a graphical representation version of the structural (inner)
part of the PLS-SEM model estimated. This command requires the installation of
the {cmd:nwcommands} suite ({stata nwinstall, all:click here to install}).

{phang}
{opt plssemplot, loadings} provides a bar plot of the loadings of indicators for their respective
latent variables.

{phang}
{opt plssemplot, outermodel} produces a visualized version of the measurement (outer) part of
the PLS-SEM model estimated. This feature will be available soon.

{phang}
{opt plssemplot, outerweights} produces a diagram showing the evolution over time of the
outer weights for the estimated  model. This diagram is useful for assessing convergence of
the estimation algorithm.
 
{phang}
{opt plssemplot, scores} provides the scatterplot matrix of the scores for the
latent variables defined in the PLS-SEM model.

{phang}
{opt plssemplot, stats(LV)} provides the scatterplot matrix for the indicators in
the block defining the latent variable {cmd:LV}.

 
{marker examples}{...}
{title:Examples}

    {hline}
{pstd}Setup{p_end}
{phang2}{cmd:. sysuse workout2, clear}{p_end}

{pstd}Model estimation{p_end}
{phang2}{cmd:. plssem (Attractive > face sexy) (Appearance > body appear attract) (Muscle > muscle strength endur) (Weight > lweight calories cweight), structural(Appearance Attractive, Muscle Appearance, Weight Appearance)}{p_end}

{pstd}Outer weights evolution{p_end}
{phang2}{cmd:. plssemplot, outerweights}{p_end}

{pstd}Scatterplot matrix for indicators related to a specific latent variable{p_end}
{phang2}{cmd:. plssemplot, stats(Appearance)}{p_end}

{pstd}Scatterplot matrix of latent variable scores{p_end}
{phang2}{cmd:. plssemplot, scores}{p_end}

{pstd}Plot of loadings and cross-loadings{p_end}
{phang2}{cmd:. plssemplot, loadings}{p_end}
{phang2}{cmd:. plssemplot, crossloadings}{p_end}

    {hline}


{marker authors}{...}
{title:Authors}

{pstd} Sergio Venturini{break}
Department of Economic and Social Sciences{break}
Università Cattolica del Sacro Cuore, Italy{break}
{browse "mailto:sergio.venturini@unicatt.it":sergio.venturini@unicatt.it}{break}

{pstd} Mehmet Mehmetoglu{break}
Department of Psychology{break}
Norwegian University of Science and Technology{break}
{browse "mailto:mehmetm@svt.ntnu.no":mehmetm@svt.ntnu.no}{break}
{p_end}
