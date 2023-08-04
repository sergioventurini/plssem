{smcl}
{* *! version 0.4.0  26Apr2023}{...}
{vieweralsosee "plssemc_mi" "help plssemc_mi"}{...}
{viewerjumpto "Syntax" "plssem_mi##syntax"}{...}
{viewerjumpto "Description" "plssem_mi##description"}{...}
{viewerjumpto "Options" "plssem_mi##options"}{...}
{viewerjumpto "Examples" "plssem_mi##examples"}{...}
{viewerjumpto "Authors" "plssem_mi##authors"}{...}
{viewerjumpto "Stored results" "plssem_mi##results"}{...}
{viewerjumpto "References" "plssem_mi##references"}{...}
{title:Title}

{p 4 17 2}
{hi:plssem_mi} {hline 2} Multiple imputation analysis of partial least
squares structural equation modelling (PLS-SEM)


{marker syntax}{...}
{title:Syntax}

{pstd}
Multiple imputation analysis of partial least squares structural equation modeling

{p 8 12 2}
{cmd:plssem_mi} (LV1 > indblock1) (LV2 > indblock2) (...) {ifin}
[{cmd:,} structural(LV2 LV1, ...) {it:{help plssem##plssemopts:plssem_options}}
mioptions({it:{help mi_estimate##options:mi_options}}) midisplay]

{p 8 12 2}
{cmd:plssemc_mi} (LV1 > indblock1) (LV2 > indblock2) (...) {ifin}
[{cmd:,} structural(LV2 LV1, ...) {it:{help plssem##plssemopts:plssem_options}}
mioptions({it:{help mi_estimate##options:mi_options}}) midisplay]

{synoptset 24 tabbed}{...}
{marker plssem_miopts}{...}
{synopthdr}
{synoptline}
{synopt:{it:{help plssem##plssemopts:plssem_options}}}options to pass to {helpb plssem} and
{helpb plssemc} for performing the primary analysis{p_end}
{synopt:{opth mio:ptions(mi_estimate##options:mi_options)}}options to pass to {helpb mi estimate}{p_end}
{synopt:{opt mid:isplay}}display the results using the {helpb mi estimate} command layout{p_end}
{synoptline}

{p 4 4 2}
The {cmd:plssem_mi} and {cmd:plssemc_mi} commands provide the completed-data analysis
(estimation) and pooling steps of multiple imputation for a PLS-SEM or PLSc-SEM analysis.


{marker description}{...}
{title:Description}

{pstd}First, the data must be {helpb mi set}. Then, a set of imputed data sets must be generated
using any of the {helpb mi impute} commands. Finally, once missing data have been imputed,
{cmd:plssem_mi} or {cmd:plssemc_mi} can be used to perform the primary analysis (i.e., PLS-SEM or
PLSc-SEM) on each imputed data set and combine them into a single multiple-imputation result.

{pstd}The syntax for the PLS-SEM or PLSc-SEM analysis mirrors exactly that of the corresponding
{helpb plssem} and {helpb plssemc} commands.


{marker options}{...}
{title:Options}

{phang}{opth mio:ptions(mi_estimate##options:mi_options)}
accepts all the options available for {helpb mi estimate}, with the exception of {opt cmdok}
which has been already included internally.

{phang}{opt mid:isplay}
shows the results using the layout implemented in {helpb mi estimate} instead of that used
by {helpb plssem} and {helpb plssemc}.


{marker examples}{...}
{title:Examples}

    {hline}
{pstd}Setup{p_end}
{phang2}{cmd:. sysuse ecsimobi, clear}{p_end}

{pstd}Generate missing values (MCAR){p_end}
{phang2}{cmd:. local prop_miss 0.92}{p_end}

{phang2}  {cmd:. foreach vn of varlist CUEX1-CUEX3 {c -(}}{p_end}
{phang2}    {cmd:2.       generate R_mcar = 1 - rbinomial(1, `prop_miss')}{p_end}
{phang2}    {cmd:3.       replace `vn' = . if R_mcar == 1}{p_end}
{phang2}    {cmd:4.       drop R_mcar}{p_end}
{phang2}    {cmd:5. {c )-}}{p_end}

{phang2}  {cmd:. foreach vn of varlist IMAG1-IMAG5 {c -(}}{p_end}
{phang2}    {cmd:2.       generate R_mcar = 1 - rbinomial(1, `prop_miss')}{p_end}
{phang2}    {cmd:3.       replace `vn' = . if R_mcar == 1}{p_end}
{phang2}    {cmd:4.       drop R_mcar}{p_end}
{phang2}    {cmd:5. {c )-}}{p_end}

{phang2}  {cmd:. foreach vn of varlist PERV1-PERV2 {c -(}}{p_end}
{phang2}    {cmd:2.       generate R_mcar = 1 - rbinomial(1, `prop_miss')}{p_end}
{phang2}    {cmd:3.       replace `vn' = . if R_mcar == 1}{p_end}
{phang2}    {cmd:4.       drop R_mcar}{p_end}
{phang2}    {cmd:5. {c )-}}{p_end}

{pstd}Set up mi{p_end}
{phang2}{cmd:. mi set wide}{p_end}
{phang2}{cmd:. mi register imputed CUEX1-CUEX3 IMAG1-IMAG5 PERV1-PERV2}{p_end}
{phang2}{cmd:. mi register passive CUSA* CUSCO CUSL* PERQ*}{p_end}

{pstd}Impute missing values{p_end}
{phang2}{cmd:. mi impute mvn CUEX* IMAG* PERV* = CUSA* CUSCO CUSL* PERQ*, add(20) dots burnin(1000) replace}{p_end}

{pstd}Perform multiple imputation{p_end}
{phang2}{cmd:. plssem_mi (Expectation > CUEX1 CUEX2 CUEX3) (Satisfaction > CUSA1 CUSA2 CUSA3) (Complaints > CUSCO) (Loyalty < CUSL1 CUSL2 CUSL3) (Image < IMAG1 IMAG2 IMAG3 IMAG4 IMAG5) (Quality > PERQ1 PERQ2 PERQ3 PERQ4 PERQ5 PERQ6 PERQ7) (Value < PERV1 PERV2), structural(Expectation Image, Quality Expectation, Value Expectation Quality, Satisfaction Value Quality Image Expectation, Complaints Satisfaction, Loyalty Complaints Satisfaction Image) wscheme("centroid") boot(200) mioptions(dots)}{p_end}

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


{marker results}{...}
{title:Stored results}

{pstd}
{cmd:plssem_mi} stores in {cmd:e()} the same results provided by {helpb mi_estimate##results:mi estimate}.


{marker references}{...}
{title:References} 

{marker Hairetal2022}{...}
{phang}
Hair, J. F., Hult, G. T. M., Ringle, C. M., and Sarstedt, M. 2022. {it:A Primer on Partial Least Squares Structural Equation Modeling (PLS-SEM)}. Third edition. Sage.

{marker Lohmoller1989}{...}
{phang}
Lohmöller, J. B. 1989. {it:Latent Variable Path Modeling with Partial Least Squares}. Heidelberg: Physica.

{marker Wold1975}{...}
{phang}
Wold, H. O. A. 1975. Path Models with Latent Variables: The NIPALS Approach.
In Blalock, H. M., Aganbegian, A., Borodkin, F. M., Boudon, R., and Cappecchi, V. (ed.), {it:Quantitative Sociology} (pp. 307-359). New York: Academic Press.
{p_end}
