cd "/Users/Sergio/Dropbox (Personal)/PLS-SEM"

/* Example 1 */
/* --------- */
use ./data/ecsimobi, clear
plssem (Expectation > CUEX1-CUEX3) (Satisfaction > CUSA1-CUSA3) ///
	(Complaints > CUSCO) (Loyalty > CUSL1-CUSL3) (Image > IMAG1-IMAG5) ///
	(Quality > PERQ1-PERQ7) (Value > PERV1-PERV2), ///
	structural(Expectation Image, Quality Expectation, Value Expectation Quality, ///
	Satisfaction Value Quality Image Expectation, Complaints Satisfaction, ///
	Loyalty Complaints Satisfaction Image) ///
	wscheme("path") digits(4) correlate(mv lv cross, cutoff(.3)) ///
	//noheader nomeastable nodiscrimtable nostructtable //boot(50)

estat indirect, effects(Loyalty Satisfaction Image, ///
	Value Quality Expectation) level(0.9) ///
	digits(5) //seed(123) boot(50) */

plssemplot, innermodel
* plssemplot, outermodel
plssemplot, outerweights

/* Example 2 */
/* --------- */
use ./data/ecsimobi, clear
plssem (Expectation > CUEX1-CUEX3) (Satisfaction > CUSA1-CUSA3) ///
	(Complaints > CUSCO) (Loyalty < CUSL1-CUSL3) (Image < IMAG1-IMAG5) ///
	(Quality > PERQ1-PERQ7) (Value < PERV1-PERV2), ///
	structural(Expectation Image, Quality Expectation, Value Expectation Quality, ///
	Satisfaction Value Quality Image Expectation, Complaints Satisfaction, ///
	Loyalty Complaints Satisfaction Image) ///
	wscheme("centroid")

/* Example 3 */
/* --------- */
use ./data/ecsimobi, clear
plssem (Expectation > CUEX1-CUEX3) (Satisfaction < CUSA1-CUSA3) ///
	(Quality > PERQ1-PERQ7), ///
	structural(Satisfaction Expectation Quality) ///
	wscheme("centroid") boot(50)

/* Example 4 */
/* --------- */
use ./data/ecsimobi, clear
plssem (Expectation > CUEX1-CUEX3) (Satisfaction > CUSA1-CUSA3) ///
	(Quality > PERQ1-PERQ7) (Image < IMAG1-IMAG5), ///
	structural(Quality Expectation, Satisfaction Expectation Quality Image, ///
	Image Expectation Quality) wscheme("path") ///
	seed(123) digits(3) //boot(50)

estat indirect, effects(Satisfaction Image Expectation, ///
	Satisfaction Quality Expectation, Satisfaction Image Quality) level(0.9) ///
	digits(5) seed(123) //boot(50)

estat vif, digits(4)

/* Example 5 */
/* --------- */
use ./data/ecsimobi, clear
plssem (Expectation > CUEX1) (Satisfaction > CUSA1-CUSA3) ///
	(Quality > PERQ1-PERQ7) (Image < IMAG1-IMAG5), ///
	structural(Quality Expectation, Satisfaction Expectation Quality Image, ///
	Image Expectation Quality) wscheme("path") ///
	seed(123) digits(3) boot(50)

/* Example 6 */
/* --------- */
use ./data/ecsimobi, clear
set seed 123
generate group = rbinomial(1, .4)
bysort group: plssem (Expectation > CUEX1-CUEX3) (Satisfaction > CUSA1-CUSA3) ///
	(Quality > PERQ1-PERQ7) (Image < IMAG1-IMAG5), ///
	structural(Quality Expectation, Satisfaction Expectation Quality Image, ///
	Image Expectation Quality) wscheme("path") ///
	seed(123) digits(5) //boot(50)

estat total, plot
	
/* Example 7 */
/* --------- */
use ./data/ecsimobi, clear
plssem (CUEX > CUEX1) (Expectation > CUEX2-CUEX3) (Satisfaction > CUSA1-CUSA3) ///
	(Complaints > CUSCO) (Loyalty > CUSL1-CUSL3) (Image > IMAG1-IMAG5) ///
	(Quality > PERQ1-PERQ7) (Value > PERV1-PERV2), ///
	structural(CUEX Expectation Image CUSA1, Quality Expectation, Value Expectation Quality, ///
	Satisfaction Value Quality Image Expectation, Complaints Satisfaction, ///
	Loyalty Complaints Satisfaction Image) ///
	wscheme("path") digits(5) stats //boot(50)

/* Example 8 */
/* --------- */
use ./data/ecsimobi, clear
plssem (Expectation > CUEX1-CUEX3) (Satisfaction > CUSA1-CUSA3) ///
	(Quality > PERQ1-PERQ7) (Image > IMAG1-IMAG5), ///
	structural(Quality Expectation, Satisfaction Expectation Quality Image) ///
	noheader nostructtable nodiscrimtable init(eigen)

factor CUEX1-CUEX3, factors(1) pcf
factor CUSA1-CUSA3, factors(1) pcf
factor PERQ1-PERQ7, factors(1) pcf
factor IMAG1-IMAG5, factors(1) pcf

/* Example 9 */
/* --------- */
use ./data/ecsimobi, clear
rename *, lower

sem (Expectation -> cuex1-cuex3) (Satisfaction -> cusa1-cusa3) ///
	(Quality -> perq1-perq7) (Quality Image <- Expectation) ///
	(Image -> imag1-imag5) (Satisfaction <- Expectation Quality Image), stand

estat teffects, nototal nodirect stand

/* Example 10 */
/* ---------- */
use ./data/ecsimobi, clear

generate group = 1
replace group = 0 in 1/125
* keep if group == 1
bysort group: plssem (Expectation > CUEX1-CUEX3) (Satisfaction > CUSA1-CUSA3) ///
	(Complaints > CUSCO) (Loyalty > CUSL1-CUSL3) (Image > IMAG1-IMAG5) ///
	(Quality > PERQ1-PERQ7) (Value > PERV1-PERV2), ///
	structural(Expectation Image, Quality Expectation, Value Expectation Quality, ///
	Satisfaction Value Quality Image Expectation, Complaints Satisfaction, ///
	Loyalty Complaints Satisfaction Image) ///
	wscheme("path") digits(4) stats //boot(50)

/* Example 11 */
/* ---------- */
clear all
set obs 1000
set seed 123
generate ITEM1 = rbinomial(6, .5)
generate ITEM2 = rbinomial(6, .5)
generate ITEM3 = rbinomial(6, .5)
/* kdensity ITEM1
kdensity ITEM2
kdensity ITEM3 */
plssem (lv1 > ITEM1) (lv2 > ITEM2-ITEM3), structural(lv1 lv2)

/* Example 12 */
/* ---------- */
use ./data/ecsimobi, clear
generate group = 1
replace group = 0 in 1/125
/* bysort group: */ plssem (Expectation > CUEX1-CUEX3) (Satisfaction > CUSA1-CUSA3) ///
	(Complaints > CUSCO) (Loyalty > CUSL1-CUSL3) (Image > IMAG1-IMAG5) ///
	(Quality > PERQ1-PERQ7) (Value > PERV1-PERV2), ///
	structural(Expectation Image, Quality Expectation, Value Expectation Quality, ///
	Satisfaction Value Quality Image Expectation, Complaints Satisfaction, ///
	Loyalty Complaints Satisfaction Image) ///
	wscheme("path") digits(5) init("eigen") //boot(50)

/* Example 13 */
/* ---------- */
use ./data/ecsimobi, clear

set seed 123
generate BIN1 = rbinomial(1, .4)
generate BIN2 = rbinomial(1, .6)
generate BIN3 = rbinomial(4, .3)
/* generate group = 1
replace group = 0 in 1/125
bysort group: */ plssem (Expectation > CUEX1-CUEX3) (Satisfaction > CUSA1-CUSA3) ///
	(Quality > PERQ1-PERQ7) (Image < IMAG1-IMAG5) (Virtual > BIN1) (Virtual2 > BIN2), ///
	structural(Quality Expectation Virtual Virtual2, Satisfaction Expectation Quality Image, ///
	Image Expectation Quality, Virtual Image, Virtual2 Quality) wscheme("path") ///
	seed(123) digits(5) binary(Virtual Virtual2) //boot(50)
	
estat unobshet

/* Example 14 */
/* ---------- */
/* [From Sanchez, G. (2013) PLS Path Modeling with R (Chapter 5)] */

// Round 1
set seed 123
use ./data/education, clear
plssem (Support > sup*) (Advising > adv*) (Tutoring > tut*) (Value > val*) ///
	(Satisfaction > sat*) (Loyalty > loy*), ///
	structural(Value Support Advising Tutoring, ///
	Satisfaction Support Advising Tutoring Value, Loyalty Satisfaction) ///
	wscheme("centroid") digits(4) tol(1e-6) //boot(50)

// Round 2
set seed 123
use ./data/education, clear
generate sup_appre = 8 - sup_under
generate loy_pleas = 8 - loy_asha
drop sup_under loy_asha
plssem (Support > sup*) (Advising > adv*) (Tutoring > tut*) (Value > val*) ///
	(Satisfaction > sat*) (Loyalty > loy*), ///
	structural(Value Support Advising Tutoring, ///
	Satisfaction Support Advising Tutoring Value, Loyalty Satisfaction) ///
	wscheme("centroid") digits(4) tol(1e-6) //boot(50)

// Round 3
set seed 123
use ./data/education, clear
drop sup_under loy_asha
plssem (Support > sup*) (Advising > adv*) (Tutoring > tut*) (Value > val*) ///
	(Satisfaction > sat*) (Loyalty > loy*), ///
	structural(Value Support Advising Tutoring, ///
	Satisfaction Support Advising Tutoring Value, Loyalty Satisfaction) ///
	wscheme("centroid") digits(4) tol(1e-6) // boot(50)

estat indirect, effects(Satisfaction Value Support, ///
	Satisfaction Value Advising, Satisfaction Value Tutoring, ///
	Loyalty Satisfaction Support) digits(5)

estat total, digits(7) plot

/* Example 15 */
/* ---------- */
/* [From Sanchez, G. (2013) PLS Path Modeling with R (Section 6.2.3)] */
use ./data/college, clear

bysort Gender: plssem (HighSchool > HS_GPA SAT_Verbal SAT_Math) ///
	(Intro > Biology1 Chemistry1 Math1 Physics1) ///
	(Medium > Biology2 Chemistry2 Math2 Physics2) (Graduation > FinalGPA), ///
	structural(Intro HighSchool, Medium Intro HighSchool, ///
	Graduation Medium Intro HighSchool) ///
	wscheme("centroid") digits(4) tol(1e-6) //boot(50)

/* Example 16 */
/* ---------- */
use ./data/ecsimobi, clear

*set seed 123
*generate group = rbinomial(1, .4)
generate group = 1
replace group = 2 in 126/250
plssem (Expectation > CUEX1-CUEX3) (Satisfaction > CUSA1-CUSA3) ///
	(Complaints > CUSCO) (Loyalty > CUSL1-CUSL3) (Image > IMAG1-IMAG5) ///
	(Quality > PERQ1-PERQ7) (Value > PERV1-PERV2), ///
	structural(Expectation Image, Quality Expectation, Value Expectation Quality, ///
	Satisfaction Value Quality Image Expectation, Complaints Satisfaction, ///
	Loyalty Complaints Satisfaction Image) ///
	wscheme("path") digits(7) tol(1e-5) ///
	group(group, reps(20) method(permutation) plot) //boot(50)

/* Example 17 */
/* ---------- */
/* [From Sanchez, G. (2013) PLS Path Modeling with R (Section 6.2.3)] */
use ./data/college, clear

set seed 123
plssem (HighSchool > HS_GPA SAT_Verbal SAT_Math) ///
	(Intro > Biology1 Chemistry1 Math1 Physics1) ///
	(Medium < Biology2 Chemistry2 Math2 Physics2) (Graduation > FinalGPA), ///
	structural(Intro HighSchool, Medium Intro HighSchool, ///
	Graduation Medium Intro HighSchool) ///
	wscheme("centroid") digits(4) tol(1e-6) correlate(mv lv cross, cutoff(.3)) ///
	group(Gender, reps(100) method(permutation) plot alpha(0.1) what(path)) ///
	//boot(50)

/* Example 18 */
/* ---------- */
use ./data/ecsimobi, clear

/* generate group = 1
replace group = 0 in 1/125
bysort group: */ plssem (Expectation > CUEX1-CUEX3) (Satisfaction > CUSA1-CUSA3) ///
	(Quality > PERQ1-PERQ7) (Image < IMAG1-IMAG5), ///
	structural(Quality Expectation, Satisfaction Expectation Quality Image, ///
	Image Expectation Quality) wscheme("path") seed(123) digits(5) //boot(50)

/* Example 19 */
/* ---------- */
use ./data/ecsimobi, clear
plssem (Expectation > CUEX1-CUEX3) (Satisfaction > CUSA1-CUSA3) ///
	(Complaints > CUSCO) (Loyalty > CUSL1-CUSL3) (Image > IMAG1-IMAG5) ///
	(Quality > PERQ1-PERQ7) (Value > PERV1-PERV2), ///
	structural(Expectation Image, Quality Expectation, Value Expectation Quality, ///
	Satisfaction Value Quality Image Expectation, Complaints Satisfaction, ///
	Loyalty Complaints Satisfaction Image) ///
	wscheme("path") digits(4) correlate(mv lv cross, cutoff(.3)) //boot(50)

estat unobshet
	
/* Example 20 */
/* ---------- */
use ./data/ecsimobi, clear
set seed 123
generate group = rbinomial(2, .5) + 1
label define group_lbl 1 "Low" 2 "Mid" 3 "High"
label value group group_lbl
plssem (Expectation > CUEX1-CUEX3) (Satisfaction > CUSA1-CUSA3) ///
	(Complaints > CUSCO) (Loyalty > CUSL1-CUSL3) (Image > IMAG1-IMAG5) ///
	(Quality > PERQ1-PERQ7) (Value > PERV1-PERV2), ///
	structural(Expectation Image, Quality Expectation, Value Expectation Quality, ///
	Satisfaction Value Quality Image Expectation, Complaints Satisfaction, ///
	Loyalty Complaints Satisfaction Image) ///
	wscheme("path") digits(4) ///
	group(group, reps(50) method(bootstrap) plot alpha(0.1) what(path)) ///
	// boot(50)

/* Example 21 */
/* ---------- */
/* [From Sanchez, G. (2013) PLS Path Modeling with R (Section 7.3.1)] */
/* Product indicator approach */
use ./data/satisfaction, clear
forvalues i = 1/3 {
	forvalues j = 1/3 {
		local ninter = (`i' - 1)*3 + `j'
		generate inter`ninter' = imag`i'*sat`j'
	}
}

plssem (Image > imag1-imag3) (Satisfaction > sat1-sat3) ///
	(Loyalty > loy1-loy3) (Inter > inter1-inter9), ///
	structural(Loyalty Image Satisfaction Inter) ///
	wscheme("centroid") digits(5) tol(1e-6) //boot(200)

plssem (Image > imag1-imag3) (Satisfaction > sat1-sat3) ///
	(Loyalty > loy1-loy3), structural(Loyalty Image*Satisfaction) ///
	wscheme("centroid") digits(5) tol(1e-6) //boot(200) ///
	// group(gender, reps(50) method(normal) plot alpha(0.1) what(loadings))

/* Two-stage path modeling approach */
use ./data/satisfaction, clear

// Stage 1
quietly plssem (Image > imag1-imag3) (Satisfaction > sat1-sat3) ///
	(Loyalty > loy1-loy3), structural(Loyalty Image Satisfaction) ///
	wscheme("centroid") digits(5) tol(1e-6) //boot(200)

// Stage 2
quietly generate ImageSatisfaction = Image*Satisfaction
plssem (Image > imag1-imag3) (Satisfaction > sat1-sat3) ///
	(Loyalty > loy1-loy3) (Inter > ImageSatisfaction), ///
	structural(Loyalty Image Satisfaction Inter) ///
	wscheme("centroid") digits(5) tol(1e-6) //boot(200)

/* Two-stage regression approach */
use ./data/satisfaction, clear
// Stage 1
quietly plssem (Image > imag1-imag3) (Satisfaction > sat1-sat3) ///
	(Loyalty > loy1-loy3), structural(Loyalty Image Satisfaction) ///
	wscheme("centroid") digits(5) tol(1e-6) //boot(200)

// Stage 2
quietly generate ImageSatisfaction = Image*Satisfaction
regress Loyalty Image ImageSatisfaction Satisfaction

/* Example 22 */
/* ---------- */
use ./data/ecsimobi, clear

set seed 123
generate group = rbinomial(1, .5) + 1
label define group_lbl 1 "Low" 2 "High"
label value group group_lbl

plssem (Expectation > CUEX1-CUEX3) (Satisfaction > CUSA1-CUSA3) ///
	(Complaints > CUSCO) (Loyalty > CUSL1-CUSL3) (Image > IMAG1-IMAG5) ///
	(Quality > PERQ1-PERQ7) (Value > PERV1-PERV2), ///
	structural(Expectation Image, Quality Expectation, Value Expectation Quality, ///
	Satisfaction Value Quality Image*Expectation, Complaints Satisfaction, ///
	Loyalty Complaints Satisfaction Image) ///
	wscheme("path") digits(4) tol(1e-6) ///
	//group(group, reps(50) method(normal) plot alpha(0.1) what(loadings)) ///
	//boot(50)

plssemplot, outerweights

/* Example 23 */
/* ---------- */
/* [From Sanchez, G. (2013) PLS Path Modeling with R (Section 8.1.3)] */
/* Repeated indicators approach */
use ./data/offense, clear

plssem (Special < FieldGoals OtherTDs) ///
	(Rushing > YardsRushAtt RushYards RushFirstDown) ///
	(Passing > YardsPassComp PassYards PassFirstDown) ///
	(Offense > YardsPassComp PassYards PassFirstDown YardsRushAtt RushYards RushFirstDown) ///
	(Scoring > PointsGame OffensTD TDGame), ///
	structural(Scoring Special Offense, Offense Passing Rushing) ///
	wscheme("centroid") digits(5) tol(1e-6) //boot(200)

// plssemplot, outerweights // currently this does not work 

/* Example 24 */
/* ---------- */
use ./data/ecsimobi, clear
generate group = rbinomial(1, .5) + 1
label define group_lbl 1 "Low" 2 "High"
label value group group_lbl
plssem (Expectation > CUEX1-CUEX3) (Satisfaction > CUSA1-CUSA3) ///
	(Complaints > CUSCO) (Loyalty > CUSL1-CUSL3) (Image > IMAG1-IMAG5) ///
	(Quality > PERQ1-PERQ7) (Value > PERV1-PERV2), digits(4) init(eigen) ///
	//group(group, reps(50) method(normal) plot alpha(0.1) what(loadings)) ///
	/* structural(Expectation Image, Quality Expectation, Value Expectation Quality, ///
	Satisfaction Value Quality Image Expectation, Complaints Satisfaction, ///
	Loyalty Complaints Satisfaction Image) wscheme("path") tol(1e-6) ///
	//boot(50) */

/* Example 25 */
/* ---------- */
/* [From Sanchez, G. (2013) PLS Path Modeling with R (Chapter 4)] */
use ./data/spainfoot, clear
generate NGCH = -1*GCH
generate NGCA = -1*GCA
plssem (Attack > GSH GSA SSH SSA) (Defense > NGCH NGCA CSH CSA) ///
	(Success > WMH WMA LWR LRWL), ///
	structural(Success Attack Defense) correlate(cross) ///
	wscheme("centroid") digits(5) tol(1e-6) //boot(200)

* plssemplot, scores
* plssemplot, crossloadings
* plssemplot, loadings
plssemplot, stats(Attack)

/* Example 26 */
/* ---------- */
/* [From Sanchez, G. (2013) PLS Path Modeling with R (Section 6.2.3)] */
use ./data/college, clear

plssem (HighSchool > HS_GPA SAT_Verbal SAT_Math) ///
	(Intro > Biology1 Chemistry1 Math1 Physics1) ///
	(Medium < Biology2 Chemistry2 Math2 Physics2) (Graduation > FinalGPA), ///
	structural(Intro HighSchool, Medium Intro HighSchool, ///
	Graduation Medium Intro HighSchool) ///
	wscheme("centroid") digits(4) tol(1e-6) correlate(mv lv cross, cutoff(.3)) ///
	group(Gender, plot what(loadings))

/* Example 27 */
/* ---------- */
/* [From Sanchez, G. (2013) PLS Path Modeling with R (Section 6.2.3)] */
use ./data/college, clear

plssem (HighSchool > HS_GPA SAT_Verbal SAT_Math) ///
	(Intro > Biology1 Chemistry1 Math1 Physics1) ///
	(Medium < Biology2 Chemistry2 Math2 Physics2) (Graduation > FinalGPA), ///
	structural(Intro HighSchool, Medium Intro HighSchool, ///
	Graduation Medium Intro HighSchool) ///
	digits(4) rawsum

/* Example 28 */
/* ---------- */
/* [From Sanchez, G. (2013) PLS Path Modeling with R (Section 7.3.1)] */
/* Product indicator approach */
use ./data/satisfaction, clear

plssem (Image > imag1-imag3) (Satisfaction > sat1-sat3) ///
	(Loyalty > loy1-loy3), structural(Loyalty Image*Satisfaction) ///
	wscheme("centroid") digits(5) tol(1e-6) ///
	group(gender, reps(10) method(permutation) groupseed(123))

/* Example 29 */
/* ---------- */
use ./data/ecsimobi, clear

/* two formative model does not work [SERGIO: NOW IT WORKS] */
plssem (Expectation < CUEX1-CUEX3) (Satisfaction < CUSA1-CUSA3), ///
	structural(Satisfaction Expectation)

/* Example 30 */
/* ---------- */
use ./data/ecsimobi, clear

set seed 123
generate BIN1 = rbinomial(1, .4)
generate BIN2 = rbinomial(1, .6)
generate BIN3 = rbinomial(4, .3)
/* generate group = 1
replace group = 0 in 1/125
bysort group: */ plssem (Expectation > CUEX1-CUEX3) (Satisfaction > CUSA1-CUSA3) ///
	(Quality > PERQ1-PERQ7) (Image < IMAG1-IMAG5) (Virtual > BIN1), ///
	structural(Quality Expectation Virtual, Satisfaction Expectation Quality Image, ///
	Image Expectation Quality, Virtual Image) wscheme("path") ///
	seed(123) digits(5) binary(Virtual) //boot(50)

/* Example 31 */
/* ---------- */
use ./data/ecsimobi, clear

plssem (Expectation > CUEX1-CUEX3) (Satisfaction > CUSA1-CUSA3) ///
	(Image > IMAG1-IMAG4), ///
	structural(Expectation Image, Satisfaction Expectation Image) ///
	wscheme("path") boot(25)

mat list e(loadings)
mat list e(loadings_bs)
mat list e(pathcoef)
mat list e(pathcoef_bs)

/* Example 32 */
/* ---------- */
/* This is the application included in the JSS paper */
use ./data/workout2, clear

plssem (Attractive > face sexy) ///
			 (Appearance > body appear attract) ///
			 (Muscle > muscle strength endur) ///
			 (Weight > lweight calories cweight), ///
			 structural(Appearance Attractive, ///
	                Muscle Appearance, ///
	                Weight Appearance) ///
	     boot(200) seed(123) stats correlate(lv)
		
estat indirect, effects(Muscle Appearance Attractive, ///
                        Weight Appearance Attractive) ///
								boot(200) seed(456)

plssemplot, loadings

/*
// SEM solution
sem (Attractive -> face sexy) ///
		(Appearance -> body appear attract) ///
		(Muscle -> muscle strength endur) ///
		(Weight -> lweight calories cweight) ///
		(Appearance <- Attractive) ///
		(Muscle <- Appearance) ///
		(Weight <- Appearance), stand

estat teffects, nototal nodirect stand

predict f*, latent

graph matrix Appearance Muscle Weight Attractive f1-f4, half
correlate Appearance Muscle Weight Attractive f1-f4
*/

/* Example 33 */
/* ---------- */
/* From Garson's book */
use ./data/jobsat, clear
plssem (SES > OccStat StdEduc) ///
			(Incentives > Incent1 Incent2) ///
			(Motivation > Motive1 Motive2), ///
			structural(Motivation Incentives SES) ///
	    //boot(25) seed(123) stats corr(lv)
		
plssemplot, loadings
plssemplot, crossloadings

plssem (SES > OccStat StdEduc) ///
			(Incentives > Incent1 Incent2) ///
			(Motivation > Motive1 Motive2), ///
			structural(Motivation Incentives SES) ///
			group(Gender, what("path")) ///
	    //boot(25) seed(123) stats corr(lv)

/* Example 34 */
/* ---------- */
use ./data/ecsimobi, clear

// Let us first make some of the items binary here using the percat
percat IMAG1-IMAG5

// Now, we have got five binary image indicators
plssem (Expectation > CUEX1-CUEX3) (Satisfaction > CUSA1-CUSA3) ///
     (Image > med_IMAG1-med_IMAG5), ///
	   structural(Expectation Image, Satisfaction Expectation Image) ///
	   rawsum digits(4)

// Let us see the range of the Image variable which is now sum of the five binary med_IMAG items
sum Image //as you see it has a range from 0 to 5 as mentioned above
corr CUEX* Expectation   
corr CUSA* Satisfaction   
corr med_IMAG* Image // the loadings in the measurement part should be equal to these corr but they are not
reg Expectation Image, beta // the coefs provided in the table are not standardised
reg Satisfaction Expectation Image, beta // the coefs in the table are not standardised

/* Example 35 */
/* ---------- */
import excel "./data/rebus_data.xlsx", firstrow clear

rename ServMøte1 Service
rename CSmat1 Food
rename CSrengjør1 Hygiene
rename ImagePsy1 Lively
rename CSlys1 Lighting
rename Lojal1 Positive_talk
rename Lojal2 Recommend
rename Kjønn Gender

keep Service Food Hygiene Lively Lighting Positive_talk Recommend Gender Rebus_gr

plssem (Tangible > Service Food Hygiene) (Atmospheric > Lively Lighting) ///
	(Loyalty > Positive_talk Recommend), structural(Loyalty Tangible Atmospheric) ///
	group(Rebus_gr, method(permutation) reps(100) plot seed(101))

// REBUS SOLUTION
plssem (Tangible > Service Food Hygiene) (Atmospheric > Lively Lighting) ///
	(Loyalty > Positive_talk Recommend), structural(Loyalty Tangible Atmospheric) ///
	wscheme(centroid) tol(1e-06) digits(7)

estat unobshet // a posteriori approach
estat unobshet, numclass(2) // a priori approach

/* Example 36 */
/* ---------- */
use ./data/workout2, clear

plssem (Attractive > face sexy) (Appearance > body appear attract) ///
	(Muscle > muscle strength endur) (Weight > lweight calories cweight), ///
	structural(Appearance Attractive, Muscle Appearance, Weight Appearance)

predict, xb residuals

/* Example 37 */
/* ---------- */
/* [From Sanchez, G. (2013) PLS Path Modeling with R (Chapter 9)] */
use ./data/futbol, clear

plssem (Attack > GSH GSA SSH SSA) (Defense > NGCH NGCA CSH CSA) ///
	(Success > WMH WMA), structural(Success Attack Defense) ///
	wscheme(centroid) tol(1e-06) digits(7)

* Cluster Analysis on Scores
/*
cluster wardslinkage Attack Defense Success, name(plssem_ward)
cluster dendrogram plssem_ward, labels(Team) xlabel(, angle(90) labsize(*.65))
cluster generate member = gr(4), name(plssem_ward)
table member
twoway (scatter Success Attack if member == 1, mlabel(Team) mlabposition(0) msize(0)) ///
	(scatter Success Attack if member == 2, mlabel(Team) mlabposition(0) msize(0)) ///
	(scatter Success Attack if member == 3, mlabel(Team) mlabposition(0) msize(0)) ///
	(scatter Success Attack if member == 4, mlabel(Team) mlabposition(0) msize(0)), ///
	legend(off)
*/

* REBUS solution
estat unobshet // a posteriori approach
// estat unobshet, numclass(5) stop(.05) dendro name("rebus_cluster")
// estat unobshet, numclass(6)
// estat unobshet, test reps(300) plot seed(123)

estat unobshet, numclass(2) // a priori approach

plssem (Attack > GSH GSA SSH SSA) (Defense > NGCH NGCA CSH CSA) ///
	(Success > WMH WMA), structural(Success Attack Defense) ///
	wscheme(centroid) tol(1e-06) ///
	group(rebus_class, reps(50) method(normal) plot)

plssem (Attack > GSH GSA SSH SSA) (Defense > NGCH NGCA CSH CSA) ///
	(Success > WMH WMA), structural(Success Attack Defense) ///
	wscheme(centroid) tol(1e-06) ///
	group(rebus_class, reps(50) method(normal) what(loadings) plot)

/* Example 38 */
/* ---------- */
use ./data/workout2, clear

plssem (Attractive > face sexy) (Appearance > body appear attract) ///
	(Muscle > muscle strength endur) ///
	(Weight > lweight calories cweight), ///
	structural(Appearance Attractive, Muscle Appearance, Weight Appearance) ///
	tol(1e-06) wscheme(centroid)

estat unobshet, test reps(100) plot

/* Example 39 */
/* ---------- */
/* Comparison of results of plssem without structural part with those of factor
	 analysis */
use ./data/workout2, clear

plssem (Attractive > face sexy) (Appearance > body appear attract) ///
	(Muscle > muscle strength endur) (Weight > lweight calories cweight), init(eigen)

// Using orthogonal rotation
factor face sexy body appear attract muscle strength endur lweight calories cweight, ///
	factors(4) pcf
rotate
predict f1 f2 f3 f4, regression

graph matrix Attractive-Weight f1-f4, half
correlate Attractive-Weight f1-f4

// Using oblique rotation
factor face sexy body appear attract muscle strength endur lweight calories cweight, ///
	factors(4) pcf
rotate, oblique oblimin
predict fo1 fo2 fo3 fo4, regression

graph matrix Attractive-Weight fo1-fo4, half
correlate Attractive-Weight fo1-fo4

twoway (scatter Attractive fo4) (function y = x, range(-4 4))
twoway (scatter Appearance fo1) (function y = x, range(-4 4))
twoway (scatter Muscle fo3) (function y = x, range(-4 4))
twoway (scatter Weight fo2) (function y = x, range(-4 4))

/* Example 40 */
/* --------- */
use ./data/ecsimobi, clear

plssem (Expectation > CUEX1-CUEX3) (Satisfaction > CUSA1-CUSA3) ///
	(Complaints > CUSCO) (Loyalty > CUSL1-CUSL3) (Image > IMAG1-IMAG5) ///
	(Quality > PERQ1-PERQ7) (Value > PERV1-PERV2), ///
	structural(Expectation Image, Quality Expectation, Value Expectation Quality, ///
	Satisfaction Value Quality Image Expectation, Complaints Satisfaction, ///
	Loyalty Complaints Satisfaction Image) ///
	wscheme("path") digits(4)

plssemplot, outerweights

/* Example 41 */
/* ---------- */
/* [From https://www.smartpls.com/documentation/sample-corporate-reputation] */
import excel "./data/Corporate Reputation Data.xlsx", ///
	sheet("Sheet1") firstrow clear
mvdecode _all, mv(-99)

plssem (LIKE > like_?) (COMP > comp_?) (CUSA > cusa) (CUSL > cusl_?) ///
			 (CSOR < csor_?) (ATTR < attr_?) (PERF < perf_?) (QUAL < qual_?), ///
			 structural(CUSA COMP, CUSA LIKE, CUSL COMP, CUSL LIKE, CUSL CUSA, ///
									LIKE QUAL, LIKE CSOR, LIKE PERF, LIKE ATTR, COMP QUAL, ///
									COMP CSOR, COMP PERF, COMP ATTR)

plssemplot, loadings
plssemplot, cross
plssemplot, innermodel
plssemplot, stats(LIKE)

estat total, plot

plssemplot, outerweights

/* Example 42 */
/* ---------- */
/* [From http://marketing-bulletin.massey.ac.nz/V26/MB_v26_T1_Wong_2016.pdf] */
import excel "./data/photocopier.xlsx", ///
	sheet("Sheet1") firstrow clear
label define type 1 "non-profit" 2 "for-profit"
label value type type

// Final model
plssem (QUALI > quali_?) (FINAN > finan_1 finan_3) ///
			 (GOVER > gover_1 gover_3) (LEADR > leadr_2 leadr_3) ///
			 (REPUT > finan_1 finan_3 gover_1 gover_3 leadr_2 leadr_3 quali_?) ///
			 (SATIS > satis_1) (PRICE > price_?) (LOYAL > loyal_?), ///
			 structural(QUALI REPUT, FINAN REPUT, GOVER REPUT, ///
									LEADR REPUT, SATIS REPUT PRICE, LOYAL REPUT SATIS PRICE) ///
			 wscheme(path) digits(5) tol(1e-06)

plssemplot, innermodel
plssemplot, loadings
plssemplot, scores

estat total, plot

predict, xb residuals

estat unobshet, test reps(20) plot

// Multigroup analysis
plssem (QUALI > quali_?) (FINAN > finan_1 finan_3) ///
			 (GOVER > gover_1 gover_3) (LEADR > leadr_2 leadr_3) ///
			 (REPUT > finan_1 finan_3 gover_1 gover_3 leadr_2 leadr_3 quali_?) ///
			 (SATIS > satis_1) (PRICE > price_?) (LOYAL > loyal_?), ///
			 structural(QUALI REPUT, FINAN REPUT, GOVER REPUT, ///
									LEADR REPUT, SATIS REPUT PRICE, LOYAL REPUT SATIS PRICE) ///
			 wscheme(path) digits(5) tol(1e-06) ///
			 group(type, method(normal) plot) //boot(50)

plssem (QUALI > quali_?) (FINAN > finan_1 finan_3) ///
			 (GOVER > gover_1 gover_3) (LEADR > leadr_2 leadr_3) ///
			 (REPUT > finan_1 finan_3 gover_1 gover_3 leadr_2 leadr_3 quali_?) ///
			 (SATIS > satis_1) (PRICE > price_?) (LOYAL > loyal_?), ///
			 structural(QUALI REPUT, FINAN REPUT, GOVER REPUT, ///
									LEADR REPUT, SATIS REPUT PRICE, LOYAL REPUT SATIS PRICE) ///
			 wscheme(path) digits(5) tol(1e-06) ///
			 group(rebus_class, method(normal) plot) //boot(50)
