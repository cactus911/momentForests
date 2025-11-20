
/*	
	Authors: Stephen P. Ryan <stephen.p.ryan@wustl.edu>, Nathan Jiang <jiang.n@wustl.edu>
*/

import delimited "table2.csv", clear

label define black_label 0 "white" 1 "black"
label values black black_label 

label define region_1966_label 1 "New England" 2 "Mid Atlantic" 3 "East North Central" 4 "West North Central" 5 "South Atlantic" 6 "East South Central" 7 "West South Central" 8 "Mountain" 9 "Pacific"
label values region_1966 region_1966_label 

*** Part 1
log using "logfile_CardStata.log", replace

momentforest lwage76 constant, ///
    z(ed76 exp76 black reg76r smsa76r region_1966 smsa66r daded momed nodaded nomomed momdad14 sinmom14) ///
	numtrees(40) ///
	seed(20251116) ///
	propstructure(0.15) ///
    discretevars(region_1966) ///
	strata(region_1966 black) ///
	testhomogeneity(true) ///
    cv(true) ///
    cvgridminleaf(50 200 400) ///
    cvgridminimp(0.1 1 10) ///
    cvgridmaxdepth(1 2 3 4 5 6 7)

momentforest lwage76 constant ed76, ///
    z(ed76 exp76 black reg76r smsa76r region_1966 smsa66r daded momed nodaded nomomed momdad14 sinmom14) ///
	numtrees(40) ///
	seed(20251116) ///
	propstructure(0.15) ///
    discretevars(region_1966) ///
	strata(region_1966 black) ///
	testhomogeneity(true) ///
    cv(true) ///
    cvgridminleaf(50 200 400) ///
    cvgridminimp(0.1 1 10) ///
    cvgridmaxdepth(1 2 3 4 5 6 7)
	
momentforest lwage76 constant ed76 exp76, ///
    z(ed76 exp76 black reg76r smsa76r region_1966 smsa66r daded momed nodaded nomomed momdad14 sinmom14) ///
	numtrees(40) ///
	seed(20251116) ///
	propstructure(0.15) ///
    discretevars(region_1966) ///
	strata(region_1966 black) ///
	testhomogeneity(true) ///
    cv(true) ///
    cvgridminleaf(50 200 400) ///
    cvgridminimp(0.1 1 10) ///
    cvgridmaxdepth(1 2 3 4 5 6 7) ///
	gen(beta)
	
save "CardStata constant ed76 exp76.dta", replace

* One-hot encoding of region_1966, omitting category 8
gen region_1966_1 = region_1966 == 1
gen region_1966_2 = region_1966 == 2
gen region_1966_3 = region_1966 == 3
gen region_1966_4 = region_1966 == 4
gen region_1966_5 = region_1966 == 5
gen region_1966_6 = region_1966 == 6
gen region_1966_7 = region_1966 == 7
gen region_1966_9 = region_1966 == 9
momentforest lwage76 constant ed76 exp76 region_1966_1 region_1966_2 region_1966_3 region_1966_4 region_1966_5 region_1966_6 region_1966_7 region_1966_9, ///
	z(ed76 exp76 black reg76r smsa76r region_1966 smsa66r daded momed nodaded nomomed momdad14 sinmom14) ///
	numtrees(40) ///
	seed(20251116) ///
	propstructure(0.30) ///
    discretevars(region_1966) ///
	strata(region_1966 black) ///
	testhomogeneity(true) ///
    cv(true) ///
    cvgridminleaf(50 200 400) ///
    cvgridminimp(0.1 1 10) ///
    cvgridmaxdepth(1 2 3 4 5 6 7)

log close	
	
****************************************
*** Plots 							 ***
****************************************

****** Partial dependence plot, education
clear all
use "CardStata constant ed76 exp76.dta", clear
rename (beta0 beta1 beta2) (beta_constant beta_ed76 beta_exp76)

tempfile pdp
save `pdp', emptyok

* Define education values to evaluate
levelsof ed76, local(edvals)

* Prepare to store results
postfile handle float ed76 double mean_black double mean_white using `pdp', replace

foreach e of local edvals {
    use "CardStata constant ed76 exp76.dta", clear
	rename (beta0 beta1 beta2) (beta_constant beta_ed76 beta_exp76)
    replace ed76 = `e'
    gen pred = beta_constant + beta_ed76*ed76 + beta_exp76*exp76

    quietly summarize pred if black == 1
    local mean_black = r(mean)

    quietly summarize pred if black == 0
    local mean_white = r(mean)

    post handle (`e') (`mean_black') (`mean_white')
}
postclose handle

use `pdp', clear
label var mean_black "Black"
label var mean_white "White"

twoway (line mean_white ed76 if inrange(ed76, 8, 18), lcolor(blue)) ///
       (line mean_black ed76 if inrange(ed76, 8, 18), lcolor(red)), ///
       legend(position(6) ring(0) col(2) order(1 "White" 2 "Black")) ///
	   title("Education and wages") ///
	   ytitle("Average predicted log wage") xtitle("Education (years)")
	
****** Partial dependence plot, experience
clear all
use "CardStata constant ed76 exp76.dta", clear
rename (beta0 beta1 beta2) (beta_constant beta_ed76 beta_exp76)

tempfile pdp
save `pdp', emptyok

* Define experience values to evaluate
levelsof exp76, local(edvals)

* Prepare to store results
postfile handle float exp76 double mean_black double mean_white using `pdp', replace

foreach e of local edvals {
    use "CardStata constant ed76 exp76.dta", clear
	rename (beta0 beta1 beta2) (beta_constant beta_ed76 beta_exp76)

    replace exp76 = `e'
    gen pred = beta_constant + beta_ed76*ed76 + beta_exp76*exp76

    quietly summarize pred if black == 1
    local mean_black = r(mean)

    quietly summarize pred if black == 0
    local mean_white = r(mean)

    post handle (`e') (`mean_black') (`mean_white')
}
postclose handle

use `pdp', clear
label var mean_black "Black"
label var mean_white "White"

twoway (line mean_white exp76 if inrange(exp76, 0, 23), lcolor(blue)) ///
       (line mean_black exp76 if inrange(exp76, 0, 23), lcolor(red)), ///
       legend(position(6) ring(0) col(2) order(1 "White" 2 "Black")) ///
	   title("Experience and wages") ///
	   ytitle("Average predicted log wage") xtitle("Experience (years)")

***** Kernel density plots
clear all
use "CardStata constant ed76 exp76.dta", clear
rename (beta0 beta1 beta2) (beta_constant beta_ed76 beta_exp76)

kdensity beta_constant, title("") xtitle("")
kdensity beta_ed76, title("") xtitle("")
kdensity beta_exp76, title("") xtitle("")
	   
** End of file **
