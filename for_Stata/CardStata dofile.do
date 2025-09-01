
/*	
	Authors: Stephen P. Ryan <stephen.p.ryan@wustl.edu>, Nathan Jiang <jiang.n@wustl.edu>
*/

import delimited "table2.csv", clear

label define black_label 0 "white" 1 "black"
label values black black_label 

label define region_1966_label 1 "New England" 2 "Mid Atlantic" 3 "East North Central" 4 "West North Central" 5 "South Atlantic" 6 "East South Central" 7 "West South Central" 8 "Mountain" 9 "Pacific"
label values region_1966 region_1966_label 

momentforest lwage76 constant ed76, ///
    z(ed76 exp76 black reg76r region_1966) ///
    discretevars(black reg76r region_1966) ///
    strata(region_1966) ///
    numtrees(1) ///
    cv(true) ///
    cvgridminleaf(25 50 100 200) ///
    cvgridminimp(0.1 0.2 0.4 0.8 1.0 1.2) ///
    cvgridmaxdepth(1 2 3 4 5 6) ///
	gen(beta)

