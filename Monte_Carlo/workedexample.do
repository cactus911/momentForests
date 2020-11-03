cd "" // set your working directory here.


/*** 1/3. Discrete Case ***/

clear all
import excel ".\1_discrete.xlsx", sheet("Sheet1") firstrow
momentforests y w x1 x2, mink(5) cv(0) msebar(0.0001) num_tree(200) bootstrap(100)

gen beta_true = x1 +10*(x2 -1)
br
scatter beta_estimated beta_true

save ".\1_discrete_result.dta",replace




/*** 2/3. Continuous Case ***/

clear all
import excel ".\2_continuous.xlsx", sheet("Sheet1") firstrow
momentforests y w x, mink(5) cv(0) msebar(0.0001) num_tree(200) bootstrap(100)

gen beta_true = sin(x)
br
twoway (scatter beta_estimated x) (scatter beta_true x)

save ".\2_continuous_result.dta",replace




/*** 3/3. Continuous Case ***/

clear all
import excel ".\3_hybrid.xlsx", sheet("Sheet1") firstrow
momentforests y w x1 x2, mink(5) cv(0) msebar(0.0001) num_tree(200) bootstrap(100)

gen beta_true = (sin(x1) * (x2 - 5))
br
scatter beta_estimated beta_true

save ".\3_hybrid_result.dta",replace


