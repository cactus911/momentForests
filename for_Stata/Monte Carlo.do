
/*	
	Authors: Stephen P. Ryan <stephen.p.ryan@wustl.edu>, Nathan Jiang <jiang.n@wustl.edu>
*/

log using "logfile_montecarlo.log", replace

local totaln = 100
local nlist 500 1000 2000 4000

*------------------------------*
* Linear DGP (Sec. 4.1)       *
* scenarios: 1=fully het, 2=one hom, 3=both hom
*------------------------------*

forvalues scenario = 1/3 {
	di as error "====================================="
	di as error "  Running simulations for scenario `scenario' "
	di as error "====================================="

	foreach n of local nlist {
		
		local n1 = `n'/20
		local n2 = `n'/10
		local n3 = `n'/5
		
		di as error "===================================="
		di as error "  Running simulations for n = `n'   "
		di as error "===================================="
		
		forvalues iteration = 1/`totaln' {		
			if mod(`iteration' - 1, 10)==0 {
				di "Starting iteration `iteration' for n = `n'"
			}
			
			* Generate data
			qui dgp_linear, n(`n') scenario(`scenario') seed(`iteration')
			
			* Moment forest (unrestricted)
			qui momentforest y x1 x2, ///
				z(z1 z2 z3) ///
				numtrees(50) /// 			
				seed(`iteration') ///
				discretevars(z3) ///
				propstructure(0.35) ///
				testhomogeneity(false) ///
				cv(true) ///
				cvgridminleaf(`n1' `n2' `n3') ///
				cvgridminimp(10 20 40) ///
				cvgridmaxdepth(3 4 5) ///
				gen("beta_mf_unrestricted_")				
			qui gen y_fitted_mf_unrestricted = beta_mf_unrestricted_0*x1 + beta_mf_unrestricted_1*x2
			
			* Moment forest (restricted)
			qui momentforest y x1 x2, ///
				z(z1 z2 z3) ///
				numtrees(50) /// 				
				seed(`iteration') ///
				discretevars(z3) ///
				propstructure(0.35) ///
				testhomogeneity(true) ///
				cv(true) ///
				cvgridminleaf(`n1' `n2' `n3') ///
				cvgridminimp(10 20 40) ///
				cvgridmaxdepth(3 4 5) ///
				gen("beta_mf_restricted_")				
			qui gen y_fitted_mf_restricted = beta_mf_restricted_0*x1 + beta_mf_restricted_1*x2
			
			qui gen iteration = `iteration'
			qui gen obs = _n
			qui order obs
			qui tempfile iteration_`iteration'
			qui save `iteration_`iteration''
		}

		* Merge results together
		clear all
		forvalues iteration = 1/`totaln' {
			append using `iteration_`iteration''
		} 

		* Check homogeneity
		bysort iteration: egen min_beta0 = min(beta_mf_restricted_0)
		bysort iteration: egen max_beta0 = max(beta_mf_restricted_0)
		bysort iteration: egen min_beta1 = min(beta_mf_restricted_1)
		bysort iteration: egen max_beta1 = max(beta_mf_restricted_1)
		bysort iteration: gen flag_beta0_homogeneous = min_beta0 == max_beta0
		bysort iteration: gen flag_beta1_homogeneous = min_beta1 == max_beta1

		tab flag_beta0_homogeneous
		tab flag_beta1_homogeneous
		
		* Compute MSE
		qui gen mse_y_mf_unrestricted = (y - y_fitted_mf_unrestricted)^2
		qui gen mse_b_mf_unrestricted = ((b1 - beta_mf_unrestricted_0)^2 + (b2-beta_mf_unrestricted_1)^2)/2
		qui gen mse_y_mf_restricted = (y - y_fitted_mf_restricted)^2
		qui gen mse_b_mf_restricted = ((b1 - beta_mf_restricted_0)^2 + (b2-beta_mf_restricted_1)^2)/2

		* Collapse and compute summary stats
		qui collapse (mean) mse_y_mf_unrestricted mse_b_mf_unrestricted beta_mf_unrestricted_0 beta_mf_unrestricted_1 mse_y_mf_restricted mse_b_mf_restricted beta_mf_restricted_0 beta_mf_restricted_1 flag_beta0_homogeneous flag_beta1_homogeneous, by(iteration)
		
		
		di "Results for moment forest (unrestricted)"
		summ mse_y_mf_unrestricted
		summ mse_b_mf_unrestricted
		
		di "Results for moment forest (restricted)"
		summ mse_y_mf_restricted
		summ mse_b_mf_restricted
		summ beta_mf_restricted_0 if flag_beta0_homogeneous == 1
		summ beta_mf_restricted_1 if flag_beta1_homogeneous == 1
	}
}

*---------------------------------------*
* Partially linear DGP (Sec. 4.2)       *
* Y = 2.5*sin(Z) + 0.25*Z^2 + 1.0*X2 + e
*---------------------------------------*

di as error "==============================================="
di as error "  Running simulations for partially linear DGP "
di as error "==============================================="

foreach n of local nlist {
	
	di as error "===================================="
	di as error "  Running simulations for n = `n'   "
	di as error "===================================="
	
	local n1 = `n'/50
	local n2 = `n'/20
	local n3 = `n'/10
	
	* Run moment forests
	forvalues iteration = 1/`totaln' {
		
		if mod(`iteration' - 1, 10)==0 {
			di "Starting iteration `iteration' for n = `n'"
		}
		
		* Generate data
		qui dgp_partiallylinear, n(`n') seed(`iteration')
		
		* Moment forest (unrestricted)
		qui momentforest y x1 x2, ///
			z(z) ///
			numtrees(50) /// 
			seed(`iteration') ///
			testhomogeneity(false) ///
			propstructure(0.35) ///
			cv(true) ///
			cvgridminleaf(`n1' `n2' `n3') ///
			cvgridminimp(10 20 40) ///
			cvgridmaxdepth(3 4 5) ///
			gen("beta_mf_unrestricted_")				
		qui gen y_fitted_mf_unrestricted = beta_mf_unrestricted_0*x1 + beta_mf_unrestricted_1*x2
		
		* Moment forest (restricted)
		qui momentforest y x1 x2, ///
			z(z) ///
			numtrees(50) /// 		
			seed(`iteration') ///
			testhomogeneity(true) ///
			propstructure(0.35) ///
			cv(true) ///
			cvgridminleaf(`n1' `n2' `n3') ///
			cvgridminimp(10 20 40) ///
			cvgridmaxdepth(3 4 5) ///
			gen("beta_mf_restricted_")				
		qui gen y_fitted_mf_restricted = beta_mf_restricted_0*x1 + beta_mf_restricted_1*x2
		
		qui gen iteration = `iteration'
		qui gen obs = _n
		qui order obs
		qui tempfile iteration_`iteration'
		qui save `iteration_`iteration''
	}

	* Merge results together
	clear all
	forvalues iteration = 1/`totaln' {
		append using `iteration_`iteration''
	} 
	
	* Save for n = 4000
	if `n' == 4000 {
		save "Partial linear model estimates n = 4000.dta", replace
	}
	
	* Check homogeneity
	bysort iteration: egen min_beta0 = min(beta_mf_restricted_0)
	bysort iteration: egen max_beta0 = max(beta_mf_restricted_0)
	bysort iteration: egen min_beta1 = min(beta_mf_restricted_1)
	bysort iteration: egen max_beta1 = max(beta_mf_restricted_1)
	bysort iteration: gen flag_beta0_homogeneous = min_beta0 == max_beta0
	bysort iteration: gen flag_beta1_homogeneous = min_beta1 == max_beta1

	tab flag_beta0_homogeneous
	tab flag_beta1_homogeneous
	
	drop if flag_beta1_homogeneous == 0
	
	* Compute MSE
	qui gen mse_y_mf_unrestricted = (y - y_fitted_mf_unrestricted)^2
	qui gen mse_b_mf_unrestricted = ((b1 - beta_mf_unrestricted_0)^2 + (b2-beta_mf_unrestricted_1)^2)/2
	qui gen mse_y_mf_restricted = (y - y_fitted_mf_restricted)^2
	qui gen mse_b_mf_restricted = ((b1 - beta_mf_restricted_0)^2 + (b2-beta_mf_restricted_1)^2)/2

	* Collapse and compute summary stats
	qui collapse (mean) mse_y_mf_unrestricted mse_b_mf_unrestricted beta_mf_unrestricted_0 beta_mf_unrestricted_1 mse_y_mf_restricted mse_b_mf_restricted beta_mf_restricted_0 beta_mf_restricted_1 flag_beta0_homogeneous flag_beta1_homogeneous, by(iteration)
	
	di "Results for moment forest (unrestricted)"
	summ mse_y_mf_unrestricted
	summ mse_b_mf_unrestricted
	
	di "Results for moment forest (restricted)"
	summ mse_y_mf_restricted
	summ mse_b_mf_restricted
	summ beta_mf_restricted_0 if flag_beta0_homogeneous == 1
	summ beta_mf_restricted_1 if flag_beta1_homogeneous == 1
}

log close

*************

use "Partial linear model estimates n = 4000.dta", clear
gen z_bin = round(z, 0.01)
collapse (mean) b1 beta_mf_mean = beta_mf_restricted_0 (sd) beta_mf_sd = beta_mf_restricted_0, by(z_bin)	
gen beta_mean_plus  = beta_mf_mean + beta_mf_sd
gen beta_mean_minus = beta_mf_mean - beta_mf_sd
twoway ///
    (scatter b1 z_bin, mcolor(blue) msymbol(o) legend(label(1 "b1"))) ///
    (rcap beta_mean_plus beta_mean_minus z_bin, lcolor(red)) ///
    (scatter beta_mean_plus z_bin, mcolor(red) msymbol(D) legend(label(2 "Mean ± SD of estimated beta")))
	
** End of file **


