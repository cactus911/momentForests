
program define momentforests

	version 15.1

	syntax varlist [if] [in] [, num_tree(numlist)] [mink(numlist)] [mink_lower(numlist)] [mink_size(numlist)] [mink_upper(numlist)] [msebar(numlist)] [msebar_lower(numlist)] [msebar_size(numlist)] [msebar_upper(numlist)] [bootstrap(numlist) ] [cv(numlist) ]
	
	
	/*** Preliminary calculations for searching parameters ***/
	quietly keep `varlist'
	token `varlist'
	
	quietly describe `varlist'
	foreach j of numlist 2/`r(k)' {
	local jj = `j'-1
	quietly order ``j'', after(``jj'')
	}

	quietly reg `1' `2'
	
	local mse = e(rmse) * e(rmse)
	local obs = e(N)
	
	
	/*** Add the Discrete vs Continuous Parameters in the last row ***/
	quietly ds, has(type byte int)
	mata : s1 = "`r(varlist)'"
	
	mata : st_view(N1=., ., "`varlist'")
	mata : st_view(N2=., ., "`r(varlist)'")
	mata : st_numscalar("n1", cols(N1))
	mata : st_numscalar("n2", cols(N2))
	
	matrix n = J(1,n1*n2,0)
	matrix parameter = J(1,n1,0)
	
	local i = 1
	foreach var_target of varlist `r(varlist)' {
		foreach var_search of varlist `varlist' {
		if "`var_target'" == "`var_search'" {
		matrix n[1,`i'] = 1
		local i = `i' + 1
		}
		else {
		local i = `i' + 1
		}
	}
	}
	
	local N : word count `varlist'
	local NN : word count `r(varlist)'
	foreach j of numlist 1/`N' {
		foreach k of numlist 1/`NN' {
			matrix parameter[1,`j'] = parameter[1,`j'] + n[1,(`k'-1)*`N' + `j']
		}
	}
	
	
	/*** Add the number of trees in the last row ***/
	matrix numtree = J(1,n1,0)
	
	if (mi("`num_tree'")) {
		 local num_tree = 200
	}
		matrix numtree[1,1] = `num_tree'

			
	/*** Add the parameters to search hyper-parameters in the last six rows ***/
	matrix search_parameter = J(6,n1,0)
		
	
	if (mi("`mink_lower'")) {
		 local mink_lower = round(0.0010 * `obs') 
	}	
		matrix search_parameter[1,1] = round(`mink_lower')
	
	
	if (mi("`mink_upper'")) {
		 local mink_upper = round(0.0220 * `obs') 
	}
		matrix search_parameter[3,1] = round(`mink_upper')
		
		
	if (mi("`mink_size'")) {
		 local mink_size = round( (`mink_upper'-`mink_lower')/10 ) 
	}
		matrix search_parameter[2,1] = round(`mink_size')
	
	
	if (mi("`msebar_lower'")) {
		//if round((1/15000)*`mse', 0.01)  > 0.01 { 
        //             local msebar_lower = round((1/15000)*`mse', 0.0001) 
        //        }
        //        else local msebar_lower = 0.01
		local msebar_lower = round((1/15000)*`mse', 0.0000001) 
	}
		matrix search_parameter[4,1] = round(`msebar_lower', 0.0000001) 
	
	
	if (mi("`msebar_upper'")) {
		 local msebar_upper = round((1/9000)*`mse', 0.0000001) 
	}
		matrix search_parameter[6,1] = round(`msebar_upper', 0.0000001) 
		
	if (mi("`msebar_size'")) {
		 local msebar_size = (`msebar_upper' - `msebar_lower')/10 
		 // local msebar_size = round( (`msebar_upper' - `msebar_lower')/50 , 0.0000001) 
	}
		matrix search_parameter[5,1] = round(`msebar_size', 0.0000001) 
	
	
	/*** Add a row for CV or not ***/
	matrix cv_parameter = J(1,n1,0)	
	local cv = 1
	if (!mi("`mink'") | !mi("`msebar'") ) {
		local cv = 0
	}
		
	matrix cv_parameter[1,1] = `cv'
	
		
	/*** Add two rows for the hyper-parameters: mink and msebar ***/	
	matrix hyper_parameter = J(2,n1,0)	
		
	if (mi("`mink'")) {
		local mink = round((`mink_lower' + `mink_upper')/2)
	}
	
	matrix hyper_parameter[1,1] = `mink'
	
	
	if (mi("`msebar'")) {
		local msebar = round( (`msebar_lower' + `msebar_upper')/2, 0.0000001) 
	}
	
	matrix hyper_parameter[2,1] = round(`msebar', 0.0000001)

	
	
	/*** Add Number of Bootstrapping : B **/
	matrix numbootstrap = J(1,n1,0)
	
	if (mi("`bootstrap'")) {
		local bootstrap = 100
	}
	
	matrix numbootstrap[1,1] = `bootstrap'
	
	
	
	// Stack all the parameters
	mata : Parameter = st_matrix("parameter")
	mata : Numtree = st_matrix("numtree")
	mata : Search_Parameter = st_matrix("search_parameter")
	mata : CV_Parameter = st_matrix("cv_parameter")
	mata : Hyper_Parameter = st_matrix("hyper_parameter")
	mata : Numbootstrap = st_matrix("numbootstrap")
		
	mata : N12 = ( N1 \ Parameter \ Numtree \ Search_Parameter \ CV_Parameter \ Hyper_Parameter \ Numbootstrap )
	
	mata : st_addobs(1)  // add one observation line for the discrete/continuous parameters
	mata : st_addobs(1)  // add one observation line for the number of trees
	mata : st_addobs(6)  // add six observation lines for the search parameters
	mata : st_addobs(1)  // add one observation line for the cv parameter
	mata : st_addobs(2)  // add two observation lines for the hyper parameters
	mata : st_addobs(1)  // add one observation line for the bootstrap parameter
	
	mata : st_store((1::rows(N12)),(1::cols(N12))', N12)	

	

	
	/* increase the capacity of stata-java plug-in*/
	quietly set java_heapmax 10g
	
	/*** Pass the Stata dataset along with the parameter rows to Java for the main estimation computation ***/
	javacall examples.SimpleRCT.StataInterface Momentforests `varlist' `if' `in',  ///		// /* if `touse' */  , args(`result')
	jars(momentforests.jar; commons-beanutils-1.9.3.jar; commons-collections4-4.1.jar; ///
	commons-io-2.6.jar; itext-1.3.jar; Jama-1.0.3.jar; jcommon-1.0.23.jar; jfreechart-1.0.19.jar; ///
	jsci-core.jar; optimization.jar; pmUtility.jar; snappy-java-1.1.7.jar) // args(estimated_beta) // `V' `N')
	
	quietly set java_heapmax
	

		
	/*** Get rid of ALL the parameter rows, so that the dataset has the same number of rows as it was at the beginning ***/		
	/*** Discrete/Continuous parameter (1) num_tree(1) search parameters (6) cv parameter (1) hyper parameters (2) Number of Bootstrapping (1) ***/
	quietly drop in L // Discrete/Continuous parameter (1)
	quietly drop in L // num_tree(1)
	quietly drop in L
	quietly drop in L
	quietly drop in L
	quietly drop in L
	quietly drop in L
	quietly drop in L // search parameters (6)
	quietly drop in L // cv parameter (1)
	quietly drop in L 
	quietly drop in L // hyper parameters (2)
	quietly drop in L // Number of Bootstrapping (1)
	
	
end


