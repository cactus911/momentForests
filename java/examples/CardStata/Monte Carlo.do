
/*	
	Authors: Stephen P. Ryan <stephen.p.ryan@wustl.edu>, Nathan Jiang <jiang.n@wustl.edu>
*/

*------------------------------*
* Linear DGP (Sec. 4.1)       *
* scenarios: 1=fully het, 2=one hom, 3=both hom
*------------------------------*

local totaln = 5

* Run moment forests
forvalues iteration = 1/`totaln' {
	dgp_linear, n(500) scenario(1) seed(`iteration')
	
	momentforest y x1 x2, ///
		z(z1 z2 z3) ///
		numtrees(100) /// 
		propstructure(0.35) ///
		seed(777) ///
		discretevars(z3) ///
		cv(true) ///
		cvgridminleaf(10) ///
		cvgridminimp(5) ///
		cvgridmaxdepth(6) ///
		gen("beta")
		
		gen iteration = `iteration'
		gen obs = _n
		order obs
		tempfile iteration_`iteration'
		save `iteration_`iteration''
}

* Merge results together
clear all
forvalues iteration = 1/`totaln' {
	append using `iteration_`iteration''
} 

* Check homogeneity
bysort iteration: egen min_beta0 = min(beta0)
bysort iteration: egen max_beta0 = max(beta0)
bysort iteration: egen min_beta1 = min(beta1)
bysort iteration: egen max_beta1 = max(beta1)
bysort iteration: gen flag_beta0_homogeneous = min_beta0 == max_beta0
bysort iteration: gen flag_beta1_homogeneous = min_beta1 == max_beta1

tab flag_beta0_homogeneous
tab flag_beta1_homogeneous

* Check estimated values for homogeneous parameters
duplicates drop
summ beta0 if flag_beta0_homogeneous == 1
summ beta1 if flag_beta1_homogeneous == 1





	
*---------------------------------------*
* Partially linear DGP (Sec. 4.2)       *
* Y = 2.5*sin(Z) + 0.25*Z^2 + 1.0*X2 + e
*---------------------------------------*

local n = 2

dgp_partiallylinear, n(4000) seed(20251003)
gen x1 = 1
rename (beta1 beta2) (beta1_true beta2_true)

momentforest y x1 x2, ///
    z(z) ///
	numtrees(10) ///
	propstructure(0.35) ///
	seed(777) ///
    cv(true) ///
    cvgridminleaf(50) ///
    cvgridminimp(10) ///
    cvgridmaxdepth(5) ///
	gen("beta`n'_")


/*
*-----------------------------------------------*
* Logit DGP (Sec. 4.3)                          *
* y ~ Bernoulli( invlogit( X'β(Z) ) ),          *
* β(Z)=(-1,2) if Z1>0 else (-1,-2)              *
*-----------------------------------------------*
dgp_logit, n(2000) seed(20251003)

momentforest y x1 x2, ///
    z(z1) ///
	numtrees(10) ///
	propstructure(0.35) ///
	seed(777) ///
    cv(true) ///
    cvgridminleaf(50) ///
    cvgridminimp(10) ///
    cvgridmaxdepth(5) ///
	gen(beta)
*/


** End of file **


