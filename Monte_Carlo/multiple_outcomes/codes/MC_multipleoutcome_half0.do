clear all
set matsize 10000


log using \\files.wustl.edu\charlie.an\RCT_MC\set2\sample_half0\CV_search

/******************************************************************************/
/******************************** DGP 1 ***************************************/
/******************************************************************************/
clear all
local n = 500
mata : st_numscalar("sample", `n')
matrix M = (J(10,1,0))
matrix v = (1,1,1,1,1,1,1,1,1,1)
matrix V = diag(v)

set seed 111
drawnorm u1 u2 u3 u4 u5 u6 u7 u8 u9 u10, n(`n') cov(V) means(M)
mkmat u1 u2 u3 u4 u5 u6 u7 u8 u9 u10
matrix u = (u1 \ u2 \ u3 \ u4 \ u5 \ u6 \ u7 \ u8 \ u9 \ u10 )
drop u*
svmat u

gen w = runiformint(0,1)
matrix beta = (J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,-30) \ J(sample,1,-10) \ J(sample,1,10) \ J(sample,1,30) \ J(sample,1,50))
svmat beta

gen y = beta1*w +u1
drop beta1 u1
order y, before(w)

matrix f1 =  (J(sample,1,1) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0))
matrix f2 =  (J(sample,1,0) \ J(sample,1,1) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0))
matrix f3 =  (J(sample,1,0) \ J(sample,1,0) \ J(sample,1,1) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0))
matrix f4 =  (J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,1) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0))
matrix f5 =  (J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,1) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0))
matrix f6 =  (J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,1) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0))
matrix f7 =  (J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,1) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0))
matrix f8 =  (J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,1) \ J(sample,1,0) \ J(sample,1,0))
matrix f9 =  (J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,1) \ J(sample,1,0))
matrix f10 = (J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,1))
svmat f1
svmat f2 
svmat f3 
svmat f4 
svmat f5 
svmat f6 
svmat f7 
svmat f8 
svmat f9 
svmat f10
rename(f11 f21 f31 f41 f51 f61 f71 f81 f91 f101) (f1 f2 f3 f4 f5 f6 f7 f8 f9 f10) 



/******************************************************************************/
/************************* Separate Regression 1 ******************************/
/******************************************************************************/

foreach i of numlist 1/10 {
eststo: reg y w if f`i' ==1
estimates store OLS`i'
predict y`i'
}

esttab OLS1 OLS2 OLS3 OLS4 OLS5 OLS6 OLS7 OLS8 OLS9 OLS10 using "\\files.wustl.edu\charlie.an\RCT_MC\set2\sample_half0\MC1_half0.csv", replace cells(b(fmt(3)) se(fmt(4) par)) keep(w) nomtitle collabels(none)


/******************************************************************************/
/*************************** Moment Forests 1 *********************************/
/******************************************************************************/

recast byte w f1 f2 f3 f4 f5 f6 f7 f8 f9 f10

momentforests y w f1 f2 f3 f4 f5 f6 f7 f8 f9 f10 
save "\\files.wustl.edu\charlie.an\RCT_MC\set2\sample_half0\MC1_half0.dta"




/******************************************************************************/
/******************************** DGP 2 ***************************************/
/******************************************************************************/
clear all
local n = 500
mata : st_numscalar("sample", `n')
matrix M = (J(10,1,0))
matrix v = (1,1,1,1,1,1,1,1,1,1)
matrix V = diag(v)

set seed 222
drawnorm u1 u2 u3 u4 u5 u6 u7 u8 u9 u10, n(`n') cov(V) means(M)
mkmat u1 u2 u3 u4 u5 u6 u7 u8 u9 u10
matrix u = (u1 \ u2 \ u3 \ u4 \ u5 \ u6 \ u7 \ u8 \ u9 \ u10 )
drop u*
svmat u

gen w = runiformint(0,1)

matrix beta = ( J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,-40) \ J(sample,1,-40) \ J(sample,1,40) \ J(sample,1,40) \ J(sample,1,40))
svmat beta

gen y = beta1*w +u1
drop beta1 u1
order y, before(w)

matrix f1 =  (J(sample,1,1) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0))
matrix f2 =  (J(sample,1,0) \ J(sample,1,1) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0))
matrix f3 =  (J(sample,1,0) \ J(sample,1,0) \ J(sample,1,1) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0))
matrix f4 =  (J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,1) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0))
matrix f5 =  (J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,1) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0))
matrix f6 =  (J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,1) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0))
matrix f7 =  (J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,1) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0))
matrix f8 =  (J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,1) \ J(sample,1,0) \ J(sample,1,0))
matrix f9 =  (J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,1) \ J(sample,1,0))
matrix f10 = (J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,1))
svmat f1
svmat f2 
svmat f3 
svmat f4 
svmat f5 
svmat f6 
svmat f7 
svmat f8 
svmat f9 
svmat f10
rename(f11 f21 f31 f41 f51 f61 f71 f81 f91 f101) (f1 f2 f3 f4 f5 f6 f7 f8 f9 f10) 


/******************************************************************************/
/************************** Separate Regression 2 *****************************/
/******************************************************************************/

foreach i of numlist 1/10 {
eststo: reg y w if f`i' ==1
estimates store OLS`i'
predict y`i'
}

esttab OLS1 OLS2 OLS3 OLS4 OLS5 OLS6 OLS7 OLS8 OLS9 OLS10 using \\files.wustl.edu\charlie.an\RCT_MC\set2\sample_half0\MC2_half0.csv, replace cells(b(fmt(3) ) se(fmt(4) par)) keep(w) nomtitle collabels(none)


/******************************************************************************/
/**************************** Moment Forests 2 ********************************/
/******************************************************************************/

recast byte w f1 f2 f3 f4 f5 f6 f7 f8 f9 f10
momentforests y w f1 f2 f3 f4 f5 f6 f7 f8 f9 f10
save "\\files.wustl.edu\charlie.an\RCT_MC\set2\sample_half0\MC2_half0.dta"





/******************************************************************************/
/******************************** DGP 3 ***************************************/
/******************************************************************************/
clear all
local n = 500
mata : st_numscalar("sample", `n')
matrix M = (J(10,1,0))
matrix v = (1,2,3,4,5,6,7,8,9,10)
matrix V = diag(v)

set seed 333
drawnorm u1 u2 u3 u4 u5 u6 u7 u8 u9 u10, n(`n') cov(V) means(M)
mkmat u1 u2 u3 u4 u5 u6 u7 u8 u9 u10
matrix u = (u1 \ u2 \ u3 \ u4 \ u5 \ u6 \ u7 \ u8 \ u9 \ u10 )
drop u*
svmat u

gen w = runiformint(0,1)

matrix beta = (J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,-30) \ J(sample,1,-10) \ J(sample,1,10) \ J(sample,1,30) \ J(sample,1,50))
svmat beta

gen y = beta1*w +u1
drop beta1 u1
order y, before(w)

matrix f1 =  (J(sample,1,1) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0))
matrix f2 =  (J(sample,1,0) \ J(sample,1,1) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0))
matrix f3 =  (J(sample,1,0) \ J(sample,1,0) \ J(sample,1,1) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0))
matrix f4 =  (J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,1) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0))
matrix f5 =  (J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,1) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0))
matrix f6 =  (J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,1) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0))
matrix f7 =  (J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,1) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0))
matrix f8 =  (J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,1) \ J(sample,1,0) \ J(sample,1,0))
matrix f9 =  (J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,1) \ J(sample,1,0))
matrix f10 = (J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,1))
svmat f1
svmat f2 
svmat f3 
svmat f4 
svmat f5 
svmat f6 
svmat f7 
svmat f8 
svmat f9 
svmat f10
rename(f11 f21 f31 f41 f51 f61 f71 f81 f91 f101) (f1 f2 f3 f4 f5 f6 f7 f8 f9 f10) 



/******************************************************************************/
/************************** Separate Regression 3 *****************************/
/******************************************************************************/

foreach i of numlist 1/10 {
eststo: reg y w if f`i' ==1
estimates store OLS`i'
predict y`i'
}

esttab OLS1 OLS2 OLS3 OLS4 OLS5 OLS6 OLS7 OLS8 OLS9 OLS10 using \\files.wustl.edu\charlie.an\RCT_MC\set2\sample_half0\MC3_half0.csv, replace cells(b(fmt(3) ) se(fmt(4) par)) keep(w) nomtitle collabels(none)


/******************************************************************************/
/**************************** Moment Forests 3 ********************************/
/******************************************************************************/

recast byte w f1 f2 f3 f4 f5 f6 f7 f8 f9 f10
momentforests y w f1 f2 f3 f4 f5 f6 f7 f8 f9 f10
save "\\files.wustl.edu\charlie.an\RCT_MC\set2\sample_half0\MC3_half0.dta"





/******************************************************************************/
/******************************** DGP 4 ***************************************/
/******************************************************************************/
clear all
local n = 500
mata : st_numscalar("sample", `n')
matrix M = (J(10,1,0))
matrix v = (1,2,3,4,5,6,7,8,9,10)
matrix V = diag(v)

set seed 444
drawnorm u1 u2 u3 u4 u5 u6 u7 u8 u9 u10, n(`n') cov(V) means(M)
mkmat u1 u2 u3 u4 u5 u6 u7 u8 u9 u10
matrix u = (u1 \ u2 \ u3 \ u4 \ u5 \ u6 \ u7 \ u8 \ u9 \ u10 )
drop u*
svmat u

gen w = runiformint(0,1)

matrix beta = ( J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,-40) \ J(sample,1,-40) \ J(sample,1,40) \ J(sample,1,40) \ J(sample,1,40))
svmat beta

gen y = beta1*w +u1
drop beta1 u1
order y, before(w)

matrix f1 =  (J(sample,1,1) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0))
matrix f2 =  (J(sample,1,0) \ J(sample,1,1) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0))
matrix f3 =  (J(sample,1,0) \ J(sample,1,0) \ J(sample,1,1) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0))
matrix f4 =  (J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,1) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0))
matrix f5 =  (J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,1) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0))
matrix f6 =  (J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,1) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0))
matrix f7 =  (J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,1) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0))
matrix f8 =  (J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,1) \ J(sample,1,0) \ J(sample,1,0))
matrix f9 =  (J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,1) \ J(sample,1,0))
matrix f10 = (J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,1))
svmat f1
svmat f2 
svmat f3 
svmat f4 
svmat f5 
svmat f6 
svmat f7 
svmat f8 
svmat f9 
svmat f10
rename(f11 f21 f31 f41 f51 f61 f71 f81 f91 f101) (f1 f2 f3 f4 f5 f6 f7 f8 f9 f10)  



/******************************************************************************/
/************************** Separate Regression 4 *****************************/
/******************************************************************************/

foreach i of numlist 1/10 {
eststo: reg y w if f`i' ==1
estimates store OLS`i'
predict y`i'
}

esttab OLS1 OLS2 OLS3 OLS4 OLS5 OLS6 OLS7 OLS8 OLS9 OLS10 using \\files.wustl.edu\charlie.an\RCT_MC\set2\sample_half0\MC4_half0.csv, replace cells(b(fmt(3) ) se(fmt(4) par)) keep(w) nomtitle collabels(none)


/******************************************************************************/
/**************************** Moment Forests 4 ********************************/
/******************************************************************************/

recast byte w f1 f2 f3 f4 f5 f6 f7 f8 f9 f10
momentforests y w f1 f2 f3 f4 f5 f6 f7 f8 f9 f10, num_tree(30) bootstrap(100)
save "\\files.wustl.edu\charlie.an\RCT_MC\set2\sample_half0\MC4_half0.dta"


log off
