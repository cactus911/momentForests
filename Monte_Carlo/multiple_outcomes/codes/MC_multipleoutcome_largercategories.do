clear all
set matsize 11000

log using D:\Data-Users\charlie.an\Documents\Stata\RCT_MC\set2\sample_largercategories\CV_search

/******************************************************************************/
/******************************** DGP 1 ***************************************/
/******************************************************************************/
clear all
local n = 100
mata : st_numscalar("sample", `n')
matrix M = (J(50,1,0))
matrix v = (1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1)
matrix V = diag(v)

set seed 111
drawnorm u1 u2 u3 u4 u5 u6 u7 u8 u9 u10 u11 u12 u13 u14 u15 u16 u17 u18 u19 u20 u21 u22 u23 u24 u25 u26 u27 u28 u29 u30 u31 u32 u33 u34 u35 u36 u37 u38 u39 u40 u41 u42 u43 u44 u45 u46 u47 u48 u49 u50, n(`n') cov(V) means(M)
mkmat u1 u2 u3 u4 u5 u6 u7 u8 u9 u10 u11 u12 u13 u14 u15 u16 u17 u18 u19 u20 u21 u22 u23 u24 u25 u26 u27 u28 u29 u30 u31 u32 u33 u34 u35 u36 u37 u38 u39 u40 u41 u42 u43 u44 u45 u46 u47 u48 u49 u50
matrix u = (u1 \ u2 \ u3 \ u4 \ u5 \ u6 \ u7 \ u8 \ u9 \ u10 \ u11 \ u12 \ u13 \ u14 \ u15 \ u16 \ u17 \ u18 \ u19 \ u20 \ u21 \ u22 \ u23 \ u24 \ u25 \ u26 \ u27 \ u28 \ u29 \ u30 \ u31 \ u32 \ u33 \ u34 \ u35 \ u36 \ u37 \ u38 \ u39 \ u40 \ u41 \ u42 \ u43 \ u44 \ u45 \ u46 \ u47 \ u48 \ u49 \ u50)
drop u*
svmat u

gen w = runiformint(0,1)

matrix beta = (J(sample,1,-30) \ J(sample,1,-30) \ J(sample,1,-30) \ J(sample,1,-30) \ J(sample,1,-30) \ J(sample,1,-20) \ J(sample,1,-20) \ J(sample,1,-20) \ J(sample,1,-20) \ J(sample,1,-20) /*
*/ \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) /*
*/ \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) /*
*/ \ J(sample,1,10) \ J(sample,1,10) \ J(sample,1,10) \ J(sample,1,10) \ J(sample,1,10) \ J(sample,1,20) \ J(sample,1,20) \ J(sample,1,20) \ J(sample,1,20) \ J(sample,1,20) /*
*/ \ J(sample,1,30) \ J(sample,1,30) \ J(sample,1,30) \ J(sample,1,30) \ J(sample,1,30))
svmat beta

gen y = beta1*w +u1
drop beta1 u1
order y, before(w)

matrix f1 =  (J(sample,1,1) \ J(sample*49,1,0) ) 
matrix f2 =  (J(sample,1,0) \ J(sample,1,1) \ J(sample*48,1,0) )
matrix f3 =  (J(sample*2,1,0) \ J(sample,1,1) \ J(sample*47,1,0) )
matrix f4 =  (J(sample*3,1,0) \ J(sample,1,1) \ J(sample*46,1,0) )
matrix f5 =  (J(sample*4,1,0) \ J(sample,1,1) \ J(sample*45,1,0) )
matrix f6 =  (J(sample*5,1,0) \ J(sample,1,1) \ J(sample*44,1,0) )
matrix f7 =  (J(sample*6,1,0) \ J(sample,1,1) \ J(sample*43,1,0) )
matrix f8 =  (J(sample*7,1,0) \ J(sample,1,1) \ J(sample*42,1,0) )
matrix f9 =  (J(sample*8,1,0) \ J(sample,1,1) \ J(sample*41,1,0) )
matrix f10 =  (J(sample*9,1,0) \ J(sample,1,1) \ J(sample*40,1,0) )
matrix f11 =  (J(sample*10,1,0) \ J(sample,1,1) \ J(sample*39,1,0) )
matrix f12 =  (J(sample*11,1,0) \ J(sample,1,1) \ J(sample*38,1,0) )
matrix f13 =  (J(sample*12,1,0) \ J(sample,1,1) \ J(sample*37,1,0) )
matrix f14 =  (J(sample*13,1,0) \ J(sample,1,1) \ J(sample*36,1,0) )
matrix f15 =  (J(sample*14,1,0) \ J(sample,1,1) \ J(sample*35,1,0) )
matrix f16 =  (J(sample*15,1,0) \ J(sample,1,1) \ J(sample*34,1,0) )
matrix f17 =  (J(sample*16,1,0) \ J(sample,1,1) \ J(sample*33,1,0) )
matrix f18 =  (J(sample*17,1,0) \ J(sample,1,1) \ J(sample*32,1,0) )
matrix f19 =  (J(sample*18,1,0) \ J(sample,1,1) \ J(sample*31,1,0) )
matrix f20 =  (J(sample*19,1,0) \ J(sample,1,1) \ J(sample*30,1,0) )
matrix f21 =  (J(sample*20,1,0) \ J(sample,1,1) \ J(sample*29,1,0) )
matrix f22 =  (J(sample*21,1,0) \ J(sample,1,1) \ J(sample*28,1,0) )
matrix f23 =  (J(sample*22,1,0) \ J(sample,1,1) \ J(sample*27,1,0) )
matrix f24 =  (J(sample*23,1,0) \ J(sample,1,1) \ J(sample*26,1,0) )
matrix f25 =  (J(sample*24,1,0) \ J(sample,1,1) \ J(sample*25,1,0) )
matrix f26 =  (J(sample*25,1,0) \ J(sample,1,1) \ J(sample*24,1,0) )
matrix f27 =  (J(sample*26,1,0) \ J(sample,1,1) \ J(sample*23,1,0) )
matrix f28 =  (J(sample*27,1,0) \ J(sample,1,1) \ J(sample*22,1,0) )
matrix f29 =  (J(sample*28,1,0) \ J(sample,1,1) \ J(sample*21,1,0) )
matrix f30 =  (J(sample*29,1,0) \ J(sample,1,1) \ J(sample*20,1,0) )
matrix f31 =  (J(sample*30,1,0) \ J(sample,1,1) \ J(sample*19,1,0) )
matrix f32 =  (J(sample*31,1,0) \ J(sample,1,1) \ J(sample*18,1,0) )
matrix f33 =  (J(sample*32,1,0) \ J(sample,1,1) \ J(sample*17,1,0) )
matrix f34 =  (J(sample*33,1,0) \ J(sample,1,1) \ J(sample*16,1,0) )
matrix f35 =  (J(sample*34,1,0) \ J(sample,1,1) \ J(sample*15,1,0) )
matrix f36 =  (J(sample*35,1,0) \ J(sample,1,1) \ J(sample*14,1,0) )
matrix f37 =  (J(sample*36,1,0) \ J(sample,1,1) \ J(sample*13,1,0) )
matrix f38 =  (J(sample*37,1,0) \ J(sample,1,1) \ J(sample*12,1,0) )
matrix f39 =  (J(sample*38,1,0) \ J(sample,1,1) \ J(sample*11,1,0) )
matrix f40 =  (J(sample*39,1,0) \ J(sample,1,1) \ J(sample*10,1,0) )
matrix f41 =  (J(sample*40,1,0) \ J(sample,1,1) \ J(sample*9,1,0) )
matrix f42 =  (J(sample*41,1,0) \ J(sample,1,1) \ J(sample*8,1,0) )
matrix f43 =  (J(sample*42,1,0) \ J(sample,1,1) \ J(sample*7,1,0) )
matrix f44 =  (J(sample*43,1,0) \ J(sample,1,1) \ J(sample*6,1,0) )
matrix f45 =  (J(sample*44,1,0) \ J(sample,1,1) \ J(sample*5,1,0) )
matrix f46 =  (J(sample*45,1,0) \ J(sample,1,1) \ J(sample*4,1,0) )
matrix f47 =  (J(sample*46,1,0) \ J(sample,1,1) \ J(sample*3,1,0) )
matrix f48 =  (J(sample*47,1,0) \ J(sample,1,1) \ J(sample*2,1,0) )
matrix f49 =  (J(sample*48,1,0) \ J(sample,1,1) \ J(sample*1,1,0) )
matrix f50 =  (J(sample*49,1,0) \ J(sample,1,1) )

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
svmat f11
svmat f12 
svmat f13 
svmat f14 
svmat f15 
svmat f16 
svmat f17 
svmat f18 
svmat f19 
svmat f20
svmat f21
svmat f22 
svmat f23 
svmat f24 
svmat f25 
svmat f26 
svmat f27 
svmat f28 
svmat f29 
svmat f30
svmat f31
svmat f32 
svmat f33 
svmat f34 
svmat f35 
svmat f36 
svmat f37 
svmat f38 
svmat f39 
svmat f40
svmat f41
svmat f42 
svmat f43 
svmat f44 
svmat f45 
svmat f46 
svmat f47 
svmat f48 
svmat f49 
svmat f50
rename(f11 f21 f31 f41 f51 f61 f71 f81 f91 f101 f111 f121 f131 f141 f151 f161 f171 f181 f191 f201 f211 f221 f231 f241 f251 f261 f271 f281 f291 f301 f311 f321 f331 f341 f351 f361 f371 f381 f391 f401 f411 f421 f431 f441 f451 f461 f471 f481 f491 f501) (f1 f2 f3 f4 f5 f6 f7 f8 f9 f10 f11 f12 f13 f14 f15 f16 f17 f18 f19 f20 f21 f22 f23 f24 f25 f26 f27 f28 f29 f30 f31 f32 f33 f34 f35 f36 f37 f38 f39 f40 f41 f42 f43 f44 f45 f46 f47 f48 f49 f50) 



/******************************************************************************/
/************************* Separate Regression 1 ******************************/
/******************************************************************************/

foreach i of numlist 1/50 {
eststo: reg y w if f`i' ==1
estimates store OLS`i'
predict y`i'
}

esttab OLS1 OLS2 OLS3 OLS4 OLS5 OLS6 OLS7 OLS8 OLS9 OLS10 using D:\Data-Users\charlie.an\Documents\Stata\RCT_MC\set2\sample_largercategories\MC1_largercategories.csv, replace cells(b(fmt(3)) se(fmt(4) par)) keep(w) nomtitle collabels(none)


/******************************************************************************/
/*************************** Moment Forests 1 *********************************/
/******************************************************************************/

recast byte w f1 f2 f3 f4 f5 f6 f7 f8 f9 f10 f11 f12 f13 f14 f15 f16 f17 f18 f19 f20 f21 f22 f23 f24 f25 f26 f27 f28 f29 f30 
momentforests y w f1 f2 f3 f4 f5 f6 f7 f8 f9 f10 f11 f12 f13 f14 f15 f16 f17 f18 f19 f20 f21 f22 f23 f24 f25 f26 f27 f28 f29 f30, num_tree(20)
save "D:\Data-Users\charlie.an\Documents\Stata\RCT_MC\set2\sample_largercategories\MC1_largercategories.dta", replace




/******************************************************************************/
/******************************** DGP 2 ***************************************/
/******************************************************************************/
clear all
local n = 200
mata : st_numscalar("sample", `n')
matrix M = (J(30,1,0))
matrix v = (1,1,1,1,1,1,1,1,1,1, 1,1,1,1,1,1,1,1,1,1, 1,1,1,1,1,1,1,1,1,1)
matrix V = diag(v)

set seed 111
drawnorm u1 u2 u3 u4 u5 u6 u7 u8 u9 u10 u11 u12 u13 u14 u15 u16 u17 u18 u19 u20 u21 u22 u23 u24 u25 u26 u27 u28 u29 u30 , n(`n') cov(V) means(M)
mkmat u1 u2 u3 u4 u5 u6 u7 u8 u9 u10 u11 u12 u13 u14 u15 u16 u17 u18 u19 u20 u21 u22 u23 u24 u25 u26 u27 u28 u29 u30 
matrix u = (u1 \ u2 \ u3 \ u4 \ u5 \ u6 \ u7 \ u8 \ u9 \ u10 \ u11 \ u12 \ u13 \ u14 \ u15 \ u16 \ u17 \ u18 \ u19 \ u20 \ u21 \ u22 \ u23 \ u24 \ u25 \ u26 \ u27 \ u28 \ u29 \ u30 )
drop u*
svmat u

gen w = runiformint(0,1)

matrix beta = (J(sample,1,-30) \ J(sample,1,-30) \ J(sample,1,-30) \ J(sample,1,-30) \ J(sample,1,-30) /*
*/ \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) /*
*/ \ J(sample,1,10) \ J(sample,1,10) \ J(sample,1,10) \ J(sample,1,10) \ J(sample,1,10)  /*
*/ \ J(sample,1,30) \ J(sample,1,30) \ J(sample,1,30) \ J(sample,1,30) \ J(sample,1,30))
svmat beta

gen y = beta1*w +u1
drop beta1 u1
order y, before(w)

matrix f1 =  (J(sample,1,1) \ J(sample*29,1,0) ) 
matrix f2 =  (J(sample,1,0) \ J(sample,1,1) \ J(sample*28,1,0) )
matrix f3 =  (J(sample*2,1,0) \ J(sample,1,1) \ J(sample*27,1,0) )
matrix f4 =  (J(sample*3,1,0) \ J(sample,1,1) \ J(sample*26,1,0) )
matrix f5 =  (J(sample*4,1,0) \ J(sample,1,1) \ J(sample*25,1,0) )
matrix f6 =  (J(sample*5,1,0) \ J(sample,1,1) \ J(sample*24,1,0) )
matrix f7 =  (J(sample*6,1,0) \ J(sample,1,1) \ J(sample*23,1,0) )
matrix f8 =  (J(sample*7,1,0) \ J(sample,1,1) \ J(sample*22,1,0) )
matrix f9 =  (J(sample*8,1,0) \ J(sample,1,1) \ J(sample*21,1,0) )
matrix f10 =  (J(sample*9,1,0) \ J(sample,1,1) \ J(sample*20,1,0) )
matrix f11 =  (J(sample*10,1,0) \ J(sample,1,1) \ J(sample*19,1,0) )
matrix f12 =  (J(sample*11,1,0) \ J(sample,1,1) \ J(sample*18,1,0) )
matrix f13 =  (J(sample*12,1,0) \ J(sample,1,1) \ J(sample*17,1,0) )
matrix f14 =  (J(sample*13,1,0) \ J(sample,1,1) \ J(sample*16,1,0) )
matrix f15 =  (J(sample*14,1,0) \ J(sample,1,1) \ J(sample*15,1,0) )
matrix f16 =  (J(sample*15,1,0) \ J(sample,1,1) \ J(sample*14,1,0) )
matrix f17 =  (J(sample*16,1,0) \ J(sample,1,1) \ J(sample*13,1,0) )
matrix f18 =  (J(sample*17,1,0) \ J(sample,1,1) \ J(sample*12,1,0) )
matrix f19 =  (J(sample*18,1,0) \ J(sample,1,1) \ J(sample*11,1,0) )
matrix f20 =  (J(sample*19,1,0) \ J(sample,1,1) \ J(sample*10,1,0) )
matrix f21 =  (J(sample*20,1,0) \ J(sample,1,1) \ J(sample*9,1,0) )
matrix f22 =  (J(sample*21,1,0) \ J(sample,1,1) \ J(sample*8,1,0) )
matrix f23 =  (J(sample*22,1,0) \ J(sample,1,1) \ J(sample*7,1,0) )
matrix f24 =  (J(sample*23,1,0) \ J(sample,1,1) \ J(sample*6,1,0) )
matrix f25 =  (J(sample*24,1,0) \ J(sample,1,1) \ J(sample*5,1,0) )
matrix f26 =  (J(sample*25,1,0) \ J(sample,1,1) \ J(sample*4,1,0) )
matrix f27 =  (J(sample*26,1,0) \ J(sample,1,1) \ J(sample*3,1,0) )
matrix f28 =  (J(sample*27,1,0) \ J(sample,1,1) \ J(sample*2,1,0) )
matrix f29 =  (J(sample*28,1,0) \ J(sample,1,1) \ J(sample*1,1,0) )
matrix f30 =  (J(sample*29,1,0) \ J(sample,1,1) )


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
svmat f11
svmat f12 
svmat f13 
svmat f14 
svmat f15 
svmat f16 
svmat f17 
svmat f18 
svmat f19 
svmat f20
svmat f21
svmat f22 
svmat f23 
svmat f24 
svmat f25 
svmat f26 
svmat f27 
svmat f28 
svmat f29 
svmat f30

rename(f11 f21 f31 f41 f51 f61 f71 f81 f91 f101 f111 f121 f131 f141 f151 f161 f171 f181 f191 f201 f211 f221 f231 f241 f251 f261 f271 f281 f291 f301 ) (f1 f2 f3 f4 f5 f6 f7 f8 f9 f10 f11 f12 f13 f14 f15 f16 f17 f18 f19 f20 f21 f22 f23 f24 f25 f26 f27 f28 f29 f30 ) 



/******************************************************************************/
/************************** Separate Regression 2 *****************************/
/******************************************************************************/

foreach i of numlist 1/30 {
eststo: reg y w if f`i' ==1
estimates store OLS`i'
predict y`i'
}

esttab OLS1 OLS2 OLS3 OLS4 OLS5 OLS6 OLS7 OLS8 OLS9 OLS10 OLS11 OLS12 OLS13 OLS14 OLS15 OLS16 OLS17 OLS18 OLS19 OLS20 OLS21 OLS22 OLS23 OLS24 OLS25 OLS26 OLS27 OLS28 OLS29 OLS30 using D:\Data-Users\charlie.an\Documents\Stata\RCT_MC\set2\sample_largercategories\MC2_largercategories.csv, replace cells(b(fmt(3) ) se(fmt(4) par)) keep(w) nomtitle collabels(none)


/******************************************************************************/
/**************************** Moment Forests 2 ********************************/
/******************************************************************************/

recast byte w f1 f2 f3 f4 f5 f6 f7 f8 f9 f10 f11 f12 f13 f14 f15 f16 f17 f18 f19 f20 f21 f22 f23 f24 f25 f26 f27 f28 f29 f30 
momentforests y w f1 f2 f3 f4 f5 f6 f7 f8 f9 f10 f11 f12 f13 f14 f15 f16 f17 f18 f19 f20 f21 f22 f23 f24 f25 f26 f27 f28 f29 f30, num_tree(20)
save "D:\Data-Users\charlie.an\Documents\Stata\RCT_MC\set2\sample_largercategories\MC2_largercategories.dta", replace





/******************************************************************************/
/******************************** DGP 3 ***************************************/
/******************************************************************************/
clear all
local n = 200
mata : st_numscalar("sample", `n')
matrix M = (J(30,1,0))
matrix v = (10,10,10,10,10,10,10,10,10,10, 10,10,10,10,10,10,10,10,10,10, 10,10,10,10,10,10,10,10,10,10)
matrix V = diag(v)

set seed 111
drawnorm u1 u2 u3 u4 u5 u6 u7 u8 u9 u10 u11 u12 u13 u14 u15 u16 u17 u18 u19 u20 u21 u22 u23 u24 u25 u26 u27 u28 u29 u30 , n(`n') cov(V) means(M)
mkmat u1 u2 u3 u4 u5 u6 u7 u8 u9 u10 u11 u12 u13 u14 u15 u16 u17 u18 u19 u20 u21 u22 u23 u24 u25 u26 u27 u28 u29 u30 
matrix u = (u1 \ u2 \ u3 \ u4 \ u5 \ u6 \ u7 \ u8 \ u9 \ u10 \ u11 \ u12 \ u13 \ u14 \ u15 \ u16 \ u17 \ u18 \ u19 \ u20 \ u21 \ u22 \ u23 \ u24 \ u25 \ u26 \ u27 \ u28 \ u29 \ u30 )
drop u*
svmat u

gen w = runiformint(0,1)

matrix beta = (J(sample,1,-30) \ J(sample,1,-30) \ J(sample,1,-30) \ J(sample,1,-30) \ J(sample,1,-30) /*
*/ \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) /*
*/ \ J(sample,1,10) \ J(sample,1,10) \ J(sample,1,10) \ J(sample,1,10) \ J(sample,1,10)  /*
*/ \ J(sample,1,30) \ J(sample,1,30) \ J(sample,1,30) \ J(sample,1,30) \ J(sample,1,30))
svmat beta

gen y = beta1*w +u1
drop beta1 u1
order y, before(w)

matrix f1 =  (J(sample,1,1) \ J(sample*29,1,0) ) 
matrix f2 =  (J(sample,1,0) \ J(sample,1,1) \ J(sample*28,1,0) )
matrix f3 =  (J(sample*2,1,0) \ J(sample,1,1) \ J(sample*27,1,0) )
matrix f4 =  (J(sample*3,1,0) \ J(sample,1,1) \ J(sample*26,1,0) )
matrix f5 =  (J(sample*4,1,0) \ J(sample,1,1) \ J(sample*25,1,0) )
matrix f6 =  (J(sample*5,1,0) \ J(sample,1,1) \ J(sample*24,1,0) )
matrix f7 =  (J(sample*6,1,0) \ J(sample,1,1) \ J(sample*23,1,0) )
matrix f8 =  (J(sample*7,1,0) \ J(sample,1,1) \ J(sample*22,1,0) )
matrix f9 =  (J(sample*8,1,0) \ J(sample,1,1) \ J(sample*21,1,0) )
matrix f10 =  (J(sample*9,1,0) \ J(sample,1,1) \ J(sample*20,1,0) )
matrix f11 =  (J(sample*10,1,0) \ J(sample,1,1) \ J(sample*19,1,0) )
matrix f12 =  (J(sample*11,1,0) \ J(sample,1,1) \ J(sample*18,1,0) )
matrix f13 =  (J(sample*12,1,0) \ J(sample,1,1) \ J(sample*17,1,0) )
matrix f14 =  (J(sample*13,1,0) \ J(sample,1,1) \ J(sample*16,1,0) )
matrix f15 =  (J(sample*14,1,0) \ J(sample,1,1) \ J(sample*15,1,0) )
matrix f16 =  (J(sample*15,1,0) \ J(sample,1,1) \ J(sample*14,1,0) )
matrix f17 =  (J(sample*16,1,0) \ J(sample,1,1) \ J(sample*13,1,0) )
matrix f18 =  (J(sample*17,1,0) \ J(sample,1,1) \ J(sample*12,1,0) )
matrix f19 =  (J(sample*18,1,0) \ J(sample,1,1) \ J(sample*11,1,0) )
matrix f20 =  (J(sample*19,1,0) \ J(sample,1,1) \ J(sample*10,1,0) )
matrix f21 =  (J(sample*20,1,0) \ J(sample,1,1) \ J(sample*9,1,0) )
matrix f22 =  (J(sample*21,1,0) \ J(sample,1,1) \ J(sample*8,1,0) )
matrix f23 =  (J(sample*22,1,0) \ J(sample,1,1) \ J(sample*7,1,0) )
matrix f24 =  (J(sample*23,1,0) \ J(sample,1,1) \ J(sample*6,1,0) )
matrix f25 =  (J(sample*24,1,0) \ J(sample,1,1) \ J(sample*5,1,0) )
matrix f26 =  (J(sample*25,1,0) \ J(sample,1,1) \ J(sample*4,1,0) )
matrix f27 =  (J(sample*26,1,0) \ J(sample,1,1) \ J(sample*3,1,0) )
matrix f28 =  (J(sample*27,1,0) \ J(sample,1,1) \ J(sample*2,1,0) )
matrix f29 =  (J(sample*28,1,0) \ J(sample,1,1) \ J(sample*1,1,0) )
matrix f30 =  (J(sample*29,1,0) \ J(sample,1,1) )


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
svmat f11
svmat f12 
svmat f13 
svmat f14 
svmat f15 
svmat f16 
svmat f17 
svmat f18 
svmat f19 
svmat f20
svmat f21
svmat f22 
svmat f23 
svmat f24 
svmat f25 
svmat f26 
svmat f27 
svmat f28 
svmat f29 
svmat f30

rename(f11 f21 f31 f41 f51 f61 f71 f81 f91 f101 f111 f121 f131 f141 f151 f161 f171 f181 f191 f201 f211 f221 f231 f241 f251 f261 f271 f281 f291 f301 ) (f1 f2 f3 f4 f5 f6 f7 f8 f9 f10 f11 f12 f13 f14 f15 f16 f17 f18 f19 f20 f21 f22 f23 f24 f25 f26 f27 f28 f29 f30 ) 



/******************************************************************************/
/************************** Separate Regression 3 *****************************/
/******************************************************************************/

foreach i of numlist 1/30 {
eststo: reg y w if f`i' ==1
estimates store OLS`i'
predict y`i'
}

esttab OLS1 OLS2 OLS3 OLS4 OLS5 OLS6 OLS7 OLS8 OLS9 OLS10 OLS11 OLS12 OLS13 OLS14 OLS15 OLS16 OLS17 OLS18 OLS19 OLS20 OLS21 OLS22 OLS23 OLS24 OLS25 OLS26 OLS27 OLS28 OLS29 OLS30 using D:\Data-Users\charlie.an\Documents\Stata\RCT_MC\set2\sample_largercategories\MC3_largercategories.csv, replace cells(b(fmt(3) ) se(fmt(4) par)) keep(w) nomtitle collabels(none)


/******************************************************************************/
/**************************** Moment Forests 3 ********************************/
/******************************************************************************/

recast byte w f1 f2 f3 f4 f5 f6 f7 f8 f9 f10 f11 f12 f13 f14 f15 f16 f17 f18 f19 f20 f21 f22 f23 f24 f25 f26 f27 f28 f29 f30 
momentforests y w f1 f2 f3 f4 f5 f6 f7 f8 f9 f10 f11 f12 f13 f14 f15 f16 f17 f18 f19 f20 f21 f22 f23 f24 f25 f26 f27 f28 f29 f30, num_tree(20)
save "D:\Data-Users\charlie.an\Documents\Stata\RCT_MC\set2\sample_largercategories\MC3_largercategories.dta", replace





/******************************************************************************/
/******************************** DGP 4 ***************************************/
/******************************************************************************/
clear all
local n = 200
mata : st_numscalar("sample", `n')
matrix M = (J(30,1,0))
matrix v = (10,10,10,10,10,20,20,20,20,20, 30,30,30,30,30,40,40,40,40,40, 50,50,50,50,50,60,60,60,60,60)
matrix V = diag(v)

set seed 111
drawnorm u1 u2 u3 u4 u5 u6 u7 u8 u9 u10 u11 u12 u13 u14 u15 u16 u17 u18 u19 u20 u21 u22 u23 u24 u25 u26 u27 u28 u29 u30 , n(`n') cov(V) means(M)
mkmat u1 u2 u3 u4 u5 u6 u7 u8 u9 u10 u11 u12 u13 u14 u15 u16 u17 u18 u19 u20 u21 u22 u23 u24 u25 u26 u27 u28 u29 u30 
matrix u = (u1 \ u2 \ u3 \ u4 \ u5 \ u6 \ u7 \ u8 \ u9 \ u10 \ u11 \ u12 \ u13 \ u14 \ u15 \ u16 \ u17 \ u18 \ u19 \ u20 \ u21 \ u22 \ u23 \ u24 \ u25 \ u26 \ u27 \ u28 \ u29 \ u30 )
drop u*
svmat u

gen w = runiformint(0,1)

matrix beta = (J(sample,1,-30) \ J(sample,1,-30) \ J(sample,1,-30) \ J(sample,1,-30) \ J(sample,1,-30) /*
*/ \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) \ J(sample,1,0) /*
*/ \ J(sample,1,10) \ J(sample,1,10) \ J(sample,1,10) \ J(sample,1,10) \ J(sample,1,10)  /*
*/ \ J(sample,1,30) \ J(sample,1,30) \ J(sample,1,30) \ J(sample,1,30) \ J(sample,1,30))
svmat beta

gen y = beta1*w +u1
drop beta1 u1
order y, before(w)

matrix f1 =  (J(sample,1,1) \ J(sample*29,1,0) ) 
matrix f2 =  (J(sample,1,0) \ J(sample,1,1) \ J(sample*28,1,0) )
matrix f3 =  (J(sample*2,1,0) \ J(sample,1,1) \ J(sample*27,1,0) )
matrix f4 =  (J(sample*3,1,0) \ J(sample,1,1) \ J(sample*26,1,0) )
matrix f5 =  (J(sample*4,1,0) \ J(sample,1,1) \ J(sample*25,1,0) )
matrix f6 =  (J(sample*5,1,0) \ J(sample,1,1) \ J(sample*24,1,0) )
matrix f7 =  (J(sample*6,1,0) \ J(sample,1,1) \ J(sample*23,1,0) )
matrix f8 =  (J(sample*7,1,0) \ J(sample,1,1) \ J(sample*22,1,0) )
matrix f9 =  (J(sample*8,1,0) \ J(sample,1,1) \ J(sample*21,1,0) )
matrix f10 =  (J(sample*9,1,0) \ J(sample,1,1) \ J(sample*20,1,0) )
matrix f11 =  (J(sample*10,1,0) \ J(sample,1,1) \ J(sample*19,1,0) )
matrix f12 =  (J(sample*11,1,0) \ J(sample,1,1) \ J(sample*18,1,0) )
matrix f13 =  (J(sample*12,1,0) \ J(sample,1,1) \ J(sample*17,1,0) )
matrix f14 =  (J(sample*13,1,0) \ J(sample,1,1) \ J(sample*16,1,0) )
matrix f15 =  (J(sample*14,1,0) \ J(sample,1,1) \ J(sample*15,1,0) )
matrix f16 =  (J(sample*15,1,0) \ J(sample,1,1) \ J(sample*14,1,0) )
matrix f17 =  (J(sample*16,1,0) \ J(sample,1,1) \ J(sample*13,1,0) )
matrix f18 =  (J(sample*17,1,0) \ J(sample,1,1) \ J(sample*12,1,0) )
matrix f19 =  (J(sample*18,1,0) \ J(sample,1,1) \ J(sample*11,1,0) )
matrix f20 =  (J(sample*19,1,0) \ J(sample,1,1) \ J(sample*10,1,0) )
matrix f21 =  (J(sample*20,1,0) \ J(sample,1,1) \ J(sample*9,1,0) )
matrix f22 =  (J(sample*21,1,0) \ J(sample,1,1) \ J(sample*8,1,0) )
matrix f23 =  (J(sample*22,1,0) \ J(sample,1,1) \ J(sample*7,1,0) )
matrix f24 =  (J(sample*23,1,0) \ J(sample,1,1) \ J(sample*6,1,0) )
matrix f25 =  (J(sample*24,1,0) \ J(sample,1,1) \ J(sample*5,1,0) )
matrix f26 =  (J(sample*25,1,0) \ J(sample,1,1) \ J(sample*4,1,0) )
matrix f27 =  (J(sample*26,1,0) \ J(sample,1,1) \ J(sample*3,1,0) )
matrix f28 =  (J(sample*27,1,0) \ J(sample,1,1) \ J(sample*2,1,0) )
matrix f29 =  (J(sample*28,1,0) \ J(sample,1,1) \ J(sample*1,1,0) )
matrix f30 =  (J(sample*29,1,0) \ J(sample,1,1) )


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
svmat f11
svmat f12 
svmat f13 
svmat f14 
svmat f15 
svmat f16 
svmat f17 
svmat f18 
svmat f19 
svmat f20
svmat f21
svmat f22 
svmat f23 
svmat f24 
svmat f25 
svmat f26 
svmat f27 
svmat f28 
svmat f29 
svmat f30

rename(f11 f21 f31 f41 f51 f61 f71 f81 f91 f101 f111 f121 f131 f141 f151 f161 f171 f181 f191 f201 f211 f221 f231 f241 f251 f261 f271 f281 f291 f301 ) (f1 f2 f3 f4 f5 f6 f7 f8 f9 f10 f11 f12 f13 f14 f15 f16 f17 f18 f19 f20 f21 f22 f23 f24 f25 f26 f27 f28 f29 f30 ) 

/******************************************************************************/
/************************** Separate Regression 4 *****************************/
/******************************************************************************/

foreach i of numlist 1/30 {
eststo: reg y w if f`i' ==1
estimates store OLS`i'
predict y`i'
}

esttab OLS1 OLS2 OLS3 OLS4 OLS5 OLS6 OLS7 OLS8 OLS9 OLS10 OLS11 OLS12 OLS13 OLS14 OLS15 OLS16 OLS17 OLS18 OLS19 OLS20 OLS21 OLS22 OLS23 OLS24 OLS25 OLS26 OLS27 OLS28 OLS29 OLS30 using D:\Data-Users\charlie.an\Documents\Stata\RCT_MC\set2\sample_largercategories\MC4_largercategories.csv, replace cells(b(fmt(3) ) se(fmt(4) par)) keep(w) nomtitle collabels(none)


/******************************************************************************/
/**************************** Moment Forests 4 ********************************/
/******************************************************************************/

recast byte w f1 f2 f3 f4 f5 f6 f7 f8 f9 f10 f11 f12 f13 f14 f15 f16 f17 f18 f19 f20 f21 f22 f23 f24 f25 f26 f27 f28 f29 f30 
momentforests y w f1 f2 f3 f4 f5 f6 f7 f8 f9 f10 f11 f12 f13 f14 f15 f16 f17 f18 f19 f20 f21 f22 f23 f24 f25 f26 f27 f28 f29 f30, num_tree(20)
save "D:\Data-Users\charlie.an\Documents\Stata\RCT_MC\set2\sample_largercategories\MC4_largercategories.dta",replace


log off
