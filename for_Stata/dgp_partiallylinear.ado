
/*	
	Authors: Stephen P. Ryan <stephen.p.ryan@wustl.edu>, Nathan Jiang <jiang.n@wustl.edu>
*/

*---------------------------------------*
* Partially linear DGP (Sec. 4.2)       *
* Y = 2.5*sin(Z) + 0.25*Z^2 + 1.0*X2 + e
*---------------------------------------*
capture program drop dgp_partiallylinear
program define dgp_partiallylinear, rclass
    version 17
    syntax , N(integer) [SEED(integer 12345)]
    set seed `seed'
    clear
    set obs `n'

    gen z   = runiform(-6, 4)        
    gen x2  = rnormal(0, sqrt(2))    
    gen eps = rnormal(0,1)
	qui summ eps
	replace eps = (eps - r(mean)) / r(sd)
	
    gen beta1 = 2.5*sin(z) + 0.25*z^2
    gen beta2 = 1.0

    gen y = beta1 + beta2*x2 + eps
	
	gen x1 = 1
	order x1, before(x2)
	rename (beta1 beta2) (b1 b2)
	
    label var y  "Outcome"
    label var z  "Z ~ N(0,1)"
    label var x2 "X2 ~ N(0,2)"
end

