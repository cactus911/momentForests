
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

    gen z   = rnormal(0,1)           // mirrored from linear section
    gen x2  = rnormal(0, sqrt(2))    // mirrored: Var=2 (change to 1 if you prefer)
    gen eps = rnormal(0,1)

    gen beta1 = 2.5*sin(z) + 0.25*z^2
    gen beta2 = 1.0

    gen y = beta1 + beta2*x2 + eps

    label var y  "Outcome"
    label var z  "Z ~ N(0,1)"
    label var x2 "X2 ~ N(0,2)"
end

