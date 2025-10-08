
/*	
	Authors: Stephen P. Ryan <stephen.p.ryan@wustl.edu>, Nathan Jiang <jiang.n@wustl.edu>
*/

*-----------------------------------------------*
* Logit DGP (Sec. 4.3)                          *
* y ~ Bernoulli( invlogit( X'β(Z) ) ),          *
* β(Z)=(-1,2) if Z1>0 else (-1,-2)              *
*-----------------------------------------------*
capture program drop dgp_logit
program define dgp_logit, rclass
    version 17
    syntax , N(integer) [SEED(integer 12345)]
    set seed `seed'
    clear
    set obs `n'

    gen z1 = rnormal(0,1)
    gen x1 = 1
    gen x2 = rnormal(0, sqrt(2))   // mirror linear: Var=2

    gen b1 = -1
    gen b2 = cond(z1>0,  2, -2)

    gen xb = x1*b1 + x2*b2
    gen p  = invlogit(xb)
    gen y  = (runiform()<p)

    label var y  "Choice: inside (1) vs outside (0)"
    label var x2 "Regressor"
    label var z1 "Heterogeneity driver"
end

