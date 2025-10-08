
/*	
	Authors: Stephen P. Ryan <stephen.p.ryan@wustl.edu>, Nathan Jiang <jiang.n@wustl.edu>
*/

*------------------------------*
* Linear DGP (Sec. 4.1)       *
* scenarios: 1=fully het, 2=one hom, 3=both hom
*------------------------------*
capture program drop dgp_linear
program define dgp_linear, rclass
    version 17
    syntax , N(integer) SCenario(integer) [SEED(integer 12345)]
    set seed `seed'
    clear
    set obs `n'

    * Instruments / covariates
    gen z1 = rnormal(0,1)
    gen z2 = runiform()
    gen z3 = cond(runiform()<0.3, 1, 2)

    gen x1 = 1
    gen x2 = rnormal(0, sqrt(2))   // variance 2 => sd = sqrt(2) per paper
    gen eps = rnormal(0,1)

    * Betas by scenario
    gen b1 = . 
    gen b2 = .
    if `scenario'==1 {
        replace b1 = cond(z1>0, -1.0, 0.33)
        //replace b2 = cond(z1>0,  1.0, 1.05)
		replace b2 = cond(z1>0,  1.0, -1.0)
    }
    else if `scenario'==2 {
        replace b1 = -1.0
        replace b2 = cond(z1>0, 1.0, 5.0)
    }
    else if `scenario'==3 {
        replace b1 = -1.0
        replace b2 =  1.0
    }
    else {
        di as err "scenario() must be 1, 2, or 3"
        exit 198
    }

    * Outcome
    gen y = x1*b1 + x2*b2 + eps

    label var y  "Outcome"
    label var x2 "Regressor X2 ~ N(0,2)"
    label var z1 "Instrument Z1 ~ N(0,1)"
    label var z2 "Instrument Z2 ~ U[0,1]"
    label var z3 "Instrument Z3 in {1,2} with P(1)=0.3"
end

