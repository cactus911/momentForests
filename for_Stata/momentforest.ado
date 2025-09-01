/*	
	Authors: Stephen P. Ryan <stephen.p.ryan@wustl.edu>, Nathan Jiang <jiang.n@wustl.edu>
*/

program define momentforest
    version 16.0

	syntax varlist(min=2) ///
    , Z(varlist) ///
      DISCRETEVARS(varlist) ///
      STRATA(name) ///
      NUMTREES(integer) ///
	  GEN(string) ///
      [ CV(string) ///
        CVGRIDMINLEAF(string) ///
        CVGRIDMINIMP(string) ///
        CVGRIDMAXDEPTH(string) ]

    * Extract y and x
	local y : word 1 of `varlist'
	local x ""
	local i = 2
	while `i' <= `: word count `varlist'' {
		local x "`x' `: word `i' of `varlist''"
		local ++i
	}

    * Quote all variable lists properly
    local x_list : subinstr local x " " " ", all
    local z_list : subinstr local z " " " ", all
    local discrete_list : subinstr local discretevars " " " ", all
	
    * Normalize cv string 
	local cvflag = lower("`cv'")

	* Cross-validation consistency checks 
	if "`cvflag'" == "false" | "`cvflag'" == "" {
		if "`cvgridminleaf'" != "" | "`cvgridminimp'" != "" | "`cvgridmaxdepth'" != "" {
			di as error "cvgrid* options are not allowed unless cv(true) is specified"
			exit 198
		}
	}
	else if "`cvflag'" == "true" {
		if "`cvgridminleaf'" == "" {
			di as error "cvgridminleaf() is required when cv(true) is specified"
			exit 198
		}
		if "`cvgridminimp'" == "" {
			di as error "cvgridminimp() is required when cv(true) is specified"
			exit 198
		}
		if "`cvgridmaxdepth'" == "" {
			di as error "cvgridmaxdepth() is required when cv(true) is specified"
			exit 198
		}
		local args `args' cv=true
	}
	else {
		di as error "cv() must be either 'true' or 'false'"
		exit 198
	}
	
    * Compose args string
    local args y=`y' x=`x_list' z=`z_list' discretevars=`discrete_list' strata=`strata' numtrees=`numtrees' gen=`gen'

    if "`cv'" == "cv" {
        local args `args' cv=true
    }

	if "`cvgridminleaf'" != "" {
		local cvgridminleaf : subinstr local cvgridminleaf " " ", ", all
		local args `args' cvgrid_minleaf=`cvgridminleaf'
	}
	if "`cvgridminimp'" != "" {
		local cvgridminimp : subinstr local cvgridminimp " " ", ", all
		local args `args' cvgrid_minimp=`cvgridminimp'
	}
	if "`cvgridmaxdepth'" != "" {
		local cvgridmaxdepth : subinstr local cvgridmaxdepth " " ", ", all
		local args `args' cvgrid_maxdepth=`cvgridmaxdepth'
	}
	
	local value_labels_clean : subinstr local value_labels `"""' `"\\""' , all
	local args `args' value_labels=`value_labels_clean'
	di as txt "value_labels_clean = `value_labels_clean'"
	*/
	
    * Call Java
	//di as text "args = `args'"
    javacall examples.CardStata.StataInterface RunCardModel, ///
        jars(momentforests.jar; sfi-api.jar; jfreechart-1.0.19.jar; utility.jar; itext-1.3.jar; jsci-core.jar) ///
        args("`args'")

end

