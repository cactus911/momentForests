
/*	
	Authors: Stephen P. Ryan <stephen.p.ryan@wustl.edu>, Nathan Jiang <jiang.n@wustl.edu>
	
	To-dos:
	1. Print everything that it would print in Java
	2. Add all the optionality from Java into the Stata javacall
	3. Create an .ado file for the javacall
	
	Working on variable names, it is inconsistently working
			
	Working on stratified random sampling, the code currently doesn't work 
	The Stata interface allows for a stratum variable outside of Z 
	The DataLens part requires that the stratum variable is in Z 

*/

import delimited "C:\Users\natha\Documents\table2.csv", clear
javacall examples.CardStata.StataInterface RunCardModel, jars(momentforests.jar; sfi-api.jar; jfreechart-1.0.19.jar; utility.jar; itext-1.3.jar; jsci-core.jar) ///
	args("y=lwage76" ///
		"x=constant ed76" ///
		"z=ed76 exp76 black reg76r region_1966"  ///
		"discretevars=black reg76r region_1966" ///
		"numtrees=1" ///
		"cv=true" ///
		"cvgrid_minleaf=25, 50, 100, 200" ///
		"cvgrid_minimp=0.1, 0.2, 0.4, 0.8, 1.0, 1.2" ///
		"cvgrid_maxdepth=1, 2, 3, 4, 5, 6" ///
		"varnames=ed76 exp76 black reg76r region_1966" /// // Want to set this up in an .ado file local allvars : list zvars
		)

