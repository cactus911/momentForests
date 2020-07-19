/*
 * The MIT License
 *
 * Copyright 2020 Stephen P. Ryan <stephen.p.ryan@wustl.edu>.
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */
package examples.SimpleRCT;

import com.stata.sfi.*;
import java.util.TreeSet;
import examples.SimpleRCT.SimpleRCTMain;

/**
 *
 * @author Stephen P. Ryan <stephen.p.ryan@wustl.edu>
 */

public class StataInterface {
    public static int Momentforests(String[] args) {

        int rc = 0;
        Jama.Matrix X;
        Jama.Matrix Y;
        int numObs;
        int numtrees;
        Jama.Matrix CVparameters;
    
        
        // 1. Load Data from Stata
        int varIndex_x, varIndex_y;
	double value_x, value_y;
	String msg_x, msg_y;
	TreeSet<Integer> exclusionTree = new TreeSet<>();    
        
        
	int nVariables = Data.getParsedVarCount(); // Get number of variables in varlist specified to javacall
	long firstObs = Data.getObsParsedIn1(); // Get first observation specified by an in restriction
	long lastObs = Data.getObsParsedIn2(); // Get last observation specified by an in restriction
	long nObs = Data.getObsTotal();
	int nObss = (int) nObs;  	   

	// find out missing values
	int counter = 0;
	for (long obs = firstObs; obs <= lastObs; obs++ ) {
	if (!Data.isParsedIfTrue(obs)) {
		exclusionTree.add(counter);
	   }
	   counter++;
	   }
        
        // should not count the last discrete/continuous parameter row
	// should not count the number of tree parameter row
	// should not count the last 6 CV parameter rows
	Y = new Jama.Matrix(counter - 1 - 1 - 6 - exclusionTree.size(), 1);
	X = new Jama.Matrix(counter - 1 - 1 - 6 - exclusionTree.size(), nVariables-1);
        
        numtrees = (int) Data.getNum(1,nObs - 1 - 5);
        
    	CVparameters = new Jama.Matrix(1,6);
	for (int cv = 5; cv >= 0; cv-- ) {
		CVparameters.set(0,5-cv, Data.getNum(1,nObs - cv));  
	}    
        
        // Loop over y to assign values 
	    for (long obs_y = 1; obs_y <= nObs - 1 - 1 - 6; obs_y++ ) {
                if (!exclusionTree.contains(obs_y-1)) {
	   	    int var_y = 1;
	            // get the real variable index for parsed variable -var-
	            varIndex_y = Data.mapParsedVarIndex(var_y);
	            value_y = Data.getNum(varIndex_y, obs_y);
			   
	 	    // Exit with error
	            if (Data.isValueMissing(value_y)) {
		        msg_y = "{err}missing values encountered" ;
			SFIToolkit.errorln(msg_y);
			}
			
	   	    int obss_y = (int) obs_y ;
	            Y.set(obss_y - 1 , 0 , value_y);
                }
	    } 		
        // Loop over x to assign values 		   
	    for(int var_x = 2; var_x <= nVariables; var_x++) {
	   	for (long obs_x = 1; obs_x<=nObs - 1 - 1 - 6; obs_x++ ) {
	            if (!exclusionTree.contains(obs_x-1)) {
		        // get the real variable index for parsed variable -var-
			varIndex_x = Data.mapParsedVarIndex(var_x);
			value_x = Data.getNum(varIndex_x, obs_x);
			   
		    // Exit with error
		    if (Data.isValueMissing(value_x)) {
			msg_x = "{err}missing values encountered" ;
		        SFIToolkit.errorln(msg_x);
			}

		    int obss_x = (int) obs_x ;
	            X.set(obss_x - 1 , var_x - 2 , value_x);
		   }
	   	}
	   }
            
            
            // 2. Variable Index for Searching over
            int[] variableSearchIndex = new int[nVariables-2];
            for (int i=0; i<variableSearchIndex.length; i++)
                {
                    variableSearchIndex[i] = i + 1;
                }    
            
            
            // 3. Discrete, Continuous Index
            double parameter;
            Boolean[] DiscreteVariables = new Boolean[nVariables - 1];
            for(int i=0; i < DiscreteVariables.length; i++){
        	parameter = Data.getNum(i+2,lastObs - 1 - 6); // We need to skip those numtree and CV parameter rows
			if (parameter >= 1) {
				DiscreteVariables[i] = true;
    		} else {
    			DiscreteVariables[i] = false;
    		}
            }

            
        SimpleRCTMain go = new SimpleRCTMain(X, Y, numtrees, CVparameters, variableSearchIndex, DiscreteVariables);
            
        Jama.Matrix estimationresults = go.EstimationResults();
        Data.addVarLong("beta_estimated");
        Data.addVarLong("se_bootstrapped");
        
            for (int i = 0; i < estimationresults.getRowDimension(); i++) {      
                Data.storeNum(Data.getParsedVarCount()+2, i+1, estimationresults.get(i,0));
                Data.storeNum(Data.getParsedVarCount()+3, i+1, estimationresults.get(i,1));
            }    
            
    return(rc);
    
    }
}
