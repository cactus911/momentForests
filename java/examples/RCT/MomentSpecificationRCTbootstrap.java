/*
 * The MIT License
 *
 * Copyright 2018 Stephen P. Ryan <stephen.p.ryan@wustl.edu>.
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
package examples.RCT;

import com.stata.sfi.*;
import Jama.Matrix;
import java.util.Random;
import JSci.maths.statistics.NormalDistribution;
import core.ContainerMoment;
import core.IntegerPartition;
import core.MomentContinuousSplitObj;
import core.MomentPartitionObj;
import core.MomentSpecification;
import core.NaiveContainer;
import java.io.BufferedReader;
import java.io.FileReader;
import java.util.TreeSet;

/**
 *
 * @author Stephen P. Ryan <stephen.p.ryan@wustl.edu>
 */
public class MomentSpecificationRCTbootstrap implements MomentSpecification {


    Jama.Matrix X;
    Jama.Matrix Xoriginal;
    Jama.Matrix Y;
    int[] searchArray;
    
    Jama.Matrix CVparameters;
    int numtrees;
   
   
    // Let it automatically detect the names of the variables from Stata dataset
    private final String[] variableNamess() {
    	
        int  numvar = Data.getParsedVarCount();
    	String[] variableNamess = new String[numvar - 1];
    
    	for (int i=1; i<numvar-1; i++)
        {
        	variableNamess[i-1] = Data.getVarName(i+1);
        }

    	return variableNamess;	
    }
    
    final private String[] variableNames = variableNamess();

        

    /**
     * Try the simplest possible thing here. Y = \tau W + \epsilon, where W is
     * randomly assigned
     *
     * So X is a vector of zeros and ones FIRST COLUMN: This is like W---the
     * treatment indicator We will skip that when evaluating splits add one
     * column to X as a set of group indicators
     */
    public MomentSpecificationRCTbootstrap() {
    }

   @Override
   public Double getPredictedY(Matrix xi, Jama.Matrix beta) {
        // constant + beta*1(W=1)
//        System.out.println("****");
//        pmUtility.prettyPrintVector(beta);
//        pmUtility.prettyPrint(xi);
        if (beta != null) {
            double yhat = xi.get(0, 0) * beta.get(0, 0);
            return yhat;
        }
        return null;
     }

  @Override
  public MomentContinuousSplitObj getFminObjective(Matrix nodeY, Matrix nodeX, int indexSplitVariable, double minProportionEachPartition, int minCountEachPartition) {
     return new MomentContinuousSplitObjRCT(indexSplitVariable, nodeX, nodeY, minProportionEachPartition, minCountEachPartition);
  }

    @Override
    public int[] getVariableIndicesToSearchOver() {
        int parsedVariables = Data.getParsedVarCount();	
        searchArray = new int[parsedVariables-2];
        for (int i=0; i<searchArray.length; i++)
        {
        	searchArray[i] = i + 1;
        }
        return searchArray;
    }

    @Override
    public MomentPartitionObj getMomentPartitionObj(Matrix nodeX, Matrix nodeY, int indexSplitVariable, IntegerPartition partition) {
        return new MomentPartitionObjRCT(partition, indexSplitVariable, nodeX, nodeY);
    }

     @Override
     public ContainerMoment computeOptimalBeta(Matrix nodeY, Matrix nodeX) {
        return new ContainerRCT(nodeX, nodeY);
     }

  //  @Override
  //  public void generateData(int numObs, Random rng, boolean addNoise) {
  //      X = new Jama.Matrix(numObs, discreteIndicator.length);
  //      Y = new Jama.Matrix(numObs, 1);

  //      double meanTreatment = 0;
  //      double meanControl = 0;
  //      int countTreatment = 0;
  //      int countControl = 0;

  //      for (int i = 0; i < numObs; i++) {
  //          double draw = rng.nextDouble();
  //          if (draw < 0.5) {
  //              X.set(i, 0, 1);
  //          }
  //          for (int j = 0; j < discreteIndicator.length - 1; j++) {
  //              if (discreteIndicator[j + 1]) {
  //                  X.set(i, j + 1, 1 + Math.floor(rng.nextDouble() * 2.0));
  //              } else {
  //                  X.set(i, j + 1, rng.nextDouble() * 2.0);
  //              }
  //          }

  //          if (addNoise) {
  //              Y.set(i, 0, Y.get(i, 0) + 0 * normal.inverse(rng.nextDouble()));
  //          }

  //          Jama.Matrix obsBeta = getBetaTruth(X.getMatrix(i, i, 0, X.getColumnDimension() - 1));

  //          if (X.get(i, 0) == 1) {
                // treatment
//                Y.set(i, 0, Y.get(i, 0) + 1);
//                if (X.get(i, 2) == 2) {
//                    Y.set(i, 0, Y.get(i, 0) + 1);
//                }
//                if (X.get(i, 3) == 1) {
//                    Y.set(i, 0, Y.get(i, 0) + 2);
//                }
//                if (X.get(i, 4) == 1) {
//                    Y.set(i, 0, Y.get(i, 0) + 3);
//                }
    //            if (X.get(i, 1) > 1.0) {
    //                Y.set(i, 0, Y.get(i, 0) + 1);
    //            }
    //            if (X.get(i, 1) < 0.50) {
    //                Y.set(i, 0, Y.get(i, 0) + 2);
    //            }
    //            if (X.get(i, 1) > 0.5 && X.get(i, 2) == 1) {
    //                Y.set(i, 0, Y.get(i, 0) - 3);
    //            }
    //            meanTreatment += Y.get(i, 0);
    //            countTreatment++;
    //        } else {
    //            meanControl += Y.get(i, 0);
    //            countControl++;
    //        }
    //    }
    //    meanTreatment /= countTreatment;
    //    meanControl /= countControl;
    //    System.out.println("meanControl: " + meanControl + " meanTreatment: " + meanTreatment + " diff: " + (meanTreatment - meanControl));
        // System.exit(0);
//        pmUtility.prettyPrint(pmUtility.concatMatrix(Y, X));
//        System.exit(0);
    // }

    @Override
    public Jama.Matrix getY() {
        return Y;
     }

    @Override
    public Jama.Matrix getX() {
        return X;
    }
    
    @Override
    public Jama.Matrix getXoriginal() {
        return Xoriginal; // In this case, Xoriginal
    }
    
    @Override
    public Jama.Matrix cvparameters() {
        return CVparameters;
    }
    
    @Override
    public int numberoftrees() {
        return numtrees;
    }

    @Override  // using the parameter row, which is the last row from Stata, to identify whether discrete or continuous
    public Boolean[] getDiscreteVector() {
    	    	
    	double parameter;
    	long obsEnd = Data.getObsParsedIn2();
        int  numvar = Data.getParsedVarCount();
    	Boolean[] discreteIndicator = new Boolean[numvar - 1];
        for(int i=0; i < discreteIndicator.length; i++){
        	parameter = Data.getNum(i+2,obsEnd - 1 - 1 - 6); // We need to skip those numtree and CV parameter rows and bootstrap parameter row
			if (parameter >= 1) {
				discreteIndicator[i] = true;
    		} else {
    			discreteIndicator[i] = false;
    		}
        }    	
        return discreteIndicator;
    }
    

    @Override
    public Matrix getBetaTruth(Matrix xi) {
        Jama.Matrix beta = new Jama.Matrix(1, 1);
        beta.set(0, 0, 0);
   //     /**
   //      * Add in some observable heterogeneity. I think this should work.
   //      */
   //     // if (xi.get(0, 1) == 1 && xi.get(0,2)==4) {
// //       if (xi.get(0, 1) == 1) {
// //           beta.set(0, 0, 20);
// //       }
// //       if (xi.get(0, 1) == 4 && xi.get(0,2)==3) {
// //           beta.set(0, 0, -10);
// //       }
        beta.set(0, 0, xi.get(0, 1) + xi.get(0, 1) * (xi.get(0, 2) - 1));
        return beta;
    }

    /* Add this part in order to Subsample with replacement: Bootstrap sampling */
    public static Jama.Matrix resample(Jama.Matrix x, long seed) {
        Random rng = new Random(seed);
        if (1 == 0) {
            return x.copy();
        }
        Jama.Matrix re = new Jama.Matrix(x.getRowDimension(), x.getColumnDimension());
        for (int i = 0; i < re.getRowDimension(); i++) {
            int index = (int) Math.floor(re.getRowDimension() * rng.nextDouble());
           
            for (int j = 0; j < re.getColumnDimension(); j++) {
                re.set(i, j, x.get(index, j));
            }
        }
        return re;
    }
    
    
    @Override   
	   public void loadData() {
  	
	   int varIndex_x, varIndex_y;
	   int rc ;
	   double value_x, value_y;
	   String msg_x, msg_y ;
	   TreeSet<Integer> exclusionTree = new TreeSet<>();
      // Get number of variables in varlist specified to javacall
	   int nVariables = Data.getParsedVarCount();
	   SFIToolkit.displayln("How many variables?  " + nVariables);
	   // Get first observation specified by an in restriction
	   long firstObs = Data.getObsParsedIn1();
	   // Get last observation specified by an in restriction
	   long lastObs = Data.getObsParsedIn2();
	   long nObs = Data.getObsTotal();
	   int nObss = (int) nObs ;
	   	   	   

	   // find out missing values
	   int counter = 0;
	   for (long obs = firstObs; obs <= lastObs; obs++ ) {
	   if (!Data.isParsedIfTrue(obs)) {
		   exclusionTree.add(counter);
	   }
	   counter++;
	   }
	   
	   // We should not count the last discrete/continuous parameter row + and the bootstrap parameter row
	   // We should not count the number of tree parameter row
	   // We should not count the last 6 CV parameter rows
           // We should not count the bootstrap parameter row
	   Y = new Jama.Matrix(counter - 1 - 1 - 1 - 6 - exclusionTree.size(), 1);
	   X = new Jama.Matrix(counter - 1 - 1 - 1 - 6 - exclusionTree.size(), nVariables-1);
	   Xoriginal = new Jama.Matrix(counter - 1 - 1 - 1 - 6 - exclusionTree.size(), nVariables-1);
	   
	   numtrees = (int) Data.getNum(1,nObs - 1 - 5 - 1);
	   
	   CVparameters = new Jama.Matrix(1,2);
	   for (int cv = 0; cv <= 1; cv++ ) {
		   CVparameters.set(0, cv, Data.getNum(1,nObs - 6 + cv)); // plug in bootstrap parameters  
	   }
	   
	   
	// Loop over variables
	   	for (long obs_y = 1; obs_y <= nObs - 1 - 1 - 1 - 6; obs_y++ ) {
                   if (!exclusionTree.contains(obs_y-1)) {
	   		   int var_y = 1;
	 		   // get the real variable index for parsed variable -var-
			   varIndex_y = Data.mapParsedVarIndex(var_y);
			   // Put value of observation obs on variable varIndex into value
			   value_y = Data.getNum(varIndex_y, obs_y);
			   
	 		   // Exit with error
			   if (Data.isValueMissing(value_y)) {
				   msg_y = "{err}missing values encountered" ;
				   SFIToolkit.errorln(msg_y);
				   // return(416);
			   }
			
	   		int obss_y = (int) obs_y ;
			Y.set(obss_y - 1 , 0 , value_y);
			// X.set(obss_y - 1 , 0 , value_y);	
                      }
		   }
	   		
	   			   
		   for(int var_x = 2; var_x <= nVariables; var_x++) {
			// Loop over observations
	   	   for (long obs_x = 1; obs_x<=nObs - 1 - 1 - 1 - 6; obs_x++ ) {
	   		if (!exclusionTree.contains(obs_x-1)) {
		   // get the real variable index for parsed variable -var-
			   varIndex_x = Data.mapParsedVarIndex(var_x);
		   // Put value of observation obs on variable varIndex into value
			   value_x = Data.getNum(varIndex_x, obs_x);
			   
		   // Exit with error
			   if (Data.isValueMissing(value_x)) {
				   msg_x = "{err}missing values encountered" ;
				   SFIToolkit.errorln(msg_x);
				   // return(4166);
			   }

			   int obss_x = (int) obs_x ;
			   X.set(obss_x - 1 , var_x - 2 , value_x);
			    
		   }
	   	   }
	   }
		  
		    SFIToolkit.displayln("Missing Observations : " + exclusionTree.size()); 
                    
		 
	        for (int i = 0; i < X.getRowDimension(); i++) {      
	            for (int j = 0; j < X.getColumnDimension(); j++) {
	                Xoriginal.set(i, j, X.get(i, j));
	            }
	        }
	        
		
                /* Add this part in order to Subsample with replacement: Bootstrap sampling */  
	        Random rng = new Random();
	        long seed = rng.nextLong();
	        Jama.Matrix XX = resample(X, seed);
            Jama.Matrix YY = resample(Y, seed);
		    
	        for (int i = 0; i < XX.getRowDimension(); i++) {      
	            for (int j = 0; j < XX.getColumnDimension(); j++) {
	                X.set(i, j, XX.get(i, j));
	            }
	        }
	        
	        for (int i = 0; i < YY.getRowDimension(); i++) {      
	            for (int j = 0; j < YY.getColumnDimension(); j++) {
	                Y.set(i, j, YY.get(i, j));
	            }
	        }
	        
          Data.addVarLong("beta_bootstrapped");
}
    

       

     @Override
     public Matrix getOutOfSampleX() {
         throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
     }

     @Override
     public NaiveContainer computeNaiveStatistics() {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
     }

    @Override
    public String getVariableName(int variableIndex) {
        return variableNames[variableIndex];
    }

    @Override
    public String getFixedEffectName(int variableIndex, int fixedEffectIndex) {
             return "var[" + variableIndex + "]=" + fixedEffectIndex;
    }

     @Override
     public String formatTreeLeafOutput(Matrix beta, Matrix variance) {
         if (beta == null) {
             return "null (shouldn't be here!)";
         }
         double b = beta.get(0, 0);
         double se = Math.sqrt(variance.get(0, 0));
         String stars = "";
        NormalDistribution normal = new NormalDistribution(0, 1);
        if (Math.abs(b / se) > Math.abs(normal.inverse(0.90))) {
            stars = "*";
        }
        if (Math.abs(b / se) > Math.abs(normal.inverse(0.95))) {
            stars = "**";
        }
        if (Math.abs(b / se) > Math.abs(normal.inverse(0.99))) {
            stars = "***";
        }
        return String.format("%.2f (%.2f) %s", b, se, stars);
     }

}
