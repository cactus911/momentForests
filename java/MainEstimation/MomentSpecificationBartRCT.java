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
package MainEstimation;


import com.stata.sfi.*;
import Jama.Matrix;
import java.util.Random;
import JSci.maths.statistics.NormalDistribution;
import MainEstimation.ContainerMoment;
import MainEstimation.IntegerPartition;
import MainEstimation.MomentContinuousSplitObj;
import MainEstimation.MomentPartitionObj;
import MainEstimation.MomentSpecification;
import MainEstimation.NaiveContainer;
import java.io.BufferedReader;
import java.io.FileReader;
import java.util.TreeSet;

/**
 *
 * @author Stephen P. Ryan <stephen.p.ryan@wustl.edu>
 */
public class MomentSpecificationBartRCT implements MomentSpecification {

    // boolean monteCarlo = false;
    
    Jama.Matrix X;
    Jama.Matrix Y;
    int[] searchArray;
    
    Jama.Matrix CVparameters;
    int numtrees;
    
    // static NormalDistribution normal = new NormalDistribution();
    
    // public int[] searchArray() {
    //    int parsedVariables = Data.getParsedVarCount();	
    //    searchArray = new int[parsedVariables];
    //    for (int i=0; i<searchArray.length-1; i++)
    //    {
    //    	searchArray[i] = i + 1;
    //    }
    //    return searchArray;
    //}
    
    // int[] searchArray = {1, 2, 3, 4, 5, 6, 7, 8, 9};
    // disceteIndicator includes the first column, which is the treatment effect?
    // okay, the way this works is that the indices above select variables from this list next
    // so this list +is+ defined across everything in the data set, but searchArray above controls
    // which of those X's are actually used when splitting the sample
    
    // Boolean[] discreteIndicator = {true, false, true, true, true, true, true, true, true, false};
    
    
    
   // final private String[] variableNames = {"treatment", "age", "symptomatic", "homosexual", "IVDrugUser", "PriorAZT", "Male", "White", "FutureDropout", "Baseline"};
   
    // Let it automatically detect the names of the variables from Stata dataset
    private final String[] variableNamess() {
    	
        int  numvar = Data.getParsedVarCount();
    	String[] variableNamess = new String[numvar - 1];
    
    	for (int i=1; i<numvar-1; i++)
        {
        	variableNamess[i-1] = Data.getVarName(i+1);
        }

    	// String cc = new Data.getVarName(1);
    	// variableNames = {"treatment", "age", "symptomatic", "homosexual", "IVDrugUser", "PriorAZT", "Male", "White", "FutureDropout", "Baseline"};
    	return variableNamess;	
    }
    
    final private String[] variableNames = variableNamess();

    
	// SFIToolkit.displayln("Variable Names : " + variableNamess[0]); 
    		
    		
    

    /**
     * Try the simplest possible thing here. Y = \tau W + \epsilon, where W is
     * randomly assigned
     *
     * So X is a vector of zeros and ones FIRST COLUMN: This is like W---the
     * treatment indicator We will skip that when evaluating splits add one
     * column to X as a set of group indicators
     */
    public MomentSpecificationBartRCT() {
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
     return new MomentContinuousSplitObjBart(indexSplitVariable, nodeX, nodeY, minProportionEachPartition, minCountEachPartition);
  }

    @Override
    public int[] getVariableIndicesToSearchOver() {
        int parsedVariables = Data.getParsedVarCount();	
        // SFIToolkit.displayln(" CHECKCHECK11 " + parsedVariables);
        searchArray = new int[parsedVariables-2];
        // SFIToolkit.displayln(" CHECKCHECK12 " + searchArray.length);
        for (int i=0; i<searchArray.length; i++)
        {
        	searchArray[i] = i + 1;
        }
        // SFIToolkit.displayln(" CHECKCHECK13 " + searchArray[0] + searchArray[1] );
        return searchArray;
    }

    @Override
    public MomentPartitionObj getMomentPartitionObj(Matrix nodeX, Matrix nodeY, int indexSplitVariable, IntegerPartition partition) {
        return new MomentPartitionObjBartRCT(partition, indexSplitVariable, nodeX, nodeY);
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
    public Jama.Matrix cvparameters() {
        return CVparameters;
    }
    
    @Override
    public int numberoftrees() {
        return numtrees;
    }
    
    @Override  // using the parameter row from Stata to identify whether discrete or continuous
    public Boolean[] getDiscreteVector() {
    	    	
    	double parameter;
    	long obsEnd = Data.getObsParsedIn2();
        int  numvar = Data.getParsedVarCount();
        // SFIToolkit.displayln(" CHECKCHECK21 " + numvar);
    	Boolean[] discreteIndicator = new Boolean[numvar - 1];
    	// SFIToolkit.displayln(" CHECKCHECK22 " + discreteIndicator.length);
        for(int i=0; i < discreteIndicator.length; i++){
        	parameter = Data.getNum(i+2,obsEnd - 1 - 6); // We need to skip those numtree and CV parameter rows
			if (parameter >= 1) {
				discreteIndicator[i] = true;
    		} else {
    			discreteIndicator[i] = false;
    		}
        }
        // SFIToolkit.displayln(" CHECKCHECK32 " + discreteIndicator[0]+ discreteIndicator[1]+ discreteIndicator[2]);
        // the below was the original predefined discreteIndicator
        // Boolean[] discreteIndicator = {true, false, true, true, true, true, true, true, true, false};
    	
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


    @Override  
	   public void loadData() {
  	
	   int varIndex_x, varIndex_y;
	   int rc ;
	   double value_x, value_y;
	   String msg_x, msg_y;
	   // rc = 0 ;
	   TreeSet<Integer> exclusionTree = new TreeSet<>();
       // Get number of variables in varlist specified to javacall
	   int nVariables = Data.getParsedVarCount();
	   SFIToolkit.displayln("How many variables?  " + nVariables); 
	   // Get first observation specified by an in restriction
	   long firstObs = Data.getObsParsedIn1();
	   // Get last observation specified by an in restriction
	    long lastObs = Data.getObsParsedIn2();
	   // create and initialize MyLong for sample size
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
	   
	   // We should not count the last discrete/continuous parameter row
	   // We should not count the number of tree parameter row
	   // We would not count the last 6 CV parameter rows
	   Y = new Jama.Matrix(counter - 1 - 1 - 6 - exclusionTree.size(), 1);
	   X = new Jama.Matrix(counter - 1 - 1 - 6 - exclusionTree.size(), nVariables-1);
	   
	   numtrees = (int) Data.getNum(1,nObs - 1 - 5);
	   
	   CVparameters = new Jama.Matrix(1,6);
	   for (int cv = 5; cv >= 0; cv-- ) {
		   CVparameters.set(0,5-cv, Data.getNum(1,nObs - cv));  
	   }
	    
	   
	// Loop over variables
	   		for (long obs_y = 1; obs_y <= nObs - 1 - 1 - 6; obs_y++ ) {
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
	   	   for (long obs_x = 1; obs_x<=nObs - 1 - 1 - 6; obs_x++ ) {
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
		  // SFIToolkit.displayln("Report!!: " + searchArray.length);  
		  // SFIToolkit.displayln("Loaded " + X.get(3, 1)); // + " observations after dropping " + exclusionTree.size() + " observations for missing data.");
		  // SFIToolkit.displayln("Loaded " + Y.get(0, 0) ); // + " observations after dropping " + exclusionTree.size() + " observations for missing data.");
	  	  // SFIToolkit.displayln("Loaded " + Y.get(3, 0) ); // + " observations after dropping " + exclusionTree.size() + " observations for missing data.");
	    // return(rc) ;  
    
  
    	
/* this is the original import data code */
//  try {
//      BufferedReader inCounter = new BufferedReader(new FileReader("simulation_er.csv"));
//      String line = inCounter.readLine(); // headers
//      int counter = 0;
//      TreeSet<Integer> exclusionTree = new TreeSet<>();
//      while (line != null) {
//          if (line != null) {
//              try {
//                  line = inCounter.readLine();
//                  int a = 0;
//                  for (int j = 0; j < 10; j++) {
//                      int b = line.indexOf(",", a);
//                      Double.parseDouble(line.substring(a, b));
//                      a = b + 1;
//                  }
//                  Double.parseDouble(line.substring(a));
//              } catch (Exception e) {
//                  System.out.println(e);
//                  System.out.println("Missing data for observation " + counter);
//                  exclusionTree.add(counter);
//              }
//              counter++;
//          }
//      }
//      inCounter.close();

//      X = new Jama.Matrix(counter - exclusionTree.size(), 10);
//      Y = new Jama.Matrix(counter - exclusionTree.size(), 1);  
 /* set up a new jama matrix so that it can be filled by the following */

//      BufferedReader in = new BufferedReader(new FileReader("simulation_6000.csv"));
//      line = in.readLine(); // headers
//      int counter2 = 0;
//      for (int i = 0; i < counter; i++) {
//          line = in.readLine();
     // System.out.println(line);
//          if (!exclusionTree.contains(i)) {
//              int a = 0;
//              int b = line.indexOf(",", a);
         // System.out.println(line.substring(a, b));
//              Y.set(counter2, 0, Double.parseDouble(line.substring(a, b)));
//              a = b + 1;
//              b = line.indexOf(",", a);
         // System.out.println(line.substring(a, b));
//              X.set(counter2, 0, Double.parseDouble(line.substring(a, b)));
//              for (int j = 0; j < 9; j++) {
//                  a = b + 1;
//                  b = line.indexOf(",", a);
             // System.out.println(line.substring(a, b));
//                  if (j < 8) {
//                      X.set(counter2, j + 1, Double.parseDouble(line.substring(a, b)));
//                  } else {
//                      X.set(counter2, j + 1, Double.parseDouble(line.substring(a)));
//                  }
//              }
//              counter2++;
//          }
//      }
//      System.out.println("Loaded " + counter2 + " observations after dropping " + exclusionTree.size() + " observations for missing data.");
//      in.close();
//  } catch (Exception e) {
//      e.printStackTrace();
//      System.exit(0);
//  }

//    if (monteCarlo) {
//        // generateData(X.getRowDimension(), new Random(), true);
//        generateData(15000, new Random(), true);
//    }
      Data.addVarLong("beta_estimated");
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
        // final private String[] variableNames = {"treatment", "age", "symptomatic", "homosexual", "IVDrugUser", "PriorAZT", "Male", "White", "FutureDropout", "Baseline"};
        // String[][] feNames = {
        		// {"V0_1", "V0_2"}, // x0
        		// {"V1"}, // x1
        		// {"1", "2", "3", "4", "5", "6", "7","8","9","10","11"}, // x2
                //{"V3_1", "V3_2"}, // x3
                //{"V4_1", "V4_2"}, // x4
                //{"V5_1", "V5_2"}, // x5
                //{"V6_1", "V6_2"}, // x6
                //{"V7_1", "V7_2"}, // x7
                //{"V8_1", "V8_2"}, // x8
                //{"V9"} // x9        	
        		
            // {"Control", "Treatment"}, // x0
            // {"Age (continuous)"}, // x1
            // {"Not Symptomatic", "Symptomatic"}, // x2
            // {"Homosexual", "Heterosexual"}, // x3
            // {"Not IV Drug User", "IV Drug User"}, // x4
            // {"No Prior AZT", "Prior AZT"}, // x5
            // {"Female", "Male"}, // x6
            // {"Not White", "White"}, // x7
            // {"Not Dropout", "Future Dropout"}, // x8
            // {"Baseline (not coded)"} // x9
        // };
        // System.out.println(variableIndex+" "+fixedEffectIndex);
        // if(monteCarlo) {
            return "var[" + variableIndex + "]=" + fixedEffectIndex;
        // }
        // return feNames[variableIndex][fixedEffectIndex]; /* naming part */
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