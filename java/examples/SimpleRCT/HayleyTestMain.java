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

import core.DataLens;
import core.MomentForest;
import core.MomentSpecification;
import core.TreeOptions;
import java.util.Arrays;
import utility.pmUtility;

/**
 *
 * @author Stephen P. Ryan <stephen.p.ryan@wustl.edu> 
 * @author Hayley Wabiszewski <hwabiszewski@wustl.edu>
 */ 

public class HayleyTestMain { 
    // Replicate table 4: RCT with saturated heterogeneity 

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        /**
         * Aiming for a simple three-step process:
         *
         * 1. Load the data 2. Specify which variables are splitting, and, of
         * those, which are discrete (this is embodied in MomentSpecification)
         * 3. Specify the objective function
         *
         * After that, we should have a default implementation that allows for:
         * 1. Running CV to obtain optimal hyper-parameters 2. Estimate a moment
         * forest 3. Estimate standard errors via bootstrapping
         */
        
        int Number_Simulations = 1; // 500; //Number of time we Monte Carlo
        int Number_Rows = 1; //Number of rows in the table I will replicate, can only do one at a time bc importing data
        Jama.Matrix Number_Leaves = new Jama.Matrix(Number_Rows,Number_Simulations);
        Jama.Matrix MSPE_Forest = new Jama.Matrix(Number_Rows,Number_Simulations);
        int Monte_Carlo_Index = -1;
        
        
        for (int n = 5000; n <= 5000; n *= 2) {  // this is how many observations to use in each monte carlo
            Monte_Carlo_Index += 1;
            int Data_Index = 0;
            for (int Loop_Number = 0; Loop_Number < Number_Simulations; Loop_Number += 1) { //Monte carlo using 500 simulations
                // System.out.println("*********** n = " + n + " (really one-tenth of that, since n here = observations*outcomes) ***********");
                MomentSpecification mySpecification = new SimpleRCTMomentSpecification(n);
                mySpecification.loadData(); // Create data using rng

                int numberTreesInForest = 1;
                // System.out.println("numTrees: " + numberTreesInForest);

                /**
                 * Initialize the moment forest
                 */
                Jama.Matrix rctX1 = mySpecification.getX();
                Jama.Matrix rctY1 = mySpecification.getY();
                Jama.Matrix balancing1 = pmUtility.getColumn(mySpecification.getX(), 0); // Treatment indicators
                // pmUtility.prettyPrint(rctX);
                
                Jama.Matrix rctX = rctX1.getMatrix(Data_Index, Data_Index + n - 1, 0, 2);
                Jama.Matrix rctY = rctY1.getMatrix(Data_Index, Data_Index + n - 1, 0, 0);
                Jama.Matrix balancing = balancing1.getMatrix(Data_Index, Data_Index + n - 1, 0, 0);
                
                DataLens forestLens = new DataLens(rctX, rctY, balancing); /* Contains X data, Y data, balancing vector (treatment indicators), and data index (just an array numbered 0 - numObs) */
                boolean verbose = true;
                MomentForest myForest = new MomentForest(mySpecification, numberTreesInForest, 314, forestLens, verbose, new TreeOptions());

                TreeOptions cvOptions = new TreeOptions(0.01, 1, 1E-3, 100); // k = 1
                myForest.setTreeOptions(cvOptions);
                /**
                 * Run a CV for the hyper-parameters and see the tree options
                 */
                /*
                boolean useCV = false; //Turn of CV for now... takes too long
                if (useCV) {
                    int numTreesCrossValidation = 1;
                    cvOptions = myForest.performCrossValidation(numTreesCrossValidation);
                    myForest.setTreeOptions(cvOptions);
                }
                */
                /**
                 * Grow the moment forest
                 */
                myForest.growForest();

                /**
                 * Compute standard errors
                 */
                
                /* Don't need standard errors for Monte Carlo    
                int numberBootstraps = 50;
                // System.out.println("Number of bootstraps: " + numberBootstraps);
                
                int numberTreesInBootForest = 10;
                BootstrapForest boot = new BootstrapForest(mySpecification, numberBootstraps, numberTreesInBootForest, 787, cvOptions);
                */

                Jama.Matrix fitX = new Jama.Matrix(10, 2);
                boolean isMonteCarlo = false; // Not sure what this does
                if (isMonteCarlo) {
                    /**
                     * Show fits for out of sample data
                     */
                    for (int i = 0; i < fitX.getRowDimension(); i++) {
                        fitX.set(i, 1, i); 
                    }
                } else {
                    //fitX = mySpecification.getX();
                    fitX = rctX;
                }

                //System.out.println("\nMoment Forest Estimates by Group");
                
                double MSPE = 0.0;       
                for (int i = 0; i < fitX.getRowDimension(); i++) {
                    Jama.Matrix xi = fitX.getMatrix(i, i, 0, mySpecification.getX().getColumnDimension()-1);
                    Jama.Matrix estimatedTreatmentEffects = myForest.getEstimatedParameters(xi);
                    //Jama.Matrix standardErrors = estimatedTreatmentEffects.times(0);
                    
                    /* Don't need standard errors for Monte Carlo
                    boolean useBoot = false; 
                    if (useBoot) {
                        standardErrors = boot.computeStandardErrors(xi);
                    }
                    */
                    
                    /* String sig = "";
                    if (Math.abs(estimatedTreatmentEffects.get(0, 0) / standardErrors.get(0, 0)) > 1.98) {
                        sig = "*";
                    // I think I can add something here to calculate the SPE for each estimate    
                    }
                    System.out.format("x = (%g,%g) tau = %g (%g) %s %n", fitX.get(i, 1), fitX.get(i, 2), estimatedTreatmentEffects.get(0, 0), standardErrors.get(0, 0), sig);
                    */
                    
                    MSPE += Math.pow(rctY.get(i, 0)-estimatedTreatmentEffects.get(0, 0)*fitX.get(i, 0),2);
                }
                MSPE /= n;
                
                double Leaves = 0.0;
                for (int i = 0; i < numberTreesInForest; i++) {
                    Leaves += myForest.getTree(i).countTerminalNodes();
                }
                Leaves /= numberTreesInForest;
                              
            Number_Leaves.set(Monte_Carlo_Index,Loop_Number,Leaves);
            MSPE_Forest.set(Monte_Carlo_Index,Loop_Number,MSPE);
            Data_Index += n;
            }            
        }
        double[][] array1 = Number_Leaves.getArray();
        double[][] array2 = MSPE_Forest.getArray();
        
        System.out.println("numLeaves " + Arrays.deepToString(array1));
        System.out.println("MSPE Forest " + Arrays.deepToString(array2));
    }

}

