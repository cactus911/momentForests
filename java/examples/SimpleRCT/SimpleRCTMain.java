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

// import com.stata.sfi.SFIToolkit;
import core.BootstrapForest;
import core.DataLens;
import core.MomentForest;
import core.MomentSpecification;
import core.TreeOptions;
import utility.pmUtility;

/**
 *
 * @author Stephen P. Ryan <stephen.p.ryan@wustl.edu>
 */
public class SimpleRCTMain {

    Jama.Matrix EstimationResults;
    
    double minProportion;
    int minCount;
    double minMSEImprovement;
    int maxDepth;
            
    public SimpleRCTMain(Jama.Matrix X, Jama.Matrix Y, int numtrees, Jama.Matrix CVparameters1, Jama.Matrix CVparameters2, boolean cv, int[] variableSearchIndex, Boolean[] DiscreteVariables, int numbootstrap) {
            

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
            
            // Write down an option here
            MomentSpecification mySpecification = new SimpleRCTMomentSpecification(X, Y, numtrees, variableSearchIndex, DiscreteVariables);
            // mySpecification.loadData();

            // int numberTreesInForest = mySpecification.numberoftrees();
            // System.out.println("numTrees: " + numberTreesInForest);

            /**
             * Initialize the moment forest
             */
            DataLens forestLens = new DataLens(mySpecification.getX(), mySpecification.getY(), pmUtility.getColumn(mySpecification.getX(), 0));
            boolean verbose = true;
            MomentForest myForest = new MomentForest(mySpecification, numtrees, 314, forestLens, verbose, new TreeOptions());
            
            TreeOptions cvOptions = new TreeOptions(1E-5, (int) CVparameters1.get(0,0), (double) CVparameters1.get(0,1), 100);
            /**
             * Run a CV for the hyper-parameters and see the tree options
             */
            boolean useCV = cv;
            if (useCV) {
                // int numTreesCrossValidation = mySpecification.numberoftrees();
                cvOptions = myForest.performCrossValidation(numtrees, CVparameters2);
            } 
            
            /* Read out hyper parameters */ 
            minProportion = cvOptions.getMinProportion();
            minCount = cvOptions.getMinCount();
            minMSEImprovement = cvOptions.getMinMSEImprovement();
            maxDepth = cvOptions.getMaxDepth();
            
            myForest.setTreeOptions(cvOptions);
            
            // SFIToolkit.displayln("test.getMincount " + cvOptions.getMinCount() + "test.MinMSEImprove " + cvOptions.getMinMSEImprovement());
            
            /**
             * Grow the moment forest
             */
            myForest.growForest();

            /**
             * Compute standard errors
             */
            int numberBootstraps = numbootstrap;
            // System.out.println("Number of bootstraps: " + numberBootstraps);
            int numberTreesInBootForest = mySpecification.numberoftrees();
            BootstrapForest boot = new BootstrapForest(mySpecification, numberBootstraps, numberTreesInBootForest, 787, cvOptions);
            
            /**
             * stack the estimation results
             */                 
            
            Jama.Matrix allX = mySpecification.getX().copy();
            EstimationResults = new Jama.Matrix(allX.getRowDimension(), 2); //Treatment effect, standard error
            
            for (int i = 0; i < allX.getRowDimension(); i++) {
                Jama.Matrix xi = allX.getMatrix(i, i, 0, mySpecification.getX().getColumnDimension() - 1);  
                Jama.Matrix estimatedTreatmentEffects = myForest.getEstimatedParameters(xi);
                Jama.Matrix standardErrors = estimatedTreatmentEffects.times(0);

                boolean useBoot = true;
                if (useBoot) {
                    standardErrors = boot.computeStandardErrors(xi);
                }
                
                EstimationResults.set(i, 0 ,estimatedTreatmentEffects.get(0, 0)); //Only coefficient is treatment effect
                EstimationResults.set(i, 1 ,standardErrors.get(0, 0));
                // Data.storeNum(resv+2, i+1, estimatedTreatmentEffects.get(0, 0));
                // Data.storeNum(resv+3, i+1, standardErrors.get(0, 0));
            }
           
    }

    
    
    
    public Jama.Matrix EstimationResults() {
        return EstimationResults;
    }

    public double minProportion() {
        return minProportion;
    }
    
    public int minCount() {
        return minCount;
    }

    public double minMSEImprovement() {
        return minMSEImprovement;
    }

    public int maxDepth() {
        return maxDepth;
    }    
        
}
