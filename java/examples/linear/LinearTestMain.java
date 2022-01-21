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
package examples.linear;

import core.BootstrapForest;
import core.DataLens;
import core.HomogeneousSearchContainer;
import core.MomentForest;
import core.MomentSpecification;
import core.TreeOptions;
import java.util.ArrayList;
import optimization.Uncmin_f77;
import optimization.Uncmin_methods;
import utility.pmUtility;

/**
 *
 * @author Stephen P. Ryan <stephen.p.ryan@wustl.edu>
 */
public class LinearTestMain {
    // Replicate table 4: RCT with saturated heterogeneity 

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        LinearTestMain go = new LinearTestMain();
    }

    public LinearTestMain() {
        // MomentSpecification mySpecification = new LinearMomentSpecification("data/airline_subset.csv");
        MomentSpecification mySpecification = new LinearMomentSpecification(2000);
        mySpecification.loadData(); // Create data using rng

        double bestMinImprovement = 0;
        int bestMinObservationsPerLeaf = 0;
        double bestAlpha = -3;

        double lowestSSE = 0;
        boolean first = true;

        /**
         * Need to figure out how to implement a pre-estimation test to see
         * which parameters are globally constant. Our first idea is to do a
         * single split and then test at that first split.
         */
        int maxTreeDepth = 1;
        /**
         * January 21, 2022: there is a question here about how to estimate the
         * set of homogeneous parameters in that should we do so after running
         * cross validation? I think the answer to that is probably yes since it
         * is primarily an issue of imposing a constraint for efficiency
         * reasons.
         *
         * Probably do this:
         *
         * 1. Run cross validation on unrestricted model 2. Run test of
         * homogeneity 3. Impose homogeneity and re-estimate
         */

        /**
         * Make sure that cross-validation is run on completely unrestricted
         * model; set all parameters to heterogeneous
         */
        mySpecification.resetHomogeneityIndex();

        int numberTreesInForest = 1;
        // System.out.println("numTrees: " + numberTreesInForest);

        /**
         * Initialize the moment forest
         *
         * Note that dataLens changes as we impose homogeneity, so need to
         * recall this after changing the homogeneity parameters in
         * MommentSpecification
         */
        DataLens forestLens = new DataLens(mySpecification.getX(), mySpecification.getY(true), mySpecification.getZ(), null);
        /* Contains X data, Y data, balancing vector (treatment indicators), and data index (just an array numbered 0 - numObs) */
        boolean verbose = true;
        boolean testParameterHomogeneity = false;

        System.out.println("************************");
        System.out.println("* Run Cross-Validation *");
        System.out.println("************************");

        for (double minImprovement = 1.0; minImprovement <= 1.0; minImprovement *= 10) {
            for (int minObservationsPerLeaf = 50; minObservationsPerLeaf <= 50; minObservationsPerLeaf *= 2) {
                MomentForest myForest = new MomentForest(mySpecification, numberTreesInForest, 314, forestLens, verbose, new TreeOptions());
                System.out.println("Minimum Improvement Threshold: " + minImprovement);
                System.out.println("Minimum Observations per Leaf: " + minObservationsPerLeaf);

                TreeOptions cvOptions = new TreeOptions(0.01, minObservationsPerLeaf, minImprovement, maxTreeDepth, testParameterHomogeneity); // k = 1
                myForest.setTreeOptions(cvOptions);
                /**
                 * Grow the moment forest
                 */
                myForest.growForest();

                myForest.getTree(0).printTree();

                /**
                 * Test vectors for assessment
                 */
//                DataLens oosDataLens = mySpecification.getOutOfSampleXYZ(100);
//                Jama.Matrix testZ = oosDataLens.getZ();
//                // pmUtility.prettyPrint(pmUtility.concatMatrix(pmUtility.concatMatrix(oosDataLens.getY(), oosDataLens.getX()), testZ));
//
//                /**
//                 * Compute out-of-sample fit for CV purposes
//                 */
//                double outOfSampleFit = 0;
//                for (int i = 0; i < testZ.getRowDimension(); i++) {
//                    Jama.Matrix zi = testZ.getMatrix(i, i, 0, testZ.getColumnDimension() - 1);
//                    Jama.Matrix b = myForest.getEstimatedParameterForest(zi);
//
//                    // System.out.print("z: " + pmUtility.stringPrettyPrint(zi) + " beta: " + pmUtility.stringPrettyPrintVector(b));
//
//                    Jama.Matrix xi = oosDataLens.getX().getMatrix(i, i, 0, mySpecification.getX().getColumnDimension() - 1);
//                    double fitY = xi.times(b).get(0, 0) + mySpecification.getHomogeneousComponent(xi);
//                    double error = fitY - (oosDataLens.getY().get(i, 0));
//                    // System.out.println(" y_i = " + oosDataLens.getY().get(i, 0) + " fitY: " + fitY + " SE: " + error*error);
//                    outOfSampleFit += error * error;
//                }
                DataLens oosDataLens = mySpecification.getOutOfSampleXYZ(100); // this should eventually be modified to come out of the data itself (or generalized in some way)
                Jama.Matrix testZ = oosDataLens.getZ();
                Jama.Matrix residualizedX = mySpecification.residualizeX(oosDataLens.getX());
                // pmUtility.prettyPrint(pmUtility.concatMatrix(pmUtility.concatMatrix(oosDataLens.getY(), oosDataLens.getX()), testZ));

                /**
                 * Compute out-of-sample fit at current homogeneous parameter
                 * vector
                 */
                double outOfSampleFit = 0;
                for (int i = 0; i < testZ.getRowDimension(); i++) {
                    Jama.Matrix zi = testZ.getMatrix(i, i, 0, testZ.getColumnDimension() - 1);
                    Jama.Matrix b = myForest.getEstimatedParameterForest(zi);
                    System.out.print("z: " + pmUtility.stringPrettyPrint(zi) + " beta: " + pmUtility.stringPrettyPrintVector(b));
                    Jama.Matrix residualizedXi = residualizedX.getMatrix(i, i, 0, residualizedX.getColumnDimension() - 1);
                    Jama.Matrix fullXi = oosDataLens.getX().getMatrix(i, i, 0, oosDataLens.getX().getColumnDimension() - 1);

                    System.out.print(" x: " + pmUtility.stringPrettyPrint(fullXi));
                    double fitY = residualizedXi.times(b).get(0, 0) + mySpecification.getHomogeneousComponent(fullXi);
                    System.out.print(" residualizedXb: " + residualizedXi.times(b).get(0, 0) + " hc: " + mySpecification.getHomogeneousComponent(fullXi));
                    double error = fitY - (oosDataLens.getY().get(i, 0));
                    System.out.println(" fitY: " + fitY + " Y: " + oosDataLens.getY().get(i, 0) + " SE: " + error * error);
                    outOfSampleFit += error * error;
                }

                double MSE = outOfSampleFit / testZ.getRowDimension();
                System.out.println("Out-of-sample MSE: " + MSE);
                if (MSE < lowestSSE || first) {
                    lowestSSE = MSE;
                    first = false;
                    bestMinImprovement = minImprovement;
                    bestMinObservationsPerLeaf = minObservationsPerLeaf;
                }
            }
        }

        System.out.println("************************");
        System.out.println("* Test for Homogeneity *");
        System.out.println("************************");
        /**
         * Step 2: determine homogeneous parameters post-CV
         */
        double minProportionInEachLeaf = 0.01;
        /**
         * Impose a max depth of 1 to test at the first split
         */
        maxTreeDepth = 1;
        testParameterHomogeneity = true;
        TreeOptions cvOptions = new TreeOptions(minProportionInEachLeaf, bestMinObservationsPerLeaf, bestMinImprovement, maxTreeDepth, testParameterHomogeneity); // k = 1
        MomentForest myForest = new MomentForest(mySpecification, numberTreesInForest, 314, forestLens, verbose, new TreeOptions());

        myForest.setTreeOptions(cvOptions);
        myForest.growForest();
        ArrayList<Integer> homogeneousParameterList = myForest.getTree(0).getIndexHomogeneousParameters(); // this is only using the first tree, is that the right way of thinking about this?
        int numParams = homogeneousParameterList.size();
        mySpecification.resetHomogeneityIndex();
        for (Integer i : homogeneousParameterList) {
            mySpecification.setHomogeneousIndex(i);
        }

        /**
         * I wonder if there should be an additional step here to recompute the
         * optimal CV parameters given the results of the homogeneity test in
         * the prior step?
         */
        System.out.println("******************************");
        System.out.println("* Estimate Constrained Model *");
        System.out.println("******************************");
        /**
         * Step 3: estimate constrained model
         */
        maxTreeDepth = 1;
        HomogeneousSearchContainer con = new HomogeneousSearchContainer(mySpecification, numberTreesInForest, verbose, bestMinImprovement, bestMinObservationsPerLeaf, maxTreeDepth, numParams);
        con.executeSearch();
        Jama.Matrix homogeneousParameters = con.getEstimatedHomogeneousParameters();
        System.out.println("Estimated homogeneous parameters: ");
        pmUtility.prettyPrintVector(homogeneousParameters);

        boolean computeSE = false;
        if (computeSE) {
            int numberBootstraps = 50;
            // System.out.println("Number of bootstraps: " + numberBootstraps);

            int numberTreesInBootForest = 10;

            BootstrapForest boot = new BootstrapForest(mySpecification, numberBootstraps, numberTreesInBootForest, 787, cvOptions);
            /**
             * TODO: have to make sure that the homogeneous parameter thing is
             * working inside the boot forest ALSO, how to compute standard
             * errors on homogeneous parameters here COME BACK TO THIS!!!
             */

            DataLens oosDataLens = mySpecification.getOutOfSampleXYZ(100);
            Jama.Matrix testZ = oosDataLens.getZ();

            for (int j = 0; j < testZ.getRowDimension(); j++) {
                Jama.Matrix zi = testZ.getMatrix(j, j, 0, testZ.getColumnDimension() - 1);
                Jama.Matrix b = myForest.getEstimatedParameterForest(zi);
                System.out.println("z: " + pmUtility.stringPrettyPrint(zi) + " beta: " + pmUtility.stringPrettyPrintVector(b) + " se: " + pmUtility.stringPrettyPrintVector(boot.computeStandardErrors(zi)));
            }
        }

        System.out.println("Lowest MSE: " + lowestSSE + " at min_N = " + bestMinObservationsPerLeaf + " min_MSE = " + bestMinImprovement + " bestAlpha: " + bestAlpha);
    }

}
