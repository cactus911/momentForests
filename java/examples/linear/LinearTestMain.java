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
import core.MomentForest;
import core.MomentSpecification;
import core.TreeOptions;
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

        // MomentSpecification mySpecification = new LinearMomentSpecification("data/airline_subset.csv");
        MomentSpecification mySpecification = new LinearMomentSpecification(1500);
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
        int maxTreeDepth = 25;

        for (double alpha = -1; alpha <= -1; alpha++) {
            mySpecification.setHomogeneousParameters(new Jama.Matrix(1, 1, alpha));

            int numberTreesInForest = 1;
            // System.out.println("numTrees: " + numberTreesInForest);

            /**
             * Initialize the moment forest
             */
            DataLens forestLens = new DataLens(mySpecification.getX(), mySpecification.getY(true), mySpecification.getZ(), null);
            /* Contains X data, Y data, balancing vector (treatment indicators), and data index (just an array numbered 0 - numObs) */
            boolean verbose = true;
            MomentForest myForest = new MomentForest(mySpecification, numberTreesInForest, 314, forestLens, verbose, new TreeOptions());

            for (double minImprovement = 1E-35; minImprovement <= 1E-35; minImprovement *= 10) {
                for (int minObservationsPerLeaf = 50; minObservationsPerLeaf <= 50; minObservationsPerLeaf *= 2) {
                    System.out.println("Alpha: "+alpha);
                    System.out.println("Minimum Improvement Threshold: " + minImprovement);
                    System.out.println("Minimum Observations per Leaf: " + minObservationsPerLeaf);
                    TreeOptions cvOptions = new TreeOptions(0.01, minObservationsPerLeaf, minImprovement, maxTreeDepth); // k = 1
                    myForest.setTreeOptions(cvOptions);
                    /**
                     * Grow the moment forest
                     */
                    myForest.growForest();

                    myForest.getTree(0).printTree();

                    /**
                     * Test vectors for assessment
                     */
                    DataLens oosDataLens = mySpecification.getOutOfSampleXYZ(100);
                    Jama.Matrix testZ = oosDataLens.getZ();
                    // pmUtility.prettyPrint(pmUtility.concatMatrix(pmUtility.concatMatrix(oosDataLens.getY(), oosDataLens.getX()), testZ));

                    boolean computeSE = false;
                    /**
                     * Compute standard errors
                     */
                    if (computeSE) {
                        int numberBootstraps = 50;
                        // System.out.println("Number of bootstraps: " + numberBootstraps);

                        int numberTreesInBootForest = 10;
                        BootstrapForest boot = new BootstrapForest(mySpecification, numberBootstraps, numberTreesInBootForest, 787, cvOptions);

                        for (int j = 0; j < testZ.getRowDimension(); j++) {
                            Jama.Matrix zi = testZ.getMatrix(j, j, 0, testZ.getColumnDimension() - 1);
                            Jama.Matrix b = myForest.getEstimatedParameterForest(zi);
                            System.out.println("z: " + pmUtility.stringPrettyPrint(zi) + " beta: " + pmUtility.stringPrettyPrintVector(b) + " se: " + pmUtility.stringPrettyPrintVector(boot.computeStandardErrors(zi)));
                        }
                    } else {
                        double inSampleFitSSE = 0;
                        for (int i = 0; i < testZ.getRowDimension(); i++) {
                            Jama.Matrix zi = testZ.getMatrix(i, i, 0, testZ.getColumnDimension() - 1);
                            Jama.Matrix b = myForest.getEstimatedParameterForest(zi);
                            // System.out.println("z: " + pmUtility.stringPrettyPrint(zi)+ " beta: " + pmUtility.stringPrettyPrintVector(b));
                            Jama.Matrix xi = oosDataLens.getX().getMatrix(i, i, 0, mySpecification.getX().getColumnDimension() - 1);
                            double fitY = xi.times(b).get(0, 0) + mySpecification.getHomogeneousComponent(xi);
                            double error = fitY - (oosDataLens.getY().get(i, 0));
                            inSampleFitSSE += error * error;
                        }
                        double MSE = inSampleFitSSE / testZ.getRowDimension();
                        System.out.println("Out-of-sample MSE: " + MSE);
                        if (MSE < lowestSSE || first) {
                            lowestSSE = MSE;
                            first = false;
                            bestMinImprovement = minImprovement;
                            bestMinObservationsPerLeaf = minObservationsPerLeaf;
                            bestAlpha = alpha;
                        }
                    }
                }
            }
        }
        System.out.println("Lowest MSE: " + lowestSSE + " at min_N = " + bestMinObservationsPerLeaf + " min_MSE = " + bestMinImprovement+" bestAlpha: "+bestAlpha);
    }

}
