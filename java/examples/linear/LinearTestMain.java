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

        int numberTreesInForest = 25;
        // System.out.println("numTrees: " + numberTreesInForest);

        /**
         * Initialize the moment forest
         */
        DataLens forestLens = new DataLens(mySpecification.getX(), mySpecification.getY(), mySpecification.getZ(), null);
        /* Contains X data, Y data, balancing vector (treatment indicators), and data index (just an array numbered 0 - numObs) */
        boolean verbose = false;
        MomentForest myForest = new MomentForest(mySpecification, numberTreesInForest, 314, forestLens, verbose, new TreeOptions());

        double lowestSSE = 0;
        boolean first = true;
        double bestMinImprovement = 0;
        int bestMinObservationsPerLeaf = 0;

        for (double minImprovement = 10; minImprovement <= 1000; minImprovement *= 2) {
            for (int minObservationsPerLeaf = 10; minObservationsPerLeaf <= 1000; minObservationsPerLeaf *= 2) {
                System.out.println("Minimum Improvement Threshold: " + minImprovement);
                System.out.println("Minimum Observations per Leaf: " + minObservationsPerLeaf);
                TreeOptions cvOptions = new TreeOptions(0.01, minObservationsPerLeaf, minImprovement, 20); // k = 1
                myForest.setTreeOptions(cvOptions);
                /**
                 * Grow the moment forest
                 */
                myForest.growForest();

                myForest.getTree(0).printTree();

                /**
                 * Test vectors for assessment
                 */
                Jama.Matrix testZ = mySpecification.getZ();

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
                    for (int j = 0; j < testZ.getRowDimension(); j++) {
                        Jama.Matrix zi = testZ.getMatrix(j, j, 0, testZ.getColumnDimension() - 1);
                        Jama.Matrix b = myForest.getEstimatedParameterForest(zi);
                        // System.out.println("z: " + pmUtility.stringPrettyPrint(zi)+ " beta: " + pmUtility.stringPrettyPrintVector(b));
                        Jama.Matrix xi = mySpecification.getX().getMatrix(j, j, 0, mySpecification.getX().getColumnDimension() - 1);
                        double fitY = xi.times(b).get(0, 0);
                        double error = fitY - mySpecification.getY().get(j, 0);
                        inSampleFitSSE += error * error;
                    }
                    double MSE = inSampleFitSSE / testZ.getRowDimension();
                    System.out.println("In-sample MSE: " + MSE);
                    if (MSE < lowestSSE || first) {
                        lowestSSE = MSE;
                        first = false;
                        bestMinImprovement = minImprovement;
                        bestMinObservationsPerLeaf = minObservationsPerLeaf;
                    }
                }
            }
        }
        System.out.println("Lowest MSE: "+lowestSSE+" at min_N = "+bestMinObservationsPerLeaf+" min_MSE = "+bestMinImprovement);
    }

}
