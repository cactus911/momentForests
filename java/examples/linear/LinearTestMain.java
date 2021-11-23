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
        MomentSpecification mySpecification = new LinearMomentSpecification(5000);
        mySpecification.loadData(); // Create data using rng

        int numberTreesInForest = 1;
        // System.out.println("numTrees: " + numberTreesInForest);

        /**
         * Initialize the moment forest
         */
        DataLens forestLens = new DataLens(mySpecification.getX(), mySpecification.getY(), mySpecification.getZ(), null);
        /* Contains X data, Y data, balancing vector (treatment indicators), and data index (just an array numbered 0 - numObs) */
        boolean verbose = true;
        MomentForest myForest = new MomentForest(mySpecification, numberTreesInForest, 314, forestLens, verbose, new TreeOptions());

        TreeOptions cvOptions = new TreeOptions(0.01, 100, 1E-1, 20); // k = 1
        myForest.setTreeOptions(cvOptions);
        /**
         * Grow the moment forest
         */
        myForest.growForest();

        /**
         * Test vectors for assessment
         */
        Jama.Matrix testX = new Jama.Matrix(4, 2);
        testX.set(0, 0, -1);
        testX.set(0, 1, -1);

        testX.set(1, 0, 1);
        testX.set(1, 1, -1);

        testX.set(2, 0, -1);
        testX.set(2, 1, 1);

        testX.set(3, 0, 1);
        testX.set(3, 1, 1);

        boolean computeSE = false;
        /**
         * Compute standard errors
         */
        if (computeSE) {
            int numberBootstraps = 50;
            // System.out.println("Number of bootstraps: " + numberBootstraps);

            int numberTreesInBootForest = 10;
            BootstrapForest boot = new BootstrapForest(mySpecification, numberBootstraps, numberTreesInBootForest, 787, cvOptions);

            for (int j = 0; j < 4; j++) {
                Jama.Matrix txj = testX.getMatrix(j, j, 0, 1);
                Jama.Matrix b = myForest.getEstimatedParameterForest(txj);
                System.out.println("z1: " + testX.get(j, 0) + " z2: " + testX.get(j, 1) + " beta: " + pmUtility.stringPrettyPrintVector(b) + " se: " + pmUtility.stringPrettyPrintVector(boot.computeStandardErrors(txj)));
            }
        } else {
            for (int j = 0; j < 4; j++) {
                Jama.Matrix txj = testX.getMatrix(j, j, 0, 1);
                Jama.Matrix b = myForest.getEstimatedParameterForest(txj);
                System.out.println("z1: " + testX.get(j, 0) + " z2: " + testX.get(j, 1) + " beta: " + pmUtility.stringPrettyPrintVector(b));
            }
        }
    }

}
