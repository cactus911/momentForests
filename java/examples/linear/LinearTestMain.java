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
        MomentSpecification mySpecification = new LinearMomentSpecification(1000);
        mySpecification.loadData(); // Create data using rng

        int numberTreesInForest = 100;
        // System.out.println("numTrees: " + numberTreesInForest);

        /**
         * Initialize the moment forest
         */
        DataLens forestLens = new DataLens(mySpecification.getX(), mySpecification.getY(), null);
        /* Contains X data, Y data, balancing vector (treatment indicators), and data index (just an array numbered 0 - numObs) */
        boolean verbose = false;
        MomentForest myForest = new MomentForest(mySpecification, numberTreesInForest, 314, forestLens, verbose, new TreeOptions());

        TreeOptions cvOptions = new TreeOptions(0.01, 100, 1E-1, 20); // k = 1
        myForest.setTreeOptions(cvOptions);
        /**
         * Grow the moment forest
         */
        myForest.growForest();

        /**
         * Compute standard errors
         */
        int numberBootstraps = 50;
        // System.out.println("Number of bootstraps: " + numberBootstraps);

        int numberTreesInBootForest = 10;
        BootstrapForest boot = new BootstrapForest(mySpecification, numberBootstraps, numberTreesInBootForest, 787, cvOptions);

        Jama.Matrix testX = new Jama.Matrix(4, 4);
        testX.set(0, 2, -1);
        testX.set(0, 3, -1);

        testX.set(1, 2, 1);
        testX.set(1, 3, -1);

        testX.set(2, 2, -1);
        testX.set(2, 3, 1);

        testX.set(3, 2, 1);
        testX.set(3, 3, 1);

        for (int j = 0; j < 4; j++) {
            Jama.Matrix txj = testX.getMatrix(j, j, 0, 3);
            Jama.Matrix b = myForest.getEstimatedParameters(txj);
            System.out.println("z1: " + testX.get(j, 2) + " z2: " + testX.get(j, 3) + " beta: " + pmUtility.stringPrettyPrintVector(b)+" se: "+pmUtility.stringPrettyPrintVector(boot.computeStandardErrors(txj)));
        }
    }

}
