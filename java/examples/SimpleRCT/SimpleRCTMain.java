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

import core.BootstrapForest;
import core.MomentForest;
import core.MomentSpecification;
import core.TreeOptions;

/**
 *
 * @author Stephen P. Ryan <stephen.p.ryan@wustl.edu>
 */
public class SimpleRCTMain {

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
        for (int n = 5000; n <= 16000; n *= 2) {
            System.out.println("*********** n = "+n+" (really one-tenth of that, since n here = observations*outcomes) ***********");
            MomentSpecification mySpecification = new SimpleRCTMomentSpecification(n);
            mySpecification.loadData();

            int numberTreesInForest = 500;
            // System.out.println("numTrees: " + numberTreesInForest);

            /**
             * Initialize the moment forest
             */
            MomentForest myForest = new MomentForest(mySpecification, numberTreesInForest, 314, mySpecification.getX(), mySpecification.getY(), false, new TreeOptions());
            /**
             * Run a CV for the hyper-parameters and see the tree options
             */
            int numTreesCrossValidation = 100;
            TreeOptions cvOptions = myForest.performCrossValidation(numTreesCrossValidation);
            myForest.setTreeOptions(cvOptions);
            /**
             * Grow the moment forest
             */
            myForest.growForest();

            /**
             * Compute standard errors
             */
            int numberBootstraps = 50;
            int numberTreesInBootForest = 50;
            BootstrapForest boot = new BootstrapForest(mySpecification, numberBootstraps, numberTreesInBootForest, 787, cvOptions);

            /**
             * Show fits for out of sample data
             */
            Jama.Matrix fitX = new Jama.Matrix(10, 2);
            for (int i = 0; i < fitX.getRowDimension(); i++) {
                fitX.set(i, 1, i);
            }

            for (int i = 0; i < fitX.getRowDimension(); i++) {
                Jama.Matrix xi = fitX.getMatrix(i, i, 0, mySpecification.getX().getColumnDimension() - 1);
                Jama.Matrix estimatedTreatmentEffects = myForest.getEstimatedParameters(xi);
                Jama.Matrix standardErrors = estimatedTreatmentEffects.times(0);
                boolean useBoot = true;
                if (useBoot) {
                    standardErrors = boot.computeStandardErrors(xi);
                }
                String sig = "";
                if (Math.abs(estimatedTreatmentEffects.get(0, 0) / standardErrors.get(0, 0)) > 1.98) {
                    sig = "*";
                }
                System.out.format("%g %g (%g) %s %n", fitX.get(i, 1), estimatedTreatmentEffects.get(0, 0), standardErrors.get(0, 0), sig);
            }

            mySpecification.computeNaiveStatistics();
        }

    }

}
