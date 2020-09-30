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
import core.ContainerMoment;
import core.DataLens;
import core.MomentForest;
import core.TreeOptions;
import java.io.BufferedWriter;
import java.io.FileWriter;
import java.util.Random;
import utility.pmUtility;

/**
 *
 * @author Stephen P. Ryan <stephen.p.ryan@wustl.edu>
 */
public class SimpleSTAR {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {

        SimpleRCTMomentSpecification mySpecification = new SimpleRCTMomentSpecification(2000);
        mySpecification.loadData();

        int numberTreesInForest = 500;
        // System.out.println("numTrees: " + numberTreesInForest);

        /**
         * Initialize the moment forest
         */
        Jama.Matrix rctX = mySpecification.getX();
        Jama.Matrix rctY = mySpecification.getY();
        Jama.Matrix balancing = pmUtility.getColumn(mySpecification.getX(), 0);
        DataLens forestLens = new DataLens(rctX, rctY, balancing);

        /**
         * I want to compare to just resampling the data, running moment on each
         * resample, and computing bootstrapped standard errors on that
         */
        boolean manualBootstrapForest = false;
        if (manualBootstrapForest) {
            int B = 500;
            Random rng = new Random(787);
            Jama.Matrix results = new Jama.Matrix(B, 1);
            for (int b = 0; b < B; b++) {
                long s = rng.nextLong();
                // System.out.println("Random seed is: "+s);
                DataLens re = forestLens.getResampledDataLens(s);
                // System.out.println(re);
                ContainerMoment c = mySpecification.computeOptimalBeta(re);
                results.set(b, 0, c.getBeta().get(0, 0));
            }
            System.out.println("Mean estimate: " + pmUtility.mean(results, 0));
            System.out.println("Standard deviation: " + pmUtility.standardDeviation(results, 0));
        }

        boolean verbose = false;
        MomentForest myForest = new MomentForest(mySpecification, numberTreesInForest, 314, forestLens, verbose, new TreeOptions());

        double MSE_threshold = 1E-4;
        double proportionInEach = 1E-5;
        int minimumObsPerLeaf = 400;
        int maxDepth = 100;
        TreeOptions cvOptions = new TreeOptions(proportionInEach, minimumObsPerLeaf, MSE_threshold, maxDepth);
        /**
         * Run a CV for the hyper-parameters and see the tree options
         */
        boolean useCV = true;
        if (useCV) {
            int numTreesCrossValidation = 100;
            cvOptions = myForest.performCrossValidation(numTreesCrossValidation);
            myForest.setTreeOptions(cvOptions);
        }
        /**
         * Grow the moment forest
         */
        System.out.print("Growing forest...");
        long t1 = System.currentTimeMillis();
        myForest.growForest();
        long t2 = System.currentTimeMillis();
        System.out.println("done [" + (t2 - t1) / 1000.0 + " s]");
        
        myForest.getTree(0).printTree();

        /**
         * Compute standard errors
         */
        int numberBootstraps = 1;
        // System.out.println("Number of bootstraps: " + numberBootstraps);
        int numberTreesInBootForest = 500;
        System.out.print("\nComputing bootstrapped standard errors...\n");
        t1 = System.currentTimeMillis();
        BootstrapForest boot = new BootstrapForest(mySpecification, numberBootstraps, numberTreesInBootForest, 787, cvOptions);
        t2 = System.currentTimeMillis();
        System.out.println("done [" + (t2 - t1) / 1000.0 + " s]");

        Jama.Matrix fitX = mySpecification.getX();
        System.out.println("\nMoment Forest Estimates In Sample");

        // we can do better than that
        // we can enumerate all the discrete combinations of student characteristics
        /**
         * dX.set(i, 0, Double.valueOf(line.substring(a, b))); // treatment
         * status dX.set(i, 1, Double.valueOf(line.substring(a, b))); // year
         * star (?) dX.set(i, 2, Double.valueOf(line.substring(a, b))); // year
         * small (?) dX.set(i, 3, Double.valueOf(line.substring(a, b))); // is
         * white dX.set(i, 4, Double.valueOf(line.substring(a, b))); // is free
         * lunch dX.set(i, 5, Double.valueOf(line.substring(a, b))); // is
         * female // dX.set(i, 5, Double.valueOf(line.substring(a, b))); //
         * school wave FE dX.set(i, 6, Double.valueOf(line.substring(a))); //
         * outcome indicator (runs from 1 to 10)
         */
        if (mySpecification.usingSTARData()) {
            if (mySpecification.isSplittingOnObservables()) {
                fitX = new Jama.Matrix(2 * 2 * 2 * 10, fitX.getColumnDimension());
                int counter = 0;
                for (int white = 0; white <= 1; white++) {
                    for (int freelunch = 0; freelunch <= 1; freelunch++) {
                        for (int female = 0; female <= 1; female++) {
                            for (int outcome = 1; outcome <= 10; outcome++) {
                                fitX.set(counter, 0, 1); // treatment
                                fitX.set(counter, 1, 0); // year (not using now)
                                fitX.set(counter, 2, 0); // year small
                                fitX.set(counter, 3, white);
                                fitX.set(counter, 4, freelunch);
                                fitX.set(counter, 5, female);
                                fitX.set(counter, 6, outcome);
                                counter++;
                            }
                        }
                    }
                }
            } else {
                fitX = new Jama.Matrix(10, fitX.getColumnDimension());
                int counter = 0;
                for (int outcome = 1; outcome <= 10; outcome++) {
                    fitX.set(counter, 0, 1); // treatment
                    fitX.set(counter, 1, 0); // year (not using now)
                    fitX.set(counter, 2, 0); // year small
                    fitX.set(counter, 3, 0);
                    fitX.set(counter, 4, 0);
                    fitX.set(counter, 5, 0);
                    fitX.set(counter, 6, outcome);
                    counter++;
                }
            }
        } else {
            fitX = new Jama.Matrix(1, 2);
            fitX.set(0, 0, 1); // treatment
            fitX.set(0, 1, 1); // group 1 (doesn't matter)
        }

        try {
            BufferedWriter out = new BufferedWriter(new FileWriter("../star_results" + mySpecification.getSpecName() + ".csv"));
            out.write("treatment,year,yearsmall,white,freelunch,female,outcome,parameter,se,t-stat\n");
            System.out.print("treatment,year,yearsmall,white,freelunch,female,outcome,parameter,se,t-stat\n");
            for (int i = 0; i < fitX.getRowDimension(); i++) {
                Jama.Matrix xi = fitX.getMatrix(i, i, 0, mySpecification.getX().getColumnDimension() - 1);
                Jama.Matrix estimatedTreatmentEffects = myForest.getEstimatedParameters(xi);
                Jama.Matrix standardErrors = estimatedTreatmentEffects.times(0);
                boolean useBoot = true;
                if (useBoot) {
                    // standardErrors = boot.computeStandardErrors(xi);
                    standardErrors = boot.computeStandardErrorsSecondWay(xi);
                }
                String sig = "";
                if (Math.abs(estimatedTreatmentEffects.get(0, 0) / standardErrors.get(0, 0)) > 1.98) {
                    sig = "*";
                }
                // System.out.format("%g %g (%g) %s %n", fitX.get(i, 6), estimatedTreatmentEffects.get(0, 0), standardErrors.get(0, 0), sig);
                for (int j = 0; j < xi.getColumnDimension(); j++) {
                    out.write(xi.get(0, j) + ",");
                    System.out.print(xi.get(0, j) + ",");
                }
                out.write(estimatedTreatmentEffects.get(0, 0) + "," + standardErrors.get(0, 0) + "," + (estimatedTreatmentEffects.get(0, 0) / standardErrors.get(0, 0)) + "\n");
                System.out.print(estimatedTreatmentEffects.get(0, 0) + "," + standardErrors.get(0, 0) + "," + (estimatedTreatmentEffects.get(0, 0) / standardErrors.get(0, 0)) + "\n");
            }
            out.close();
        } catch (Exception e) {
            e.printStackTrace();
        }

        mySpecification.computeNaiveStatistics();
    }
}
