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
package examples.logitVSL;

import core.OutOfSampleStatisticsContainer;
import Jama.Matrix;
import core.DataLens;
import core.HomogeneousSearchContainer;
import core.MomentForest;
import core.MomentSpecification;

import core.TreeOptions;
import java.util.ArrayList;
import java.util.Random;
import utility.pmUtility;

/**
 *
 * @author Stephen P. Ryan <stephen.p.ryan@wustl.edu>
 */
public class LogitVSLMain {

    private ArrayList<Integer> homogeneousParameterList;
    private final long rngSeed;
    private final int numObs;
    private double outOfSampleYMSE;
    private double estimatedBetaVersusTruthMSE;
    private Jama.Matrix estimatedHomogeneousParameters;
    private final boolean detectHomogeneity;

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        LogitVSLMain go = new LogitVSLMain();
        go.execute();
    }

    public LogitVSLMain() {
        this.rngSeed = 0;
        this.numObs = 0;
        this.detectHomogeneity = true;
    }

    private void execute() {
        Random rng = new Random(rngSeed);

        MomentSpecification mySpecification = new LogitVSLMomentSpecification("d://dropbox/vsl/Replication_Materials/MomentTreesVSL/data/data_full_hazard_ny_fake_FE.csv");
        mySpecification.loadData(rng.nextLong()); // Create data using rng

        double bestMinImprovement = 0.8;
        int bestMinObservationsPerLeaf = 25;
        int bestMaxDepth = 5;

        double lowestSSE = 0;
        boolean first = true;

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
        /* Contains X data, Y data, balancing vector (treatment indicators), and data index (just an array numbered 0 - numObs) */
        boolean verbose = true;
        boolean testParameterHomogeneity;

        long rngBaseSeedMomentForest = rng.nextLong();
        long rngBaseSeedOutOfSample = rng.nextLong();

        boolean runCV = false;
        if (runCV) {
            if (verbose) {
                System.out.println("************************");
                System.out.println("* Run Cross-Validation *");
                System.out.println("************************");
            }

            for (int minObservationsPerLeaf = 50; minObservationsPerLeaf <= 200; minObservationsPerLeaf *= 2) {
                for (double minImprovement = 0.1; minImprovement <= 12.1; minImprovement *= 10) {
                    for (int maxDepth = 4; maxDepth >= 0; maxDepth--) {
                        // for (int maxDepth = 1; maxDepth >= 1; maxDepth--) {
                        verbose = true;
                        OutOfSampleStatisticsContainer cvResults = mySpecification.computeOutOfSampleStatistics(numberTreesInForest, rngBaseSeedMomentForest, verbose, minObservationsPerLeaf, minImprovement, maxDepth, rngBaseSeedOutOfSample);
                        double combinationMSE = cvResults.getOutOfSampleMeasureY();
                        String star = "";
                        if (combinationMSE <= lowestSSE || first) {
                            lowestSSE = combinationMSE;
                            first = false;
                            bestMinImprovement = minImprovement;
                            bestMinObservationsPerLeaf = minObservationsPerLeaf;
                            bestMaxDepth = maxDepth;
                            star = "(*)";
                        }
                        System.out.println("detectHomogeneity: " + detectHomogeneity + " minMSE: " + minImprovement + " minObs: " + minObservationsPerLeaf + " maxDepth: " + maxDepth + " Out-of-sample MSE: " + combinationMSE + " " + star);
                    }
                }
            }

            System.out.println("Lowest MSE: " + lowestSSE + " at min_N = " + bestMinObservationsPerLeaf + " min_MSE = " + bestMinImprovement + " maxDepth: " + bestMaxDepth);
        } else {
            bestMinObservationsPerLeaf = 200;
            bestMinImprovement = 0.1;
            bestMaxDepth = 1;
        }

        numberTreesInForest = 1;
        mySpecification.resetHomogeneityIndex();
        
        if (detectHomogeneity) {
            if (verbose) {
                System.out.println("************************");
                System.out.println("* Test for Homogeneity *");
                System.out.println("************************");
            }
            /**
             * Step 2: determine homogeneous parameters post-CV
             */
            double minProportionInEachLeaf = 0.01;

            DataLens forestLens = new DataLens(mySpecification.getX(), mySpecification.getY(), mySpecification.getZ(), null);

            testParameterHomogeneity = true;
            TreeOptions cvOptions = new TreeOptions(minProportionInEachLeaf, bestMinObservationsPerLeaf, bestMinImprovement, bestMaxDepth, testParameterHomogeneity); // k = 1
            MomentForest myForest = new MomentForest(mySpecification, numberTreesInForest, rngBaseSeedMomentForest, forestLens, verbose, new TreeOptions());

            myForest.setTreeOptions(cvOptions);
            myForest.growForest();

            ArrayList<Integer> hpl = myForest.applyHomogeneityVotes(true);
            setHomogeneousParameterList(hpl);
            // System.out.println("After setting hp list");

            /**
             * Estimate values of those homogeneous parameters
             */
            if (!hpl.isEmpty()) {
                System.out.println("Initializing search container");
                // numberTreesInForest = 1;
                HomogeneousSearchContainer con = new HomogeneousSearchContainer(mySpecification, numberTreesInForest, verbose, bestMinImprovement, bestMinObservationsPerLeaf, bestMaxDepth,
                        getHomogeneousParameterList(), rngBaseSeedMomentForest, rngBaseSeedOutOfSample);

                System.out.println("Calling execute search");
                con.executeSearch();
                System.out.println("Post search");
                Jama.Matrix homogeneousParameters = con.getEstimatedHomogeneousParameters();
                System.out.print("Post-HomogeneousSearchContainer Estimated homogeneous parameters: ");
                pmUtility.prettyPrintVector(homogeneousParameters); // this is a compact vector of parameters

                int K = mySpecification.getHomogeneousParameterVector().getRowDimension();
                Jama.Matrix expandedHomogeneousParameterVector = new Jama.Matrix(K, 1);
                int counter = 0;
                for (int k = 0; k < K; k++) {
                    if (mySpecification.getHomogeneousIndex()[k]) {
                        expandedHomogeneousParameterVector.set(k, 0, homogeneousParameters.get(counter, 0));
                        counter++;
                    }
                }
                setEstimatedHomogeneousParameters(expandedHomogeneousParameterVector);
            }
        }

        /**
         * Compute out-of-sample measures of fit (against Y, and true beta)
         */
        // numberTreesInForest = 1;
        System.out.println("Building final tree");
        // OutOfSampleStatisticsContainer results = mySpecification.computeOutOfSampleStatistics(numberTreesInForest, rngBaseSeedMomentForest, verbose, bestMinObservationsPerLeaf,                bestMinImprovement, bestMaxDepth, rngBaseSeedOutOfSample);

        testParameterHomogeneity = false;
        DataLens forestLens = new DataLens(mySpecification.getX(), mySpecification.getY(), mySpecification.getZ(), null);
        TreeOptions cvOptions = new TreeOptions(0.01, bestMinObservationsPerLeaf, bestMinImprovement, bestMaxDepth, testParameterHomogeneity); // k = 1
        MomentForest myForest = new MomentForest(mySpecification, numberTreesInForest, rngBaseSeedMomentForest, forestLens, verbose, new TreeOptions());

        myForest.setTreeOptions(cvOptions);
        myForest.growForest();
        
        myForest.getTree(0).printTree();

    }

    public double getEstimatedBetaVersusTruthMSE() {
        return estimatedBetaVersusTruthMSE;
    }

    public void setEstimatedBetaVersusTruthMSE(double estimatedBetaVersusTruthMSE) {
        this.estimatedBetaVersusTruthMSE = estimatedBetaVersusTruthMSE;
    }

    private void setHomogeneousParameterList(ArrayList<Integer> homogeneousParameterList) {
        this.homogeneousParameterList = homogeneousParameterList;
    }

    public ArrayList<Integer> getHomogeneousParameterList() {
        return homogeneousParameterList;
    }

    /**
     * @return the MSE
     */
    public double getOutOfSampleYMSE() {
        // if (outOfSampleYMSE > 3) {
        // jt.append("MSPE: " + outOfSampleYMSE + " seed: " + rngSeed + "\n");
        // }
        return outOfSampleYMSE;
    }

    /**
     * @param outOfSampleYMSE the MSE to set
     */
    public void setOutOfSampleYMSE(double outOfSampleYMSE) {
        this.outOfSampleYMSE = outOfSampleYMSE;
    }

    /**
     * @return the estimateHomogeneousParameter
     */
    public Jama.Matrix getEstimatedHomogeneousParameters() {
        return estimatedHomogeneousParameters;
    }

    public void setEstimatedHomogeneousParameters(Matrix estimatedHomogeneousParameters) {
        this.estimatedHomogeneousParameters = estimatedHomogeneousParameters;
    }

}
