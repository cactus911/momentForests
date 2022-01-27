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
import java.awt.BorderLayout;
import java.util.ArrayList;
import java.util.Random;
import javax.swing.JFrame;
import utility.JTextAreaAutoscroll;
import utility.pmUtility;

/**
 *
 * @author Stephen P. Ryan <stephen.p.ryan@wustl.edu>
 */
public class LinearTestMain {

    private ArrayList<Integer> homogeneousParameterList;
    private final long rngSeed;
    private final int numObs;
    private double outOfSampleYMSE;
    private double estimatedBetaVersusTruthMSE;
    private double estimateHomogeneousParameter;
    private final boolean detectHomogeneity;

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
//        JFrame f = new JFrame("Monte Carlo");
//        f.setBounds(300, 300, 1500, 500);
//        f.getContentPane().setLayout(new BorderLayout());
        JTextAreaAutoscroll jt = new JTextAreaAutoscroll();
//        f.getContentPane().add(jt, BorderLayout.CENTER);
//        f.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
//        f.setVisible(true);

        /**
         * It appears that we get classification basically correct What do we
         * want to put into the Monte Carlo?
         *
         * Different levels of observations (n) With and without parameter
         * homogeneity testing Show how often we get it right with it on MSE in
         * out-of-sample prediction across two methods The parameter value when
         * it is homogeneous?
         *
         * Should have some focus on parameter values, since the MSE converges
         * to idiosyncratic levels (Which is cool but not what we are directly
         * interested in)
         *
         * How to implement that? We know the truth here, so query the
         * specification for its guess of the parameter vector for a bunch of
         * X,Z combinations, run l2-norm on that? Done that, seems to be working
         * really nicely.
         */
        boolean detectHomogeneity = true;
        for (int numObs = 1000; numObs <= 4000; numObs *= 2) {
            Random rng = new Random(22);

            double[] homogeneousClassificationRate = new double[2];
            double Y_MSE = 0;
            double Y_MSE_var = 0;
            double beta_MSE = 0;
            double beta_MSE_var = 0;
            double meanHomogeneousParameter = 0; // this is going to be a little weird since this may be a vector, and that vector may change size across runs!
            double varHomogeneousParameter = 0;

            int numMonteCarlos = 16;

            ArrayList<LinearTestMain> parallelLTM = new ArrayList<>();

            for (int m = 0; m < numMonteCarlos; m++) {
                LinearTestMain go = new LinearTestMain(rng.nextLong(), numObs, detectHomogeneity);
                parallelLTM.add(go);
            }
            // parallelLTM.parallelStream().forEach(e -> {
            parallelLTM.stream().forEach(e -> {
                e.execute();
            });

            // for (int m = 0; m < numMonteCarlos; m++) {
            for (LinearTestMain m : parallelLTM) {
                if (detectHomogeneity) {
                    ArrayList<Integer> hList = m.getHomogeneousParameterList();
                    for (Integer h : hList) {
                        homogeneousClassificationRate[h] = homogeneousClassificationRate[h] + 1.0;
                    }
                    // we'll figure out mean homogeneous parameter here when i get back (how to to deal with it sometimes not showing up?)
                }
                Y_MSE += m.getOutOfSampleYMSE();
                beta_MSE += m.getEstimatedBetaVersusTruthMSE();
            }
            Y_MSE /= parallelLTM.size();
            beta_MSE /= parallelLTM.size();
            for (LinearTestMain m : parallelLTM) {
                Y_MSE_var += Math.pow(Y_MSE - m.getOutOfSampleYMSE(), 2);
                beta_MSE_var += Math.pow(beta_MSE - m.getEstimatedBetaVersusTruthMSE(), 2);
            }
            jt.append("---------------------------------------------------------\n");
            jt.append("n = " + numObs + " Classification Rate by Parameter\n");
            for (int i = 0; i < homogeneousClassificationRate.length; i++) {
                jt.append(i + ". " + (homogeneousClassificationRate[i] / (0.0 + numMonteCarlos)) + "\n");
            }
            jt.append("Y_MSE: " + Y_MSE + " (" + Y_MSE_var + ")\n");
            jt.append("beta_MSE: " + beta_MSE + " (" + beta_MSE_var + ")\n");
            if (detectHomogeneity) {
                jt.append("mean homogeneous parameter: " + meanHomogeneousParameter + " (" + varHomogeneousParameter + ")\n");
            }
            jt.append("---------------------------------------------------------\n");
        }
    }

    public LinearTestMain(long rngSeed, int numObs, boolean detectHomogeneity) {
        this.rngSeed = rngSeed;
        this.numObs = numObs;
        this.detectHomogeneity = detectHomogeneity;
    }

    private void execute() {
        Random rng = new Random(rngSeed);

        // MomentSpecification mySpecification = new LinearMomentSpecification("data/airline_subset.csv");
        MomentSpecification mySpecification = new LinearMomentSpecification(numObs);
        mySpecification.loadData(rng.nextLong()); // Create data using rng

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
        boolean verbose = false;
        boolean testParameterHomogeneity = false;

        long rngBaseSeedMomentForest = rng.nextLong();
        long rngBaseSeedOutOfSample = rng.nextLong();

        boolean runCV = false;
        if (runCV) {
            if (verbose) {
                System.out.println("************************");
                System.out.println("* Run Cross-Validation *");
                System.out.println("************************");
            }

            for (double minImprovement = 1.0; minImprovement <= 1.0; minImprovement *= 10) {
                for (int minObservationsPerLeaf = 50; minObservationsPerLeaf <= 50; minObservationsPerLeaf *= 2) {
                    double combinationMSE = computeOutOfSampleMSE(mySpecification, numberTreesInForest, rngBaseSeedMomentForest, forestLens, verbose, minObservationsPerLeaf, minImprovement, maxTreeDepth, rngBaseSeedOutOfSample);
                    System.out.println("Out-of-sample MSE: " + combinationMSE);
                    if (combinationMSE < lowestSSE || first) {
                        lowestSSE = combinationMSE;
                        first = false;
                        bestMinImprovement = minImprovement;
                        bestMinObservationsPerLeaf = minObservationsPerLeaf;
                    }
                }
            }

            System.out.println("Lowest MSE: " + lowestSSE + " at min_N = " + bestMinObservationsPerLeaf + " min_MSE = " + bestMinImprovement + " bestAlpha: " + bestAlpha);
        }

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
            /**
             * Impose a max depth of 1 to test at the first split
             */
            maxTreeDepth = 1;
            testParameterHomogeneity = true;
            TreeOptions cvOptions = new TreeOptions(minProportionInEachLeaf, bestMinObservationsPerLeaf, bestMinImprovement, maxTreeDepth, testParameterHomogeneity); // k = 1
            MomentForest myForest = new MomentForest(mySpecification, numberTreesInForest, rngBaseSeedMomentForest, forestLens, verbose, new TreeOptions());

            myForest.setTreeOptions(cvOptions);
            myForest.growForest();
            System.out.println("Done with growforest");
            ArrayList<Integer> hpl = myForest.getTree(0).getIndexHomogeneousParameters(); // this is only using the first tree, is that the right way of thinking about this?
            System.out.println("Post get homogeneous parameters");

            for (Integer i : hpl) {
                mySpecification.setHomogeneousIndex(i);
            }
            setHomogeneousParameterList(hpl);
            System.out.println("After setting hp list");

            /**
             * Estimate values of those homogeneous parameters
             */
            maxTreeDepth = 1;
            System.out.println("Initializing search container");
            HomogeneousSearchContainer con = new HomogeneousSearchContainer(mySpecification, numberTreesInForest, verbose, bestMinImprovement, bestMinObservationsPerLeaf, maxTreeDepth,
                    getHomogeneousParameterList().size(), rngBaseSeedMomentForest, rngBaseSeedOutOfSample);
            System.out.println("Calling execute search");
            con.executeSearch();
            System.out.println("Post search");
            Jama.Matrix homogeneousParameters = con.getEstimatedHomogeneousParameters();
            System.out.println("Estimated homogeneous parameters: ");
            pmUtility.prettyPrintVector(homogeneousParameters);
            // outOfSampleYMSE = con.computeOutOfSampleMSE();
        }

        /**
         * Compute out-of-sample measures of fit (against Y, and true beta)
         */
        outOfSampleYMSE = computeOutOfSampleMSE(mySpecification, numberTreesInForest, rngBaseSeedMomentForest, forestLens, verbose, bestMinObservationsPerLeaf,
                bestMinImprovement, maxTreeDepth, rngBaseSeedOutOfSample);
        setEstimatedBetaVersusTruthMSE(computeOutOfSampleMSEInParameterSpace(mySpecification, numberTreesInForest, rngBaseSeedMomentForest, forestLens, verbose, bestMinObservationsPerLeaf,
                bestMinImprovement, maxTreeDepth, testParameterHomogeneity, rngBaseSeedOutOfSample));
    }

    public double getEstimatedBetaVersusTruthMSE() {
        return estimatedBetaVersusTruthMSE;
    }

    public void setEstimatedBetaVersusTruthMSE(double estimatedBetaVersusTruthMSE) {
        this.estimatedBetaVersusTruthMSE = estimatedBetaVersusTruthMSE;
    }

    private double computeOutOfSampleMSE(MomentSpecification mySpecification, int numberTreesInForest, long rngBaseSeedMomentForest, DataLens forestLens, boolean verbose,
            int minObservationsPerLeaf, double minImprovement, int maxTreeDepth, long rngBaseSeedOutOfSample) {

        DataLens homogenizedForestLens = new DataLens(mySpecification.getX(), mySpecification.getY(true), mySpecification.getZ(), null);

        boolean testParameterHomogeneity = false;
        MomentForest myForest = new MomentForest(mySpecification, numberTreesInForest, rngBaseSeedMomentForest, homogenizedForestLens, verbose, new TreeOptions());
        TreeOptions cvOptions = new TreeOptions(0.01, minObservationsPerLeaf, minImprovement, maxTreeDepth, testParameterHomogeneity); // k = 1
        myForest.setTreeOptions(cvOptions);
        /**
         * Grow the moment forest
         */
        myForest.growForest();

        /**
         * Test vectors for assessment
         */
        DataLens oosDataLens = mySpecification.getOutOfSampleXYZ(2000, rngBaseSeedOutOfSample); // this should eventually be modified to come out of the data itself (or generalized in some way)
        Jama.Matrix testZ = oosDataLens.getZ();
        Jama.Matrix residualizedX = mySpecification.residualizeX(oosDataLens.getX());
        // pmUtility.prettyPrint(pmUtility.concatMatrix(pmUtility.concatMatrix(oosDataLens.getY(), oosDataLens.getX()), testZ));

        /**
         * Compute out-of-sample fit at current homogeneous parameter vector
         */
        double outOfSampleFit = 0;
        for (int i = 0; i < testZ.getRowDimension(); i++) {
            Jama.Matrix zi = testZ.getMatrix(i, i, 0, testZ.getColumnDimension() - 1);
            Jama.Matrix b = myForest.getEstimatedParameterForest(zi);

            Jama.Matrix residualizedXi = residualizedX.getMatrix(i, i, 0, residualizedX.getColumnDimension() - 1);
            Jama.Matrix fullXi = oosDataLens.getX().getMatrix(i, i, 0, oosDataLens.getX().getColumnDimension() - 1);

            // note, xi is part of x that is not homogeneous
            // for the homogeneous component below, need to find the other parts of x not contained in xi
            // not only that, the way that that method is specified is that it takes the ENTIRE xi to get the homogeneous part, where I have pulled out the subset of X already in xi
            // so how to get the whole row?
            // i am going to do two things here: one, remove residualization from getOutOfSampleXYZ
            // added a new method to residualize an X matrix
            double fitY = residualizedXi.times(b).get(0, 0) + mySpecification.getHomogeneousComponent(fullXi);
            double error = fitY - (oosDataLens.getY().get(i, 0));

            outOfSampleFit += error * error;

            boolean outputFits = false;
            if (outputFits) {
                System.out.print("z: " + pmUtility.stringPrettyPrint(zi) + " beta: " + pmUtility.stringPrettyPrintVector(b));
                System.out.print(" x: " + pmUtility.stringPrettyPrint(fullXi));
                System.out.print(" residualizedXb: " + residualizedXi.times(b).get(0, 0) + " hc: " + mySpecification.getHomogeneousComponent(fullXi));
                System.out.println(" fitY: " + fitY + " Y: " + oosDataLens.getY().get(i, 0) + " SE: " + error * error);
            }
        }
        return outOfSampleFit / testZ.getRowDimension();
    }

    private double computeOutOfSampleMSEInParameterSpace(MomentSpecification mySpecification, int numberTreesInForest, long rngBaseSeedMomentForest, DataLens forestLens, boolean verbose,
            int minObservationsPerLeaf, double minImprovement, int maxTreeDepth, boolean testParameterHomogeneity, long rngBaseSeedOutOfSample) {
        MomentForest myForest = new MomentForest(mySpecification, numberTreesInForest, rngBaseSeedMomentForest, forestLens, verbose, new TreeOptions());
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
        DataLens oosDataLens = mySpecification.getOutOfSampleXYZ(100, rngBaseSeedOutOfSample); // this should eventually be modified to come out of the data itself (or generalized in some way)
        Jama.Matrix testZ = oosDataLens.getZ();

        /**
         * Compute out-of-sample fit at current homogeneous parameter vector
         */
        double outOfSampleFit = 0;
        for (int i = 0; i < testZ.getRowDimension(); i++) {
            Jama.Matrix zi = testZ.getMatrix(i, i, 0, testZ.getColumnDimension() - 1);
            Jama.Matrix b = myForest.getEstimatedParameterForest(zi);

            // going to compare directly to the true parameter vector in this method instead of using fit of Y
            Jama.Matrix bTruth = mySpecification.getBetaTruth(zi);
            outOfSampleFit += (b.minus(bTruth)).norm2();
        }

        return outOfSampleFit / testZ.getRowDimension(); // mse
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
    public double getEstimateHomogeneousParameter() {
        return estimateHomogeneousParameter;
    }

}
