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

import Jama.Matrix;
import core.DataLens;
import core.HomogeneousSearchContainer;
import core.MomentForest;
import core.MomentSpecification;
import core.TreeOptions;
import java.awt.BorderLayout;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Random;
import javax.swing.JFrame;
import javax.swing.JScrollPane;
import javax.swing.JTextArea;
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
    private Jama.Matrix estimatedHomogeneousParameters;
    private final boolean detectHomogeneity;
    JTextArea jt;

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        JFrame f = new JFrame("Monte Carlo");
        f.setBounds(300, 300, 1500, 500);
        f.getContentPane().setLayout(new BorderLayout());
        JTextArea jt = new JTextAreaAutoscroll();
        f.getContentPane().add(new JScrollPane(jt), BorderLayout.CENTER);
        f.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        f.setVisible(true);

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
        boolean[] d = {false, true};
        // boolean[] d = {!true};
        for (boolean detectHomogeneity : d) {
            // boolean detectHomogeneity = !true;
            if (detectHomogeneity) {
                jt.append("***** ESTIMATING HOMOGENEOUS PARAMETERS *****\n");
            } else {
                jt.append("***** UNRESTRICTED MODEL *****\n");
            }

            for (int numObs = 500; numObs <= 4000; numObs *= 2) {
                Random rng = new Random(22);

                int numParameters = 2;

                double[] homogeneousClassificationRate = new double[numParameters];
                double[] avgParameter = new double[numParameters];
                double Y_MSE = 0;
                double Y_MSE_var = 0;
                double beta_MSE = 0;
                double beta_MSE_var = 0;

                int numMonteCarlos = 100;

                ArrayList<LinearTestMain> parallelLTM = new ArrayList<>();

                for (int m = 0; m < numMonteCarlos; m++) {
                    LinearTestMain go = new LinearTestMain(rng.nextLong(), numObs, detectHomogeneity, jt);
                    parallelLTM.add(go);
                }
                parallelLTM.parallelStream().forEach(e -> {
                    //    parallelLTM.stream().forEach(e -> {
                    e.execute();
                });

                // for (int m = 0; m < numMonteCarlos; m++) {
                int[] counts = new int[numParameters];
                double[] standardErrors = new double[numParameters];
                for (LinearTestMain m : parallelLTM) {
                    if (detectHomogeneity) {
                        ArrayList<Integer> hList = m.getHomogeneousParameterList();
                        if (hList.size() > 0) {
                            for (Integer h : hList) {
                                // jt.append("Detected homogeneity on parameter index " + h + " hList.size = " + hList.size());
                                // jt.append(" estimatedParemeterSize: " + m.getEstimatedHomogeneousParameters().getRowDimension() + "\n");
                                homogeneousClassificationRate[h] = homogeneousClassificationRate[h] + 1.0;
                                for (int i = 0; i < numParameters; i++) {
                                    if (h == i) {
                                        counts[i]++; // this is counting how many times this parameter is chosen as homogeneous
                                        avgParameter[i] += m.getEstimatedHomogeneousParameters().get(i, 0);
                                        // jt.append("Estimated Homogeneous Parameter: "+pmUtility.stringPrettyPrintVector(m.getEstimatedHomogeneousParameters())+"\n");
                                        // jt.append("Average parameter sum for "+i+" is: "+avgParameter[i]+" Running count: "+counts[i]+"\n");
                                    }
                                }
                            }
                        }
                    }
                    Y_MSE += m.getOutOfSampleYMSE();
                    beta_MSE += m.getEstimatedBetaVersusTruthMSE();
                }

                if (detectHomogeneity) {
                    for (int i = 0; i < numParameters; i++) {
                        if (counts[i] > 0) {
                            avgParameter[i] = avgParameter[i] / counts[i];
                        }
                    }
                }

                // needed average to compute standard errors (kind of clunky, but whatever)
                for (LinearTestMain m : parallelLTM) {
                    if (detectHomogeneity) {
                        ArrayList<Integer> hList = m.getHomogeneousParameterList();
                        if (hList.size() > 0) {
                            for (Integer h : hList) {
                                for (int i = 0; i < numParameters; i++) {
                                    if (h == i) {
                                        standardErrors[i] += Math.pow(m.getEstimatedHomogeneousParameters().get(i, 0) - avgParameter[i], 2);
                                    }
                                }
                            }
                        }
                    }
                }
                if (detectHomogeneity) {
                    for (int i = 0; i < numParameters; i++) {
                        if (counts[i] > 0) {
                            standardErrors[i] = Math.sqrt(standardErrors[i] / counts[i]);
                        }
                    }
                }

                Y_MSE /= parallelLTM.size();
                beta_MSE /= parallelLTM.size();
                for (LinearTestMain m : parallelLTM) {
                    Y_MSE_var += Math.pow(Y_MSE - m.getOutOfSampleYMSE(), 2);
                    beta_MSE_var += Math.pow(beta_MSE - m.getEstimatedBetaVersusTruthMSE(), 2);
                }
                jt.append("---------------------------------------------------------\n");
                jt.append("n = " + numObs + " Classification Rate by Parameter\n");
                if (detectHomogeneity) {
                    for (int i = 0; i < homogeneousClassificationRate.length; i++) {
                        jt.append(i + ". " + (homogeneousClassificationRate[i] / (0.0 + numMonteCarlos)) + "\n");
                    }
                }
                jt.append("Y_MSE: " + Y_MSE + " (" + Y_MSE_var + ")\n");
                jt.append("beta_MSE: " + beta_MSE + " (" + beta_MSE_var + ")\n");
                if (detectHomogeneity) {
                    for (int i = 0; i < numParameters; i++) {
                        String s = String.format("x[%d]: mean: %g se: %g [%d]\n", i, avgParameter[i], standardErrors[i], counts[i]);
                        // jt.append("Parameter " + i + " mean: " + avgParameter[i] + " se: " + standardErrors[i] + " [" + counts[i] + "]\n");
                        jt.append(s);
                    }
                }
                jt.append("---------------------------------------------------------\n");
            }
        }
    }

    public LinearTestMain(long rngSeed, int numObs, boolean detectHomogeneity, JTextArea jt) {
        this.rngSeed = rngSeed;
        this.jt = jt;
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
        /* Contains X data, Y data, balancing vector (treatment indicators), and data index (just an array numbered 0 - numObs) */
        boolean verbose = false;
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

            for (double minImprovement = 1.0; minImprovement <= 1.0; minImprovement *= 10) {
                for (int minObservationsPerLeaf = 50; minObservationsPerLeaf <= 50; minObservationsPerLeaf *= 2) {
                    double combinationMSE = computeOutOfSampleMSE(mySpecification, numberTreesInForest, rngBaseSeedMomentForest, verbose, minObservationsPerLeaf, minImprovement, maxTreeDepth, rngBaseSeedOutOfSample);
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

            DataLens forestLens = new DataLens(mySpecification.getX(), mySpecification.getY(true), mySpecification.getZ(), null);

            testParameterHomogeneity = true;
            TreeOptions cvOptions = new TreeOptions(minProportionInEachLeaf, bestMinObservationsPerLeaf, bestMinImprovement, maxTreeDepth, testParameterHomogeneity); // k = 1
            MomentForest myForest = new MomentForest(mySpecification, numberTreesInForest, rngBaseSeedMomentForest, forestLens, verbose, new TreeOptions());

            myForest.setTreeOptions(cvOptions);
            myForest.growForest();
            // System.out.println("Done with growforest");
            ArrayList<Integer> hpl = myForest.getTree(0).getIndexHomogeneousParameters(); // this is only using the first tree, is that the right way of thinking about this?
            // System.out.println("Post get homogeneous parameters");

            // tell the specification that these parameters have been determined to be homogeneous
            Collections.sort(hpl); // ensure that indices are ascending (this can cause some weird problems elsewhere due to my bad coding skills if not)
            for (Integer i : hpl) {
                mySpecification.setHomogeneousIndex(i);
                // System.out.println(i);
            }
            setHomogeneousParameterList(hpl);
            // System.out.println("After setting hp list");

            /**
             * Estimate values of those homogeneous parameters
             */
            maxTreeDepth = 1;
            boolean cheatToVerifyWorking = false;
            if (cheatToVerifyWorking) {
                // this is here to verify that the code is working in that we should get a lower OOS MSE when the truth is imposed (it works)
                hpl.clear();
                hpl.add(0);
                hpl.add(1);
                mySpecification.resetHomogeneityIndex();
                mySpecification.setHomogeneousIndex(0);
                mySpecification.setHomogeneousParameter(0, -1.0);
                mySpecification.setHomogeneousIndex(1);
                mySpecification.setHomogeneousParameter(1, 1.0);
                setEstimatedHomogeneousParameters(mySpecification.getHomogeneousParameterVector());

            } else {
                if (hpl.size() > 0) {
                    // System.out.println("Initializing search container");
                    HomogeneousSearchContainer con = new HomogeneousSearchContainer(mySpecification, numberTreesInForest, verbose, bestMinImprovement, bestMinObservationsPerLeaf, maxTreeDepth,
                            getHomogeneousParameterList(), rngBaseSeedMomentForest, rngBaseSeedOutOfSample);
                    // System.out.println("Calling execute search");
                    con.executeSearch();
                    // System.out.println("Post search");
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
        }

        /**
         * Compute out-of-sample measures of fit (against Y, and true beta)
         */
        outOfSampleYMSE = computeOutOfSampleMSE(mySpecification, numberTreesInForest, rngBaseSeedMomentForest, verbose, bestMinObservationsPerLeaf,
                bestMinImprovement, maxTreeDepth, rngBaseSeedOutOfSample);
        setEstimatedBetaVersusTruthMSE(computeOutOfSampleMSEInParameterSpace(mySpecification, numberTreesInForest, rngBaseSeedMomentForest, verbose, bestMinObservationsPerLeaf,
                bestMinImprovement, maxTreeDepth, rngBaseSeedOutOfSample));
    }

    public double getEstimatedBetaVersusTruthMSE() {
        return estimatedBetaVersusTruthMSE;
    }

    public void setEstimatedBetaVersusTruthMSE(double estimatedBetaVersusTruthMSE) {
        this.estimatedBetaVersusTruthMSE = estimatedBetaVersusTruthMSE;
    }

    private double computeOutOfSampleMSE(MomentSpecification mySpecification, int numberTreesInForest, long rngBaseSeedMomentForest, boolean verbose,
            int minObservationsPerLeaf, double minImprovement, int maxTreeDepth, long rngBaseSeedOutOfSample) {

        boolean allParametersHomogeneous = true;
        for (boolean b : mySpecification.getHomogeneousIndex()) {
            if (!b) {
                allParametersHomogeneous = false;
            }
        }

        MomentForest myForest = null;
        if (!allParametersHomogeneous) {
            DataLens homogenizedForestLens = new DataLens(mySpecification.getX(), mySpecification.getY(true), mySpecification.getZ(), null);

            // System.out.println("\nComputing OOS MSE\n");
            // System.out.println("Homogeneous parameter length in spec: "+mySpecification.getHomogeneousIndex().length);
            boolean testParameterHomogeneity = false;
            myForest = new MomentForest(mySpecification, numberTreesInForest, rngBaseSeedMomentForest, homogenizedForestLens, verbose, new TreeOptions());
            TreeOptions cvOptions = new TreeOptions(0.01, minObservationsPerLeaf, minImprovement, maxTreeDepth, testParameterHomogeneity); // k = 1
            myForest.setTreeOptions(cvOptions);
            /**
             * Grow the moment forest
             */
            myForest.growForest();
        }

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

            Jama.Matrix fullXi = oosDataLens.getX().getMatrix(i, i, 0, oosDataLens.getX().getColumnDimension() - 1);
            Jama.Matrix zi = testZ.getMatrix(i, i, 0, testZ.getColumnDimension() - 1);

            // note, xi is part of x that is not homogeneous
            // for the homogeneous component below, need to find the other parts of x not contained in xi
            // not only that, the way that that method is specified is that it takes the ENTIRE xi to get the homogeneous part, where I have pulled out the subset of X already in xi
            // so how to get the whole row?
            // i am going to do two things here: one, remove residualization from getOutOfSampleXYZ
            // added a new method to residualize an X matrix
            double fitY = mySpecification.getHomogeneousComponent(fullXi);
            if (!allParametersHomogeneous) {
                Jama.Matrix b = myForest.getEstimatedParameterForest(zi);
                Jama.Matrix residualizedXi = residualizedX.getMatrix(i, i, 0, residualizedX.getColumnDimension() - 1);
                fitY += residualizedXi.times(b).get(0, 0);
            }
            double error = fitY - (oosDataLens.getY().get(i, 0));

            outOfSampleFit += error * error;

            boolean outputFits = false;
            if (outputFits) {
                Jama.Matrix bTruth = mySpecification.getBetaTruth(zi);
                // have to reconstruct a composite beta from homogeneous and heterogeneous parameters
                int heterogeneousCounter = 0;
                Jama.Matrix compositeEstimatedBeta = new Jama.Matrix(bTruth.getRowDimension(), 1);
                for (int k = 0; k < bTruth.getRowDimension(); k++) {
                    if (mySpecification.getHomogeneousIndex()[k]) {
                        compositeEstimatedBeta.set(k, 0, mySpecification.getHomogeneousParameter(k));
                    } else {
                        Jama.Matrix b = myForest.getEstimatedParameterForest(zi);
                        compositeEstimatedBeta.set(k, 0, b.get(heterogeneousCounter, 0));
                        heterogeneousCounter++;
                    }
                }

                System.out.print("z: " + pmUtility.stringPrettyPrint(zi) + " beta: " + pmUtility.stringPrettyPrintVector(compositeEstimatedBeta) + " trueBeta: " + pmUtility.stringPrettyPrintVector(bTruth));
                System.out.print(" x: " + pmUtility.stringPrettyPrint(fullXi));
                if (!allParametersHomogeneous) {
                    Jama.Matrix residualizedXi = residualizedX.getMatrix(i, i, 0, residualizedX.getColumnDimension() - 1);
                    Jama.Matrix b = myForest.getEstimatedParameterForest(zi);
                    System.out.print(" residualizedX'b: " + residualizedXi.times(b).get(0, 0));
                }
                System.out.print(" hc: " + mySpecification.getHomogeneousComponent(fullXi));
                System.out.println(" fitY: " + fitY + " Y: " + oosDataLens.getY().get(i, 0) + " SE: " + error * error);
            }
        }
        return outOfSampleFit / testZ.getRowDimension();

    }

    private double computeOutOfSampleMSEInParameterSpace(MomentSpecification mySpecification, int numberTreesInForest, long rngBaseSeedMomentForest, boolean verbose,
            int minObservationsPerLeaf, double minImprovement, int maxTreeDepth, long rngBaseSeedOutOfSample) {
        // System.out.println("\nComputing OOS In Parameter Space\n");
        // System.out.println("Homogeneous parameter length in spec: "+mySpecification.getHomogeneousIndex().length);

        boolean allParametersHomogeneous = true;
        for (boolean b : mySpecification.getHomogeneousIndex()) {
            if (!b) {
                allParametersHomogeneous = false;
            }
        }

        MomentForest myForest = null;
        if (!allParametersHomogeneous) {
            DataLens homogenizedForestLens = new DataLens(mySpecification.getX(), mySpecification.getY(true), mySpecification.getZ(), null);

            myForest = new MomentForest(mySpecification, numberTreesInForest, rngBaseSeedMomentForest, homogenizedForestLens, verbose, new TreeOptions());
            TreeOptions cvOptions = new TreeOptions(0.01, minObservationsPerLeaf, minImprovement, maxTreeDepth, false); // k = 1
            myForest.setTreeOptions(cvOptions);
            /**
             * Grow the moment forest
             */
            myForest.growForest();

            myForest.getTree(0).printTree();
        }

        /**
         * Test vectors for assessment
         */
        DataLens oosDataLens = mySpecification.getOutOfSampleXYZ(2000, rngBaseSeedOutOfSample); // this should eventually be modified to come out of the data itself (or generalized in some way)
        Jama.Matrix testZ = oosDataLens.getZ();

        /**
         * Compute out-of-sample fit at current homogeneous parameter vector
         */
        double outOfSampleFit = 0;
        for (int i = 0; i < testZ.getRowDimension(); i++) {
            Jama.Matrix zi = testZ.getMatrix(i, i, 0, testZ.getColumnDimension() - 1);

            // going to compare directly to the true parameter vector in this method instead of using fit of Y
            Jama.Matrix bTruth = mySpecification.getBetaTruth(zi);

            // have to reconstruct a composite beta from homogeneous and heterogeneous parameters
            int heterogeneousCounter = 0;
            Jama.Matrix compositeEstimatedBeta = new Jama.Matrix(bTruth.getRowDimension(), 1);

            for (int k = 0; k < bTruth.getRowDimension(); k++) {
                if (mySpecification.getHomogeneousIndex()[k]) {
                    compositeEstimatedBeta.set(k, 0, mySpecification.getHomogeneousParameter(k));
                } else {
                    Jama.Matrix b = myForest.getEstimatedParameterForest(zi);
                    compositeEstimatedBeta.set(k, 0, b.get(heterogeneousCounter, 0));
                    heterogeneousCounter++;
                }
            }
            if (i == 0) {
                String hString = "[ ";
                for (int k = 0; k < bTruth.getRowDimension(); k++) {
                    if (mySpecification.getHomogeneousIndex()[k]) {
                        hString = hString + "X ";
                    } else {
                        hString = hString + "O ";
                    }
                }
                hString = hString + "]";
                // jt.append("Composite estimated beta: " + pmUtility.stringPrettyPrintVector(compositeEstimatedBeta) + " " + hString + "\n");
            }
            //pmUtility.prettyPrintVector(compositeEstimatedBeta);

            outOfSampleFit += (compositeEstimatedBeta.minus(bTruth)).norm2();
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
    public Jama.Matrix getEstimatedHomogeneousParameters() {
        return estimatedHomogeneousParameters;
    }

    public void setEstimatedHomogeneousParameters(Matrix estimatedHomogeneousParameters) {
        this.estimatedHomogeneousParameters = estimatedHomogeneousParameters;
    }

}
