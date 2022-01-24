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

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        JFrame f =new JFrame("Monte Carlo");
        f.setBounds(300,300,1500,500);
        f.getContentPane().setLayout(new BorderLayout());
        JTextAreaAutoscroll jt = new JTextAreaAutoscroll();
        f.getContentPane().add(jt, BorderLayout.CENTER);
        f.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        f.setVisible(true);
        
        for (int numObs = 250; numObs < 500 * 2 * 2 * 2 * 2 * 2 * 2; numObs *= 2) {
            Random rng = new Random(22);
            
            double[] homogeneousClassificationRate = new double[2];
            int numMonteCarlos = 500;

            ArrayList<LinearTestMain> parallelLTM = new ArrayList<>();

            for (int m = 0; m < numMonteCarlos; m++) {
                LinearTestMain go = new LinearTestMain(rng.nextLong(), numObs);
                parallelLTM.add(go);
            }
            parallelLTM.parallelStream().forEach(e -> {
                e.execute();
            });

            for (int m = 0; m < numMonteCarlos; m++) {
                ArrayList<Integer> hList = parallelLTM.get(m).getHomogeneousParameterList();
                for (Integer h : hList) {
                    homogeneousClassificationRate[h] = homogeneousClassificationRate[h] + 1.0;
                }
            }
            jt.append("---------------------------------------------------------\n");
            jt.append("n = " + numObs + " Classification Rate by Parameter\n");
            for (int i = 0; i < homogeneousClassificationRate.length; i++) {
                jt.append(i + ". " + (homogeneousClassificationRate[i] / (0.0 + numMonteCarlos))+"\n");
            }
            jt.append("---------------------------------------------------------\n");
        }
    }

    public LinearTestMain(long rngSeed, int numObs) {
        this.rngSeed = rngSeed;
        this.numObs = numObs;
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
                    MomentForest myForest = new MomentForest(mySpecification, numberTreesInForest, rngBaseSeedMomentForest, forestLens, verbose, new TreeOptions());
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
                    DataLens oosDataLens = mySpecification.getOutOfSampleXYZ(100, rngBaseSeedOutOfSample); // this should eventually be modified to come out of the data itself (or generalized in some way)
                    Jama.Matrix testZ = oosDataLens.getZ();
                    Jama.Matrix residualizedX = mySpecification.residualizeX(oosDataLens.getX());
                    // pmUtility.prettyPrint(pmUtility.concatMatrix(pmUtility.concatMatrix(oosDataLens.getY(), oosDataLens.getX()), testZ));

                    /**
                     * Compute out-of-sample fit at current homogeneous
                     * parameter vector
                     */
                    double outOfSampleFit = 0;
                    for (int i = 0; i < testZ.getRowDimension(); i++) {
                        Jama.Matrix zi = testZ.getMatrix(i, i, 0, testZ.getColumnDimension() - 1);
                        Jama.Matrix b = myForest.getEstimatedParameterForest(zi);

                        Jama.Matrix residualizedXi = residualizedX.getMatrix(i, i, 0, residualizedX.getColumnDimension() - 1);
                        Jama.Matrix fullXi = oosDataLens.getX().getMatrix(i, i, 0, oosDataLens.getX().getColumnDimension() - 1);

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

            System.out.println("Lowest MSE: " + lowestSSE + " at min_N = " + bestMinObservationsPerLeaf + " min_MSE = " + bestMinImprovement + " bestAlpha: " + bestAlpha);
        }

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
        ArrayList<Integer> hpl = myForest.getTree(0).getIndexHomogeneousParameters(); // this is only using the first tree, is that the right way of thinking about this?
        int numParams = hpl.size();
        mySpecification.resetHomogeneityIndex();
        for (Integer i : hpl) {
            mySpecification.setHomogeneousIndex(i);
        }
        setHomogeneousParameterList(hpl);

        boolean estimateConstrainedModel = false;
        if (estimateConstrainedModel) {
            /**
             * I wonder if there should be an additional step here to recompute
             * the optimal CV parameters given the results of the homogeneity
             * test in the prior step?
             */
            if (verbose) {
                System.out.println("******************************");
                System.out.println("* Estimate Constrained Model *");
                System.out.println("******************************");
            }
            /**
             * Step 3: estimate constrained model
             */
            maxTreeDepth = 1;
            HomogeneousSearchContainer con = new HomogeneousSearchContainer(mySpecification, numberTreesInForest, verbose, bestMinImprovement, bestMinObservationsPerLeaf, maxTreeDepth, numParams, rngBaseSeedMomentForest, rngBaseSeedOutOfSample);
            con.executeSearch();
            Jama.Matrix homogeneousParameters = con.getEstimatedHomogeneousParameters();
            System.out.println("Estimated homogeneous parameters: ");
            pmUtility.prettyPrintVector(homogeneousParameters);

            boolean computeSE = false;
            if (computeSE) {
                int numberBootstraps = 50;
                // System.out.println("Number of bootstraps: " + numberBootstraps);

                int numberTreesInBootForest = 10;

                System.out.println("Bootstrap forest not working; see below");
                System.exit(0);
                // have to figure out how the random number generator fits into this
                BootstrapForest boot = new BootstrapForest(mySpecification, numberBootstraps, numberTreesInBootForest, 787, cvOptions);
                /**
                 * TODO: have to make sure that the homogeneous parameter thing
                 * is working inside the boot forest ALSO, how to compute
                 * standard errors on homogeneous parameters here COME BACK TO
                 * THIS!!!
                 */

                DataLens oosDataLens = mySpecification.getOutOfSampleXYZ(100, rngBaseSeedOutOfSample);
                Jama.Matrix testZ = oosDataLens.getZ();

                for (int j = 0; j < testZ.getRowDimension(); j++) {
                    Jama.Matrix zi = testZ.getMatrix(j, j, 0, testZ.getColumnDimension() - 1);
                    Jama.Matrix b = myForest.getEstimatedParameterForest(zi);
                    System.out.println("z: " + pmUtility.stringPrettyPrint(zi) + " beta: " + pmUtility.stringPrettyPrintVector(b) + " se: " + pmUtility.stringPrettyPrintVector(boot.computeStandardErrors(zi)));
                }
            }
        }
    }

    private void setHomogeneousParameterList(ArrayList<Integer> homogeneousParameterList) {
        this.homogeneousParameterList = homogeneousParameterList;
    }

    public ArrayList<Integer> getHomogeneousParameterList() {
        return homogeneousParameterList;
    }

}
