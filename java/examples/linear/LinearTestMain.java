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
import examples.PartialLinearGasDemand.HomogeneousParameterSorter;
import java.awt.GridLayout;
import java.util.ArrayList;
import java.util.Random;
import java.util.concurrent.atomic.AtomicInteger;
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
    private JTextArea jt;
    private final int dimX;

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        JFrame f = new JFrame("Monte Carlo");
        f.setBounds(100, 100, 1100, 500);
        f.getContentPane().setLayout(new GridLayout(1, 2));
        JTextAreaAutoscroll jt1 = new JTextAreaAutoscroll();
        f.getContentPane().add(new JScrollPane(jt1));
        JTextAreaAutoscroll jt2 = new JTextAreaAutoscroll();
        f.getContentPane().add(new JScrollPane(jt2));
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
        for (int dimX = 2; dimX < 10; dimX++) {
            for (int numObs = 1000; numObs <= 1000; numObs *= 2) {

                double YMSE_unrestricted = 0;
                double YMSE_SD_unrestricted = 0;
                double YMSE_restricted = 0;
                double YMSE_SD_restricted = 0;

                double betaMSE_unrestricted = 0;
                double betaMSE_restricted = 0;
                double betaMSE_SD_unrestricted = 0;
                double betaMSE_SD_restricted = 0;

                double beta1_mean = 0;
                double beta1_SD = 0;
                double beta2_mean = 0;
                double beta2_SD = 0;

                double classificationRate1 = 0;
                double classificationRate2 = 0;

                JTextAreaAutoscroll jam = new JTextAreaAutoscroll();

                boolean[] d = {false, true};
                // boolean[] d = {!true};
                for (boolean detectHomogeneity : d) {
                    // boolean detectHomogeneity = !true;
                    if (detectHomogeneity) {
                        jam = jt1;
                        // jam.append("***** ESTIMATING HOMOGENEOUS PARAMETERS *****\n");
                    } else {
                        jam = jt2;
                        // jam.append("***** UNRESTRICTED MODEL *****\n");
                    }
                    jam.append("\\\\ \n");

                    jam.append("Dimensionality of X: "+dimX+"\n");
                    
                    Random rng = new Random(22);

                    // int dimX = 3;
                    double[] homogeneousClassificationRate = new double[dimX];
                    double[] avgParameter = new double[dimX];
                    double Y_MSE = 0;
                    double Y_MSE_var = 0;
                    double beta_MSE = 0;
                    double beta_MSE_var = 0;

                    int numMonteCarlos = 48;

                    ArrayList<LinearTestMain> parallelLTM = new ArrayList<>();

                    for (int m = 0; m < numMonteCarlos; m++) {
                        LinearTestMain go;
                        if (numMonteCarlos == 1) {
                            // 8621193992485539638L
                            go = new LinearTestMain(8621193992485539638L, numObs, detectHomogeneity, jam, dimX);
                        } else {
                            go = new LinearTestMain(rng.nextLong(), numObs, detectHomogeneity, jam, dimX);
                        }
                        parallelLTM.add(go);
                    }

                    AtomicInteger bomb = new AtomicInteger();

                    parallelLTM.parallelStream().forEach(e -> {
                        //    parallelLTM.stream().forEach(e -> {
                        e.execute();
                        bomb.incrementAndGet();
                        System.out.println("Finished " + bomb.get() + " iterations.");
                    });

                    // for (int m = 0; m < numMonteCarlos; m++) {
                    int[] counts = new int[dimX];
                    double[] standardErrors = new double[dimX];
                    for (LinearTestMain m : parallelLTM) {
                        if (detectHomogeneity) {
                            ArrayList<Integer> hList = m.getHomogeneousParameterList();
                            if (!hList.isEmpty()) {
                                for (Integer h : hList) {
                                    // jt.append("Detected homogeneity on parameter index " + h + " hList.size = " + hList.size());
                                    // jt.append(" estimatedParemeterSize: " + m.getEstimatedHomogeneousParameters().getRowDimension() + "\n");
                                    homogeneousClassificationRate[h] = homogeneousClassificationRate[h] + 1.0;
                                    for (int i = 0; i < dimX; i++) {
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
                        for (int i = 0; i < dimX; i++) {
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
                                    for (int i = 0; i < dimX; i++) {
                                        if (h == i) {
                                            standardErrors[i] += Math.pow(m.getEstimatedHomogeneousParameters().get(i, 0) - avgParameter[i], 2);
                                        }
                                    }
                                }
                            }
                        }
                    }
                    if (detectHomogeneity) {
                        for (int i = 0; i < dimX; i++) {
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
                    Y_MSE_var /= parallelLTM.size();
                    beta_MSE_var /= parallelLTM.size();

                    // jt.append("---------------------------------------------------------\n");
                    // jt.append("n = " + numObs + " Classification Rate by Parameter\n");
                    if (detectHomogeneity) {
                        for (int i = 0; i < homogeneousClassificationRate.length; i++) {
                            // jt.append(i + ". " + (homogeneousClassificationRate[i] / (0.0 + numMonteCarlos)) + "\n");
                        }
                    }
                    // jt.append("Y_MSE: " + Y_MSE + " (" + Y_MSE_var + ")\n");
                    // jt.append("beta_MSE: " + beta_MSE + " (" + beta_MSE_var + ")\n");
                    if (detectHomogeneity) {
                        for (int i = 0; i < dimX; i++) {
                            String s = String.format("x[%d]: mean: %g se: %g [%d]\n", i, avgParameter[i], standardErrors[i], counts[i]);
                            // jt.append("Parameter " + i + " mean: " + avgParameter[i] + " se: " + standardErrors[i] + " [" + counts[i] + "]\n");
                            // jt.append(s);
                        }
                    }
                    // jt.append("---------------------------------------------------------\n");

                    if (detectHomogeneity) {
                        YMSE_SD_restricted = Y_MSE_var;
                        YMSE_restricted = Y_MSE;
                        betaMSE_SD_restricted = beta_MSE_var;
                        betaMSE_restricted = beta_MSE;
                        beta1_mean = avgParameter[0];
                        beta2_mean = avgParameter[1];
                        beta1_SD = standardErrors[0];
                        beta2_SD = standardErrors[1];
                        classificationRate1 = homogeneousClassificationRate[0] / numMonteCarlos;
                        classificationRate2 = homogeneousClassificationRate[1] / numMonteCarlos;
                    } else {
                        YMSE_SD_unrestricted = Y_MSE_var;
                        YMSE_unrestricted = Y_MSE;
                        betaMSE_SD_unrestricted = beta_MSE_var;
                        betaMSE_unrestricted = beta_MSE;
                    }
                }
                LinearMonteCarloTable tf = new LinearMonteCarloTable(numObs, YMSE_unrestricted, YMSE_SD_unrestricted, YMSE_restricted, YMSE_SD_restricted, betaMSE_unrestricted, betaMSE_restricted, betaMSE_SD_unrestricted, betaMSE_SD_restricted, beta1_mean, beta1_SD, beta2_mean, beta2_SD, classificationRate1, classificationRate2);
                jam.append(tf.toString());
            }
        }
        System.out.println("Execution finished.");
    }

    public LinearTestMain(long rngSeed, int numObs, boolean detectHomogeneity, JTextArea jt, int dimX) {
        this.rngSeed = rngSeed;
        this.jt = jt;
        this.numObs = numObs;
        this.detectHomogeneity = detectHomogeneity;
        this.dimX = dimX;
    }

    private void execute() {
        System.out.println("**************** rngSeed = " + rngSeed + " ****************");
        Random rng = new Random(rngSeed);

        // MomentSpecification mySpecification = new LinearMomentSpecification("data/airline_subset.csv");
        MomentSpecification mySpecification = new LinearMomentSpecification(numObs, dimX);
        mySpecification.loadData(rng.nextLong()); // Create data using rng

//        Jama.Matrix dx = mySpecification.getX();
//        Jama.Matrix dxp = dx.transpose();
//        Jama.Matrix xpx = dxp.times(dx);
//        xpx.timesEquals(1.0/numObs);
//        jt.append(pmUtility.stringPrettyPrint(xpx)+"\n");
        double bestMinImprovement = 4.0;
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
        boolean verbose = false;
        if (numberTreesInForest == 1) {
            verbose = !true;
        }
        boolean testParameterHomogeneity;

        long rngBaseSeedMomentForest = rng.nextLong();
        long rngBaseSeedOutOfSample = rng.nextLong();

        boolean runCV = !true;
        if (runCV) {
            if (verbose) {
                System.out.println("************************");
                System.out.println("* Run Cross-Validation *");
                System.out.println("************************");
            }

            // for (int minObservationsPerLeaf = numObs/10; minObservationsPerLeaf <= 4 * numObs/10; minObservationsPerLeaf *= 2) {
            for (int minObservationsPerLeaf = 25; minObservationsPerLeaf <= 100; minObservationsPerLeaf *= 2) {
                // for (double minImprovement = 0.1; minImprovement <= 10.0; minImprovement *= 10) {
                double[] improvementLevels = {1, 5, 10}; // {1, 5, 10, 20}; // , 10.0, 20.0};
                for (double minImprovement : improvementLevels) {
                    for (int maxDepth = Math.min(3, numObs / (2 * minObservationsPerLeaf)); maxDepth >= 0; maxDepth--) {
                        computeOutOfSampleMSEInParameterSpace(mySpecification, numberTreesInForest, rngBaseSeedMomentForest, verbose, minObservationsPerLeaf, minImprovement, maxDepth, rngBaseSeedOutOfSample);
                        double combinationMSE = outOfSampleYMSE;
                        String star = "";
                        if (combinationMSE <= lowestSSE || first) {
                            lowestSSE = combinationMSE;
                            first = false;
                            bestMinImprovement = minImprovement;
                            bestMinObservationsPerLeaf = minObservationsPerLeaf;
                            bestMaxDepth = maxDepth;
                            star = "(*)";
                        }
                        if (star.contains("*") || numberTreesInForest == 1) {
                            System.out.println("minMSE: " + minImprovement + " minObs: " + minObservationsPerLeaf + " maxDepth: " + maxDepth + " Out-of-sample MSE: " + combinationMSE + " " + star);
                        }
                    }
                }
            }

            System.out.println("Lowest MSE: " + lowestSSE + " at min_N = " + bestMinObservationsPerLeaf + " min_MSE = " + bestMinImprovement + " maxDepth: " + bestMaxDepth);
            if (numberTreesInForest == 1) {
                // jt.append("Lowest MSE: " + lowestSSE + " at min_N = " + bestMinObservationsPerLeaf + " min_MSE = " + bestMinImprovement + " maxDepth: " + bestMaxDepth + "\n");
            }
        } else {
            bestMinObservationsPerLeaf = 50;
            bestMinImprovement = 10;
            if (numObs == 500) {
                bestMaxDepth = 2;
            }
            if (numObs == 1000) {
                bestMaxDepth = 3;
            }
            if (numObs == 2000) {
                bestMaxDepth = 4;
            }
            if (numObs == 4000) {
                bestMaxDepth = 5;
            }
        }
//        bestMinObservationsPerLeaf = 150;
//        bestMinImprovement = 1.0;
        // bestMaxDepth = 1;

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

            numberTreesInForest = 20;

            DataLens forestLens = new DataLens(mySpecification.getX(), mySpecification.getY(), mySpecification.getZ(), null);

            testParameterHomogeneity = false;
            TreeOptions cvOptions = new TreeOptions(minProportionInEachLeaf, bestMinObservationsPerLeaf, bestMinImprovement, bestMaxDepth, testParameterHomogeneity); // k = 1
            MomentForest myForest = new MomentForest(mySpecification, numberTreesInForest, rngBaseSeedMomentForest, forestLens, verbose, new TreeOptions());

            myForest.setTreeOptions(cvOptions);
            myForest.growForest();
            myForest.testHomogeneity();
            // TreeMoment loblolly = myForest.getTree(0);
            // loblolly.testHomogeneity();

            // should I implement a voting scheme here across all the trees for discerning homogeneity?
            // idea is to take the majority vote across trees
            // then take the average value from those trees that voted yes
            // let's try it and see what happens to classification rates
            ArrayList<Integer> hpl = new ArrayList<>();
            ArrayList<Double> hplStartingValues = new ArrayList<>();
            boolean[] voteIndexHomogeneity = myForest.getHomogeneityVotes(jt);
            double[] startingValues = myForest.getHomogeneityStartingValues();
            for (int i = 0; i < voteIndexHomogeneity.length; i++) {
                if (voteIndexHomogeneity[i]) {
                    System.out.println("Adding index " + i + " to homogeneous list with starting value: " + startingValues[i]);
                    hpl.add(i);
                    hplStartingValues.add(startingValues[i]);
                }
//                if(i==0 && !voteIndexHomogeneity[i]) {
//                    jt.append("BAD SEED: "+rngSeed+"\n");
//                }
            }

            // System.out.println("Done with growforest");
            // ArrayList<Integer> hpl = myForest.getTree(0).getIndexHomogeneousParameters(); // this is only using the first tree, is that the right way of thinking about this?
            // ArrayList<Double> hplStartingValues = myForest.getTree(0).getValueHomogeneousParameters();
            // System.out.println("Post get homogeneous parameters");
            // tell the specification that these parameters have been determined to be homogeneous
//            System.out.println("Unsorted");            
//                System.out.print(hpl+" ");
//                System.out.println(hplStartingValues);
            // need to sort together
            HomogeneousParameterSorter sorter = new HomogeneousParameterSorter();
            sorter.sort(hpl, hplStartingValues);
//            System.out.println("Sorted");
//            System.out.print(hpl+" ");
//            System.out.println(hplStartingValues);

            mySpecification.resetHomogeneityIndex();
            for (int i = 0; i < hpl.size(); i++) {
                mySpecification.setHomogeneousIndex(hpl.get(i));
                mySpecification.setHomogeneousParameter(hpl.get(i), hplStartingValues.get(i));
                // System.out.println(i);
            }
            setHomogeneousParameterList(hpl);
            // System.out.println("After setting hp list");

            // this seems to be working
            boolean testPostClassificationConvergence = false;
            if (testPostClassificationConvergence) {
                hpl.clear();
                hpl.add(0);
                mySpecification.resetHomogeneityIndex();
                mySpecification.setHomogeneousIndex(0);
                mySpecification.setHomogeneousParameter(0, -1.0);
            }

            /**
             * Estimate values of those homogeneous parameters
             */
            boolean cheatToVerifyWorking = false;
            if (cheatToVerifyWorking) {
                // this is here to verify that the code is working in that we should get a lower OOS MSE when the truth is imposed (it works)
                hpl.clear();
                hpl.add(1);
                // hpl.add(1);
                mySpecification.resetHomogeneityIndex();
                mySpecification.setHomogeneousIndex(0);
                mySpecification.setHomogeneousParameter(0, -1.0);
                // mySpecification.setHomogeneousIndex(1);
                // mySpecification.setHomogeneousParameter(1, 1.0);
                setEstimatedHomogeneousParameters(mySpecification.getHomogeneousParameterVector());
            } else {
                if (!hpl.isEmpty()) {
                    // System.out.println("Initializing search container");
                    numberTreesInForest = 1;
                    HomogeneousSearchContainer con = new HomogeneousSearchContainer(mySpecification, numberTreesInForest, verbose, bestMinImprovement, bestMinObservationsPerLeaf, bestMaxDepth,
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
        numberTreesInForest = 1; // using larger numbers here really improves fit!
        setEstimatedBetaVersusTruthMSE(computeOutOfSampleMSEInParameterSpace(mySpecification, numberTreesInForest, rngBaseSeedMomentForest, verbose, bestMinObservationsPerLeaf,
                bestMinImprovement, bestMaxDepth, rngBaseSeedOutOfSample));
    }

    public double getEstimatedBetaVersusTruthMSE() {
        return estimatedBetaVersusTruthMSE;
    }

    public void setEstimatedBetaVersusTruthMSE(double estimatedBetaVersusTruthMSE) {
        this.estimatedBetaVersusTruthMSE = estimatedBetaVersusTruthMSE;
    }

    private double computeOutOfSampleMSEInParameterSpace(MomentSpecification mySpecification, int numberTreesInForest, long rngBaseSeedMomentForest, boolean verbose,
            int minObservationsPerLeaf, double minImprovement, int maxTreeDepth, long rngBaseSeedOutOfSample) {
        // System.out.println("\nComputing OOS In Parameter Space\n");
        // System.out.println("Homogeneous parameter length in spec: "+mySpecification.getHomogeneousIndex().length);

        MomentForest myForest;

        double minProportionPerLeaf = 0.01;

        DataLens homogenizedForestLens = new DataLens(mySpecification.getX(), mySpecification.getY(), mySpecification.getZ(), null);

        myForest = new MomentForest(mySpecification, numberTreesInForest, rngBaseSeedMomentForest, homogenizedForestLens, verbose, new TreeOptions());
        TreeOptions cvOptions = new TreeOptions(minProportionPerLeaf, minObservationsPerLeaf, minImprovement, maxTreeDepth, false); // k = 1
        myForest.setTreeOptions(cvOptions);
        /**
         * Grow the moment forest
         */
        myForest.growForest();

//        System.out.println("First tree:");
//        myForest.getTree(0).printTree();
//        System.exit(0);
        /**
         * Test vectors for assessment
         */
        DataLens oosDataLens = mySpecification.getOutOfSampleXYZ(numObs, rngBaseSeedOutOfSample); // this should eventually be modified to come out of the data itself (or generalized in some way)
        Jama.Matrix testZ = oosDataLens.getZ();
        Jama.Matrix testX = oosDataLens.getX();
        Jama.Matrix testY = oosDataLens.getY();

        Random rng = new Random(rngBaseSeedOutOfSample - 3);

        double outOfSampleFit = 0;
        outOfSampleYMSE = 0;
        for (int i = 0; i < testZ.getRowDimension(); i++) {
            Jama.Matrix zi = testZ.getMatrix(i, i, 0, testZ.getColumnDimension() - 1);
            Jama.Matrix xi = testX.getMatrix(i, i, 0, testX.getColumnDimension() - 1);
            double yi = testY.get(i, 0);

            // going to compare directly to the true parameter vector in this method instead of using fit of Y
            Jama.Matrix bTruth = mySpecification.getBetaTruth(zi, rng);

            // have to reconstruct a composite beta from homogeneous and heterogeneous parameters
            Jama.Matrix compositeEstimatedBeta = myForest.getEstimatedParameterForest(zi);
            outOfSampleYMSE += mySpecification.getGoodnessOfFit(yi, xi, compositeEstimatedBeta);

            if (i < 10) {
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

            outOfSampleFit += pmUtility.sumSquaredElements((compositeEstimatedBeta.minus(bTruth)));
        }

        // ChartGenerator.makeXYScatter(xyc, "Fit Beta", "zi", "beta");
        // jt.append("betaMSE: " + (outOfSampleFit / testZ.getRowDimension()) + " \t [" + rngSeed + "]\n");
        outOfSampleYMSE /= testZ.getRowDimension();
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
        if (outOfSampleYMSE > 3) {
            jt.append("MSPE: " + outOfSampleYMSE + " seed: " + rngSeed + "\n");
        }
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
