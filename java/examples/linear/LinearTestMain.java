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

import JSci.maths.statistics.NormalDistribution;
import Jama.Matrix;
import com.google.common.math.Quantiles;
import core.ChartGenerator;
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
import javax.swing.SwingUtilities;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;
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
    private double outOfSampleFitNWEstimator = 0;
    private double outOfSampleFitParametricEstimator = 0;
    private double estimatedBetaVersusTruthMSE;
    private double betaMSE_parametric = 0;
    private Jama.Matrix estimatedHomogeneousParameters;
    private Jama.Matrix parametricParameters;
    private final boolean detectHomogeneity;
    private JTextArea jt;
    private final int dimX;
    private final int monteCarloIndex;

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
         * X,Z combinations, run Frobenius norm on that? Done that, seems to be
         * working really nicely.
         */
        /**
         * Number of Monte Carlos to run
         */
        int numMonteCarlos = 1;

        for (int dimX = 2; dimX <= 2; dimX++) {
            for (int numObs = 1000; numObs <= 100000; numObs *= 2) {

                double YMSE_unrestricted = 0;
                double YMSE_SD_unrestricted = 0;
                double YMSE_restricted = 0;
                double YMSE_SD_restricted = 0;

                double Y_nonparametric_MSE = 0;
                double Y_parametric_MSE = 0;

                double betaMSE_unrestricted = 0;
                double betaMSE_restricted = 0;
                double betaMSE_parametric = 0;
                double betaMSE_SD_unrestricted = 0;
                double betaMSE_SD_restricted = 0;

                double[] beta_mean = new double[dimX];
                double[] beta_mean_parametric = new double[dimX];
                double[] beta_SD = new double[dimX];
                double[] classificationRate = new double[dimX];

                JTextAreaAutoscroll jam = new JTextAreaAutoscroll();

                // boolean[] d = {false};
                boolean[] d = {false, true};

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

                    // jam.append("Dimensionality of X: " + dimX + "\n");
                    Random rng = new Random(22);

                    // int dimX = 3;
                    double[] homogeneousClassificationRate = new double[dimX];
                    double[] avgParameter = new double[dimX];
                    double Y_MSE = 0;

                    double Y_MSE_var = 0;
                    double beta_MSE = 0;
                    double beta_MSE_var = 0;

                    ArrayList<LinearTestMain> parallelLTM = new ArrayList<>();

                    for (int m = 0; m < numMonteCarlos; m++) {
                        LinearTestMain go;
                        if (numMonteCarlos == 1) {
                            // 8621193992485539638L
                            go = new LinearTestMain(m, 8621193992485539638L, numObs, detectHomogeneity, jam, dimX);
                        } else {
                            go = new LinearTestMain(m, rng.nextLong(), numObs, detectHomogeneity, jam, dimX);
                        }
                        parallelLTM.add(go);
                    }

                    AtomicInteger bomb = new AtomicInteger();

                    // parallelLTM.parallelStream().forEach(e -> {
                        parallelLTM.stream().forEach(e -> {
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
                                            if (m.getEstimatedHomogeneousParameters() != null) {
                                                avgParameter[i] += m.getEstimatedHomogeneousParameters().get(i, 0);
                                            }
                                            // jt.append("Estimated Homogeneous Parameter: "+pmUtility.stringPrettyPrintVector(m.getEstimatedHomogeneousParameters())+"\n");
                                            // jt.append("Average parameter sum for "+i+" is: "+avgParameter[i]+" Running count: "+counts[i]+"\n");
                                        }
                                    }
                                }
                            }
                        }
                        Y_MSE += m.getOutOfSampleYMSE();
                        if (!detectHomogeneity) {
                            Y_nonparametric_MSE += m.getOutOfSampleFitNWEstimator();
                            Y_parametric_MSE += m.getOutOfSampleFitParametricEstimator();
                            betaMSE_parametric += m.getBetaMSE_parametric();
                        }
                        // jt2.append("Rolling: "+Y_nonparametric_MSE+" Marginal: "+m.getOutOfSampleFitNWEstimator()+"\n");
                        beta_MSE += m.getEstimatedBetaVersusTruthMSE();
                        // jt2.append("forest_beta_MSE: "+beta_MSE+" marginal: "+m.getEstimatedBetaVersusTruthMSE()+"\n");
                        // jt2.append("OLS_beta_MSE: "+betaMSE_parametric+" marginal: "+m.getBetaMSE_parametric()+"\n");
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
                            if (!hList.isEmpty()) {
                                for (Integer h : hList) {
                                    for (int i = 0; i < dimX; i++) {
                                        if (h == i) {
                                            if (m.getEstimatedHomogeneousParameters() != null) {
                                                standardErrors[i] += Math.pow(m.getEstimatedHomogeneousParameters().get(i, 0) - avgParameter[i], 2);
                                            }
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
                    if (!detectHomogeneity) {
                        Y_nonparametric_MSE /= parallelLTM.size();
                        Y_parametric_MSE /= parallelLTM.size();
                        betaMSE_parametric /= parallelLTM.size();
                    }
                    // jt2.append("Average NP_YMSE: "+Y_nonparametric_MSE+"\n");
                    beta_MSE /= parallelLTM.size();

                    for (LinearTestMain m : parallelLTM) {
                        Y_MSE_var += Math.pow(Y_MSE - m.getOutOfSampleYMSE(), 2);
                        beta_MSE_var += Math.pow(beta_MSE - m.getEstimatedBetaVersusTruthMSE(), 2);
                    }
                    Y_MSE_var /= parallelLTM.size();
                    beta_MSE_var /= parallelLTM.size();

                    if (detectHomogeneity) {
                        YMSE_SD_restricted = Y_MSE_var;
                        YMSE_restricted = Y_MSE;
                        betaMSE_SD_restricted = beta_MSE_var;
                        betaMSE_restricted = beta_MSE;
                        for (int k = 0; k < dimX; k++) {
                            beta_mean[k] = avgParameter[k];
                            beta_SD[k] = standardErrors[k];
                            classificationRate[k] = homogeneousClassificationRate[k] / numMonteCarlos;
                        }
                    } else {
                        YMSE_SD_unrestricted = Y_MSE_var;
                        YMSE_unrestricted = Y_MSE;
                        betaMSE_SD_unrestricted = beta_MSE_var;
                        betaMSE_unrestricted = beta_MSE;

                        Jama.Matrix avgParametricParameters = new Jama.Matrix(dimX, 1);
                        for (LinearTestMain m : parallelLTM) {
                            avgParametricParameters.plusEquals(m.getParametricParameters()); 
                        }
                        avgParametricParameters.timesEquals(1.0 / parallelLTM.size());
                        for (int i = 0; i < avgParametricParameters.getRowDimension(); i++) {
                            beta_mean_parametric[i] = avgParametricParameters.get(i, 0);
                        }

                    }
                }

                LinearMonteCarloTable tf = new LinearMonteCarloTable(numObs, YMSE_unrestricted, YMSE_SD_unrestricted, YMSE_restricted, YMSE_SD_restricted,
                        betaMSE_unrestricted, betaMSE_restricted, betaMSE_SD_unrestricted, betaMSE_SD_restricted, beta_mean, beta_SD, classificationRate, Y_nonparametric_MSE, Y_parametric_MSE,
                        beta_mean_parametric, betaMSE_parametric);
                jam.append(tf.toString());
            }
        }
        System.out.println("Execution finished.");
    }

    public LinearTestMain(int monteCarloIndex, long rngSeed, int numObs, boolean detectHomogeneity, JTextArea jt, int dimX) {
        this.rngSeed = rngSeed;
        this.monteCarloIndex = monteCarloIndex;
        this.jt = jt;
        this.numObs = numObs;
        this.detectHomogeneity = detectHomogeneity;
        this.dimX = dimX;
        parametricParameters = new Jama.Matrix(dimX, 1);
    }

    private void publish(String s, JTextArea jelly) {
        SwingUtilities.invokeLater(() -> {
            jelly.append(s);
        });
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

        int numberTreesInForest = 50;
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

        boolean computeAlternatives = false;
        if (computeAlternatives) {
            System.out.println("Computing nonparametric fits...");
            computeOutOfSampleNonparametricFits(mySpecification, rngBaseSeedOutOfSample);
            System.out.println("Computing completely parametric fits...");
            computeOutOfSampleParametricFits(mySpecification, rngBaseSeedOutOfSample);
        }

        // partial linear model results from CV:
        // 500: 50, 10, 2
        // 
        boolean runCV = !false;
        if (runCV) {
            if (verbose) {
                System.out.println("************************");
                System.out.println("* Run Cross-Validation *");
                System.out.println("************************");
            }

            // for (int minObservationsPerLeaf = numObs/10; minObservationsPerLeaf <= 4 * numObs/10; minObservationsPerLeaf *= 2) {
            for (int minObservationsPerLeaf = 25; minObservationsPerLeaf <= 25; minObservationsPerLeaf *= 2) {
                // for (double minImprovement = 0.1; minImprovement <= 10.0; minImprovement *= 10) {
                double[] improvementLevels = {2}; // {1,2,5,10}; // , 10, 20, 50}; // 0.1, 0.5, 1.0}; // {1, 5, 10, 20}; // , 10.0, 20.0};
                for (double minImprovement : improvementLevels) {
                    // for (int maxDepth = Math.min(5, numObs / (2 * minObservationsPerLeaf)); maxDepth >= 0; maxDepth--) {
                    for (int maxDepth = 8; maxDepth >= 0; maxDepth--) {
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
                        if (star.contains("*") || 1 == 1) {
                            System.out.println(monteCarloIndex + ". minMSE: " + minImprovement + " minObs: " + minObservationsPerLeaf + " maxDepth: " + maxDepth + " Out-of-sample MSE: " + combinationMSE + " " + star);
                        }
                    }
                }
            }

            System.out.println("Lowest MSE: " + lowestSSE + " at min_N = " + bestMinObservationsPerLeaf + " min_MSE = " + bestMinImprovement + " maxDepth: " + bestMaxDepth);
            // jt.append("Lowest MSE: " + lowestSSE + " at min_N = " + bestMinObservationsPerLeaf + " min_MSE = " + bestMinImprovement + " maxDepth: " + bestMaxDepth+"\n");
            if (numberTreesInForest == 10) {
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

            bestMinObservationsPerLeaf = 10;
            bestMinImprovement = 5.0;
            bestMaxDepth = 6;
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

            numberTreesInForest = 100;

            DataLens forestLens = new DataLens(mySpecification.getX(), mySpecification.getY(), mySpecification.getZ(), null);

            testParameterHomogeneity = false;
            TreeOptions cvOptions = new TreeOptions(minProportionInEachLeaf, bestMinObservationsPerLeaf, bestMinImprovement, bestMaxDepth, testParameterHomogeneity); // k = 1
            MomentForest myForest = new MomentForest(mySpecification, numberTreesInForest, rngBaseSeedMomentForest, forestLens, verbose, new TreeOptions());
            myForest.setTreeOptions(cvOptions);
            System.out.println("Growing forest for homogeneity testing...");
            myForest.growForest();
            System.out.println("----- Call to testHomogeneity -----");
            myForest.testHomogeneity();

            ArrayList<Integer> hpl = new ArrayList<>();
            ArrayList<Double> hplStartingValues = new ArrayList<>();

            boolean verboseVoting = true;
            boolean[] voteIndexHomogeneity = myForest.getHomogeneityVotes(jt, verboseVoting);

            double[] startingValues = myForest.getHomogeneityStartingValues();
            for (int i = 0; i < voteIndexHomogeneity.length; i++) {
                if (voteIndexHomogeneity[i]) {
                    if (verbose) {
                        System.out.println("Adding index " + i + " to homogeneous list with starting value: " + startingValues[i]);
                        jt.append("Adding index " + i + " to homogeneous list with starting value: " + startingValues[i] + "\n");
                    }
                    hpl.add(i);
                    hplStartingValues.add(startingValues[i]);
                }
//                if(i==0 && !voteIndexHomogeneity[i]) {
//                    jt.append("BAD SEED: "+rngSeed+"\n");
//                }
            }

            HomogeneousParameterSorter sorter = new HomogeneousParameterSorter();
            sorter.sort(hpl, hplStartingValues);

            mySpecification.resetHomogeneityIndex();
            for (int i = 0; i < hpl.size(); i++) {
                mySpecification.setHomogeneousIndex(hpl.get(i));
                mySpecification.setHomogeneousParameter(hpl.get(i), hplStartingValues.get(i));
                // System.out.println(i);
            }
            setHomogeneousParameterList(hpl);

            boolean executeSearch = false;

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
            } else if (executeSearch) {
                if (!hpl.isEmpty()) {
                    // System.out.println("Initializing search container");
                    numberTreesInForest = 10;
                    HomogeneousSearchContainer con = new HomogeneousSearchContainer(mySpecification, numberTreesInForest, verbose, bestMinImprovement, bestMinObservationsPerLeaf, bestMaxDepth,
                            getHomogeneousParameterList(), rngBaseSeedMomentForest, rngBaseSeedOutOfSample);
                    // System.out.println("Calling execute search");

                    con.executeSearch();

                    // System.out.println("Post search");
                    Jama.Matrix homogeneousParameters = con.getEstimatedHomogeneousParameters();
                    if (verbose || 1 == 2) {
                        System.out.print("Post-HomogeneousSearchContainer Estimated homogeneous parameters: ");
                        jt.append("Post search: " + pmUtility.stringPrettyPrintVector(homogeneousParameters) + "\n");
                        pmUtility.prettyPrintVector(homogeneousParameters); // this is a compact vector of parameters
                    }

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
        System.out.println("Computing final trees...");
        numberTreesInForest = 10; // using larger numbers here really improves fit!
        verbose = false;
        setEstimatedBetaVersusTruthMSE(computeOutOfSampleMSEInParameterSpace(mySpecification, numberTreesInForest, rngBaseSeedMomentForest, verbose, bestMinObservationsPerLeaf,
                bestMinImprovement, bestMaxDepth, rngBaseSeedOutOfSample));
        System.out.println("Finished execute.");
    }

    public double getEstimatedBetaVersusTruthMSE() {
        return estimatedBetaVersusTruthMSE;
    }

    public void setEstimatedBetaVersusTruthMSE(double estimatedBetaVersusTruthMSE) {
        this.estimatedBetaVersusTruthMSE = estimatedBetaVersusTruthMSE;
    }

    private double computeOutOfSampleNonparametricFits(MomentSpecification mySpecification, long rngBaseSeedOutOfSample) {
        DataLens oosDataLens = mySpecification.getOutOfSampleXYZ(numObs, rngBaseSeedOutOfSample); // this should eventually be modified to come out of the data itself (or generalized in some way)
        Jama.Matrix testZ = oosDataLens.getZ();
        Jama.Matrix testX = oosDataLens.getX();
        Jama.Matrix testY = oosDataLens.getY();

        // do CV for optimal NW bandwidth parameter
        double optimalBandwidthNWValue = 1.0;
        double optimalBandwidthNW = 1.0;
        boolean firstNWBandwidth = true;
        for (double bandwidth = 3.0; bandwidth >= 1E-1; bandwidth /= 2.0) {
            double bVal = 0;
            for (int i = 0; i < mySpecification.getX().getRowDimension(); i++) {
                Jama.Matrix xi = testX.getMatrix(i, i, 0, mySpecification.getX().getColumnDimension() - 1);
                Jama.Matrix zi = testZ.getMatrix(i, i, 0, mySpecification.getZ().getColumnDimension() - 1);
                double nwEstimatorY = getNWEstimatorY(pmUtility.concatMatrix(xi, zi), mySpecification.getY(), mySpecification.getX(), mySpecification.getZ(), bandwidth);
                if (i < -10) {
                    System.out.println("y_i: " + testY.get(i, 0) + " yhat_i: " + nwEstimatorY + " x_i: " + pmUtility.stringPrettyPrint(xi) + " zi: " + pmUtility.stringPrettyPrint(zi));
                }
                bVal += Math.pow(testY.get(i, 0) - nwEstimatorY, 2);
            }
            bVal /= mySpecification.getX().getRowDimension();
            if (bVal < optimalBandwidthNWValue || firstNWBandwidth) {
                firstNWBandwidth = false;
                optimalBandwidthNWValue = bVal;
                optimalBandwidthNW = bandwidth;
            }
            // System.out.println(bandwidth + " " + bVal);
        }
        // System.out.println("Optimal NW bandwidth: " + optimalBandwidthNW + " MSE(Y): " + optimalBandwidthNWValue);
        // System.exit(0);

        outOfSampleFitNWEstimator = 0;

        for (int i = 0; i < testZ.getRowDimension(); i++) {
            Jama.Matrix zi = testZ.getMatrix(i, i, 0, testZ.getColumnDimension() - 1);
            Jama.Matrix xi = testX.getMatrix(i, i, 0, testX.getColumnDimension() - 1);

            double nwEstimatorY = getNWEstimatorY(pmUtility.concatMatrix(xi, zi), mySpecification.getY(), mySpecification.getX(), mySpecification.getZ(), optimalBandwidthNW);
            outOfSampleFitNWEstimator += Math.pow(testY.get(i, 0) - nwEstimatorY, 2);
        }
        outOfSampleFitNWEstimator /= testZ.getRowDimension();

        // System.out.println("NW Estimator MSE(Y): " + outOfSampleFitNWEstimator);
        return outOfSampleFitNWEstimator;
    }

    private double computeOutOfSampleParametricFits(MomentSpecification mySpecification, long rngBaseSeedOutOfSample) {
        DataLens oosDataLens = mySpecification.getOutOfSampleXYZ(numObs, rngBaseSeedOutOfSample); // this should eventually be modified to come out of the data itself (or generalized in some way)
        Jama.Matrix testX = oosDataLens.getX();
        Jama.Matrix testY = oosDataLens.getY();
        Jama.Matrix testZ = oosDataLens.getZ();

        outOfSampleFitParametricEstimator = 0;

        boolean[] fillerIndex = new boolean[mySpecification.getHomogeneousIndex().length];
        for (int i = 0; i < fillerIndex.length; i++) {
            fillerIndex[i] = false;
        }
        ContainerLinear cl = new ContainerLinear(new DataLens(mySpecification.getX(), mySpecification.getY(), mySpecification.getZ(), null), fillerIndex, null, false);
        cl.computeBetaAndErrors();
        Jama.Matrix beta = cl.getBeta();
        setParametricParameters(beta);

        Random rng = new Random(rngBaseSeedOutOfSample - 3);
        betaMSE_parametric = 0;

        for (int i = 0; i < testX.getRowDimension(); i++) {
            Jama.Matrix xi = testX.getMatrix(i, i, 0, testX.getColumnDimension() - 1);
            outOfSampleFitParametricEstimator += mySpecification.getGoodnessOfFit(testY.get(i, 0), xi, beta);

            Jama.Matrix zi = testZ.getMatrix(i, i, 0, testZ.getColumnDimension() - 1);

            // going to compare directly to the true parameter vector in this method instead of using fit of Y
            Jama.Matrix bTruth = mySpecification.getBetaTruth(zi, rng);
            betaMSE_parametric += pmUtility.sumSquaredElements((beta.minus(bTruth)));
            if (i < -10) {
                System.out.print("OLS: ");
                pmUtility.prettyPrintVector(beta);
            }
        }
        outOfSampleFitParametricEstimator /= testX.getRowDimension();
        betaMSE_parametric /= testX.getRowDimension();

        // System.out.println("Parametric MSE(Y): " + outOfSampleFitParametricEstimator+" MSE(beta): "+betaMSE_parametric);
        return outOfSampleFitParametricEstimator;
    }

    private double computeOutOfSampleMSEInParameterSpace(MomentSpecification mySpecification, int numberTreesInForest, long rngBaseSeedMomentForest, boolean verbose,
            int minObservationsPerLeaf, double minImprovement, int maxTreeDepth, long rngBaseSeedOutOfSample) {
        // System.out.println("\nComputing OOS In Parameter Space\n");
        // System.out.println("Homogeneous parameter length in spec: "+mySpecification.getHomogeneousIndex().length);

        MomentForest myForest;

        double minProportionPerLeaf = 0.01;

        DataLens homogenizedForestLens = new DataLens(mySpecification.getX(), mySpecification.getY(), mySpecification.getZ(), null);

//        for(Boolean b : mySpecification.getHomogeneousIndex()) {
//            System.out.println(b);
//        }
//        for(int i=0;i<mySpecification.getHomogeneousParameterVector().getRowDimension();i++) {
//            System.out.println(mySpecification.getHomogeneousParameterVector().get(i,0));
//        }
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
        
        XYSeries betaTruthXY = new XYSeries("True Beta");
        XYSeries betaEstimateXY = new XYSeries("Estimated Beta");
        XYSeries betaEstimateXY5 = new XYSeries("Estimated Beta 5th Percentile");
        XYSeries betaEstimateXY95 = new XYSeries("Estimated Beta 95th Percentile");

        outOfSampleYMSE = 0;
        for (int i = 0; i < testZ.getRowDimension(); i++) {
            Jama.Matrix zi = testZ.getMatrix(i, i, 0, testZ.getColumnDimension() - 1);
            
            // zi.set(0, 0, -1.5 + i * 1.0 / (testZ.getRowDimension()));
            
            Jama.Matrix xi = testX.getMatrix(i, i, 0, testX.getColumnDimension() - 1);
            double yi = testY.get(i, 0);

            Jama.Matrix bTruth = mySpecification.getBetaTruth(zi, rng);
            betaTruthXY.add(zi.get(0,0), bTruth.get(0,0));
            
            // have to reconstruct a composite beta from homogeneous and heterogeneous parameters
            Jama.Matrix compositeEstimatedBeta = myForest.getEstimatedParameterForest(zi);
            ArrayList<Jama.Matrix> allParametersInForest = myForest.getAllEstimatedParametersFromForest(zi);
            
            ArrayList<Double> betaZeroList = new ArrayList<>();
            for(int k=0;k<allParametersInForest.size();k++) {
                betaZeroList.add(allParametersInForest.get(k).get(0,0));
            }
            betaEstimateXY5.add(zi.get(0,0), Quantiles.percentiles().index(5).compute(betaZeroList));
            betaEstimateXY95.add(zi.get(0,0), Quantiles.percentiles().index(95).compute(betaZeroList));
            
            betaEstimateXY.add(zi.get(0,0), compositeEstimatedBeta.get(0,0));
            if (i == 0 && 1 == 2) {
                jt.append("Composite beta: " + pmUtility.stringPrettyPrintVector(compositeEstimatedBeta) + "\n");
            }
            outOfSampleYMSE += mySpecification.getGoodnessOfFit(yi, xi, compositeEstimatedBeta);

            if (i < -10) {
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

                System.out.print("tree: ");
                pmUtility.prettyPrintVector(compositeEstimatedBeta);

            }
            //pmUtility.prettyPrintVector(compositeEstimatedBeta);
            outOfSampleFit += pmUtility.sumSquaredElements((compositeEstimatedBeta.minus(bTruth)));
        }

        XYSeriesCollection xyc = new XYSeriesCollection(betaTruthXY);
        xyc.addSeries(betaEstimateXY);
        // xyc.addSeries(betaEstimateXY5);
        // xyc.addSeries(betaEstimateXY95);
        
        ChartGenerator.makeXYScatter(xyc, "Fit Beta", "zi", "beta");
        // jt.append("betaMSE: " + (outOfSampleFit / testZ.getRowDimension()) + " \t [" + rngSeed + "]\n");
        outOfSampleYMSE /= testZ.getRowDimension();
        // System.out.println("Forest MSE(Y): "+outOfSampleYMSE+" MSE(beta): "+outOfSampleFit/testZ.getRowDimension());

        if (1 == 2) {
            jt.append("MSE(beta): " + outOfSampleFit / testZ.getRowDimension() + "\n");
        }
        return outOfSampleFit / testZ.getRowDimension(); // mse
    }

    private void setHomogeneousParameterList(ArrayList<Integer> homogeneousParameterList) {
        this.homogeneousParameterList = homogeneousParameterList;
    }

    public ArrayList<Integer> getHomogeneousParameterList() {
        return homogeneousParameterList;
    }

    private double getNWEstimatorY(Jama.Matrix pointX, Jama.Matrix allY, Jama.Matrix allX, Jama.Matrix allZ, double h) {
        double nwTop = 0;
        double nwBottom = 0;
        for (int i = 0; i < allX.getRowDimension(); i++) {
            Jama.Matrix rowXi = allX.getMatrix(i, i, 0, allX.getColumnDimension() - 1);
            Jama.Matrix rowZi = allZ.getMatrix(i, i, 0, allZ.getColumnDimension() - 1);
            double Kh = Kh(pointX.minus(pmUtility.concatMatrix(rowXi, rowZi)), h);
            if (i < -10) {
                System.out.println("h = " + h + " Kh: " + Kh + " xi: " + pmUtility.stringPrettyPrint(pmUtility.concatMatrix(rowXi, rowZi)) + " pointXi: " + pmUtility.stringPrettyPrint(pointX));
            }
            nwTop += Kh * allY.get(i, 0);
            nwBottom += Kh;
        }
        double nw = nwTop / nwBottom;
        if (Double.isNaN(nw)) {
            // System.out.println("nwTop: "+nwTop+" nwBottom: "+nwBottom);
            return 0; // i'm not sure this is right!
        }
        return nwTop / nwBottom;
    }

    private double Kh(Jama.Matrix pointX, double h) {
        NormalDistribution normal = new NormalDistribution();
        double v = 1;
        for (int i = 0; i < pointX.getColumnDimension(); i++) {
            v *= (1.0 / h) * normal.probability(pointX.get(0, i) / h);
        }
        return v;
    }

    /**
     * @return the MSE
     */
    public double getOutOfSampleYMSE() {
//        if (outOfSampleYMSE > 3) {
//            jt.append("MSPE: " + outOfSampleYMSE + " seed: " + rngSeed + "\n");
//        }
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

    /**
     * @return the outOfSampleFitNWEstimator
     */
    public double getOutOfSampleFitNWEstimator() {
        return outOfSampleFitNWEstimator;
    }

    /**
     * @return the outOfSampleFitParametricEstimator
     */
    public double getOutOfSampleFitParametricEstimator() {
        return outOfSampleFitParametricEstimator;
    }

    /**
     * @return the parametricParameters
     */
    public Jama.Matrix getParametricParameters() {
        return parametricParameters;
    }

    /**
     * @param parametricParameters the parametricParameters to set
     */
    public void setParametricParameters(Jama.Matrix parametricParameters) {
        this.parametricParameters = parametricParameters;
    }

    /**
     * @return the betaMSE_parametric
     */
    public double getBetaMSE_parametric() {
        return betaMSE_parametric;
    }

    /**
     * @return the monteCarloIndex
     */
    public int getMonteCarloIndex() {
        return monteCarloIndex;
    }

}
