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
package examples.logit;

import core.OutOfSampleStatisticsContainer;
import Jama.Matrix;
import core.DataLens;
import core.HomogeneousSearchContainer;
import core.MomentForest;
import core.MomentSpecification;
import core.TreeMoment;
import core.TreeOptions;
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
public class LogitTestMain {

    private ArrayList<Integer> homogeneousParameterList;
    private final long rngSeed;
    private final int numObs;
    private double outOfSampleYMSE;
    private double estimatedBetaVersusTruthMSE;
    private Jama.Matrix estimatedHomogeneousParameters;
    private final boolean detectHomogeneity;
    private JTextArea jt;
    static private JTextArea debugOutputArea = new JTextArea();

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        JFrame f = new JFrame("Monte Carlo");
        f.setBounds(100, 100, 1500, 500);
        f.getContentPane().setLayout(new GridLayout(1, 2));
        JTextAreaAutoscroll jt1 = new JTextAreaAutoscroll();
        f.getContentPane().add(new JScrollPane(jt1));
        JTextAreaAutoscroll jt2 = new JTextAreaAutoscroll();
        f.getContentPane().add(new JScrollPane(jt2));
        f.getContentPane().add(new JScrollPane(debugOutputArea));
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
        for (int numObs = 500; numObs <= 4000; numObs *= 2) {

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
                    //jam.append("***** UNRESTRICTED MODEL *****\n");
                }
                jam.append("\\\\ \n");

                Random rng = new Random(22);

                int numParameters = 2;

                double[] homogeneousClassificationRate = new double[numParameters];
                double[] avgParameter = new double[numParameters];
                double Y_MSE = 0;
                double Y_MSE_var = 0;
                double beta_MSE = 0;
                double beta_MSE_var = 0;

                int numMonteCarlos = 50;

                ArrayList<LogitTestMain> parallelLTM = new ArrayList<>();

                for (int m = 0; m < numMonteCarlos; m++) {
                    LogitTestMain go;
                    if (numMonteCarlos == 1) {
                        go = new LogitTestMain(4606544446801080638L, numObs, detectHomogeneity, jam);
                    } else {
                        go = new LogitTestMain(rng.nextLong(), numObs, detectHomogeneity, jam);
                    }
                    parallelLTM.add(go);
                }

                AtomicInteger bomb = new AtomicInteger();

                parallelLTM.parallelStream().forEach(e -> {
                    // parallelLTM.stream().forEach(e -> {
                    e.execute();
                    bomb.incrementAndGet();
                    System.out.println("Finished " + bomb.get() + " iterations.");
                });

                // for (int m = 0; m < numMonteCarlos; m++) {
                int[] counts = new int[numParameters];
                double[] standardErrors = new double[numParameters];
                for (LogitTestMain m : parallelLTM) {
                    if (detectHomogeneity) {
                        ArrayList<Integer> hList = m.getHomogeneousParameterList();
                        if (!hList.isEmpty()) {
                            for (Integer h : hList) {
                                // jt1.append("Detected homogeneity on parameter index " + h + " hList.size = " + hList.size());
                                // jt1.append(" estimatedParameterSize: " + m.getEstimatedHomogeneousParameters().getRowDimension() + "\n");
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
                for (LogitTestMain m : parallelLTM) {
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
                for (LogitTestMain m : parallelLTM) {
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
                    for (int i = 0; i < numParameters; i++) {
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
            LogitMonteCarloTable tf = new LogitMonteCarloTable(numObs, YMSE_unrestricted, YMSE_SD_unrestricted, YMSE_restricted, YMSE_SD_restricted, betaMSE_unrestricted, betaMSE_restricted,
                    betaMSE_SD_unrestricted, betaMSE_SD_restricted, beta1_mean, beta1_SD, beta2_mean, beta2_SD, classificationRate1, classificationRate2);
            jam.append(tf.toString());
        }
    }

    public LogitTestMain(long rngSeed, int numObs, boolean detectHomogeneity, JTextArea jt) {
        this.rngSeed = rngSeed;
        this.jt = jt;
        this.numObs = numObs;
        this.detectHomogeneity = detectHomogeneity;
    }

    private void execute() {
        System.out.println("**************** rngSeed = " + rngSeed + " ****************");
        Random rng = new Random(rngSeed);

        MomentSpecification mySpecification = new LogitMomentSpecification(numObs);
        mySpecification.loadData(rng.nextLong()); // Create data using rng

        double bestMinImprovement = 0.8;
        int bestMinObservationsPerLeaf = 25;
        int bestMaxDepth = 2;

        double lowestSSE = 0;
        boolean first = true;

        /**
         * Make sure that cross-validation is run on completely unrestricted
         * model; set all parameters to heterogeneous
         */
        mySpecification.resetHomogeneityIndex();

        int numberTreesInForest = 100;
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

        boolean runCV = !false;
        if (runCV) {
            if (verbose) {
                System.out.println("************************");
                System.out.println("* Run Cross-Validation *");
                System.out.println("************************");
            }

            for (int minObservationsPerLeaf = 50; minObservationsPerLeaf <= 200; minObservationsPerLeaf *= 2) {
                for (double minImprovement = 0.1; minImprovement <= 0.1; minImprovement *= 2) {
                    for (int maxDepth = Math.min(5, numObs / (2 * minObservationsPerLeaf)); maxDepth >= 0; maxDepth--) {
                        // for (int maxDepth = 1; maxDepth >= 1; maxDepth--) {
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
                        debugOutputArea.append("detectHomogeneity: " + detectHomogeneity + " minMSE: " + minImprovement + " minObs: " + minObservationsPerLeaf + " maxDepth: " + maxDepth + " Out-of-sample MSE: " + combinationMSE + " " + star + "\n");

                    }
                }
            }

            System.out.println("Lowest MSE: " + lowestSSE + " at min_N = " + bestMinObservationsPerLeaf + " min_MSE = " + bestMinImprovement + " maxDepth: " + bestMaxDepth);
            jt.append("detectHomogeneity: " + detectHomogeneity + " CV Lowest MSE: " + lowestSSE + " at min_N = " + bestMinObservationsPerLeaf + " min_MSE = " + bestMinImprovement + " maxDepth: " + bestMaxDepth + "\n");
        } else {
            bestMinObservationsPerLeaf = 20;
            bestMinImprovement = 0.1;
            if (numObs == 500) {
                bestMaxDepth = 1;
                bestMinObservationsPerLeaf = 50;
            }
            if (numObs == 1000) {
                bestMaxDepth = 1;
                bestMinObservationsPerLeaf = 100;
            }
            if (numObs == 2000) {
                bestMaxDepth = 1;
                bestMinObservationsPerLeaf = 100;
            }
            if (numObs == 4000) {
                bestMaxDepth = 1;
                bestMinObservationsPerLeaf = 100;
            }
        }
        
        bestMaxDepth = 2;

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

            testParameterHomogeneity = false;
            TreeOptions cvOptions = new TreeOptions(minProportionInEachLeaf, bestMinObservationsPerLeaf, bestMinImprovement, bestMaxDepth, testParameterHomogeneity); // k = 1
            MomentForest myForest = new MomentForest(mySpecification, numberTreesInForest, rngBaseSeedMomentForest, forestLens, verbose, new TreeOptions());

            myForest.setTreeOptions(cvOptions);
            myForest.growForest();

            ArrayList<Integer> hpl = new ArrayList<>();

            /**
             * Have each tree in the forest vote for homogeneity; take majority vote as classification
             */
            for (int k = 0; k < mySpecification.getNumParams(); k++) {
                double voteShareHomogeneousParameterK = 0;
                for (int jelly = 0; jelly < myForest.getForestSize(); jelly++) {
                    myForest.getTree(jelly).testHomogeneity();
                    ArrayList<Integer> hpl_jelly = myForest.getTree(jelly).getIndexHomogeneousParameters();
                    if (hpl_jelly.contains(k)) {
                        voteShareHomogeneousParameterK++;
                    }
                }
                voteShareHomogeneousParameterK /= myForest.getForestSize();
                System.out.println("Proportion of trees that voted for parameter "+k+" as homogeneous: " + voteShareHomogeneousParameterK);
                if (voteShareHomogeneousParameterK >= 0.5) {
                    hpl.add(k);
                }
            }

            /**
             * For each tree that voted for homogeneity, average homogeneous parameter values
             */
            // ArrayList<Double> hplStartingValues = new ArrayList<>();
            for(Integer k : hpl) {
                double counter = 0;
                double avgValueParameterK = 0;
                for(int m=0;m<myForest.getForestSize();m++) {
                    TreeMoment mt = myForest.getTree(m);
                    if(mt.getIndexHomogeneousParameters().contains(k)) {
                        counter++;
                        avgValueParameterK += mt.getValueHomogeneousParameters().get(mt.getIndexHomogeneousParameters().indexOf(k));
                    }
                }
                avgValueParameterK /= counter;
                mySpecification.setHomogeneousIndex(k);
                mySpecification.setHomogeneousParameter(k, avgValueParameterK);
                System.out.println("Setting parameter "+k+" to homogeneous with value of "+avgValueParameterK);
            }
            

//            for (int i = 0; i < hpl.size(); i++) {
//                mySpecification.setHomogeneousIndex(hpl.get(i));
//                mySpecification.setHomogeneousParameter(hpl.get(i), hplStartingValues.get(i));
//                // System.out.println(i);
//            }
            setHomogeneousParameterList(hpl);
            // System.out.println("After setting hp list");

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
                mySpecification.setHomogeneousIndex(1);
                mySpecification.setHomogeneousParameter(1, 1.0);
                // mySpecification.setHomogeneousIndex(1);
                // mySpecification.setHomogeneousParameter(1, 1.0);
                setEstimatedHomogeneousParameters(mySpecification.getHomogeneousParameterVector());
            } else {
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
        }

        /**
         * Compute out-of-sample measures of fit (against Y, and true beta)
         */
        // numberTreesInForest = 1;
        OutOfSampleStatisticsContainer results = mySpecification.computeOutOfSampleStatistics(numberTreesInForest, rngBaseSeedMomentForest, verbose, bestMinObservationsPerLeaf,
                bestMinImprovement, bestMaxDepth, rngBaseSeedOutOfSample);
        setEstimatedBetaVersusTruthMSE(results.getMeanSquaredErrorParameters());
        setOutOfSampleYMSE(results.getOutOfSampleMeasureY());
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
