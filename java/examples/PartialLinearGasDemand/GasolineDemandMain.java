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
package examples.PartialLinearGasDemand;

import JSci.maths.statistics.NormalDistribution;
import Jama.Matrix;
import core.ChartGenerator;
import core.DataLens;
import core.HomogeneousSearchContainer;
import core.MomentForest;
import core.MomentSpecification;
import core.TreeMoment;
import core.TreeOptions;
import java.awt.GridLayout;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Random;
import javax.swing.JFrame;
import javax.swing.JScrollPane;
import javax.swing.JTextArea;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;
import utility.JTextAreaAutoscroll;
import utility.pmUtility;

/**
 *
 * @author Stephen P. Ryan <stephen.p.ryan@wustl.edu>
 */
public class GasolineDemandMain {

    private ArrayList<Integer> homogeneousParameterList;
    private Jama.Matrix estimatedHomogeneousParameters;
    private final boolean detectHomogeneity;
    final private JTextArea jt;

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        JFrame f = new JFrame("Gasoline Demand Partial Linear Model Application");
        f.setBounds(100, 100, 1500, 500);
        f.getContentPane().setLayout(new GridLayout(1, 2));
        JTextAreaAutoscroll jt1 = new JTextAreaAutoscroll();
        f.getContentPane().add(new JScrollPane(jt1));
        JTextAreaAutoscroll jt2 = new JTextAreaAutoscroll();
        f.getContentPane().add(new JScrollPane(jt2));
        f.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        f.setVisible(true);

        GasolineDemandMain go = new GasolineDemandMain(jt1);
        go.execute();

    }

    public GasolineDemandMain(JTextArea jt) {
        this.jt = jt;
        detectHomogeneity = true;
    }

    private void execute() {
        Random rng = new Random(777);
        MomentSpecification mySpecification = new GasolineSpecification("gasolineData.csv");

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
         */
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

            ArrayList<computeStuff> cvList = new ArrayList<>();
            for (int minObservationsPerLeaf = 10; minObservationsPerLeaf <= 60; minObservationsPerLeaf += 5) {
                for (double minImprovement = 16; minImprovement <= 16; minImprovement *= 2) {
                    for (int maxDepth = 3; maxDepth <= 4; maxDepth++) {
                        cvList.add(new computeStuff(mySpecification, numberTreesInForest, rngBaseSeedMomentForest, verbose, minObservationsPerLeaf, minImprovement, maxDepth, rngBaseSeedOutOfSample, false));
                    }
                }
            }

            cvList.parallelStream().forEach(s -> {
                s.computeOutOfSampleMSEInParameterSpace();
            });

            for (computeStuff s : cvList) {
                double combinationMSE = s.getMSE();
                String star = "";
                if (combinationMSE <= lowestSSE || first) {
                    lowestSSE = combinationMSE;
                    first = false;
                    bestMinImprovement = s.getMinImprovement();
                    bestMinObservationsPerLeaf = s.getMinObservationsPerLeaf();
                    bestMaxDepth = s.getMaxTreeDepth();
                    star = "(*)";
                }
                System.out.println("minMSE: " + s.getMinImprovement() + " minObs: " + s.getMinObservationsPerLeaf() + " maxDepth: " + s.getMaxTreeDepth() + " Out-of-sample MSE: " + combinationMSE + " " + star);
                jt.append("minMSE: " + s.getMinImprovement() + " minObs: " + s.getMinObservationsPerLeaf() + " maxDepth: " + s.getMaxTreeDepth() + " Out-of-sample MSE: " + combinationMSE + " " + star + "\n");
            }

            System.out.println("Lowest MSE: " + lowestSSE + " at min_N = " + bestMinObservationsPerLeaf + " min_MSE = " + bestMinImprovement + " maxDepth: " + bestMaxDepth);
            jt.append("Lowest MSE: " + lowestSSE + " at min_N = " + bestMinObservationsPerLeaf + " min_MSE = " + bestMinImprovement + " maxDepth: " + bestMaxDepth + "\n");
        } else {
            bestMinObservationsPerLeaf = 55;
            bestMinImprovement = 16.0;
            bestMaxDepth = 3;
        }

        bestMaxDepth = 1;

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
            MomentForest myForest = new MomentForest(mySpecification, 1, rngBaseSeedMomentForest, forestLens, verbose, new TreeOptions());

            myForest.setTreeOptions(cvOptions);
            myForest.growForest();

            TreeMoment loblolly = myForest.getTree(0);
            loblolly.testHomogeneity();

            // System.out.println("Done with growforest");
            ArrayList<Integer> hpl = myForest.getTree(0).getIndexHomogeneousParameters(); // this is only using the first tree, is that the right way of thinking about this?
            ArrayList<Double> hplStartingValues = myForest.getTree(0).getValueHomogeneousParameters();
            // System.out.println("Post get homogeneous parameters");

            // need to sort together
            HomogeneousParameterSorter sorter = new HomogeneousParameterSorter();
            sorter.sort(hpl, hplStartingValues);

            // Collections.sort(hpl); // ensure that indices are ascending (this can cause some weird problems elsewhere due to my bad coding skills if not)
            // Collections.sort(hplStartingValues);
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
                mySpecification.setHomogeneousIndex(1);
                mySpecification.setHomogeneousParameter(1, 1.0);
                // mySpecification.setHomogeneousIndex(1);
                // mySpecification.setHomogeneousParameter(1, 1.0);
                setEstimatedHomogeneousParameters(mySpecification.getHomogeneousParameterVector());
            } else {
                if (!hpl.isEmpty()) {
                    // System.out.println("Initializing search container");
                    numberTreesInForest = 10;
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
        verbose = false;
        numberTreesInForest = 48;
        computeStuff pain = new computeStuff(mySpecification, numberTreesInForest, rngBaseSeedMomentForest, verbose, bestMinObservationsPerLeaf,
                bestMinImprovement, bestMaxDepth, rngBaseSeedOutOfSample, true);
        pain.computeOutOfSampleMSEInParameterSpace();
        double outOfSampleFit = pain.getMSE();

        System.out.println("Out of sample SSE: " + outOfSampleFit);
    }

    private class computeStuff {

        MomentSpecification mySpecification;
        int numberTreesInForest;
        long rngBaseSeedMomentForest;
        boolean verbose;
        int minObservationsPerLeaf;
        double minImprovement;
        int maxTreeDepth;
        long rngBaseSeedOutOfSample;
        boolean generatePlots;

        double MSE;

        public computeStuff(MomentSpecification mySpecification, int numberTreesInForest, long rngBaseSeedMomentForest, boolean verbose, int minObservationsPerLeaf, double minImprovement, int maxTreeDepth, long rngBaseSeedOutOfSample, boolean generatePlots) {
            this.mySpecification = mySpecification;
            this.numberTreesInForest = numberTreesInForest;
            this.rngBaseSeedMomentForest = rngBaseSeedMomentForest;
            this.verbose = verbose;
            this.minObservationsPerLeaf = minObservationsPerLeaf;
            this.minImprovement = minImprovement;
            this.maxTreeDepth = maxTreeDepth;
            this.rngBaseSeedOutOfSample = rngBaseSeedOutOfSample;
            this.generatePlots = generatePlots;
        }

        public double getMSE() {
            return MSE;
        }

        public int getMaxTreeDepth() {
            return maxTreeDepth;
        }

        public double getMinImprovement() {
            return minImprovement;
        }

        public int getMinObservationsPerLeaf() {
            return minObservationsPerLeaf;
        }

        public void computeOutOfSampleMSEInParameterSpace() {
            // System.out.println("\nComputing OOS In Parameter Space\n");
            // System.out.println("Homogeneous parameter length in spec: "+mySpecification.getHomogeneousIndex().length);

            MomentForest myForest;

            DataLens homogenizedForestLens = new DataLens(mySpecification.getX(), mySpecification.getY(), mySpecification.getZ(), null);

            myForest = new MomentForest(mySpecification, numberTreesInForest, rngBaseSeedMomentForest, homogenizedForestLens, verbose, new TreeOptions());
            TreeOptions cvOptions = new TreeOptions(0.01, minObservationsPerLeaf, minImprovement, maxTreeDepth, false); // k = 1
            myForest.setTreeOptions(cvOptions);
            /**
             * Grow the moment forest
             */
            myForest.growForest();

            myForest.getTree(0).printTree();
            /**
             * Test vectors for assessment
             */
            DataLens oosDataLens = mySpecification.getOutOfSampleXYZ(2000, rngBaseSeedOutOfSample); // this should eventually be modified to come out of the data itself (or generalized in some way)
            Jama.Matrix testZ = oosDataLens.getZ();
            Jama.Matrix testX = oosDataLens.getX();
            Jama.Matrix testY = oosDataLens.getY();

            double outOfSampleFit = 0;
            for (int i = 0; i < testZ.getRowDimension(); i++) {
                Jama.Matrix zi = testZ.getMatrix(i, i, 0, testZ.getColumnDimension() - 1);
                Jama.Matrix xi = testX.getMatrix(i, i, 0, testX.getColumnDimension() - 1);
                double yi = testY.get(i, 0);

                // have to reconstruct a composite beta from homogeneous and heterogeneous parameters
                Jama.Matrix compositeEstimatedBeta = myForest.getEstimatedParameterForest(zi);
                outOfSampleFit += mySpecification.getGoodnessOfFit(yi, xi, compositeEstimatedBeta);
            }

            if (generatePlots) {
                /**
                 * This would be a fine place to make the graph ala Schmalensee
                 * and Stoker
                 *
                 * Nonparametric plot of log total gallons (Y axis) against log
                 * age (X axis)
                 *
                 */

                Jama.Matrix Z = mySpecification.getZ();
                Jama.Matrix X = mySpecification.getX();
                Jama.Matrix Y = mySpecification.getY();
                boolean useOutOfSample = !true;
                if (useOutOfSample) {
                    Z = testZ;
                    X = testX;
                    Y = testY;
                }

                XYSeriesCollection xyc = new XYSeriesCollection();
                XYSeriesCollection xycData = new XYSeriesCollection();
                XYSeriesCollection kernelAge = new XYSeriesCollection();

                NormalDistribution normal = new NormalDistribution(0, 0.05);

                int[] incomeBinArray = {1, 4, 6, 9};
                String[] incomeString = {"Less Than $10,000", "$25,000 to $34,999", "$50,000 to $74,999", "$100,000 to $149,999"};

                // for (int incomeBin : incomeBinArray) {
                for (int jk = 0; jk < incomeBinArray.length; jk++) {
                    int incomeBin = incomeBinArray[jk];
                    XYSeries xy = new XYSeries(incomeString[jk]);
                    XYSeries xy2 = new XYSeries(incomeString[jk]);
                    XYSeries xyKernel = new XYSeries(incomeString[jk]);

                    Random rng = new Random(rngBaseSeedOutOfSample);
                    for (double age = 16; age <= 92; age++) {
                        double avgLogGallonsFit = 0;
                        double avgLogGallonsData = 0;
                        double sumNPWeights = 0;
                        double counter = 0;
                        for (int i = 0; i < Z.getRowDimension(); i++) {
                            Jama.Matrix zi = Z.getMatrix(i, i, 0, Z.getColumnDimension() - 1);
                            Jama.Matrix xi = X.getMatrix(i, i, 0, X.getColumnDimension() - 1);

                            // try a smoothed estimator of the average in this area
                            if ((int) zi.get(0, 5) == incomeBin) {
                                double weight = normal.probability(zi.get(0, 3) - Math.log(age));
                                avgLogGallonsData += Y.get(i, 0) * weight;
                                sumNPWeights += weight;

                                // we are going to integrate out all the other characteristics, just setting age
                                zi.set(0, 3, Math.log(age));
                                Jama.Matrix compositeEstimatedBeta = myForest.getEstimatedParameterForest(zi);
                                avgLogGallonsFit += mySpecification.getPredictedY(xi, compositeEstimatedBeta, rng);
                                counter++;
                            }
                        }
                        avgLogGallonsFit /= counter;
                        avgLogGallonsData /= sumNPWeights;
                        xy.add(Math.log(age), avgLogGallonsFit);
                        xy2.add(Math.log(age), avgLogGallonsData);
                        xyKernel.add(Math.log(age), sumNPWeights / counter);
                        // System.out.println("age: "+age+" data: "+avgLogGallonsData+" fit: "+avgLogGallonsFit);
                    }
                    xyc.addSeries(xy);
                    xycData.addSeries(xy2);
                    kernelAge.addSeries(xyKernel);
                }
                ChartGenerator.makeXYLine(xyc, "Fitted Gasoline Consumption", "Log Age", "Log Gallons Gasoline");
                ChartGenerator.makeXYLine(xycData, "Actual Gasoline Consumption", "Log Age", "Log Gallons Gasoline");
                ChartGenerator.makeXYLine(kernelAge, "Kernel Density", "Log Age", "Log Age Density");
            }

            MSE = outOfSampleFit / testZ.getRowDimension(); // mse
        }
    }

    private void setHomogeneousParameterList(ArrayList<Integer> homogeneousParameterList) {
        this.homogeneousParameterList = homogeneousParameterList;
    }

    public ArrayList<Integer> getHomogeneousParameterList() {
        return homogeneousParameterList;
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
