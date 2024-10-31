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
package examples.CardIV;

import Jama.Matrix;
import core.DataLens;
import core.HomogeneousSearchContainer;
import core.MomentForest;
import core.TreeOptions;
import examples.PartialLinearGasDemand.HomogeneousParameterSorter;
import java.awt.GridLayout;
import java.util.ArrayList;
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
public class MainCardIV {

    private ArrayList<Integer> homogeneousParameterList;
    private Jama.Matrix estimatedHomogeneousParameters;
    private final boolean detectHomogeneity;
    final private JTextArea jt;

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        JFrame f = new JFrame("Card IV Partial Linear Application");
        f.setBounds(100, 100, 1500, 500);
        f.getContentPane().setLayout(new GridLayout(1, 2));
        JTextAreaAutoscroll jt1 = new JTextAreaAutoscroll();
        f.getContentPane().add(new JScrollPane(jt1));
        JTextAreaAutoscroll jt2 = new JTextAreaAutoscroll();
        f.getContentPane().add(new JScrollPane(jt2));
        f.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        f.setVisible(true);

        MainCardIV go = new MainCardIV(jt1);
        go.execute();

    }

    public MainCardIV(JTextArea jt) {
        this.jt = jt;
        detectHomogeneity = true;
    }

    private void execute() {
        Random rng = new Random(777);
        // MomentSpecificationCardIV mySpecification = new MomentSpecificationCardIV("d:/git/momentforests/java/examples/cardiv/IV test.csv");
        MomentSpecificationCardIV mySpecification = new MomentSpecificationCardIV("d:/git/momentforests/java/examples/cardiv/table3.csv");

        double bestMinImprovement = 4.0;
        int bestMinObservationsPerLeaf = 25;
        int bestMaxDepth = 5;

        double minOutOfSampleFit = 0;
        double minInSampleFit = 0;
        boolean first = true;

        /**
         * Make sure that cross-validation is run on completely unrestricted
         * model; set all parameters to heterogeneous
         */
        mySpecification.resetHomogeneityIndex();

        int numberTreesInForest = 1;
        // System.out.println("numTrees: " + numberTreesInForest);

        /*
         * Initialize the moment forest
         */
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

            // NEED TO UPDATE
            ArrayList<computeFitStatistics> cvList = new ArrayList<>();
            for (int minObservationsPerLeaf = 50; minObservationsPerLeaf <= 50; minObservationsPerLeaf *= 2) {
                for (double minImprovement = 1.0; minImprovement <= 1.0; minImprovement *= 5) {
                    for (int maxDepth = 1; maxDepth >= 0; maxDepth--) {
                        cvList.add(new computeFitStatistics(mySpecification, numberTreesInForest, rngBaseSeedMomentForest, verbose, minObservationsPerLeaf, minImprovement, maxDepth, rngBaseSeedOutOfSample, false));
                    }
                }
            }
            cvList.parallelStream().forEach(s -> {
                s.computeOutOfSampleMSE();
            });

            for (computeFitStatistics s : cvList) {
                double combinationMSE = s.getMSE();
                String star = "";
                String starIn = "";
                if (combinationMSE <= minOutOfSampleFit || first) {
                    minOutOfSampleFit = combinationMSE;
                    bestMinImprovement = s.getMinImprovement();
                    bestMinObservationsPerLeaf = s.getMinObservationsPerLeaf();
                    bestMaxDepth = s.getMaxTreeDepth();
                    star = "(*)";
                }
                System.out.print("best: " + minInSampleFit + " this: " + s.getInSampleFit() + " " + (minInSampleFit > s.getInSampleFit()));
                if (s.getInSampleFit() < minInSampleFit || first) {
                    System.out.println("detected lower");
                    minInSampleFit = s.getInSampleFit();
                    starIn = "(**)";
                }

                if (first) {
                    first = false;
                }
                System.out.println(starIn);
                System.out.println("minMSE: " + s.getMinImprovement() + " minObs: " + s.getMinObservationsPerLeaf() + " maxDepth: " + s.getMaxTreeDepth() + " Out-of-sample MSE: " + combinationMSE + " " + star);
                jt.append("minMSE: " + s.getMinImprovement() + " minObs: " + s.getMinObservationsPerLeaf() + " maxDepth: " + s.getMaxTreeDepth() + " Out-of-sample MSE: " + combinationMSE + " " + star + " In-Sample Fit: " + s.getInSampleFit() + " " + starIn + "\n");
            }

            System.out.println("Lowest MSE: " + minOutOfSampleFit + " at min_N = " + bestMinObservationsPerLeaf + " min_MSE = " + bestMinImprovement + " maxDepth: " + bestMaxDepth);
            jt.append("Lowest MSE: " + minOutOfSampleFit + " at min_N = " + bestMinObservationsPerLeaf + " min_MSE = " + bestMinImprovement + " maxDepth: " + bestMaxDepth + "\n");
            jt.append("Best in-sample fit: " + minInSampleFit + "\n");
        } else {
            // this is what CV gave me Oct 24 2024
            bestMinObservationsPerLeaf = 50;
            bestMinImprovement = 5.0;
            bestMaxDepth = 0; // 3;
        }

        mySpecification.resetHomogeneityIndex();
        if (detectHomogeneity && 1 == 1) {
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
            myForest.testHomogeneity();

            ArrayList<Integer> hpl = new ArrayList<>();
            ArrayList<Double> hplStartingValues = new ArrayList<>();
            boolean[] voteIndexHomogeneity = myForest.getHomogeneityVotes(jt, true);
            double[] startingValues = myForest.getHomogeneityStartingValues();
            for (int i = 0; i < voteIndexHomogeneity.length; i++) {
                if (voteIndexHomogeneity[i]) {
                    System.out.println("Adding index " + i + " to homogeneous list with starting value: " + startingValues[i]);
                    hpl.add(i);
                    hplStartingValues.add(startingValues[i]);
                }
            }

            HomogeneousParameterSorter sorter = new HomogeneousParameterSorter();
            sorter.sort(hpl, hplStartingValues);
            mySpecification.resetHomogeneityIndex();
            for (int i = 0; i < hpl.size(); i++) {
                mySpecification.setHomogeneousIndex(hpl.get(i));
                mySpecification.setHomogeneousParameter(hpl.get(i), hplStartingValues.get(i));
            }
            setHomogeneousParameterList(hpl);
            // System.out.println("After setting hp list");
            // this seems to be working
            
            /**
             * OK, key insight: if all the parameters are homogeneous, don't have to do anything
             * Just sent max depth to 0 and then grow the forest as usual.
             */
            boolean allParametersHomogeneous = false;
            if(hpl.size()==mySpecification.getNumParams()) {
                allParametersHomogeneous = true;
                bestMaxDepth = 0;
            }

            /*
             * Estimate values of those homogeneous parameters
             */
            if (!hpl.isEmpty() && !allParametersHomogeneous) {
                System.out.println("Initializing search container");
                numberTreesInForest = 1; // 10
                HomogeneousSearchContainer con = new HomogeneousSearchContainer(mySpecification, numberTreesInForest, verbose, bestMinImprovement, bestMinObservationsPerLeaf, bestMaxDepth,
                        getHomogeneousParameterList(), rngBaseSeedMomentForest, rngBaseSeedOutOfSample);
                System.out.println("Calling execute search");
                con.executeSearch();
                System.out.println("Post search");
                Jama.Matrix homogeneousParameters = con.getEstimatedHomogeneousParameters();
                System.out.print("Post-HomogeneousSearchContainer Estimated homogeneous parameters: ");
                pmUtility.prettyPrintVector(homogeneousParameters); // this is a compact vector of parameters

                int K = mySpecification.getHomogeneousParameterVector().getRowDimension();
                System.out.print("Post-HomogeneousSearchContainer Length of homogeneous parameter vector: " + K);
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
        verbose = false;
        numberTreesInForest = 1000;

        computeFitStatistics fitStats = new computeFitStatistics(mySpecification, numberTreesInForest, rngBaseSeedMomentForest, verbose, bestMinObservationsPerLeaf,
                bestMinImprovement, bestMaxDepth, rngBaseSeedOutOfSample, false);
        fitStats.computeOutOfSampleMSE();
        double outOfSampleFit = fitStats.getMSE();

        System.out.println("Best minimum observations in leaf: " + bestMinObservationsPerLeaf);
        System.out.println("Best minimum improvement: " + bestMinImprovement);
        System.out.println("Best maximum depth: " + bestMaxDepth);

        System.out.println("Out of sample SSE: " + outOfSampleFit);
        System.out.println("Number of valid trees in forest: " + fitStats.myForest.getForestSize());

        double[] countVariableSplitsInForest = fitStats.getSplitVariables();
        System.out.println("Number of split variables: " + countVariableSplitsInForest.length);

        System.out.println("Forest split on the following variables:");
        for (int i = 0; i < countVariableSplitsInForest.length; i++) {
            if (countVariableSplitsInForest[i] > 0) {
                System.out.format("%20s [%.2f%%] %n", mySpecification.getVariableName(i), 100.0 * countVariableSplitsInForest[i] / numberTreesInForest);
            }
        }
        
    }

    private class computeFitStatistics {

        MomentSpecificationCardIV mySpecification;
        int numberTreesInForest;
        long rngBaseSeedMomentForest;
        boolean verbose;
        int minObservationsPerLeaf;
        double minImprovement;
        int maxTreeDepth;
        long rngBaseSeedOutOfSample;
        boolean generatePlots;

        double MSE;
        double inSampleFit;

        MomentForest myForest;

        public computeFitStatistics(MomentSpecificationCardIV mySpecification, int numberTreesInForest, long rngBaseSeedMomentForest, boolean verbose, int minObservationsPerLeaf, double minImprovement, int maxTreeDepth, long rngBaseSeedOutOfSample, boolean generatePlots) {
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

        public double[] getSplitVariables() {
            return myForest.getNumberTimesTreesInForestSplitOnAGivenVariableIndex();
        }

        public void computeOutOfSampleMSE() {
            // System.out.println("\nComputing OOS In Parameter Space\n");
            // System.out.println("Homogeneous parameter length in spec: "+mySpecification.getHomogeneousIndex().length);

            DataLens overallLens = new DataLens(mySpecification.getX(), mySpecification.getY(), mySpecification.getZ(), null);
            DataLens[] split = overallLens.randomlySplitSample(0.9, 383);
            DataLens estimatingLens = split[0];
            DataLens oosDataLens = split[1];
            
            // System.out.println(estimatingLens.getNumObs()+" "+oosDataLens.getNumObs());

            myForest = new MomentForest(mySpecification, numberTreesInForest, rngBaseSeedMomentForest, estimatingLens, verbose, new TreeOptions());

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

            Jama.Matrix testZ = oosDataLens.getZ();
            Jama.Matrix testY = oosDataLens.getY();
            Jama.Matrix testX = oosDataLens.getX();

            double outOfSampleFit = 0;
            
            for (int i = 0; i < testZ.getRowDimension(); i++) {
                Jama.Matrix zi = testZ.getMatrix(i, i, 0, testZ.getColumnDimension() - 1);
                Jama.Matrix xi = testX.getMatrix(i, i, 0, testX.getColumnDimension() - 1);
                double yi = testY.get(i, 0);

                // have to reconstruct a composite beta from homogeneous and heterogeneous parameters
                Jama.Matrix compositeEstimatedBeta = myForest.getEstimatedParameterForest(zi);
                outOfSampleFit += mySpecification.getGoodnessOfFit(yi, xi, compositeEstimatedBeta);               
                // I think I just need to evaluate the moments at each observation for the beta that the forest has given
                if(i<10) {
                    System.out.println(pmUtility.stringPrettyPrintVector(compositeEstimatedBeta)+"\t"+pmUtility.stringPrettyPrint(zi));
                }
            }
            
            inSampleFit = 0;

            // System.out.println(mySpecification.getNumMoments());

            // System.out.println("Starting in-sample Z");
            for (int i = 0; i < estimatingLens.getZ().getRowDimension(); i++) {
                Jama.Matrix zi = estimatingLens.getZ().getMatrix(i, i, 0, estimatingLens.getZ().getColumnDimension() - 1);
                Jama.Matrix xi = estimatingLens.getX().getMatrix(i, i, 0, estimatingLens.getX().getColumnDimension() - 1);
                double yi = estimatingLens.getY().get(i, 0);
                // have to reconstruct a composite beta from homogeneous and heterogeneous parameters
                Jama.Matrix compositeEstimatedBeta = myForest.getEstimatedParameterForest(zi);
                
                if(i==0 && 1==1) {
                    for(int f=0;f<myForest.getForestSize();f++) {
                        pmUtility.prettyPrintVector(myForest.getTree(f).getEstimatedBeta(zi));
                    }
                }

                inSampleFit += mySpecification.getGoodnessOfFit(yi, xi, compositeEstimatedBeta);
                // inMoments.plusEquals(civIn.getGi(compositeEstimatedBeta, i));
            }
            // inMoments.timesEquals(1.0 / mySpecification.getZ().getRowDimension());
            // inSampleFit = ((inMoments.transpose()).times(inMoments)).get(0, 0);
            inSampleFit /= mySpecification.getZ().getRowDimension();
            System.out.println("In-sample SSE: "+inSampleFit);

            MSE = outOfSampleFit / testZ.getRowDimension(); // mse
        }

        private double getInSampleFit() {
            return inSampleFit;
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
