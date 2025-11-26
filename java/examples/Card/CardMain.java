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
package examples.Card;

import Jama.Matrix;
import core.DataLens;
import core.HomogeneousParameterSorter;
import core.HomogeneousSearchContainer;
import core.MomentForest;
import core.MomentSpecification;
import core.TreeMoment;
import core.TreeOptions;
import java.awt.BorderLayout;

import java.awt.GridLayout;
import java.io.BufferedWriter;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.Random;
import javax.swing.JFrame;
import javax.swing.JScrollPane;
import javax.swing.JTextArea;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;
import utility.JTextAreaAutoscroll;
import utility.pmUtility;

/**
 *
 * @author Stephen P. Ryan <stephen.p.ryan@wustl.edu>
 */
public class CardMain {

    private ArrayList<Integer> homogeneousParameterList;
    private Jama.Matrix estimatedHomogeneousParameters;
    private final boolean detectHomogeneity;
    final private JTextArea jt;

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        JFrame f = new JFrame("Card Partial Linear Model Application");
        f.setBounds(100, 100, 1500, 500);
        f.getContentPane().setLayout(new GridLayout(1, 2));
        JTextAreaAutoscroll jt1 = new JTextAreaAutoscroll();
        f.getContentPane().add(new JScrollPane(jt1));
        JTextAreaAutoscroll jt2 = new JTextAreaAutoscroll();
        f.getContentPane().add(new JScrollPane(jt2));
        f.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        f.setVisible(true);

        CardMain go = new CardMain(jt1);
        go.execute();

    }

    public CardMain(JTextArea jt) {
        this.jt = jt;
        detectHomogeneity = true;
    }

    private void execute() {
        Random rng = new Random(777);
        MomentSpecification mySpecification = new CardSpecification("C:/Users/nathan/Documents/GitHub/momentForests/java/examples/Card/table2.csv");

        double bestMinImprovement = 4.0;
        int bestMinObservationsPerLeaf = 25;
        int bestMaxDepth = 0;

        double minOutOfSampleFit = 0;
        double minInSampleFit = 0;
        boolean first = true;

        /**
         * Make sure that cross-validation is run on completely unrestricted
         * model; set all parameters to heterogeneous
         */
        mySpecification.resetHomogeneityIndex();

        int numberTreesInForest = 40;
        // System.out.println("numTrees: " + numberTreesInForest);

        /*
         * Initialize the moment forest
         */
        boolean verbose = true;
        boolean testParameterHomogeneity;

        long rngBaseSeedMomentForest = rng.nextLong();
        long rngBaseSeedOutOfSample = rng.nextLong();

        boolean runCV = true;
        if (runCV) {
            if (verbose) {
                System.out.println("************************");
                System.out.println("* Run Cross-Validation *");
                System.out.println("************************");
            }

            // NEED TO UPDATE
            ArrayList<computeFitStatistics> cvList = new ArrayList<>();
            for (int minObservationsPerLeaf = 50; minObservationsPerLeaf <= 200; minObservationsPerLeaf *= 2) {
                for (double minImprovement = 0.1; minImprovement <= 10; minImprovement *= 10) {
                    for (int maxDepth = 7; maxDepth >= 5; maxDepth--) {
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
            bestMinObservationsPerLeaf = 25;
            bestMinImprovement = 0.1;
            bestMaxDepth = 5;
        }

        mySpecification.resetHomogeneityIndex();
        if (detectHomogeneity && 1 == 1 && bestMaxDepth > 0) {
            if (verbose) {
                System.out.println("************************");
                System.out.println("* Test for Homogeneity *");
                System.out.println("************************");
            }
            /**
             * Step 2: determine homogeneous parameters post-CV
             */
            double minProportionInEachLeaf = 0.01;
            
            DataLens forestLens = new DataLens(mySpecification.getX(), mySpecification.getY(), mySpecification.getZ(), null, mySpecification.getStratificationIndex());

            testParameterHomogeneity = true;
            TreeOptions cvOptions = new TreeOptions(minProportionInEachLeaf, bestMinObservationsPerLeaf, bestMinImprovement, bestMaxDepth, testParameterHomogeneity); // k = 1
            MomentForest myForest = new MomentForest(mySpecification, numberTreesInForest, rngBaseSeedMomentForest, forestLens, verbose, new TreeOptions());

            myForest.setTreeOptions(cvOptions);
            myForest.growForest();

            if (verbose) {
                TreeMoment loblolly = myForest.getTree(0);
                System.out.println("***************************************************");
                System.out.println("* Before homogeneity testing: printing first tree *");
                System.out.println("***************************************************");
                loblolly.printTree();
            }

            myForest.testHomogeneity(verbose);

            // loblolly.testHomogeneity(verbose);
            // should I implement a voting scheme here across all the trees for discerning homogeneity?
            // idea is to take the majority vote across trees
            // then take the average value from those trees that voted yes
            // let's try it and see what happens to classification rates
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
            boolean testPostClassificationConvergence = !true;
            if (testPostClassificationConvergence) {
                hpl.clear();
                hpl.add(0);
                mySpecification.resetHomogeneityIndex();
                mySpecification.setHomogeneousIndex(0);
                mySpecification.setHomogeneousParameter(0, -1.0);
            }


            /*
             * Estimate values of those homogeneous parameters
             */
            if (!hpl.isEmpty()) {
                if (hpl.size() == mySpecification.getNumParams()) {
                    // all homogenous, don't need to optimize, just set to stump and let it run
                    bestMaxDepth = 0;
                } else {
                    System.out.println("Initializing search container");
                    numberTreesInForest = 40; // 10
                    HomogeneousSearchContainer con = new HomogeneousSearchContainer(mySpecification, numberTreesInForest, verbose, bestMinImprovement, bestMinObservationsPerLeaf, bestMaxDepth,
                            getHomogeneousParameterList(), rngBaseSeedMomentForest, rngBaseSeedOutOfSample);
                    System.out.println("Calling execute search");
                    con.executeSearch();
                    System.out.println("Post search");
                    Jama.Matrix homogeneousParameters = con.getEstimatedHomogeneousParameters();
                    System.out.print("Post-HomogeneousSearchContainer Estimated homogeneous parameters: ");
                    pmUtility.prettyPrintVector(homogeneousParameters); // this is a compact vector of parameters

                    // pmUtility.prettyPrintVector(mySpecification.getHomogeneousParameterVector());
                    // holy smokes we aren't setting the homogeneous parameters in mySpec to be what we found here???
                    // doing so below, why wasn't this done before?!?!?
                    int K = mySpecification.getHomogeneousParameterVector().getRowDimension();
                    System.out.println("Post-HomogeneousSearchContainer Length of homogeneous parameter vector: " + K);
                    Jama.Matrix expandedHomogeneousParameterVector = new Jama.Matrix(K, 1);
                    int counter = 0;
                    for (int k = 0; k < K; k++) {
                        if (mySpecification.getHomogeneousIndex()[k]) {
                            expandedHomogeneousParameterVector.set(k, 0, homogeneousParameters.get(counter, 0));
                            mySpecification.setHomogeneousParameter(k, homogeneousParameters.get(counter, 0));
                            counter++;
                        }
                    }
                    System.out.println("Specification homogeneous parameter vector: ");
                    pmUtility.prettyPrintVector(mySpecification.getHomogeneousParameterVector());
                    // System.exit(0);
                    setEstimatedHomogeneousParameters(expandedHomogeneousParameterVector);
                }
            }
        }

        /**
         * Compute out-of-sample measures of fit (against Y, and true beta)
         */
        verbose = true;
        numberTreesInForest = 40; // 50
        computeFitStatistics fitStats = new computeFitStatistics(mySpecification, numberTreesInForest, rngBaseSeedMomentForest, verbose, bestMinObservationsPerLeaf,
                bestMinImprovement, bestMaxDepth, rngBaseSeedOutOfSample, true);
        fitStats.computeOutOfSampleMSE();
        double outOfSampleFit = fitStats.getMSE();

        System.out.println("Best minimum observations in leaf: " + bestMinObservationsPerLeaf);
        System.out.println("Best minimum improvement: " + bestMinImprovement);
        System.out.println("Best maximum depth: " + bestMaxDepth);

        System.out.println("Out of sample SSE: " + outOfSampleFit);
        System.out.println("Number of trees in forest: " + numberTreesInForest);

        double[] countVariableSplitsInForest = fitStats.getSplitVariables();
        System.out.println("Number of split variables: " + countVariableSplitsInForest.length);

        System.out.println("Forest split on the following variables:");
        for (int i = 0; i < countVariableSplitsInForest.length; i++) {
            if (countVariableSplitsInForest[i] > 0) {
                System.out.format("%20s [%.2f%%] %n", mySpecification.getVariableName(i), 100.0 * countVariableSplitsInForest[i] / numberTreesInForest);
            }
        }
        
        fitStats.outputInSampleFits();
    }

    private class computeFitStatistics {

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
        double inSampleFit;

        MomentForest myForest;

        public computeFitStatistics(MomentSpecification mySpecification, int numberTreesInForest, long rngBaseSeedMomentForest, boolean verbose, int minObservationsPerLeaf, double minImprovement, int maxTreeDepth, long rngBaseSeedOutOfSample, boolean generatePlots) {
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
            DataLens overallLens = new DataLens(mySpecification.getX(), mySpecification.getY(), mySpecification.getZ(), null, mySpecification.getStratificationIndex());
            //DataLens[] split = overallLens.randomlySplitSample(0.9, 383);
            DataLens[] split = overallLens.randomlySplitSampleByStrata(0.9, rngBaseSeedMomentForest);
            DataLens estimatingLens = split[0];
            DataLens oosDataLens = split[1];

            myForest = new MomentForest(mySpecification, numberTreesInForest, rngBaseSeedMomentForest, estimatingLens, verbose, new TreeOptions());
            TreeOptions cvOptions = new TreeOptions(0.01, minObservationsPerLeaf, minImprovement, maxTreeDepth, false); // k = 1
            myForest.setTreeOptions(cvOptions);
            /**
             * Grow the moment forest
             */
            myForest.growForest();
            if (verbose) {
                System.out.println("First tree in forest estimated as:");
                myForest.getTree(0).printTree();
            }
            /**
             * Test vectors for assessment
             */
            // DataLens oosDataLens = mySpecification.getOutOfSampleXYZ(1500, rngBaseSeedOutOfSample); // this should eventually be modified to come out of the data itself (or generalized in some way)
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

            inSampleFit = 0;
            for (int i = 0; i < mySpecification.getZ().getRowDimension(); i++) {
                Jama.Matrix zi = mySpecification.getZ().getMatrix(i, i, 0, mySpecification.getZ().getColumnDimension() - 1);
                Jama.Matrix xi = mySpecification.getX().getMatrix(i, i, 0, mySpecification.getX().getColumnDimension() - 1);
                double yi = mySpecification.getY().get(i, 0);
                // have to reconstruct a composite beta from homogeneous and heterogeneous parameters
                Jama.Matrix compositeEstimatedBeta = myForest.getEstimatedParameterForest(zi);

                if (i == 0 && 1 == 2) {
                    for (int f = 0; f < myForest.getForestSize(); f++) {
                        pmUtility.prettyPrintVector(myForest.getTree(f).getEstimatedBeta(zi));
                    }
                }

                inSampleFit += mySpecification.getGoodnessOfFit(yi, xi, compositeEstimatedBeta);
            }
            inSampleFit /= mySpecification.getZ().getRowDimension();

            MSE = outOfSampleFit / testZ.getRowDimension(); // mse

            if (generatePlots && 1 == 2) {
                JFrame f = new JFrame("Plots");
                f.setBounds(100, 100, 800, 800);
                f.getContentPane().setLayout(new BorderLayout());
                XYSeries expSeriesBlack = new XYSeries("Experience -- Black");
                XYSeries eduSeriesBlack = new XYSeries("Education -- Black");
                XYSeries expSeriesWhite = new XYSeries("Experience -- White");
                XYSeries eduSeriesWhite = new XYSeries("Education -- White");

                XYSeriesCollection expCollection = new XYSeriesCollection();
                for (int experience = 0; experience <= 23; experience++) {
                    double black = 1;
                    double[][] demographics = {{1, 12, experience, 0, black, 1, 0, 1, 0, 12, 12, 0, 0, 0, 1, 0}};
                    Jama.Matrix z_guy = new Jama.Matrix(demographics, 1, mySpecification.getZ().getColumnDimension());
                    // pmUtility.prettyPrint(z_guy);
                    Jama.Matrix x_guy = testX.getMatrix(0, 0, 0, testX.getColumnDimension() - 1);
                    if (x_guy.getColumnDimension() > 1) {
                        x_guy.set(0, 1, 12);
                    }
                    if (x_guy.getColumnDimension() > 2) {
                        x_guy.set(0, 2, experience);
                    }

                    Jama.Matrix compositeEstimatedBetaGuy = myForest.getEstimatedParameterForest(z_guy);
                    double estimatedLogWage = (x_guy.times(compositeEstimatedBetaGuy)).get(0, 0);
                    expSeriesBlack.add(experience, (estimatedLogWage));

                    z_guy.set(0, 4, 0.0); // white
                    compositeEstimatedBetaGuy = myForest.getEstimatedParameterForest(z_guy);
                    estimatedLogWage = (x_guy.times(compositeEstimatedBetaGuy)).get(0, 0);
                    expSeriesWhite.add(experience, (estimatedLogWage));
                }
                expCollection.addSeries(expSeriesBlack);
                expCollection.addSeries(expSeriesWhite);

                XYSeriesCollection eduCollection = new XYSeriesCollection();
                for (int education = 8; education <= 18; education++) {
                    double black = 1;
                    double[][] demographics = {{1, education, 8, 0, black, 1, 0, 1, 0, 12, 12, 0, 0, 0, 1, 0}};
                    Jama.Matrix z_guy = new Jama.Matrix(demographics, 1, mySpecification.getZ().getColumnDimension());
                    // pmUtility.prettyPrint(z_guy);
                    Jama.Matrix x_guy = testX.getMatrix(0, 0, 0, testX.getColumnDimension() - 1);
                    if (x_guy.getColumnDimension() > 1) {
                        x_guy.set(0, 1, education);
                    }
                    if (x_guy.getColumnDimension() > 2) {
                        x_guy.set(0, 2, 8);
                    }
                    Jama.Matrix compositeEstimatedBetaGuy = myForest.getEstimatedParameterForest(z_guy);
                    double estimatedLogWage = (x_guy.times(compositeEstimatedBetaGuy)).get(0, 0);
                    eduSeriesBlack.add(education, (estimatedLogWage));

                    z_guy.set(0, 4, 0.0); // white
                    compositeEstimatedBetaGuy = myForest.getEstimatedParameterForest(z_guy);
                    estimatedLogWage = (x_guy.times(compositeEstimatedBetaGuy)).get(0, 0);
                    eduSeriesWhite.add(education, (estimatedLogWage));
                }
                eduCollection.addSeries(eduSeriesBlack);
                eduCollection.addSeries(eduSeriesWhite);

                JFreeChart chart = ChartFactory.createXYLineChart("Experience and Wages", "Experience", "Wages", expCollection);
                NumberAxis axis = (NumberAxis) chart.getXYPlot().getRangeAxis();
                axis.setAutoRangeIncludesZero(false);
                ChartPanel panel = new ChartPanel(chart);
                f.getContentPane().add(panel, BorderLayout.CENTER);
                f.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
                f.setVisible(true);

                JFrame f2 = new JFrame("Plots");
                f2.setBounds(900, 100, 800, 800);
                f2.getContentPane().setLayout(new BorderLayout());
                JFreeChart chart2 = ChartFactory.createXYLineChart("Education and Wages", "Education", "Wages", eduCollection);
                NumberAxis axis2 = (NumberAxis) chart2.getXYPlot().getRangeAxis();
                axis2.setAutoRangeIncludesZero(false);
                ChartPanel panel2 = new ChartPanel(chart2);
                f2.getContentPane().add(panel2, BorderLayout.CENTER);
                f2.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
                f2.setVisible(true);
            }
        }

        public void outputInSampleFits() {
            boolean verboseInSample = false;
            DataLens overallLens = new DataLens(mySpecification.getX(), mySpecification.getY(), mySpecification.getZ(), null, mySpecification.getStratificationIndex());

            myForest = new MomentForest(mySpecification, numberTreesInForest, rngBaseSeedMomentForest, overallLens, verboseInSample, new TreeOptions());
            TreeOptions cvOptions = new TreeOptions(0.01, minObservationsPerLeaf, minImprovement, maxTreeDepth, false); // k = 1
            myForest.setTreeOptions(cvOptions);
            /**
             * Grow the moment forest
             */
            myForest.growForest();
            if (verbose) {
                System.out.println("First tree in forest estimated as:");
                myForest.getTree(0).printTree();
                System.out.println("Second tree in forest estimated as:");
                myForest.getTree(1).printTree();
            }

            try {
                BufferedWriter out = new BufferedWriter(new FileWriter("estimatedParametersByObservation.csv"));
                inSampleFit = 0;
                for (int i = 0; i < mySpecification.getZ().getRowDimension(); i++) {
                    Jama.Matrix zi = mySpecification.getZ().getMatrix(i, i, 0, mySpecification.getZ().getColumnDimension() - 1);
                    Jama.Matrix xi = mySpecification.getX().getMatrix(i, i, 0, mySpecification.getX().getColumnDimension() - 1);
                    
                    //System.out.println("zi");
                    //pmUtility.prettyPrint(zi);
                                        
                    double yi = mySpecification.getY().get(i, 0);
                    // have to reconstruct a composite beta from homogeneous and heterogeneous parameters
                    Jama.Matrix compositeEstimatedBeta = myForest.getEstimatedParameterForest(zi);
                    
                    // Counterfactual: white
                    Jama.Matrix zi_white = zi.copy();
                    zi_white.set(0, 4, 0.0);      
                    Jama.Matrix compositeEstimatedBeta_white = myForest.getEstimatedParameterForest(zi_white);
                    
                    // Counterfactual: black
                    Jama.Matrix zi_black = zi.copy();
                    zi_black.set(0, 4, 1.0);      
                    Jama.Matrix compositeEstimatedBeta_black = myForest.getEstimatedParameterForest(zi_black);
          
                    // Counterfactual: not South
                    Jama.Matrix zi_not_south = zi.copy();
                    zi_not_south.set(0, 5, 0.0);      
                    Jama.Matrix compositeEstimatedBeta_not_south = myForest.getEstimatedParameterForest(zi_not_south);
                                        
                    // Counterfactual: South
                    Jama.Matrix zi_south = zi.copy();
                    zi_south.set(0, 5, 1.0);      
                    Jama.Matrix compositeEstimatedBeta_south = myForest.getEstimatedParameterForest(zi_south);
                    
                    // Counterfactual: not SMSA
                    Jama.Matrix zi_not_smsa = zi.copy();
                    zi_not_smsa.set(0, 6, 0.0);      
                    Jama.Matrix compositeEstimatedBeta_not_smsa = myForest.getEstimatedParameterForest(zi_not_smsa);
                                        
                    // Counterfactual: SMSA
                    Jama.Matrix zi_smsa = zi.copy();
                    zi_smsa.set(0, 6, 1.0);      
                    Jama.Matrix compositeEstimatedBeta_smsa = myForest.getEstimatedParameterForest(zi_smsa);
                    
                    for(int j=0;j<zi.getColumnDimension();j++) {
                        out.write(zi.get(0,j)+",");
                    }
                    for(int j=0;j<compositeEstimatedBeta.getRowDimension();j++) {
                        out.write(compositeEstimatedBeta.get(j,0)+",");
                    }
                    for(int j=0;j<compositeEstimatedBeta_white.getRowDimension();j++) {
                        out.write(compositeEstimatedBeta_white.get(j,0)+",");
                    }  
                    for(int j=0;j<compositeEstimatedBeta_black.getRowDimension();j++) {
                        out.write(compositeEstimatedBeta_black.get(j,0)+",");
                    }
                    for(int j=0;j<compositeEstimatedBeta_not_south.getRowDimension();j++) {
                        out.write(compositeEstimatedBeta_not_south.get(j,0)+",");
                    }  
                    for(int j=0;j<compositeEstimatedBeta_south.getRowDimension();j++) {
                        out.write(compositeEstimatedBeta_south.get(j,0)+",");
                    }
                    for(int j=0;j<compositeEstimatedBeta_not_smsa.getRowDimension();j++) {
                        out.write(compositeEstimatedBeta_not_smsa.get(j,0)+",");
                    }  
                    for(int j=0;j<compositeEstimatedBeta_smsa.getRowDimension();j++) {
                        out.write(compositeEstimatedBeta_smsa.get(j,0)+",");
                    }
 
                    out.write("\n");
                                  
                    inSampleFit += mySpecification.getGoodnessOfFit(yi, xi, compositeEstimatedBeta);
                }
                inSampleFit /= mySpecification.getZ().getRowDimension();
                System.out.println("In-sample fit: " + inSampleFit);
                out.close();
            } catch (Exception e) {
                e.printStackTrace();
                System.exit(0);
            }

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
