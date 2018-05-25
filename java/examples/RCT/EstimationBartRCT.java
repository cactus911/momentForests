/*
 * The MIT License
 *
 * Copyright 2018 Stephen P. Ryan <stephen.p.ryan@wustl.edu>.
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
package examples.RCT;

import JSci.maths.statistics.NormalDistribution;
import core.TreeMoment;
import java.io.BufferedWriter;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.Random;
import java.util.TreeSet;
import utility.pmUtility;

/**
 *
 * @author stephen.p.ryan
 */
public class EstimationBartRCT {

    public static void main(String[] args) {
        EstimationBartRCT go = new EstimationBartRCT();
    }

    public EstimationBartRCT() {
        MomentSpecificationBartRCT bart = new MomentSpecificationBartRCT();
        bart.loadData();

        int bestK = 0;
        double bestAlpha = 0;
        double bestMSEBar = 0;
        int bestDepth = 0;
        boolean first = true;
        double bestMSPE = 0;
        double bestProportion = 0;

        int numObs = bart.getX().getRowDimension();
        Jama.Matrix allX = bart.getX().copy();

        // try implementing a true OOS prediction error here?
        for (double proportion = 0.5; proportion <= 0.5; proportion += 0.05) {
            int numObsEstimateLeafValues = (int) Math.floor(numObs * proportion);
            int numPredictObs = 10;
            int numObsGrowTreeStructure = numObs - numObsEstimateLeafValues - numPredictObs; // halfObs / 2;
            System.out.println("numObs: " + numObs + " halfObs: " + numObsEstimateLeafValues + " predictObs: " + numPredictObs);

            Jama.Matrix treeX = new Jama.Matrix(numObsGrowTreeStructure, bart.getX().getColumnDimension());
            Jama.Matrix treeY = new Jama.Matrix(numObsGrowTreeStructure, 1);
            Jama.Matrix honestX = new Jama.Matrix(numObsEstimateLeafValues, bart.getX().getColumnDimension());
            Jama.Matrix honestY = new Jama.Matrix(numObsEstimateLeafValues, 1);
            Jama.Matrix predictX = new Jama.Matrix(numPredictObs, bart.getX().getColumnDimension());
            Jama.Matrix predictY = new Jama.Matrix(numPredictObs, 1);

            TreeSet<Integer> growTreeSet = new TreeSet<>();
            int count = 0;
            while (count < numObsGrowTreeStructure) {
                if (count < numObsGrowTreeStructure) {
                    int index = (int) Math.floor(Math.random() * numObs);
                    if (!growTreeSet.contains(index)) {
                        growTreeSet.add(index);
                        count++;
                    }
                }
            }
            TreeSet<Integer> predictSet = new TreeSet<>();
            count = 0;
            while (count < numPredictObs) {
                if (count < numPredictObs) {
                    int index = (int) Math.floor(Math.random() * numObs);
                    if (!predictSet.contains(index) && !growTreeSet.contains(index)) {
                        predictSet.add(index);
                        count++;
                    }
                }
            }
            System.out.println("honestX row: " + honestX.getRowDimension() + " col: " + honestX.getColumnDimension());
            System.out.println("treeSet.size(): " + growTreeSet.size());

            int countTree = 0;
            int countPredict = 0;
            int countHonest = 0;
            for (int i = 0; i < numObs; i++) {
                if (growTreeSet.contains(i)) {
                    for (int j = 0; j < bart.getX().getColumnDimension(); j++) {
                        treeX.set(countTree, j, bart.getX().get(i, j));
                    }
                    treeY.set(countTree, 0, bart.getY().get(i, 0));
                    countTree++;
                } else if (predictSet.contains(i)) {
                    for (int j = 0; j < bart.getX().getColumnDimension(); j++) {
                        predictX.set(countPredict, j, bart.getX().get(i, j));
                    }
                    predictY.set(countPredict, 0, bart.getY().get(i, 0));
                    countPredict++;
                } else {
                    for (int j = 0; j < bart.getX().getColumnDimension(); j++) {
                        honestX.set(countHonest, j, bart.getX().get(i, j));
                    }
                    honestY.set(countHonest, 0, bart.getY().get(i, 0));
                    countHonest++;
                }
                // System.out.println(i+" "+treeSet.contains(i)+" "+countTree+" "+countHonest);
            }

            double minProportionEachPartition = 0.01;
            int minCountEachPartition = 25;
            double improvementThreshold = 1E-3;
            int maxDepth = 100;
            boolean verbose = true;

            // for (maxDepth = 1; maxDepth <= 20; maxDepth++) {
            for (improvementThreshold = 1E-7; improvementThreshold <= 1E-7; improvementThreshold *= 10) {
                for (minCountEachPartition = 25; minCountEachPartition <= 25; minCountEachPartition += 5) {
                    System.out.println("MSE_bar: " + improvementThreshold + " k: " + minCountEachPartition);
                    TreeMoment momentTree = new TreeMoment(null, bart, treeX, treeY, bart.getDiscreteVector(), verbose, minProportionEachPartition, minCountEachPartition, improvementThreshold, true, maxDepth);
                    momentTree.determineSplit();
                    momentTree.printTree();

                    double MSE = 0;
                    double MSPE = 0;
                    double counter = 0;
                    int nullCounter = 0;
                    for (int i = 0; i < treeX.getRowDimension(); i++) {
                        Double predictedY = momentTree.getPredictedY(treeX.getMatrix(i, i, 0, treeX.getColumnDimension() - 1));
                        if (predictedY != null) {
                            MSE += Math.pow(treeY.get(i, 0) - predictedY, 2);
                            counter++;
                        } else {
                            pmUtility.prettyPrint(treeX.getMatrix(i, i, 0, treeX.getColumnDimension() - 1));
                            nullCounter++;
                        }
                    }
                    MSE /= counter;

                    boolean computeHonestTree = true;
                    if (computeHonestTree) {
                        boolean TIMING_DEBUG = true;

                        long t1 = System.currentTimeMillis();
                        for (int i = 0; i < honestX.getRowDimension(); i++) {
                            Jama.Matrix xi = honestX.getMatrix(i, i, 0, honestX.getColumnDimension() - 1);
                            Jama.Matrix yi = honestY.getMatrix(i, i, 0, honestY.getColumnDimension() - 1);
                            momentTree.sortXToCorrectLeafs(yi, xi);
                        }
                        momentTree.consolidateHonestData();
                        long t2 = System.currentTimeMillis();
                        if (TIMING_DEBUG) {
                            System.out.println("Sorted observations into tree in " + (t2 - t1) + " ms.");
                        }

                        t1 = System.currentTimeMillis();
                        momentTree.estimateHonestTree();
                        t2 = System.currentTimeMillis();
                        if (TIMING_DEBUG) {
                            System.out.println("Recomputed estimates in " + (t2 - t1) + " ms.");
                        }
                    }

                    counter = 0;
                    int nullCounterMSPE = 0;
                    for (int i = 0; i < predictX.getRowDimension(); i++) {
                        Double predictedY = momentTree.getPredictedY(predictX.getMatrix(i, i, 0, predictX.getColumnDimension() - 1));
                        if (predictedY != null) {
                            MSPE += Math.pow(predictY.get(i, 0) - predictedY, 2);
                            counter++;
                        } else {
                            pmUtility.prettyPrint(predictX.getMatrix(i, i, 0, honestX.getColumnDimension() - 1));
                            nullCounterMSPE++;
                        }
                    }
                    MSPE /= counter;

                    System.out.print("maxDepth: " + maxDepth + " " + " k: " + minCountEachPartition + " " + " alpha: " + minProportionEachPartition + " " + " mse_bar: " + improvementThreshold + " ");
                    System.out.print("MSE: " + MSE + " MSPE: " + MSPE + " nulls: " + nullCounter + " ");
                    // System.out.println("MSPE: " + MSPE + " nulls: " + nullCounterMSPE);
                    if (MSPE < bestMSPE || first) {
                        System.out.println(" * ");
                        bestK = minCountEachPartition;
                        bestAlpha = minProportionEachPartition;
                        bestMSEBar = improvementThreshold;
                        bestDepth = maxDepth;
                        bestMSPE = MSPE;
                        bestProportion = proportion;
                        first = false;
                    } else {
                        System.out.println();
                    }
                }
            }
        }
        System.out.println("CV parameters:");
        System.out.println("\t k = " + bestK);
        System.out.println("\t alpha = " + bestAlpha);
        System.out.println("\t mse_bar = " + bestMSEBar);
        System.out.println("\t depth = " + bestDepth);
        System.out.println("\t proportion used to estimate tree = " + bestProportion);
        System.out.println("Best MSPE: " + bestMSPE);

        int numObsToEstimateTreeStructure = (int) Math.floor(numObs * bestProportion);
        System.out.println("numObs: " + numObs + " halfObs: " + numObsToEstimateTreeStructure);
        Jama.Matrix treeX = new Jama.Matrix(numObsToEstimateTreeStructure, bart.getX().getColumnDimension());
        Jama.Matrix treeY = new Jama.Matrix(numObsToEstimateTreeStructure, 1);
        Jama.Matrix honestX = new Jama.Matrix(numObs - numObsToEstimateTreeStructure, bart.getX().getColumnDimension());
        Jama.Matrix honestY = new Jama.Matrix(numObs - numObsToEstimateTreeStructure, 1);
        TreeSet<Integer> treeSet = new TreeSet<>();
        int count = 0;

        Random honestRNG = new Random(787);

        while (count < numObsToEstimateTreeStructure) {
            if (count < numObsToEstimateTreeStructure) {
                int index = (int) Math.floor(honestRNG.nextDouble() * numObs);
                if (!treeSet.contains(index)) {
                    treeSet.add(index);
                    count++;
                }
            }
        }
        System.out.println("honestX row: " + honestX.getRowDimension() + " col: " + honestX.getColumnDimension());
        System.out.println("treeSet.size(): " + treeSet.size());

        int countTree = 0;
        int countHonest = 0;
        for (int i = 0; i < numObs; i++) {
            if (treeSet.contains(i)) {
                for (int j = 0; j < bart.getX().getColumnDimension(); j++) {
                    treeX.set(countTree, j, bart.getX().get(i, j));
                }
                treeY.set(countTree, 0, bart.getY().get(i, 0));
                countTree++;

            } else {
                for (int j = 0; j < bart.getX().getColumnDimension(); j++) {
                    honestX.set(countHonest, j, bart.getX().get(i, j));
                }
                honestY.set(countHonest, 0, bart.getY().get(i, 0));
                countHonest++;
            }
            // System.out.println(i+" "+treeSet.contains(i)+" "+countTree+" "+countHonest);
        }

//        System.out.println("treeX:");
//        pmUtility.prettyPrint(treeX);
//
//        System.out.println("honestX:");
//        pmUtility.prettyPrint(honestX);
        // this would be the place to introduce resampling of the treeX/Y and the honestX/Y to produce forest estimates
        ArrayList<TreeMoment> forest = new ArrayList<>();
        int numTrees = 500;

        boolean verbose = false;
        boolean printTrees= false;

        Random rng = new Random();
        for (int fi = 0; fi < numTrees; fi++) {
            long seed = rng.nextLong();
            TreeMoment momentTree = new TreeMoment(null, bart, resample(treeX, seed), resample(treeY, seed), bart.getDiscreteVector(), verbose, bestAlpha, bestK, bestMSEBar, true, bestDepth);
            if (verbose) {
                System.out.println("----------------------- mi: " + fi + " -----------------------");
            }
            momentTree.determineSplit();

            seed = rng.nextLong();
            Jama.Matrix thisHonestX = resample(honestX, seed);
            Jama.Matrix thisHonestY = resample(honestY, seed);

//            Jama.Matrix doublecheck = pmUtility.OLS(pmUtility.getColumn(thisHonestX, 0), thisHonestY, true);
//            pmUtility.prettyPrint(pmUtility.concatMatrix(thisHonestY, thisHonestX));
//            System.out.println("Manual: " + doublecheck.get(1, 0));
            for (int i = 0; i < thisHonestX.getRowDimension(); i++) {
                Jama.Matrix xi = thisHonestX.getMatrix(i, i, 0, thisHonestX.getColumnDimension() - 1);
                Jama.Matrix yi = thisHonestY.getMatrix(i, i, 0, thisHonestY.getColumnDimension() - 1);
                momentTree.sortXToCorrectLeafs(yi, xi);
            }
            momentTree.consolidateHonestData();
            momentTree.estimateHonestTree();
            if (verbose || printTrees) {
                System.out.println("********* tree "+fi+" *********");
                momentTree.printTree();
            }
            forest.add(momentTree);
        }

        boolean writeOutput = true;
        if (writeOutput) {
            try {
                BufferedWriter bartOut = new BufferedWriter(new FileWriter("bartTau.csv"));
                TreeSet<Double> uniqueTau = new TreeSet<>();
                TreeSet<Double> significantUniqueTau = new TreeSet<>();
                for (int i = 0; i < allX.getRowDimension(); i++) {
                    // for (int i = 0; i < 25; i++) {
                    Jama.Matrix xi = allX.getMatrix(i, i, 0, allX.getColumnDimension() - 1);
                    xi.set(0, 0, 1); // set the treatment to one to see what the treatment effect is
                    Jama.Matrix tau = new Jama.Matrix(forest.size(), 1);
                    for (int mi = 0; mi < forest.size(); mi++) {
                        TreeMoment m = forest.get(mi);
                        Jama.Matrix beta = m.getEstimatedBeta(xi);
                        // pmUtility.prettyPrint(beta);
                        // tau.set(mi, 0, m.getPredictedY(xi));
                        if (verbose) {
                            System.out.println("mi: " + mi);
                        }
                        tau.set(mi, 0, beta.get(0, 0));
                    }
                    uniqueTau.add(pmUtility.mean(tau, 0));
                    String stars = "";
                    NormalDistribution normal = new NormalDistribution();
//                 System.out.println(normal.inverse(0.05));
//                 System.out.println(normal.inverse(0.025));
//                 System.out.println(normal.inverse(0.005));
//                 System.exit(0);

                    boolean useZStat = false;
                    if (useZStat) {
                        if (Math.abs(pmUtility.mean(tau, 0) / pmUtility.standardDeviation(tau, 0)) > Math.abs(normal.inverse(0.05))) {
                            stars = "*";
                        }
                        if (Math.abs(pmUtility.mean(tau, 0) / pmUtility.standardDeviation(tau, 0)) > Math.abs(normal.inverse(0.025))) {
                            stars = "**";
                        }
                        if (Math.abs(pmUtility.mean(tau, 0) / pmUtility.standardDeviation(tau, 0)) > Math.abs(normal.inverse(0.005))) {
                            stars = "***";
                        }
                    } else {
                        if (pmUtility.percentile(tau, 0, 0.05) * pmUtility.percentile(tau, 0, 0.95) > 0) {
                            stars = "*";
                        }
                        if (pmUtility.percentile(tau, 0, 0.025) * pmUtility.percentile(tau, 0, 0.975) > 0) {
                            stars = "**";
                        }
                        if (pmUtility.percentile(tau, 0, 0.005) * pmUtility.percentile(tau, 0, 0.995) > 0) {
                            stars = "***";
                        }
                    }
                    if (!stars.equals("")) {
                        significantUniqueTau.add(pmUtility.mean(tau, 0));
                    }
                    System.out.format("%.3f \t %s \t %.3f \t [%.3f, %.3f] ", pmUtility.mean(tau, 0), stars, pmUtility.standardDeviation(tau, 0), pmUtility.percentile(tau, 0, 0.025), pmUtility.percentile(tau, 0, 0.975));
                    pmUtility.prettyPrint(xi);
                    bartOut.write(pmUtility.mean(tau, 0) + "," + pmUtility.standardDeviation(tau, 0) + ",");
                    for (int j = 0; j < xi.getColumnDimension(); j++) {
                        bartOut.write(xi.get(0, j) + ",");
                    }
                    bartOut.write("\n");
                }
                System.out.println("Estimating " + uniqueTau.size() + " unique treatment effects.");
                System.out.println("Estimating " + significantUniqueTau.size() + " statistically significant unique treatment effects.");
                bartOut.close();
            } catch (Exception e) {
                e.printStackTrace();
            }
        }
    }

    public static Jama.Matrix resample(Jama.Matrix x, long seed) {
        Random rng = new Random(seed);
        if (1 == 0) {
            return x.copy();
        }
        Jama.Matrix re = new Jama.Matrix(x.getRowDimension(), x.getColumnDimension());
        for (int i = 0; i < re.getRowDimension(); i++) {
            int index = (int) Math.floor(re.getRowDimension() * rng.nextDouble());
            for (int j = 0; j < re.getColumnDimension(); j++) {
                re.set(i, j, x.get(index, j));
            }
        }
        return re;
    }

}
