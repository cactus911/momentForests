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
package core;

import java.util.ArrayList;
import java.util.Random;
import java.util.TreeSet;
import javax.swing.JTextArea;

/**
 * Class containing a collection of trees and utility classes for interacting
 * with it
 *
 * @author Stephen P. Ryan <stephen.p.ryan@wustl.edu>
 */

/* TO DO:
    1. getEstimatedParameters method prints the estimated parameters and appears to assume beta is 1x1 and not a vector. Need to allow for beta to be a vector.
 */
public class MomentForest {

    ArrayList<TreeMoment> forest = new ArrayList<>();
    int numberTreesInForest;
    TreeOptions treeOptions;
    DataLens forestLens;
    long forestSeed;
    boolean verbose;

    int numObs;
    MomentSpecification spec;

    // TreeMoment momentTree = new TreeMoment(null, spec, treeX, treeY, spec.getDiscreteVector(), verbose,
    // minProportionEachPartition, minCountEachPartition, improvementThreshold, true, maxDepth, null, null);
    public MomentForest(MomentSpecification spec, int numberTreesInForest, long forestSeed, DataLens forestLens,
            boolean verbose, TreeOptions options) {
        this.spec = spec;
        this.numberTreesInForest = numberTreesInForest;
        this.verbose = verbose;
        numObs = forestLens.getNumObs();
        this.forestLens = forestLens;
        treeOptions = options;
        this.forestSeed = forestSeed;
    }

    public TreeMoment getTree(int i) {
        return forest.get(i);
    }

    public void growForest() {
        /**
         * Reinitialize the forest in case we are doing a CV (or whatever). Do
         * not want to be stacking a bunch of trees of various options into the
         * forest each time we call growForest!
         */
        forest = new ArrayList<>();
        double proportionObservationsToEstimateTreeStructure = 0.35;

        Random rng = new Random(forestSeed);

        for (int i = 0; i < numberTreesInForest; i++) {
            // resample the forestLens, then split it
            DataLens resampled;
            DataLens[] split;
            if (forestLens.balancingVector == null) {
                resampled = forestLens.getResampledDataLens(rng.nextLong());
                split = resampled.randomlySplitSample(proportionObservationsToEstimateTreeStructure, rng.nextLong());
            } else {
                resampled = forestLens.getResampledDataLensWithBalance(rng.nextLong());
                split = resampled.randomlySplitSampleWithBalance(proportionObservationsToEstimateTreeStructure, rng.nextLong());
            }

            DataLens lensGrow = split[0];
            DataLens lensHonest = split[1];

            forest.add(new TreeMoment(null, spec, lensGrow,
                    spec.getDiscreteVector(), verbose, treeOptions.getMinProportion(), treeOptions.getMinCount(), treeOptions.getMinMSEImprovement(), true, treeOptions.getMaxDepth(),
                    lensHonest, treeOptions.isTestParameterHomogeneity()));
        }

        boolean useParallel = true;

        if (!useParallel) {
            for (int i = 0; i < numberTreesInForest; i++) {
                forest.get(i).determineSplit();
                forest.get(i).estimateHonestTree();
            }
        } else {
            forest.parallelStream().forEach((tree) -> {
                tree.determineSplit();
                tree.estimateHonestTree();
            });
        }
        // System.out.format("Memory usage: %,d bytes %n", (Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory()));
    }

    /**
     * Get parameters associated with a given vector of observables zi
     *
     * @param zi Observable vector
     * @return
     */
    public Jama.Matrix getEstimatedParameterForest(Jama.Matrix zi) {
        Jama.Matrix estimatedParameters = forest.get(0).getEstimatedBeta(zi);
//        String s = "[ " + forest.get(0).getEstimatedBeta(zi).get(0, 0) + " "; //Assuming beta is 1 by 1?
        for (int i = 1; i < forest.size(); i++) {
            estimatedParameters = estimatedParameters.plus(forest.get(i).getEstimatedBeta(zi));
//            s = s.concat(forest.get(i).getEstimatedBeta(zi).get(0, 0) + " ");
        }
//        s = s.concat("]");
//        System.out.println(s);
//        System.exit(0);
        estimatedParameters.timesEquals(1.0 / forest.size());
        return estimatedParameters;
    }
    
    public int getForestSize() {
        return forest.size();
    }

    public void setTreeOptions(TreeOptions cvOptions) {
        this.treeOptions = cvOptions;
    }

    public double[] getCountSplitVariables() {

        double countEachVariableSplit[] = new double[spec.getDiscreteVector().length];

        for (TreeMoment tree : forest) {
            TreeSet<Integer> splitTree = new TreeSet<>();
            tree.getIndexSplitVariables(splitTree);

            for (int i : splitTree) {
                countEachVariableSplit[i] = countEachVariableSplit[i] + 1.0;
            }
        }
        return countEachVariableSplit;
    }

    public boolean[] getHomogeneityVotes(JTextArea jt) {
        int[] voteCounts = new int[spec.getHomogeneousIndex().length];
        for (int i = 0; i < numberTreesInForest; i++) {
            ArrayList<Integer> hpl = getTree(i).getIndexHomogeneousParameters();
            // ArrayList<Double> hplStartingValues = getTree(i).getValueHomogeneousParameters();
            for (Integer h : hpl) {
                voteCounts[h] = voteCounts[h] + 1;
            }
        }
        boolean[] votes = new boolean[voteCounts.length];
        for (int i = 0; i < votes.length; i++) {
            votes[i] = voteCounts[i] > Math.floorDiv(numberTreesInForest, 2);
        }
        if (verbose) {
            System.out.print("votes: ");
        }
        for (int i = 0; i < voteCounts.length; i++) {
            // System.out.print(voteCounts[i]+"/"+votes[i]+" ");
            if (verbose) {
                System.out.format("%.2f%% ", (100.0 * voteCounts[i] / numberTreesInForest));
            }
            double pct = 100.0 * voteCounts[i] / numberTreesInForest;
//            jt.append(i+". votes: "+voteCounts[i]+" out of "+numberTreesInForest+" ("+pct+")\n");
            if (voteCounts[i] < numberTreesInForest) {
                // System.out.println("Detected variance in voting on parameter "+i+": "+voteCounts[i]);
                // System.exit(0);
            }
        }
        if(verbose) {
            System.out.println("");
        }

        return votes;
    }

    public double[] getHomogeneityStartingValues() {
        double[] startingValues = new double[spec.getHomogeneousIndex().length];
        double[] voteCounts = new double[spec.getHomogeneousIndex().length];
        for (int i = 0; i < numberTreesInForest; i++) {
            ArrayList<Integer> hpl = getTree(i).getIndexHomogeneousParameters();
            ArrayList<Double> hplStartingValues = getTree(i).getValueHomogeneousParameters();
//            for(Double d : hplStartingValues) {
//                System.out.print(d+" ");
//            }
            // System.out.println();
            for (int j = 0; j < hpl.size(); j++) {
                int index = hpl.get(j);
                voteCounts[index] = voteCounts[index] + 1;
                startingValues[index] = startingValues[index] + hplStartingValues.get(j);
            }
        }
        for (int i = 0; i < startingValues.length; i++) {
            if (voteCounts[i] > 0) {
                startingValues[i] = startingValues[i] / voteCounts[i];
            }
            if (Double.isNaN(startingValues[i])) {
                System.out.println("Detected starting value of NaN:");
                System.out.println("voteCounts: " + voteCounts[i]);

                for (int t = 0; t < numberTreesInForest; t++) {
                    System.out.println("Tree " + t + ": ");
                    ArrayList<Integer> hpl = getTree(t).getIndexHomogeneousParameters();
                    ArrayList<Double> hplStartingValues = getTree(t).getValueHomogeneousParameters();
                    for (Double d : hplStartingValues) {
                        System.out.print(d + " ");
                    }
                    System.out.println("");
                    System.out.print("Index homogeneous parameters: ");
                    for (Integer h : hpl) {
                        System.out.print(h + " ");
                    }

                    System.out.println("");
                }
                System.exit(0);
            }
        }
        return startingValues;
    }

    public void testHomogeneity() {
        boolean useParallel = false;

        if (!useParallel) {
            for (int i = 0; i < numberTreesInForest; i++) {
                forest.get(i).testHomogeneity();
            }
        } else {
            forest.parallelStream().forEach((tree) -> {
                tree.testHomogeneity();
            });
        }
    }

}
