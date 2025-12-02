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
import java.util.stream.Collectors;
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
    double proportionObservationsToEstimateTreeStructure;
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
        this.proportionObservationsToEstimateTreeStructure = spec.getProportionObservationsToEstimateTreeStructure();
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

        Random rng = new Random(forestSeed);

        for (int i = 0; i < numberTreesInForest; i++) {
            // resample the forestLens, then split it
            DataLens resampled;
            DataLens[] split;
            if (forestLens.balancingVector != null) {
                //System.out.println("Tree " + i + ": Using balancing vector");
                resampled = forestLens.getResampledDataLensWithBalance(rng.nextLong());
                split = resampled.randomlySplitSampleWithBalance(proportionObservationsToEstimateTreeStructure, rng.nextLong());

            } else if (forestLens.strataColumnIndex != null) {
                //System.out.println("Tree " + i + ": Using strata column " + forestLens.strataColumnIndex);
                resampled = forestLens.getResampledDataLens(rng.nextLong());
                split = resampled.randomlySplitSampleByStrata(proportionObservationsToEstimateTreeStructure, rng.nextLong());

            } else {
                //System.out.println("Tree " + i + ": Using simple random sampling");
                resampled = forestLens.getResampledDataLens(rng.nextLong());
                split = resampled.randomlySplitSample(proportionObservationsToEstimateTreeStructure, rng.nextLong());
            }

            DataLens lensGrow = split[0];
            DataLens lensHonest = split[1];

            forest.add(new TreeMoment(null, spec, lensGrow,
                    spec.getDiscreteVector(), verbose, treeOptions.getMinProportion(), treeOptions.getMinCount(), treeOptions.getMinMSEImprovement(), true, treeOptions.getMaxDepth(),
                    lensHonest, treeOptions.isTestParameterHomogeneity(), rng.nextLong()));
        }

        boolean useParallel = true;

        if (!useParallel) {
            for (int i = 0; i < numberTreesInForest; i++) {
                forest.get(i).determineSplit();
                forest.get(i).estimateHonestTree();
            }
        } else {
            forest.parallelStream().forEach((tree) -> {
                // System.out.println("Split");
                tree.determineSplit();
                // System.out.println("Honest");
                tree.estimateHonestTree();
            });
        }
        // System.out.format("Memory usage: %,d bytes %n", (Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory()));

        // prune out invalid trees? these would be ones where the minimizer failed?
        ArrayList<TreeMoment> validForest = new ArrayList<>();
        for (int i = 0; i < numberTreesInForest; i++) {
            if (forest.get(i).isValidTree()) {
                validForest.add(forest.get(i));
            } else {
                System.out.println("Removing tree " + i + " due to failed estimation.");
            }
        }
        forest = validForest;
        numberTreesInForest = forest.size();
    }

    /**
     * Get parameters associated with a given vector of observables zi
     *
     * @param zi Observable vector
     * @return
     */
    public Jama.Matrix getEstimatedParameterForest(Jama.Matrix zi) {
        Jama.Matrix estimatedParameters = forest.get(0).getEstimatedBeta(zi);
        // String s = "[ " + forest.get(0).getEstimatedBeta(zi).get(0, 0) + " "; //Assuming beta is 1 by 1?
        for (int i = 1; i < forest.size(); i++) {
            estimatedParameters = estimatedParameters.plus(forest.get(i).getEstimatedBeta(zi));
            // s = s.concat(forest.get(i).getEstimatedBeta(zi).get(0, 0) + " ");
        }

        estimatedParameters.timesEquals(1.0 / forest.size());
        return estimatedParameters;
    }

    public ArrayList<Jama.Matrix> getAllEstimatedParametersFromForest(Jama.Matrix zi) {
        ArrayList<Jama.Matrix> parameterList = new ArrayList<>();
        for (int i = 0; i < forest.size(); i++) {
            parameterList.add(forest.get(i).getEstimatedBeta(zi));
        }
        return parameterList;
    }

    public int getForestSize() {
        return forest.size();
    }

    public void setTreeOptions(TreeOptions cvOptions) {
        this.treeOptions = cvOptions;
    }

    public double[] getNumberTimesTreesInForestSplitOnAGivenVariableIndex() {

        double countEachVariableSplit[] = new double[spec.getDiscreteVector().length];

        for (TreeMoment tree : forest) {
            TreeSet<Integer> splitTree = new TreeSet<>();
            tree.getEnumerationOfAllSplitVariablesInThisTree(splitTree);

            for (int i : splitTree) {
                countEachVariableSplit[i] = countEachVariableSplit[i] + 1.0;
            }
        }
        return countEachVariableSplit;
    }

    public boolean[] getHomogeneityVotes(JTextArea jt, boolean verboseVoting) {

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
            // majority classification rule (greater than 50%)
            // votes[i] = voteCounts[i] > Math.floorDiv(numberTreesInForest, 2);
            // greater than some percentage
            votes[i] = voteCounts[i] > Math.floor(0.5 * numberTreesInForest);
        }
        if (verboseVoting) {
            System.out.print("votes: ");
        }
        for (int i = 0; i < voteCounts.length; i++) {
            // System.out.print(voteCounts[i]+"/"+votes[i]+" ");
            if (verboseVoting) {
                System.out.format("%.2f%% ", (100.0 * voteCounts[i] / numberTreesInForest));
            }
            double pct = 100.0 * voteCounts[i] / numberTreesInForest;
            if (verboseVoting) {
                jt.append(i + ". votes: " + voteCounts[i] + " out of " + numberTreesInForest + " (" + pct + "%): " + votes[i] + "\n");
            }
            if (voteCounts[i] < numberTreesInForest) {
                // System.out.println("Detected variance in voting on parameter "+i+": "+voteCounts[i]);
                // System.exit(0);
            }
        }
        if (verboseVoting) {
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
    
    public ArrayList<Double> getTestStatistics() {
        return forest.parallelStream()
            .map(e -> e.getTestStatistic())
            .collect(Collectors.toCollection(ArrayList::new));
    }
    
    public ArrayList<Double> getRestrictedThetas() {
        return forest.parallelStream()
            .map(e -> e.getRestrictedTheta())
            .collect(Collectors.toCollection(ArrayList::new));
    }

    public void testHomogeneity(boolean verbose) {
        boolean useParallel = true;

        double[] averageTestValues = new double[spec.getNumParams()];

        if (!useParallel) {
            for (int i = 0; i < numberTreesInForest; i++) {
                // System.out.println("========== Tree " + i + " ==========");
                // System.out.println("Testing homogeneity");
                forest.get(i).testHomogeneity();
            }
        } else {
            forest.parallelStream().forEach((tree) -> {
                tree.testHomogeneity();
                // System.out.println("Finished testing homogeneity in a tree");
            });
        }
        // System.out.println("Finished testing homogeneity in all the trees");

        // Prune out invalid 
        // System.out.println("Prune out invalid trees");
        ArrayList<TreeMoment> validForest = new ArrayList<>();
        for (int i = 0; i < forest.size(); i++) {
            TreeMoment t = forest.get(i);
            if (t.isValidTree()) {
                validForest.add(t);
            } else {
                System.out.println("Removing tree " + i + " due to failed homogeneity test.");
            }
        }
        forest = validForest;
        numberTreesInForest = forest.size();

        if (verbose) {
            System.out.println("Number of valid trees after homogeneity testing: " + numberTreesInForest);
        }

        // Compute average test values over remaining valid trees
        for (TreeMoment tree : forest) {
            double[] testValues_i = tree.getTestValues();
            for (int k = 0; k < testValues_i.length; k++) {
                if (!Double.isNaN(testValues_i[k]) && !Double.isInfinite(testValues_i[k])) {
                    averageTestValues[k] += testValues_i[k];
                } else {
                    System.out.println("Skipping invalid test value (NaN or Inf) for parameter " + k);
                }
            }
        }
        for (int k = 0; k < averageTestValues.length; k++) {
            averageTestValues[k] = averageTestValues[k] / numberTreesInForest;
            if (verbose) {
                System.out.println("parameter " + k + ": average DM test statistic: " + averageTestValues[k]);
            }
        }
    }
}
