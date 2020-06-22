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

import Jama.Matrix;
import java.util.ArrayList;
import java.util.Random;
import java.util.TreeSet;
import utility.pmUtility;

/**
 * Class containing a collection of trees and utility classes for interacting
 * with it
 *
 * @author Stephen P. Ryan <stephen.p.ryan@wustl.edu>
 */
public class MomentForest {

    ArrayList<TreeMoment> forest;
    int numberTreesInForest;
    TreeOptions treeOptions;
    Jama.Matrix forestX;
    Jama.Matrix forestY;
    long forestSeed;
    boolean verbose;

    int numObs;
    MomentSpecification spec;

    // TreeMoment momentTree = new TreeMoment(null, spec, treeX, treeY, spec.getDiscreteVector(), verbose,
    // minProportionEachPartition, minCountEachPartition, improvementThreshold, true, maxDepth, null, null);
    public MomentForest(MomentSpecification spec, int numberTreesInForest, long forestSeed, Jama.Matrix forestX, Jama.Matrix forestY,
            boolean verbose, TreeOptions options) {
        this.spec = spec;
        this.numberTreesInForest = numberTreesInForest;
        this.verbose = verbose;
        numObs = forestX.getRowDimension();
        this.forestX = forestX;
        this.forestY = forestY;
        treeOptions = options;
        this.forestSeed = forestSeed;
    }

    public void growForest() {
        int numObsToEstimateTreeStructure = (int) Math.floor(numObs * 0.25);

        Jama.Matrix treeX = new Jama.Matrix(numObsToEstimateTreeStructure, forestX.getColumnDimension());
        Jama.Matrix treeY = new Jama.Matrix(numObsToEstimateTreeStructure, 1);
        Jama.Matrix honestX = new Jama.Matrix(numObs - numObsToEstimateTreeStructure, forestX.getColumnDimension());
        Jama.Matrix honestY = new Jama.Matrix(numObs - numObsToEstimateTreeStructure, 1);
        TreeSet<Integer> treeSet = new TreeSet<>();
        int count = 0;

        Random honestRNG = new Random(forestSeed);

        while (count < numObsToEstimateTreeStructure) {
            if (count < numObsToEstimateTreeStructure) {
                int index = (int) Math.floor(honestRNG.nextDouble() * numObs);
                if (!treeSet.contains(index)) {
                    treeSet.add(index);
                    count++;
                }
            }
        }

        /**
         * Populate data used for determining splits (treeX, treeY) and
         * estimating values conditional on those splits (honestX, honestY)
         */
        int countTree = 0;
        int countHonest = 0;
        for (int q = 0; q < numObs; q++) {
            if (treeSet.contains(q)) {
                for (int j = 0; j < forestX.getColumnDimension(); j++) {
                    treeX.set(countTree, j, forestX.get(q, j));
                }
                treeY.set(countTree, 0, forestY.get(q, 0));
                countTree++;

            } else {
                for (int j = 0; j < forestX.getColumnDimension(); j++) {
                    honestX.set(countHonest, j, forestX.get(q, j));
                }
                honestY.set(countHonest, 0, forestY.get(q, 0));
                countHonest++;
            }
        }

        Random rng = new Random(honestRNG.nextLong());
        forest = new ArrayList<>();

        for (int i = 0; i < numberTreesInForest; i++) {
            long seed = rng.nextLong();
            long seedHonest = rng.nextLong();
            forest.add(new TreeMoment(null, spec, resample(treeX, seed), resample(treeY, seed),
                    spec.getDiscreteVector(), verbose, treeOptions.getMinProportion(), treeOptions.getMinCount(), treeOptions.getMinMSEImprovement(), true, treeOptions.getMaxDepth(),
                    resample(honestX, seedHonest), resample(honestY, seedHonest)));
        }

        forest.parallelStream().forEach((tree) -> {
            tree.determineSplit();
            Jama.Matrix thisHonestX = tree.getHonestXtemp();
            Jama.Matrix thisHonestY = tree.getHonestYtemp();

            for (int i = 0; i < thisHonestX.getRowDimension(); i++) {
                Jama.Matrix xi = thisHonestX.getMatrix(i, i, 0, thisHonestX.getColumnDimension() - 1);
                Jama.Matrix yi = thisHonestY.getMatrix(i, i, 0, thisHonestY.getColumnDimension() - 1);
                tree.sortXToCorrectLeafs(yi, xi);
            }
            tree.consolidateHonestData();
            tree.estimateHonestTree();
        });
    }

    public Matrix getEstimatedParameters(Matrix x) {
        Jama.Matrix estimatedParameters = forest.get(0).getEstimatedBeta(x);
        String s = "[ " + forest.get(0).getEstimatedBeta(x).get(0, 0) + " ";
        for (int i = 1; i < forest.size(); i++) {
            estimatedParameters = estimatedParameters.plus(forest.get(i).getEstimatedBeta(x));
            s = s.concat(forest.get(i).getEstimatedBeta(x).get(0, 0) + " ");
        }
        s = s.concat("]");
//        System.out.println(s);
        estimatedParameters.timesEquals(1.0 / forest.size());
        return estimatedParameters;
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

    public TreeOptions performCrossValidation() {
        double proportion = 0.33;
        int numObsEstimateLeafValues = (int) Math.floor(numObs * proportion);
        int numPredictObs = (int) Math.floor(numObs * 0.33); // it was originally 10;
        int numObsGrowTreeStructure = numObs - numObsEstimateLeafValues - numPredictObs; // halfObs / 2;
        // SFIToolkit.displayln("numObs: " + numObs + " halfObs: " + numObsEstimateLeafValues + " predictObs: " + numPredictObs);

        Jama.Matrix treeX = new Jama.Matrix(numObsGrowTreeStructure, spec.getX().getColumnDimension());
        Jama.Matrix treeY = new Jama.Matrix(numObsGrowTreeStructure, 1);
        Jama.Matrix honestX = new Jama.Matrix(numObsEstimateLeafValues, spec.getX().getColumnDimension());
        Jama.Matrix honestY = new Jama.Matrix(numObsEstimateLeafValues, 1);
        Jama.Matrix predictX = new Jama.Matrix(numPredictObs, spec.getX().getColumnDimension());
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

        int countTree = 0;
        int countPredict = 0;
        int countHonest = 0;
        for (int i = 0; i < numObs; i++) {
            if (growTreeSet.contains(i)) {
                for (int j = 0; j < spec.getX().getColumnDimension(); j++) {
                    treeX.set(countTree, j, spec.getX().get(i, j));
                }
                treeY.set(countTree, 0, spec.getY().get(i, 0));
                countTree++;
            } else if (predictSet.contains(i)) {
                for (int j = 0; j < spec.getX().getColumnDimension(); j++) {
                    predictX.set(countPredict, j, spec.getX().get(i, j));
                }
                predictY.set(countPredict, 0, spec.getY().get(i, 0));
                countPredict++;
            } else {
                for (int j = 0; j < spec.getX().getColumnDimension(); j++) {
                    honestX.set(countHonest, j, spec.getX().get(i, j));
                }
                honestY.set(countHonest, 0, spec.getY().get(i, 0));
                countHonest++;
            }
            // System.out.println(i+" "+treeSet.contains(i)+" "+countTree+" "+countHonest);
        }

        boolean verboseTreeConstruction = true;
        double minProportionEachPartition = 1E-5;
        double improvementThreshold = 50.0;
        double bestMSPE = 0;
        double bestMSEBar = 0;
        boolean first = true;
        int bestK = 10;

        TreeOptions options = new TreeOptions();
        options.setMaxDepth(100);
        options.setMinMSEImprovement(1E-10);

        for (improvementThreshold = 0.20; improvementThreshold >= 0.001; improvementThreshold *= 0.9) {
            for (int minCountEachPartition = 10; minCountEachPartition <= 10; minCountEachPartition += 10) {
                Random rng = new Random(667);
                options.setMinCount(minCountEachPartition);
                options.setMinMSEImprovement(improvementThreshold);
                // TreeMoment momentTree = new TreeMoment(null, spec, treeX, treeY, spec.getDiscreteVector(), verbose,                        minProportionEachPartition, minCountEachPartition, improvementThreshold, true, maxDepth, null, null);
                // momentTree.determineSplit();
                // momentTree.printTree();
                MomentForest momentForest = new MomentForest(spec, 1, rng.nextLong(), treeX, treeY, verboseTreeConstruction, options);
                momentForest.growForest();

                double MSPE = 0;

                int counter = 0;
                int nullCounterMSPE = 0;
                for (int i = 0; i < predictX.getRowDimension(); i++) {
                    Jama.Matrix xi = predictX.getMatrix(i, i, 0, predictX.getColumnDimension() - 1);
                    Double predictedY = spec.getPredictedY(xi, momentForest.getEstimatedParameters(xi));
                    if (predictedY != null) {
                        MSPE += Math.pow(predictY.get(i, 0) - predictedY, 2);
                        counter++;
                    } else {
                        pmUtility.prettyPrint(predictX.getMatrix(i, i, 0, honestX.getColumnDimension() - 1));
                        nullCounterMSPE++;
                    }
                }
                MSPE /= counter;

                // SFIToolkit.displayln("CV_index: " + CV_Index + " " );
                // System.out.print("maxDepth: " + options.getMaxDepth() + " " + " k: " + options.getMinCount() + " " + " alpha: " + options.getMinProportion()
                // + " " + " mse_bar: " + options.getMinMSEImprovement() + " ");
                // System.out.print("MSPE: " + MSPE + " nulls: " + nullCounterMSPE + " ");
                String s = "";
                if (MSPE < bestMSPE || first) {
                    s = " * ";
                    bestK = minCountEachPartition;
                    // best = minProportionEachPartition;
                    bestMSEBar = improvementThreshold;
                    // bestDepth = maxDepth;
                    bestMSPE = MSPE;
                    // bestProportion = proportion;
                    first = false;
                }
                System.out.format("maxDepth: %d k: %d mse_bar: %g mspe: %g nulls: %d %s %n", options.getMaxDepth(), options.getMinCount(), options.getMinMSEImprovement(),
                         MSPE, nullCounterMSPE, s);
                // System.out.println("MSPE: " + MSPE + " nulls: " + nullCounterMSPE);

            }
        }
        System.out.println("Optimal minimum number of observations in each leaf: " + bestK + " MSPE: " + bestMSPE);
        System.out.println("Optimal improvement threshold: " + bestMSEBar);

        options.setMinCount(bestK);
        options.setMinMSEImprovement(bestMSEBar);

        return options;
    }

    public void setTreeOptions(TreeOptions cvOptions) {
        this.treeOptions = cvOptions;
    }

}
