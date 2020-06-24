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

    public TreeMoment getTree(int i) {
        return forest.get(i);
    }

    public void growForest() {
        double proportionObservationsToEstimateTreeStructure = 0.15;
        int numObsToEstimateTreeStructure = (int) Math.floor(numObs * proportionObservationsToEstimateTreeStructure);
        // System.out.println("numObs Estimate Tree Structure: "+numObsToEstimateTreeStructure+" out of "+numObs);
        double ratioTreatment = pmUtility.mean(forestX, 0);
        int goalTreatmentObsTreeX = (int) Math.round(ratioTreatment * numObsToEstimateTreeStructure);
        int countTreatmentTreeX = 0;

        /**
         * Want to figure out how to balance these treeX/Y and honestX/Y
         * matrices on treatment status.
         */
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
                if (countTreatmentTreeX < goalTreatmentObsTreeX) {
                    while (forestX.get(index, 0) == 0) {
                        index = (int) Math.floor(honestRNG.nextDouble() * numObs);
                    }
                    countTreatmentTreeX++;
                } else {
                    while (forestX.get(index, 0) == 1) {
                        index = (int) Math.floor(honestRNG.nextDouble() * numObs);
                    }
                }
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

        // System.out.println("Mean treeX: "+pmUtility.mean(treeX, 0)+" Mean HonestX: "+pmUtility.mean(honestX, 0));
        Random rng = new Random(honestRNG.nextLong());
        forest = new ArrayList<>();

        Jama.Matrix balancingTreeX = pmUtility.getColumn(treeX, 0);
        Jama.Matrix balancingHonestX = pmUtility.getColumn(honestX, 0);

        for (int i = 0; i < numberTreesInForest; i++) {
            long seed = rng.nextLong();
            long seedHonest = rng.nextLong();

            boolean useResampleMatrix = true;
            if (!useResampleMatrix) {
                forest.add(new TreeMoment(null, spec, resample(treeX, seed, balancingTreeX), resample(treeY, seed, balancingTreeX),
                        spec.getDiscreteVector(), verbose, treeOptions.getMinProportion(), treeOptions.getMinCount(), treeOptions.getMinMSEImprovement(), true, treeOptions.getMaxDepth(),
                        resample(honestX, seedHonest, balancingHonestX), resample(honestY, seedHonest, balancingHonestX)));
            } else {
                ResamplingLens reTreeX = new ResamplingLens(treeX, seed, balancingTreeX);
                ResamplingLens reTreeY = new ResamplingLens(treeY, seed, balancingTreeX);
                ResamplingLens reHonestX = new ResamplingLens(honestX, seedHonest, balancingHonestX);
                ResamplingLens reHonestY = new ResamplingLens(honestY, seedHonest, balancingHonestX);
                forest.add(new TreeMoment(null, spec, reTreeX, reTreeY,
                        spec.getDiscreteVector(), verbose, treeOptions.getMinProportion(), treeOptions.getMinCount(), treeOptions.getMinMSEImprovement(), true, treeOptions.getMaxDepth(),
                        reHonestX, reHonestY));
            }
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
            
            tree.clearEstimationData();
            tree.clearHonestyData();
        });
        // System.out.format("Memory usage: %,d bytes %n", (Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory()));
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

    public static Jama.Matrix resample(Jama.Matrix x, long seed, Jama.Matrix balancingVector) {
        /**
         * Think about balancing on treatment/control here for better
         * small-sample splitting performance?
         *
         * How to do that? Do an acceptance/rejection method to balance at
         * half/half? Need to make sure that it doesn't get stuck Maybe ensure
         * that the ratio of treatment and control is held (close to) constant
         *
         * Since we use this method to resample X and Y, need to pass a
         * balancing vector for the case when we are resampling Y
         */

        double ratio = pmUtility.mean(balancingVector, 0);
        int goalTreatment = (int) Math.round(ratio * balancingVector.getRowDimension());

        Random rng = new Random(seed);

        int countTreatment = 0;
        Jama.Matrix re = new Jama.Matrix(x.getRowDimension(), x.getColumnDimension());
        for (int i = 0; i < re.getRowDimension(); i++) {
            int index = (int) Math.floor(re.getRowDimension() * rng.nextDouble());
            double treatmentIndicator = balancingVector.get(index, 0);
            if (countTreatment < goalTreatment) {
                while (treatmentIndicator == 0) {
                    index = (int) Math.floor(re.getRowDimension() * rng.nextDouble());
                    treatmentIndicator = balancingVector.get(index, 0);
                }
                countTreatment++;
            } else {
                while (treatmentIndicator == 1) {
                    index = (int) Math.floor(re.getRowDimension() * rng.nextDouble());
                    treatmentIndicator = balancingVector.get(index, 0);
                }
            }
            for (int j = 0; j < re.getColumnDimension(); j++) {
                re.set(i, j, x.get(index, j));
            }
        }
        return re;
    }

    public TreeOptions performCrossValidation(int numTrees) {
        int numPredictObs = (int) Math.floor(numObs * 0.5);
        int numObsEstimateTree = numObs - numPredictObs; // halfObs / 2;
        // SFIToolkit.displayln("numObs: " + numObs + " halfObs: " + numObsEstimateLeafValues + " predictObs: " + numPredictObs);

        double ratioTreatment = pmUtility.mean(spec.getX(), 0);
        int goalTreatmentObsTreeX = (int) Math.round(ratioTreatment * numObsEstimateTree);
        int countTreatmentTreeX = 0;
        // System.out.println("ratio: " + ratioTreatment + " goalTreeX: " + goalTreatmentObsTreeX + " goalPredictX: " + goalTreatmentObsPredictX + " (out of " + numPredictObs + ")");

        Jama.Matrix treeX = new Jama.Matrix(numObsEstimateTree, spec.getX().getColumnDimension());
        Jama.Matrix treeY = new Jama.Matrix(numObsEstimateTree, 1);
        Jama.Matrix predictX = new Jama.Matrix(numPredictObs, spec.getX().getColumnDimension());
        Jama.Matrix predictY = new Jama.Matrix(numPredictObs, 1);

        /**
         * Balancing split of sample into two parts via treatment status
         */
        TreeSet<Integer> growTreeSet = new TreeSet<>();
        int count = 0;
        while (count < numObsEstimateTree) {
            if (count < numObsEstimateTree) {
                int index = (int) Math.floor(Math.random() * numObs);
                if (countTreatmentTreeX < goalTreatmentObsTreeX) {
                    while (spec.getX().get(index, 0) == 0) {
                        index = (int) Math.floor(Math.random() * numObs);
                    }

                } else {
                    while (spec.getX().get(index, 0) == 1) {
                        index = (int) Math.floor(Math.random() * numObs);
                    }
                }
                if (!growTreeSet.contains(index)) {
                    growTreeSet.add(index);
                    if (spec.getX().get(index, 0) == 1) {
                        countTreatmentTreeX++;
                    }
                    count++;
                }
            }
        }
        // System.out.println("countTreatmentTreeX: " + countTreatmentTreeX + " (out of " + numObsGrowTreeStructure + "); ratio = " + ((countTreatmentTreeX + 0.0) / (numObsGrowTreeStructure + 0.0)));

        /**
         * Obviously the other thing to consider is stratifying by group
         */
        int countTree = 0;
        int countPredict = 0;
        for (int i = 0; i < numObs; i++) {
            if (growTreeSet.contains(i)) {
                for (int j = 0; j < spec.getX().getColumnDimension(); j++) {
                    treeX.set(countTree, j, spec.getX().get(i, j));
                }
                treeY.set(countTree, 0, spec.getY().get(i, 0));
                countTree++;
            } else {
                for (int j = 0; j < spec.getX().getColumnDimension(); j++) {
                    predictX.set(countPredict, j, spec.getX().get(i, j));
                }
                predictY.set(countPredict, 0, spec.getY().get(i, 0));
                countPredict++;
            }
            // System.out.println(i+" "+treeSet.contains(i)+" "+countTree+" "+countHonest);
        }

        // System.out.println("sum treeX: " + pmUtility.sum(treeX, 0) + " mean treeX: " + pmUtility.mean(treeX, 0) + " mean honestX: " + pmUtility.mean(honestX, 0) + " mean predictX: " + pmUtility.mean(predictX, 0));
        // System.exit(0);
        boolean verboseTreeConstruction = false;
        double minProportionEachPartition = 1E-5;
        double improvementThreshold = 50.0;
        double bestMSPE = 0;
        double bestMSEBar = 0;
        boolean first = true;
        int bestK = 10;

        TreeOptions options = new TreeOptions();
        options.setMaxDepth(100);
        options.setMinMSEImprovement(1E-10);

        for (improvementThreshold = 5E-5; improvementThreshold <= 0.5; improvementThreshold *= 2.5) {
            for (int minCountEachPartition = 11; minCountEachPartition >= 1; minCountEachPartition -= 1) {
                Random rng = new Random(667);
                options.setMinCount(minCountEachPartition);
                options.setMinMSEImprovement(improvementThreshold);
                
                MomentForest momentForest = new MomentForest(spec, numTrees, rng.nextLong(), treeX, treeY, verboseTreeConstruction, options);
                momentForest.growForest();
                /**
                 * We can easily implement cross-fitting here by swapping the
                 * prediction and estimation data sets and adding another MSPE
                 * component
                 */
                MomentForest momentForestSwitch = new MomentForest(spec, numTrees, rng.nextLong(), predictX, predictY, verboseTreeConstruction, options);
                momentForestSwitch.growForest();

//                System.out.println("-------");
//                System.out.println("Tree A:");
//                momentForest.getTree(0).printTree();
//                System.out.println("Tree B:");
//                momentForestSwitch.getTree(0).printTree();
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
                        pmUtility.prettyPrint(predictX.getMatrix(i, i, 0, predictX.getColumnDimension() - 1));
                        nullCounterMSPE++;
                    }
                }
                /**
                 * Cross-fitting part (reverse estimation data with prediction
                 * data)
                 */
                for (int i = 0; i < treeX.getRowDimension(); i++) {
                    Jama.Matrix xi = treeX.getMatrix(i, i, 0, treeX.getColumnDimension() - 1);
                    Double predictedY = spec.getPredictedY(xi, momentForestSwitch.getEstimatedParameters(xi));
                    if (predictedY != null) {
                        MSPE += Math.pow(treeY.get(i, 0) - predictedY, 2);
                        counter++;
                    } else {
                        pmUtility.prettyPrint(treeX.getMatrix(i, i, 0, treeX.getColumnDimension() - 1));
                        nullCounterMSPE++;
                    }
                }

                MSPE /= counter;

                // SFIToolkit.displayln("CV_index: " + CV_Index + " " );
                // System.out.print("maxDepth: " + options.getMaxDepth() + " " + " k: " + options.getMinCount() + " " + " alpha: " + options.getMinProportion()
                // + " " + " mse_bar: " + options.getMinMSEImprovement() + " ");
                // System.out.print("MSPE: " + MSPE + " nulls: " + nullCounterMSPE + " ");
                String s = "";
                if (MSPE == bestMSPE) {
                    s = " - ";
                }
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
                // System.out.format("maxDepth: %d k: %d mse_bar: %g mspe: %g nulls: %d %s %n", options.getMaxDepth(), options.getMinCount(), options.getMinMSEImprovement(), MSPE, nullCounterMSPE, s);
                // System.out.println("MSPE: " + MSPE + " nulls: " + nullCounterMSPE);

            }
        }
        System.out.print("Optimal minimum number of observations in each leaf: " + bestK + " MSPE: " + bestMSPE);
        System.out.println(" Optimal improvement threshold: " + bestMSEBar);

        options.setMinCount(bestK);
        options.setMinMSEImprovement(bestMSEBar);

        return options;
    }

    public void setTreeOptions(TreeOptions cvOptions) {
        this.treeOptions = cvOptions;
    }

}
