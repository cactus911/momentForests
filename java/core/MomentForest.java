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

import com.stata.sfi.*;
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
        double proportionObservationsToEstimateTreeStructure = 0.5;

        Random rng = new Random(forestSeed);

        for (int i = 0; i < numberTreesInForest; i++) {
            // resample the forestLens, then split it
            DataLens resampled = forestLens.getResampledDataLensWithBalance(rng.nextLong());
            DataLens[] split = resampled.randomlySplitSample(proportionObservationsToEstimateTreeStructure, rng.nextLong());
            DataLens lensGrow = split[0];
            DataLens lensHonest = split[1];

            forest.add(new TreeMoment(null, spec, lensGrow,
                    spec.getDiscreteVector(), verbose, treeOptions.getMinProportion(), treeOptions.getMinCount(), treeOptions.getMinMSEImprovement(), true, treeOptions.getMaxDepth(),
                    lensHonest));
        }
        SFIToolkit.displayln("test.getMincount " + treeOptions.getMinCount() + "test.MinMSEImprove " + treeOptions.getMinMSEImprovement() );
        forest.parallelStream().forEach((tree) -> {
            tree.determineSplit();
            tree.estimateHonestTree();
        });
        // System.out.format("Memory usage: %,d bytes %n", (Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory()));
    }

    public Jama.Matrix getEstimatedParameters(Jama.Matrix x) {
        Jama.Matrix estimatedParameters = forest.get(0).getEstimatedBeta(x);
        String s = "[ " + forest.get(0).getEstimatedBeta(x).get(0, 0) + " ";
        for (int i = 1; i < forest.size(); i++) {
            estimatedParameters = estimatedParameters.plus(forest.get(i).getEstimatedBeta(x));
            s = s.concat(forest.get(i).getEstimatedBeta(x).get(0, 0) + " ");
        }
        s = s.concat("]");
        System.out.println(s);
        estimatedParameters.timesEquals(1.0 / forest.size());
        return estimatedParameters;
    }

    public TreeOptions performCrossValidation(int numTrees, Jama.Matrix CVparameters2) {
        Random rng = new Random(forestSeed);
        DataLens[] split = forestLens.randomlySplitSample(0.5, rng.nextLong());
        DataLens growLens = split[0];
        DataLens predictLens = split[1];

        // System.out.println(growLens.getNumObs()+" "+predictLens.getNumObs());
        // System.out.println("sum treeX: " + pmUtility.sum(treeX, 0) + " mean treeX: " + pmUtility.mean(treeX, 0) + " mean honestX: " + pmUtility.mean(honestX, 0) + " mean predictX: " + pmUtility.mean(predictX, 0));
        // System.exit(0);
        boolean verboseTreeConstruction = false;
        double bestMSPE = 0;
        double bestMSEBar = 0;
        boolean first = true;
        int bestK = 10;

        TreeOptions options = new TreeOptions();
        options.setMaxDepth(100);
        options.setMinMSEImprovement(1E-10);
        
        // for (double improvementThreshold = (double) CVparameters2.get(0,3); improvementThreshold <= 0.5; improvementThreshold += 0.05) {
           // for (int minCountEachPartition = (int) CVparameters2.get(0,0); minCountEachPartition <= (int) CVparameters2.get(0,2) / 2; minCountEachPartition += (int) CVparameters2.get(0,1) ) {
        // for (double improvementThreshold = 0.01; improvementThreshold <= 0.5; improvementThreshold += 0.05) {
                 // for (int minCountEachPartition = 10; minCountEachPartition <= growLens.getNumObs() / 2; minCountEachPartition *= 2) {
         for (double improvementThreshold = (double) CVparameters2.get(0,3); improvementThreshold <= (double) CVparameters2.get(0,5); improvementThreshold += (double) CVparameters2.get(0,4)) {
            for (int minCountEachPartition = (int) CVparameters2.get(0,0); minCountEachPartition <= (int) CVparameters2.get(0,2) /*growLens.getNumObs() / 2*/; minCountEachPartition += (int) CVparameters2.get(0,1) /**= 2*/) {
                rng = new Random(667);
                options.setMinCount(minCountEachPartition);
                options.setMinMSEImprovement(improvementThreshold);

                MomentForest momentForest = new MomentForest(spec, numTrees, rng.nextLong(), growLens, verboseTreeConstruction, options);
                momentForest.growForest();
                /**
                 * We can easily implement cross-fitting here by swapping the
                 * prediction and estimation data sets and adding another MSPE
                 * component
                 */
                MomentForest momentForestSwitch = new MomentForest(spec, numTrees, rng.nextLong(), predictLens, verboseTreeConstruction, options);
                momentForestSwitch.growForest();

                double MSPE = 0;

                int counter = 0;
                int nullCounterMSPE = 0;
                for (int i = 0; i < predictLens.getNumObs(); i++) {
                    Jama.Matrix xi = predictLens.getRowAsJamaMatrix(i);
                    Double predictedY = spec.getPredictedY(xi, momentForest.getEstimatedParameters(xi));
                    if (predictedY != null) {
                        MSPE += Math.pow(predictLens.getY(i) - predictedY, 2);
                        counter++;
                    } else {
                        pmUtility.prettyPrint(xi);
                        nullCounterMSPE++;
                    }
                }
                /**
                 * Cross-fitting part (reverse estimation data with prediction
                 * data)
                 */
                for (int i = 0; i < growLens.getNumObs(); i++) {
                    Jama.Matrix xi = growLens.getRowAsJamaMatrix(i);
                    Double predictedY = spec.getPredictedY(xi, momentForestSwitch.getEstimatedParameters(xi));
                    if (predictedY != null) {
                        MSPE += Math.pow(growLens.getY(i) - predictedY, 2);
                        counter++;
                    } else {
                        pmUtility.prettyPrint(xi);
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
                if (MSPE <= bestMSPE || first) {
                    s = " * ";
                    bestK = minCountEachPartition;
                    // best = minProportionEachPartition;
                    bestMSEBar = improvementThreshold;
                    // bestDepth = maxDepth;
                    bestMSPE = MSPE;
                    // bestProportion = proportion;
                    first = false;

                    System.out.format("maxDepth: %d k: %d mse_bar: %g mspe: %g nulls: %d %s %n", options.getMaxDepth(), options.getMinCount(), options.getMinMSEImprovement(), MSPE, nullCounterMSPE, s);
                    System.out.println("MSPE: " + MSPE + " nulls: " + nullCounterMSPE);
                    System.out.println("Example tree:");
                    momentForest.getTree(0).printTree();
                    System.out.println("---------------");
                }

            }
        }
        SFIToolkit.displayln("Optimal minimum number of observations in each leaf: " + bestK + " MSPE: " + bestMSPE);
        SFIToolkit.displayln(" Optimal improvement threshold: " + bestMSEBar);
        System.out.print("Optimal minimum number of observations in each leaf: " + bestK + " MSPE: " + bestMSPE);
        System.out.println(" Optimal improvement threshold: " + bestMSEBar);

        options.setMinCount(bestK);
        options.setMinMSEImprovement(bestMSEBar);

        return options;
    }
    
    public static double roundAvoid(double value, int places) {
    double scale = Math.pow(10, places);
    return Math.round(value * scale) / scale;
    }
    
    public void setTreeOptions(TreeOptions cvOptions) {
        this.treeOptions = cvOptions;
    }

}
