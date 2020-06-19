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

/**
 * Class containing a collection of trees and utility classes for interacting
 * with it
 *
 * @author Stephen P. Ryan <stephen.p.ryan@wustl.edu>
 */
public class MomentForest {

    ArrayList<TreeMoment> forest;

    public MomentForest(MomentSpecification spec, int numberTreesInForest, long forestSeed, Jama.Matrix forestX, Jama.Matrix forestY) {
        double minProportion = 0.0001;
        int minCount = 50;
        double minMSEImprovement = 1E-5;
        int maxDepth = 10;

        int numObs = forestX.getRowDimension();
        int numObsToEstimateTreeStructure = (int) Math.floor(numObs * 0.5);

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
        boolean verbose = false;

        for (int i = 0; i < numberTreesInForest; i++) {
            long seed = rng.nextLong();
            long seedHonest = rng.nextLong();
            forest.add(new TreeMoment(null, spec, resample(treeX, seed), resample(treeY, seed),
                    spec.getDiscreteVector(), verbose, minProportion, minCount, minMSEImprovement, true, maxDepth,
                    resample(honestX, seedHonest), resample(honestY, seedHonest)));
        }
    }

    public void growForest() {
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
        for (int i = 1; i < forest.size(); i++) {
            estimatedParameters.plusEquals(forest.get(i).getEstimatedBeta(x));
        }
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

}
