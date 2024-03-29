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
import utility.pmUtility;

/**
 *
 * @author Stephen P. Ryan <stephen.p.ryan@wustl.edu>
 */

/*
    TO DO:
        1. BootstrapForest takes a balancing vector as input and uses it for each forest in the bootstrap. See lines 45 and 55.
*/
public class BootstrapForest {

    ArrayList<MomentForest> forestList = new ArrayList<>();

    public BootstrapForest(MomentSpecification spec, int numberBootstraps, int numberTreesInForest, long randomSeed, TreeOptions options) {
        // DataLens originalLens = new DataLens(spec.getX(), spec.getY(), pmUtility.getColumn(spec.getX(), 0));
        DataLens originalLens = new DataLens(spec.getX(), spec.getY(), spec.getBalancingVector());
        
        Random rng = new Random(randomSeed);
        for (int i = 0; i < numberBootstraps; i++) {
            long seed = rng.nextLong();
//            forestList.add(new MomentForest(spec, numberTreesInForest, rng.nextLong(), 
//                    MomentForest.resample(spec.getX(), seed, pmUtility.getColumn(spec.getX(), 0)), 
//                    MomentForest.resample(spec.getY(), seed, pmUtility.getColumn(spec.getX(), 0)), 
//                    false, options));
            forestList.add(new MomentForest(spec, numberTreesInForest, rng.nextLong(),
                    originalLens.getResampledDataLensWithBalance(seed),
                    false, options));
        }
        forestList.parallelStream().forEach((forest) -> forest.growForest());
    }

    public Jama.Matrix computeStandardErrors(Matrix xi) {
        Jama.Matrix results = new Jama.Matrix(forestList.size(), forestList.get(0).getEstimatedParameters(xi).getRowDimension());
        for (int r = 0; r < forestList.size(); r++) {
            Jama.Matrix estimatedParameter = forestList.get(r).getEstimatedParameters(xi);
            for (int i = 0; i < estimatedParameter.getRowDimension(); i++) {
                results.set(r, i, estimatedParameter.get(i, 0));
            }
        }

        Jama.Matrix standardErrors = new Jama.Matrix(results.getColumnDimension(), 1);
        for (int i = 0; i < standardErrors.getRowDimension(); i++) {
            standardErrors.set(i, 0, pmUtility.standardDeviation(results, i));
        }
        // pmUtility.prettyPrintVector(results);
        return standardErrors;
    }

}
