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
public class BootstrapForest {

    ArrayList<MomentForest> forestList = new ArrayList<>();

    public BootstrapForest(MomentSpecification spec, int numberBootstraps, int numberTreesInForest, long randomSeed, TreeOptions options) {
        DataLens originalLens = new DataLens(spec.getX(), spec.getY(), pmUtility.getColumn(spec.getX(), 0));
        
        Random rng = new Random(randomSeed);
        for (int i = 0; i < numberBootstraps; i++) {
            long seed = rng.nextLong();
            forestList.add(new MomentForest(spec, numberTreesInForest, rng.nextLong(),
                    // originalLens.getResampledDataLensWithBalance(seed),
                    originalLens,
                    false, options));
        }
        forestList.parallelStream().forEach((forest) -> forest.growForest());
    }

    /**
     * This method returns a Jama.Matrix of standard errors. 
     * 
     * @param xi Vector of observables that we will compute standard errors at
     * @return
     */
    public Jama.Matrix computeStandardErrors(Matrix xi) {
        // Jama.Matrix beta = forestList.get(0).getEstimatedParameters(xi);
        // pmUtility.prettyPrint(beta);
        // System.out.println(beta.getRowDimension()+" by " +beta.getColumnDimension());
        
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
        pmUtility.prettyPrintVector(results);
        pmUtility.prettyPrintVector(standardErrors);
        
        return standardErrors;
    }
    
    public Jama.Matrix computeStandardErrorsSecondWay(Matrix xi) {
        // Jama.Matrix beta = forestList.get(0).getEstimatedParameters(xi);
        // pmUtility.prettyPrint(beta);
        // System.out.println(beta.getRowDimension()+" by " +beta.getColumnDimension());
        
        // take a single forest, look at how the parameter vector varies across trees within that forest to get a within-forest SE
        // maybe repeat with all forests, average the SE across forests
        
        Jama.Matrix results = new Jama.Matrix(forestList.size(), forestList.get(0).getEstimatedParameters(xi).getRowDimension());
        for (int r = 0; r < forestList.size(); r++) {
            Jama.Matrix estimatedParameter = forestList.get(r).getEstimatedStandardDeviation(xi);
            for (int i = 0; i < estimatedParameter.getRowDimension(); i++) {
                results.set(r, i, estimatedParameter.get(i, 0));
            }
        }

        Jama.Matrix standardErrors = new Jama.Matrix(results.getColumnDimension(), 1);
        for (int i = 0; i < standardErrors.getRowDimension(); i++) {
            standardErrors.set(i, 0, pmUtility.mean(results, i));
        }
        // pmUtility.prettyPrintVector(results);
        // pmUtility.prettyPrintVector(standardErrors);
        
        return standardErrors;
    }

}
