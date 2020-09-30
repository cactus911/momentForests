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

import Jama.Matrix;
import core.ContainerMoment;
import core.DataLens;
import utility.pmUtility;

/**
 *
 * @author Stephen P. Ryan <stephen.p.ryan@wustl.edu>
 */
public class ContainerRCT extends ContainerMoment {

    double mse;
    Jama.Matrix beta;
    Jama.Matrix variance;

    public ContainerRCT(DataLens lens) {
        computeBetaAndErrors(lens);
    }

    private void computeBetaAndErrors(DataLens lens) {
        double meanTreatment = 0;
        double meanControl = 0;

        int countTreatment = 0;
        int countControl = 0;
        for (int i = 0; i < lens.getNumObs(); i++) {
            // System.out.println(lens.getY(i)+" "+lens.getX(i,0)+" "+lens.getX(i,1));
            if (lens.getX(i, 0) == 0) {
                meanControl += lens.getY(i);
                countControl++;
            } else {
                meanTreatment += lens.getY(i);
                countTreatment++;
            }
        }
        
        meanControl /= countControl;
        meanTreatment /= countTreatment;

        if (countControl == 0 || countTreatment == 0) {
            // System.out.println("Setting beta to NULL");
            // System.out.println("numObs_control: "+countControl+" numObs_treatment: "+countTreatment);
            beta = null;
            variance = null;
            mse = Double.POSITIVE_INFINITY;
        } else {
            // try using OLS to get beta and SSE
            Jama.Matrix w = new Jama.Matrix(lens.getNumObs(), 1);
            Jama.Matrix Y = new Jama.Matrix(lens.getNumObs(), 1);
            Jama.Matrix X = new Jama.Matrix(lens.getNumObs(), 1);
            for (int i = 0; i < lens.getNumObs(); i++) {
                w.set(i, 0, lens.getX(i, 0));
                Y.set(i, 0, lens.getY(i));
                X.set(i, 0, lens.getX(i,1));
            }
            
            // System.out.println("data:");
            // pmUtility.prettyPrintFirstTen(pmUtility.concatMatrix(pmUtility.concatMatrix(Y, w), X));
            
            Jama.Matrix olsBeta = pmUtility.OLSsvd(w, Y, true);
            // System.out.print("beta: ");
            // pmUtility.prettyPrintVector(olsBeta);
            // System.out.println("-------------");

            // beta = new Jama.Matrix(1, 1, meanTreatment - meanControl);
            beta = new Jama.Matrix(1, 1, olsBeta.get(1, 0));

            variance = beta.times(0.0);

            // pmUtility.prettyPrintVector(Y);
            double sse = 0;
            for (int i = 0; i < lens.getNumObs(); i++) {
                if (lens.getX(i, 0) == 0) {
                    // sse += Math.pow(lens.getY(i) - meanControl, 2);
                    sse += Math.pow(lens.getY(i) - olsBeta.get(0, 0), 2);
                } else {
                    // sse += Math.pow(lens.getY(i) - meanTreatment, 2);
                    sse += Math.pow(lens.getY(i) - olsBeta.get(0, 0) - olsBeta.get(1, 0), 2);
                }
            }
            mse = sse; // why am i not dividing this by sample size? because you need to add errors of different sample sizes across splits to see where the overall improvement is, and you cannot do that if they are in means
        }

        // turns out that OLS is more efficient here (probably since we impose that the intercepts are the same across both treatment and control?)
    }

    @Override
    public Matrix getBeta() {
        return beta;
    }

    @Override
    public double getMSE() {
        return mse;
    }

    @Override
    public Matrix getVariance() {
        return variance;
    }

}
