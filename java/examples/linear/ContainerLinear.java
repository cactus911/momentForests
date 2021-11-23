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
package examples.linear;

import Jama.Matrix;
import core.ContainerMoment;
import core.DataLens;
import utility.pmUtility;

/**
 *
 * @author Stephen P. Ryan <stephen.p.ryan@wustl.edu>
 */
public class ContainerLinear extends ContainerMoment {

    double mse;
    Jama.Matrix beta;
    Jama.Matrix variance;
    boolean debugVerbose = false;

    public ContainerLinear(DataLens lens) {
        computeBetaAndErrors(lens);
    }

    private void computeBetaAndErrors(DataLens lens) {
        Jama.Matrix X = lens.getX();
        Jama.Matrix Y = lens.getY();

//        pmUtility.prettyPrint(pmUtility.concatMatrix(Y, X));
//        beta = pmUtility.OLS(X, Y, false);
//        System.out.print("beta OLS: ");
//        pmUtility.prettyPrintVector(beta);
        // System.exit(0);
        
        /**
         * Need to implement filter here to separate X and Z for putting in the
         * OLS that avoids the whole issue of splitting these matrices over and
         * over
         *
         * Done that now. Just have to get it working properly again.
         */
        try {
            beta = pmUtility.OLSsvd(X, Y, false);

            Jama.Matrix fittedY = X.times(beta);
            double sse = 0;
            for (int i = 0; i < fittedY.getRowDimension(); i++) {
                sse += Math.pow(Y.get(i, 0) - fittedY.get(i, 0), 2);
            }
            mse = sse; // cannot divide by n since we have unbalanced samples in each leaf
            if (debugVerbose) {
                System.out.format("sse: %g ", +sse);
                pmUtility.prettyPrintVector(beta);
            }
        } catch (Exception e) {
            // e.printStackTrace();
            if (debugVerbose) {
                System.out.println("Matrix not invertible");
            }
            beta = null;
            mse = Double.POSITIVE_INFINITY;
        }
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
