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
package examples.linearIV;

import Jama.Matrix;
import core.ContainerMoment;
import core.DataLens;
import optimization.Uncmin_f77;
import optimization.Uncmin_methods;
import utility.pmUtility;

/**
 * Container for IV (Instrumental Variables) regression estimation.
 *
 * Key difference from ContainerLinear (OLS):
 * - Moment conditions use instruments (Z'e) rather than regressors (X'e)
 * - First column of X is the endogenous regressor
 * - Last column of X is the instrument
 * - Middle columns are exogenous regressors that serve as their own instruments
 *
 * @author Stephen P. Ryan <stephen.p.ryan@wustl.edu>
 */
public class ContainerLinearIV extends ContainerMoment implements Uncmin_methods {

    double goodnessOfFit = -666; // make sure that we give back a crazy number if it is not called
    Jama.Matrix beta;
    boolean debugVerbose = false;
    Jama.Matrix X;
    Jama.Matrix Y;

    boolean failedEstimation = false;

    boolean[] homogeneityIndex;
    Jama.Matrix homogeneityParameters;
    boolean allParametersHomogeneous;
    LinearIVMomentSpecification spec;

    public ContainerLinearIV(DataLens lens, boolean[] homogeneityIndex, Jama.Matrix homogeneityParameters, boolean allParametersHomogeneous, LinearIVMomentSpecification spec) {
        Y = lens.getY();
        X = lens.getX();

        this.spec = spec;
        this.homogeneityIndex = homogeneityIndex;
        this.homogeneityParameters = homogeneityParameters;
        this.allParametersHomogeneous = allParametersHomogeneous;
    }

    @Override
    public void computeBetaAndErrors() {
        if (Y.getRowDimension() < 30) {
            beta = null;
            goodnessOfFit = Double.POSITIVE_INFINITY;
        } else {
            try {
                int numParamsToOptimize = spec.getNumParams();
                if (!allParametersHomogeneous) {
                    for (boolean b : homogeneityIndex) {
                        if (b) {
                            numParamsToOptimize = numParamsToOptimize - 1; // homogeneous parameter imposed externally
                        }
                    }
                }

                Uncmin_f77 minimizer = new Uncmin_f77(false);
                double[] guess = new double[numParamsToOptimize + 1];

                double[] xpls = new double[numParamsToOptimize + 1];
                double[] fpls = new double[2];
                double[] gpls = new double[numParamsToOptimize + 1];
                int[] itrmcd = new int[2];
                double[][] a = new double[numParamsToOptimize + 1][numParamsToOptimize + 1];
                double[] udiag = new double[numParamsToOptimize + 1];
                double[] typsiz = new double[numParamsToOptimize + 1];
                for (int i = 1; i <= numParamsToOptimize; i++) {
                    typsiz[i] = 1.0;
                }

                double[] fscale = {0, 1.0E-8};
                int[] method = {0, 1};
                int[] iexp = {0, 0};
                int[] msg = {0, 1};
                int[] ndigit = {0, 8};
                int[] itnlim = {0, 150};
                int[] iagflg = {0, 0};
                int[] iahflg = {0, 0};
                double[] dlt = {0, 1};
                double[] gradtl = {0, 1E-8};
                double[] stepmx = {0, 1E8};
                double[] steptl = {0, 1E-8};

                minimizer.optif9_f77(numParamsToOptimize, guess, this, typsiz, fscale, method, iexp, msg, ndigit, itnlim, iagflg, iahflg, dlt, gradtl, stepmx, steptl, xpls, fpls, gpls, itrmcd, a, udiag);

                // put in something here about a failed estimation
                if (itrmcd[1] == 4 || itrmcd[1] == 5) {
                    failedEstimation = true;
                }

                // check that optimizer didn't shoot off into extremes
                for (int i = 0; i < xpls.length; i++) {
                    if (xpls[i] < -10 || xpls[i] > 10) {
                        failedEstimation = true;
                    }
                }

                // objective should be close to zero in these exactly-identified cases
                if (f_to_minimize(xpls) > 10) {
                    failedEstimation = true;
                }

                Jama.Matrix betaUncmin = new Jama.Matrix(spec.getNumParams(), 1);
                int counter = 0;

                for (int i = 0; i < spec.getNumParams(); i++) {
                    if (homogeneityIndex[i] && !allParametersHomogeneous) {
                        betaUncmin.set(i, 0, homogeneityParameters.get(i, 0));
                    } else {
                        betaUncmin.set(i, 0, xpls[counter + 1]);
                        counter++;
                    }
                }

                beta = betaUncmin.copy();
                double sse = 0;

                Jama.Matrix fit = (X.getMatrix(0, X.getRowDimension() - 1, 0, X.getColumnDimension() - 2)).times(beta);
                for (int i = 0; i < Y.getRowDimension(); i++) {
                    sse += Math.pow(Y.get(i, 0) - fit.get(i, 0), 2);
                }
                goodnessOfFit = sse;

                if (failedEstimation) {
                    beta = null;
                    goodnessOfFit = Double.POSITIVE_INFINITY;
                }

            } catch (Exception e) {
                e.printStackTrace();
                beta = null;
                goodnessOfFit = Double.POSITIVE_INFINITY;
            }
        }
    }

    /**
     * Accumulate moment conditions using instruments (skipping endogenous regressor).
     *
     * @param g Matrix of moments to augment
     * @param error_i Error of ith observation at the current guess of parameters
     * @param i Observation to compute the moment
     */
    public void addGi(Jama.Matrix g, double error_i, int i) {
        /**
         * We are going to skip the first column of X as that is the endogenous
         * regressor
         */
        for (int k = 1; k < X.getColumnDimension(); k++) {
            g.set(k - 1, 0, g.get(k - 1, 0) + error_i * X.get(i, k));
        }
    }

    @Override
    public Matrix getGi(Matrix beta, int i) {
        double fit = (X.getMatrix(i, i, 0, X.getColumnDimension() - 2).times(beta)).get(0, 0);
        double error = fit - Y.get(i, 0);

        Jama.Matrix gi = new Jama.Matrix(X.getColumnDimension() - 1, 1);

        for (int k = 1; k < X.getColumnDimension(); k++) {
            gi.set(k - 1, 0, error * X.get(i, k));
        }

        return gi;
    }

    @Override
    public Matrix getJacobianNoDivision(Matrix beta) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public Jama.Matrix getMomentGWithoutDivision(Jama.Matrix beta) {
        int numMoments = spec.getNumMoments();
        Jama.Matrix g = new Jama.Matrix(numMoments, 1);
        Jama.Matrix fittedY = (X.getMatrix(0, X.getRowDimension() - 1, 0, X.getColumnDimension() - 2)).times(beta);
        Jama.Matrix e = fittedY.minus(Y);

        for (int i = 0; i < X.getRowDimension(); i++) {
            addGi(g, e.get(i, 0), i);
        }
        return g;
    }

    @Override
    public double computeMeasureOfFit(Jama.Matrix beta) {
        Jama.Matrix fittedY = (X.getMatrix(0, X.getRowDimension() - 1, 0, X.getColumnDimension() - 2)).times(beta);
        Jama.Matrix e = fittedY.minus(Y);

        return e.normF();
    }

    private double getMomentObjectiveFunction(Jama.Matrix beta, boolean debugMoment) {
        int numMoments = spec.getNumMoments();
        Jama.Matrix g = new Jama.Matrix(numMoments, 1);

        Jama.Matrix fittedY = (X.getMatrix(0, X.getRowDimension() - 1, 0, X.getColumnDimension() - 2)).times(beta);

        Jama.Matrix e = fittedY.minus(Y);

        for (int i = 0; i < X.getRowDimension(); i++) {
            addGi(g, e.get(i, 0), i);
        }
        g.timesEquals(1.0 / Y.getRowDimension());
        double q = ((g.transpose()).times(g)).get(0, 0); // gmm with identity weighting matrix

        return q;
    }

    @Override
    public Matrix getBeta() {
        return beta;
    }

    @Override
    public double getGoodnessOfFit() {
        return goodnessOfFit;
    }

    @Override
    public Matrix getVariance(Matrix beta) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public double getMomentFunctionImposingHomogeneity(int k, double value) {
        Jama.Matrix betaHomogeneous = beta.copy();
        betaHomogeneous.set(k, 0, value);
        return getMomentObjectiveFunction(betaHomogeneous, false);
    }

    @Override
    public double f_to_minimize(double[] x) {

        int counter = 0;
        Jama.Matrix b = new Jama.Matrix(spec.getNumParams(), 1);
        for (int i = 0; i < spec.getNumParams(); i++) {
            if (homogeneityIndex[i] && !allParametersHomogeneous) {
                b.set(i, 0, homogeneityParameters.get(i, 0));
            } else {
                b.set(i, 0, x[counter + 1]);
                counter++;
            }
        }

        return getMomentObjectiveFunction(b, false);
    }

    @Override
    public void gradient(double[] x, double[] g) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public void hessian(double[] x, double[][] h) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public double getMomentFunctionValue(Jama.Matrix b) {
        return getMomentObjectiveFunction(b, false);
    }

    boolean didEstimatorFail() {
        return failedEstimation;
    }

}
