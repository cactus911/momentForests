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
package examples.CardIV;

import Jama.Matrix;
import core.ContainerMoment;
import core.DataLens;
import optimization.Uncmin_f77;
import optimization.Uncmin_methods;
import utility.pmUtility;

/**
 *
 * @author Stephen P. Ryan <stephen.p.ryan@wustl.edu>
 */
public class ContainerCardIV extends ContainerMoment implements Uncmin_methods {

    double goodnessOfFit = -666; // make sure that we give back a crazy number if it is not called
    Jama.Matrix beta;
    boolean debugVerbose = false;
    Jama.Matrix X;
    Jama.Matrix Y;

    boolean failedEstimation = false;

    boolean[] homogeneityIndex;
    Jama.Matrix homogeneityParameters;
    boolean allParametersHomogeneous;
    MomentSpecificationCardIV spec;

    public ContainerCardIV(DataLens lens, boolean[] homogeneityIndex, Jama.Matrix homogeneityParameters, boolean allParametersHomogeneous, MomentSpecificationCardIV spec) {
        Y = lens.getY();
        X = lens.getX();

        this.spec = spec;
        this.homogeneityIndex = homogeneityIndex;
        this.homogeneityParameters = homogeneityParameters;
        this.allParametersHomogeneous = allParametersHomogeneous;
    }

    @Override
    public void computeBetaAndErrors() {
        // System.out.println("In here");

        if (Y.getRowDimension() < 30) {
            // System.out.println("Too few observations");

            beta = null;
            goodnessOfFit = Double.POSITIVE_INFINITY;
        } else {
            // System.out.println("Minimizer");
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

                // System.out.print("Composite beta inside ContainerCardIV: ");
                // pmUtility.prettyPrintVector(beta);
                // System.out.println(X.getRowDimension()+" "+X.getColumnDimension());
                Jama.Matrix fit = (X.getMatrix(0, X.getRowDimension() - 1, 0, X.getColumnDimension() - 2)).times(beta);
                for (int i = 0; i < Y.getRowDimension(); i++) {
                    sse += Math.pow(Y.get(i, 0) - fit.get(i, 0), 2);
                }
                goodnessOfFit = sse;
                // goodnessOfFit = f_to_minimize(xpls);
//                    System.out.print("\t\tFound uncmin beta: ");
//                    pmUtility.prettyPrintVector(betaUncmin);
//                    Jama.Matrix betaOLS = pmUtility.OLSsvd(X, Y, false);
//                    System.out.print("\t\tCompared to betaOLS: ");
//                    pmUtility.prettyPrintVector(betaOLS);

                if (failedEstimation) {
                    System.out.print("Optimizer flagged as failing [" + itrmcd[1] + "]: ");
                    pmUtility.prettyPrint(new Jama.Matrix(xpls, 1));
                    System.out.println("f_failed: " + f_to_minimize(xpls));
                    beta = null;
                    goodnessOfFit = Double.POSITIVE_INFINITY;
                }

            } catch (Exception e) {
                e.printStackTrace();
                System.exit(0);
                if (allParametersHomogeneous) {
                    e.printStackTrace();
                }
                if (debugVerbose) {
                    System.out.println("Matrix not invertible");
                    //System.exit(0);
                }
                beta = null;
                goodnessOfFit = Double.POSITIVE_INFINITY;
            }
        }
    }

    /**
     * This is the much faster way of constructing moments, since I am computing
     * errors once and then passing them. Could still make this faster, I'm
     * sure, but re-using matrices and so forth.
     *
     * @param g Matrix of moments to augment
     * @param error_i Error of ith observation at the current guess of
     * parameters
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

    /**
     * This is supremely slow, so use the other approach (passing errors), but
     * this is what I need to generalize how the covariance matrix is
     * constructed.
     *
     * @param beta Parameter vector
     * @param i Observation to evaluate
     * @return
     */
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

    /**
     *
     * NOT IMPORTANT: only used in Wald test for forest, which is no longer used
     * (should get rid of these, actually)
     *
     * @param beta
     * @return
     */
    @Override
    public Matrix getJacobianNoDivision(Matrix beta) {
        throw new UnsupportedOperationException("Not supported yet."); // Generated from nbfs://nbhost/SystemFileSystem/Templates/Classes/Code/GeneratedMethodBody
    }

    @Override
    public Jama.Matrix getMomentGWithoutDivision(Jama.Matrix beta) {
        /*
         * This version implements the moment-based GMM for IV regression.
         * The moment conditions are based on instruments'e, where instruments are the instruments, and e is the residual.
         */
        int numMoments = spec.getNumMoments(); // Number of instruments (columns of instruments)
        Jama.Matrix g = new Jama.Matrix(numMoments, 1); // x'e, one row for each x
        // (X.getMatrix(i, i, 0, X.getColumnDimension() - 2).times(beta)).get(0, 0);
        Jama.Matrix fittedY = (X.getMatrix(0, X.getRowDimension() - 1, 0, X.getColumnDimension() - 2)).times(beta);
        Jama.Matrix e = fittedY.minus(Y);

        for (int i = 0; i < X.getRowDimension(); i++) {
            // Instead of X'e (as in OLS), use instruments'e (in IV)
            addGi(g, e.get(i, 0), i);
        }
        return g;
    }

    @Override
    public double computeMeasureOfFit(Jama.Matrix beta) {
        // just return the GMM objective function (I think this works, no?)
        // this isn't anything critical, so don't worry about it

        // actually, i think that it is better to use SSE
        Jama.Matrix fittedY = (X.getMatrix(0, X.getRowDimension() - 1, 0, X.getColumnDimension() - 2)).times(beta);
        Jama.Matrix e = fittedY.minus(Y);

        return e.normF();
    }

    private double getMomentObjectiveFunction(Jama.Matrix beta, boolean debugMoment) {
        // Jama.Matrix runningTotal = new Jama.Matrix(X.getRowDimension(), X.getColumnDimension());
        int numMoments = spec.getNumMoments();
        Jama.Matrix g = new Jama.Matrix(numMoments, 1); // x'e, one row for each x

        // Jama.Matrix fittedY = X.times(beta);
        Jama.Matrix fittedY = (X.getMatrix(0, X.getRowDimension() - 1, 0, X.getColumnDimension() - 2)).times(beta);

        Jama.Matrix e = fittedY.minus(Y);
        // Jama.Matrix omega = new Jama.Matrix(numMoments, numMoments);

        // turns out using gi method here is crazy slow!        
        for (int i = 0; i < X.getRowDimension(); i++) {
            // Jama.Matrix gi = getGi(e.get(i, 0), i);
//            pmUtility.prettyPrintVector(gi);
            // g.plusEquals(gi);
            addGi(g, e.get(i, 0), i);
            // omega.plusEquals(gi.times(gi.transpose()));
        }
        // omega.timesEquals(1.0 / Y.getRowDimension());
        g.timesEquals(1.0 / Y.getRowDimension());
        // double q = (((g.transpose()).times(omega.inverse())).times(g)).get(0, 0); // this is the continuous updating estimator (CUE)
        // that appears to sometimes generate perverse decreases in fit when splitting (probably due to some numerical instability in inversion, and the confounding of fits versus variance)
        double q = ((g.transpose()).times(g)).get(0, 0); // this is gmm with identity weighting matrix

        // pmUtility.prettyPrintVector(beta);
        // pmUtility.prettyPrint(g);
        // key point here is that we estimate the model using GMM
        // but we report goodness-of-fit when searching over the splits
//        System.exit(0);
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

    /**
     * SHOULDN'T BE NEEDED
     *
     * @param beta
     * @return
     */
    @Override
    public Matrix getVariance(Matrix beta) {
        throw new UnsupportedOperationException("Not supported yet."); // Generated from nbfs://nbhost/SystemFileSystem/Templates/Classes/Code/GeneratedMethodBody
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

//        System.out.print("Compose beta in f_to_minimize ContainerCardIV: ");
//        pmUtility.prettyPrintVector(b);
        return getMomentObjectiveFunction(b, false);
    }

    @Override
    public void gradient(double[] x, double[] g) {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    public void hessian(double[] x, double[][] h) {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    public double getMomentFunctionValue(Jama.Matrix b) {
        return getMomentObjectiveFunction(b, false);
    }

    boolean didEstimatorFail() {
        return failedEstimation;
    }

}
