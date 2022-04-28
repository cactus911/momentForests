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
package examples.logit;

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
public class ContainerLogit extends ContainerMoment implements Uncmin_methods {

    double goodnessOfFit = -666; // make sure that we give back a crazy number if it is not called
    Jama.Matrix containerBeta;
    Jama.Matrix containerVariance;
    boolean debugVerbose = false;
    DataLens lens;
    Jama.Matrix X;
    Jama.Matrix Y;
    boolean[] homogeneityIndex;
    Jama.Matrix homogeneityParameters;
    boolean allParametersHomogeneous;

    public ContainerLogit(DataLens lens, boolean[] homogeneityIndex, Jama.Matrix homogeneityParameters, boolean allParametersHomogeneous) {
        this.lens = lens;
        X = lens.getX();
        Y = lens.getY();
        this.homogeneityIndex = homogeneityIndex;
        this.homogeneityParameters = homogeneityParameters;
        this.allParametersHomogeneous = allParametersHomogeneous;
        // computeBetaAndErrors();
    }

    public void computeBetaAndErrors() {
        if (Y.getRowDimension() < 30) {
            // System.out.println("Too few observations");
            containerBeta = null;
            goodnessOfFit = Double.POSITIVE_INFINITY;
        } else {
            try {
                Uncmin_f77 minimizer = new Uncmin_f77(false);

                int numParams = X.getColumnDimension();
                if (!allParametersHomogeneous) {
                    for (boolean b : homogeneityIndex) {
                        if (b) {
                            numParams = numParams - 1; // homogeneous parameter imposed externally
                        }
                    }
                }
                // System.out.println("numParams: " + numParams);

                double[] guess = new double[numParams + 1];

                double[] xpls = new double[numParams + 1];
                double[] fpls = new double[2];
                double[] gpls = new double[numParams + 1];
                int[] itrmcd = new int[2];
                double[][] a = new double[numParams + 1][numParams + 1];
                double[] udiag = new double[numParams + 1];
                double[] typsiz = new double[numParams + 1];
                for (int i = 1; i <= numParams; i++) {
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
                double[] stepmx = {0, 1.0};
                double[] steptl = {0, 1E-8};
                // System.out.println("In uncmin");
                minimizer.optif9_f77(numParams, guess, this, typsiz, fscale, method, iexp, msg, ndigit, itnlim, iagflg, iahflg, dlt, gradtl, stepmx, steptl, xpls, fpls, gpls, itrmcd, a, udiag);
                // System.out.println("Out of uncmin");

                Jama.Matrix betaUncmin = new Jama.Matrix(X.getColumnDimension(), 1);

                int counter = 0;

                for (int i = 0; i < X.getColumnDimension(); i++) {
                    if (homogeneityIndex[i]) {
                        betaUncmin.set(i, 0, homogeneityParameters.get(i, 0));
                    } else {
                        betaUncmin.set(i, 0, xpls[counter + 1]);
                        counter++;
                    }
                }

                containerBeta = betaUncmin.copy();

                /**
                 * This is super important: the objective function used to
                 * estimate the parameters and the goodness of fit may be VERY
                 * different objects! We need something here that is additive
                 * across subsplits so that we get an apples-to-apples
                 * comparison of the fit when thinking about where to split on
                 * Z!
                 *
                 * In containerlinear, for example, this is the sum of squared
                 * errors. Here, maybe use the LLH, even when the parameters are
                 * estimated using moments? This is ONLY to be used for the
                 * splitting criterion.
                 */
                // goodnessOfFit = f_to_minimize(xpls);
                goodnessOfFit = computeLLH(betaUncmin);

                if (debugVerbose) {
                    System.out.format("ContainerLogit.computeBetaAndErrors SSE: %g ", +goodnessOfFit);
                    pmUtility.prettyPrintVector(containerBeta);
                    System.exit(0);
                }
            } catch (Exception e) {
                e.printStackTrace();
                if (debugVerbose) {
                    System.out.println("Matrix not invertible");
                }
                containerBeta = null;
                goodnessOfFit = Double.POSITIVE_INFINITY;
            }
        }
    }

    @Override
    public Jama.Matrix getGi(Jama.Matrix beta, int i) {
        int numMoments = X.getColumnDimension();

        // try using the derivatives with respect to beta as the moments here (two X's)
        if (numMoments != 2) {
            System.out.println("Hardwired moments for two parameters in ContainerLogit.java");
            System.exit(0);
        }

        double utility = X.get(i, 0) * beta.get(0, 0) + X.get(i, 1) * beta.get(1, 0);
        double shareInside = Math.exp(utility) / (1.0 + Math.exp(utility));

        Jama.Matrix gi = new Jama.Matrix(numMoments, 1);

        boolean useScores = false;
        if (useScores) {
            if (Y.get(i, 0) == 1) {
                gi.set(0, 0, X.get(i, 0) * (1.0 - shareInside));
                gi.set(1, 0, X.get(i, 1) * (1.0 - shareInside));
            } else {
                gi.set(0, 0, -X.get(i, 0) * shareInside);
                gi.set(1, 0, -X.get(i, 1) * shareInside);
            }
        } else {
            double ei = Y.get(i, 0) - shareInside;
            gi.set(0, 0, X.get(i, 0) * ei);
            gi.set(1, 0, X.get(i, 1) * ei);
        }

        return gi;
    }

    @Override
    public Jama.Matrix getMomentGWithoutDivision(Jama.Matrix beta) {
        // Jama.Matrix runningTotal = new Jama.Matrix(X.getRowDimension(), X.getColumnDimension());
        int numMoments = X.getColumnDimension();
        Jama.Matrix g = new Jama.Matrix(numMoments, 1); // x'e, one row for each x

        // turns out using gi method here is crazy slow!        
        for (int i = 0; i < X.getRowDimension(); i++) {
            Jama.Matrix gi = getGi(beta, i);
            g.plusEquals(gi);
        }
        // cannot have this here since we divide by different n in different places!
        // g.timesEquals(1.0 / Y.getRowDimension());
        return g;
    }

    private double getMomentObjectiveFunctionValue(Jama.Matrix beta) {
        // LLH is much much faster
        boolean useGMM = false;
        if (useGMM) {
            int numMoments = X.getColumnDimension();
            Jama.Matrix g = new Jama.Matrix(numMoments, 1); // x'e, one row for each x
            Jama.Matrix omega = new Jama.Matrix(numMoments, numMoments);

            // turns out using gi method here is crazy slow!        
            for (int i = 0; i < X.getRowDimension(); i++) {
                Jama.Matrix gi = getGi(beta, i);
                g.plusEquals(gi);
                omega.plusEquals(gi.times(gi.transpose()));
            }
            omega.timesEquals(1.0 / Y.getRowDimension());
            g.timesEquals(1.0 / Y.getRowDimension());
            // double q = (((g.transpose()).times(omega.inverse())).times(g)).get(0, 0); // this is the continuous updating estimator (CUE)
            // that appears to sometimes generate perverse decreases in fit when splitting (probably due to some numerical instability in inversion, and the confounding of fits versus variance)
            double q = ((g.transpose()).times(g)).get(0, 0); // this is gmm with identity weighting matrix

            return q;
        }

        return computeLLH(beta);
    }

    private double computeLLH(Jama.Matrix beta) {
        double llh = 0;
        for (int i = 0; i < X.getRowDimension(); i++) {
            double u = beta.get(0, 0) * X.get(i, 0) + beta.get(1, 0) * X.get(i, 1);
            double insideShare = Math.exp(u) / (1.0 + Math.exp(u));
            if (Y.get(i, 0) == 1) {
                llh += Math.log(insideShare);
            } else {
                llh += Math.log(1.0 - insideShare);
            }
        }
        return -llh;
    }

    public static double computeLLHi(double y, Jama.Matrix xi, Jama.Matrix beta) {
        double llh = 0;
        double u = beta.get(0, 0) * xi.get(0, 0) + beta.get(1, 0) * xi.get(0, 1);
        double insideShare = Math.exp(u) / (1.0 + Math.exp(u));
        if (y == 1) {
            llh += Math.log(insideShare);
        } else {
            llh += Math.log(1.0 - insideShare);
        }

        return -llh;
    }

    @Override
    public Matrix getBeta() {
        return containerBeta;
    }

    @Override
    public double getGoodnessOfFit() {
        return goodnessOfFit;
    }

    @Override
    public double getMomentFunctionImposingHomogeneity(int k, double value) {
        Jama.Matrix betaHomogeneous = containerBeta.copy();
        betaHomogeneous.set(k, 0, value);
        return getMomentObjectiveFunctionValue(betaHomogeneous);
    }

    @Override
    public double f_to_minimize(double[] x) {
        /**
         * Is this the place to impose homogeneity? Just read off the x[] vector
         * one at a time (we do this elsewhere when searching over the
         * homogeneous parameters)
         */
        int counter = 0;
        Jama.Matrix b = new Jama.Matrix(X.getColumnDimension(), 1);
        for (int i = 0; i < X.getColumnDimension(); i++) {
            if (homogeneityIndex[i]) {
                b.set(i, 0, homogeneityParameters.get(i, 0));
            } else {
                b.set(i, 0, x[counter + 1]);
                counter++;
            }
        }

        boolean debug = true;
        if (debug) {
            boolean imposedHomogeneity = false;
            for (boolean bp : homogeneityIndex) {
                if (bp) {
                    imposedHomogeneity = true;
                }
            }
            if (imposedHomogeneity) {
                // System.out.println("-----");
                for (boolean bp : homogeneityIndex) {
//                     System.out.print(bp + " ");
                }
//                 System.out.println("");
//                 pmUtility.prettyPrintVector(homogeneityParameters);
//                 pmUtility.prettyPrint(new Jama.Matrix(x, 1));
                // pmUtility.prettyPrintVector(b);
                // System.out.println("Computing LLH");
                // double llh = getMomentObjectiveFunctionValue(b);
                // System.out.println("Back from that");
                // System.exit(0);
            }
        }

        return getMomentObjectiveFunctionValue(b);
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
        return getMomentObjectiveFunctionValue(b);
    }

    @Override
    public double computeMeasureOfFit(Matrix beta) {
        return computeLLH(beta);
    }

    @Override
    public Matrix getVariance(Matrix beta) {
        throw new UnsupportedOperationException("Not supported yet."); // Generated from nbfs://nbhost/SystemFileSystem/Templates/Classes/Code/GeneratedMethodBody
    }

    @Override
    public Matrix getJacobianNoDivision(Matrix beta) {
        throw new UnsupportedOperationException("Not supported yet."); // Generated from nbfs://nbhost/SystemFileSystem/Templates/Classes/Code/GeneratedMethodBody
    }

}
