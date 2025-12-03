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
import core.MomentSpecification;

import java.awt.BorderLayout;
import javax.swing.JFrame;
import optimization.Uncmin_f77;
import optimization.Uncmin_methods;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;
import utility.pmUtility;

/**
 *
 * @author Stephen P. Ryan <stephen.p.ryan@wustl.edu>
 */
public class ContainerLinear extends ContainerMoment implements Uncmin_methods {

    double goodnessOfFit = -666; // make sure that we give back a crazy number if it is not called
    Jama.Matrix beta;
    Jama.Matrix variance;
    boolean debugVerbose = false;
    DataLens lens;
    Jama.Matrix X;
    Jama.Matrix Y;
    boolean[] homogeneityIndex;
    Jama.Matrix homogeneityParameters;
    boolean allParametersHomogeneous;
    boolean failedEstimation = false;
    MomentSpecification spec;
    
    public ContainerLinear(DataLens lens, boolean[] homogeneityIndex, Jama.Matrix homogeneityParameters, boolean allParametersHomogeneous, MomentSpecification spec) {
        this.lens = lens;
        // computeBetaAndErrors();
        X = lens.getX();
        Y = lens.getY();
        this.homogeneityIndex = homogeneityIndex;
        this.homogeneityParameters = homogeneityParameters;
        this.allParametersHomogeneous = allParametersHomogeneous;
        this.spec = spec;
    }

    @Override
    public void computeBetaAndErrors() {
    	// System.out.println("In compute beta and errors");
        if (Y.getRowDimension() < Math.max(30, X.getColumnDimension())) {
            //System.out.println("Too few observations, dim(Y)  = "+Y.getRowDimension() );
            beta = null;
            goodnessOfFit = Double.POSITIVE_INFINITY;
        } else {
            try {
                int numParams = spec.getNumParams();
                if (!allParametersHomogeneous) {
                    for (boolean b : homogeneityIndex) {
                        if (b) {
                            numParams = numParams - 1; // homogeneous parameter imposed externally
                        }
                    }
                }

                Uncmin_f77 minimizer = new Uncmin_f77(false);
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
                int[] ndigit = {0, 12};
                int[] itnlim = {0, 1500};
                int[] iagflg = {0, 0};
                int[] iahflg = {0, 0};
                double[] dlt = {0, 1};
                double[] gradtl = {0, 1E-15};
                double[] stepmx = {0, 1E8};
                double[] steptl = {0, 1E-15};

                if (numParams == spec.getNumParams() && 1==1) {
                    // System.out.println("no homogeneous parameters, just use OLS");
                    // Jama.Matrix olsBeta = pmUtility.OLSsvd(X, Y, false);
                    Jama.Matrix olsBeta = pmUtility.OLS(X, Y, false);
                    
                    // System.out.println("Average of Y is :"+pmUtility.mean(Y, 0));
                    for (int i = 0; i < olsBeta.getRowDimension(); i++) {
                        xpls[i + 1] = olsBeta.get(i, 0);
                    }
                    
                    // System.out.println("Y, X");
//                    pmUtility.prettyPrint(pmUtility.concatMatrix(Y, X));
//                    
//                    pmUtility.prettyPrintVector(olsBeta);
//                    
//                    System.exit(0);
                } else {
                    boolean tryResidualizing = false;
                    if (tryResidualizing) {
                        Jama.Matrix Yres = Y.copy();
                        Jama.Matrix Xres = null;
                        boolean first = true;
                        for (int k = 0; k < X.getColumnDimension(); k++) {
                            if (homogeneityIndex[k]) {
                                for (int i = 0; i < Y.getRowDimension(); i++) {
                                    Yres.set(i, 0, Yres.get(i, 0) - X.get(i, k) * homogeneityParameters.get(k, 0));
                                }
                            } else {
                                if (first) {
                                    first = false;
                                    Xres = pmUtility.getColumn(X, k);
                                } else {
                                    Xres = pmUtility.concatMatrix(Xres, pmUtility.getColumn(X, k));
                                }
                            }
                        }

                        Jama.Matrix olsBeta = pmUtility.OLSsvd(Xres, Yres, false);

                        // System.out.println("Average of Y is :"+pmUtility.mean(Y, 0));
                        for (int i = 0; i < olsBeta.getRowDimension(); i++) {
                            xpls[i + 1] = olsBeta.get(i, 0);
                        }
                    } else {
                    	//System.out.println("ContainerLinear running optimization");
                        minimizer.optif9_f77(numParams, guess, this, typsiz, fscale, method, iexp, msg, ndigit, itnlim, iagflg, iahflg, dlt, gradtl, stepmx, steptl, xpls, fpls, gpls, itrmcd, a, udiag);
                    }
                }
                
                boolean testFailure = false;
                if (testFailure) {
                    // put in something here about a failed estimation
                    if (itrmcd[1] == 4 || itrmcd[1] == 5) {
                        failedEstimation = true;
                        System.out.println("failure: uncmin failed to converge");
                    }

                    // check that optimizer didn't shoot off into extremes
                    for (int i = 0; i < xpls.length; i++) {
                        if (xpls[i] < -10 || xpls[i] > 10) {
                            failedEstimation = true;
                            System.out.println("failure: parameter magnitude too large");
                        }
                    }

                    // objective should be close to zero in these exactly-identified cases
                    if (f_to_minimize(xpls) > 10) {
                        System.out.println("failure: fmin too large");
                        failedEstimation = true;
                    }
                }
                
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
                
                beta = betaUncmin.copy();
                double sse = 0;
                Jama.Matrix fit = X.times(beta);
                for (int i = 0; i < Y.getRowDimension(); i++) {
                    sse += Math.pow(Y.get(i, 0) - fit.get(i, 0), 2);
                }
                goodnessOfFit = sse;
                
                //pmUtility.prettyPrint(beta);
                //System.out.println(goodnessOfFit);
                
                if (failedEstimation) {
                    beta = null;
                    goodnessOfFit = Double.POSITIVE_INFINITY;
                }
                
                if (debugVerbose) {
                    System.out.format("ContainerLinear.computeBetaAndErrors SSE: %g ", +goodnessOfFit);
                    pmUtility.prettyPrintVector(beta);

                    // ok i want to check that things are working here just for a single split of the data
                    // compute the variance using the OLS formula, compute the variance using the GMM formula
                    // make sure they are the same thing, and if not, figure out why not
                    Jama.Matrix olsVariance = getVariance(beta);
                    System.out.println("OLS variance:");
                    pmUtility.prettyPrint(olsVariance);

                    // compute variance formula, B = G'omega-1G, then B-1, compare to above
                    // get G first
                    Jama.Matrix G = getJacobianNoDivision(beta);
                    G.timesEquals(1.0 / Y.getRowDimension());
                    Jama.Matrix omega = new Jama.Matrix(X.getColumnDimension(), X.getColumnDimension());
                    Jama.Matrix fits = X.times(beta);
                    Jama.Matrix errors = fits.minus(Y);
                    System.out.println("sigma2: " + pmUtility.sumSquaredElements(errors));
                    for (int i = 0; i < Y.getRowDimension(); i++) {
                        Jama.Matrix gi = getGi(beta, i);
                        omega.plusEquals(gi.times(gi.transpose()));
                    }
                    omega.timesEquals(1.0 / Y.getRowDimension());
                    Jama.Matrix B = (G.transpose()).times(omega.inverse()).times(G);

                    // that line above why the numbers are not identical? we multiple X'X by X'Xinv by X'X. Maybe that cancellation is numerically imperfect?
                    System.out.println("G:");
                    pmUtility.prettyPrint(G);

                    System.out.println("G'*omega.inverse:");
                    pmUtility.prettyPrint(G.transpose().times(omega.inverse()));

                    System.out.println("G'*omega.inverse*G:");
                    pmUtility.prettyPrint(G.transpose().times(omega.inverse()).times(G));

                    System.out.println("omega:");
                    pmUtility.prettyPrint(omega);

                    // i think that inversion, multiplication combination is why they are not identical
                    // as the sample size grows, the covariance matrices get SUPER close (to the 6th decimal place) to each other
                    System.out.println("B:");
                    pmUtility.prettyPrint(B);
                    System.out.println("B inverse:");
                    pmUtility.prettyPrint(B.inverse());
                    System.out.println("Divided by n:");
                    pmUtility.prettyPrint(B.inverse().times(1.0 / Y.getRowDimension()));

                    System.exit(0);
                } 
            } catch (Exception e) {
                if (allParametersHomogeneous) {
                    e.printStackTrace();
                }
                if (debugVerbose) {
                    e.printStackTrace();
                    System.out.println("Something went wrong...Matrix not invertible?");
                    System.out.println("X:");           
                    pmUtility.prettyPrint(X);
                    System.out.println("Y:");           
                    pmUtility.prettyPrint(Y);                                      
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
        // Jama.Matrix gi = new Jama.Matrix(X.getColumnDimension(), 1);
        for (int k = 0; k < X.getColumnDimension(); k++) {
            g.set(k, 0, g.get(k, 0) + error_i * X.get(i, k));
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
        double fit = (X.getMatrix(i, i, 0, X.getColumnDimension() - 1).times(beta)).get(0, 0);
        double error = fit - Y.get(i, 0);
        Jama.Matrix gi = new Jama.Matrix(X.getColumnDimension(), 1);
        for (int k = 0; k < X.getColumnDimension(); k++) {
            gi.set(k, 0, error * X.get(i, k));
        }

        return gi;
    }

    @Override
    public Matrix getJacobianNoDivision(Jama.Matrix beta) {
        // in OLS case, super simple
        // just X'X
        return (X.transpose()).times(X);

        // in general it is this
        /**
         * int numMoments = X.getColumnDimension(); int numParameters =
         * beta.getRowDimension(); Jama.Matrix Jacobian = new
         * Jama.Matrix(numMoments, numParameters);
         *
         * for (int obs = 0; obs < X.getRowDimension(); obs++) { for (int i = 0;
         * i < numMoments; i++) { for (int j = 0; j < numParameters; j++) {
         * Jacobian.set(i, j, Jacobian.get(i, j) - X.get(obs, i) * X.get(obs,
         * j)); } } }
         *
         * pmUtility.prettyPrint(Jacobian);
         * pmUtility.prettyPrint((X.transpose()).times(X)); System.exit(0);
         *
         * return Jacobian;
         */
    }

    @Override
    public Jama.Matrix getMomentGWithoutDivision(Jama.Matrix beta) {
        /**
         * Let's implement the moment-based version of OLS here (need this for a
         * variety of reasons, also will extend nicely to other models more
         * directly this way)
         */
        // Jama.Matrix runningTotal = new Jama.Matrix(X.getRowDimension(), X.getColumnDimension());
        int numMoments = spec.getNumMoments();
        Jama.Matrix g = new Jama.Matrix(numMoments, 1); // x'e, one row for each x
        Jama.Matrix fittedY = X.times(beta);
        Jama.Matrix e = fittedY.minus(Y);

        // pmUtility.prettyPrint(pmUtility.concatMatrix(Y, pmUtility.concatMatrix(X, e)), 20);
        
        // turns out using gi method here is crazy slow!
        // because i was recalculating vectors for individual observations each time; totally unnecessary
        for (int i = 0; i < X.getRowDimension(); i++) {
            // Jama.Matrix gi = getGi(e.get(i, 0), i);
            // g.plusEquals(gi);
            addGi(g, e.get(i, 0), i);
        }
        // cannot have this here since we divide by different n in different places!
        // g.timesEquals(1.0 / Y.getRowDimension());
        return g;
    }

    @Override
    public double computeMeasureOfFit(Jama.Matrix beta) {
        Jama.Matrix fittedY = X.times(beta);
        Jama.Matrix e = fittedY.minus(Y);

        return ((e.transpose()).times(e)).get(0, 0);
    }

    private double getMomentObjectiveFunction(Jama.Matrix beta, boolean debugMoment) {
        /**
         * Let's implement the moment-based version of OLS here (need this for a
         * variety of reasons, also will extend nicely to other models more
         * directly this way)
         */
        // Jama.Matrix runningTotal = new Jama.Matrix(X.getRowDimension(), X.getColumnDimension());
        int numMoments = spec.getNumMoments();
        Jama.Matrix g = new Jama.Matrix(numMoments, 1); // x'e, one row for each x
        Jama.Matrix fittedY = X.times(beta);
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
        
        if (debugMoment) {
            System.out.println("beta:");
            pmUtility.prettyPrint(beta);
            System.out.println("e:");
            pmUtility.prettyPrintVector(e);

            XYSeries xy = new XYSeries("X1");

            for (int k = 0; k < X.getColumnDimension(); k++) {
                System.out.println("x" + k + ":");
                pmUtility.prettyPrintVector(pmUtility.getColumn(X, k));
                // System.out.println("running total" + k + ":");
                // pmUtility.prettyPrintVector(pmUtility.getColumn(runningTotal, k));
            }

            for (int i = 0; i < X.getRowDimension(); i++) {
                xy.add(X.get(i, 1), e.get(i, 0));
            }

            System.out.println("g:");
            pmUtility.prettyPrintVector(g);
            // System.out.println("omega inverse:");
            // pmUtility.prettyPrint(omega.inverse());
            System.out.println("q: " + q + " SSE: " + pmUtility.sumSquaredElements(e));

            // let's graph the fit to see what is going on here geometrically
            XYSeriesCollection xyc = new XYSeriesCollection();

            JFreeChart chart = ChartFactory.createScatterPlot("Fit", "X", "Y", xyc);
            JFrame f = new JFrame("Fit");
            f.setBounds(100, 100, 1000, 1000);
            // f.setVisible(true);
            f.getContentPane().setLayout(new BorderLayout());
            f.getContentPane().add(new ChartPanel(chart));
        }

        // key point here is that we estimate the model using GMM
        // but we report goodness-of-fit when searching over the splits
//        System.exit(0);
        // System.out.print(q + " ");
        // pmUtility.prettyPrintVector(beta);
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
    public Matrix getVariance(Jama.Matrix b) {
        // System.out.print("Computing variance for OLS with beta = ");
        // pmUtility.prettyPrintVector(b);
        double sse = 0;
        Jama.Matrix fit = X.times(b);
        for (int i = 0; i < Y.getRowDimension(); i++) {
            sse += Math.pow(Y.get(i, 0) - fit.get(i, 0), 2);
        }
        System.out.println("sse in OLS: " + sse);
        Jama.Matrix xpx = (X.transpose()).times(X);
        Jama.Matrix xpxInv = xpx.inverse();
        // double sigma2 = sse / (Y.getRowDimension() - X.getColumnDimension());
        double sigma2 = sse / (Y.getRowDimension());
        System.out.println("sigma2 in getVariance: " + sigma2 + " n = " + Y.getRowDimension());

        boolean debugHere = false;
        if (debugHere) {
            System.out.println("X'X:");
            pmUtility.prettyPrint(xpx);
            System.out.println("X'X times sigma2");
            pmUtility.prettyPrint(xpx.times(sigma2));
        }

        variance = xpxInv.times(sigma2);
        return variance;
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
            if (homogeneityIndex[i]) {
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

    public boolean didEstimatorFail() {
        return failedEstimation;
    }

}
