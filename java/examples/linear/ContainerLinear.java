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

    public ContainerLinear(DataLens lens) {
        this.lens = lens;
        // computeBetaAndErrors();
    }

    public void computeBetaAndErrors() {
        Jama.Matrix X = lens.getX();
        Jama.Matrix Y = lens.getY();

        if (Y.getRowDimension() < 30) {
            // System.out.println("Too few observations");
            beta = null;
            goodnessOfFit = Double.POSITIVE_INFINITY;
        } else {
            try {
                boolean useUncmin = false;
                if (useUncmin) {
                    Uncmin_f77 minimizer = new Uncmin_f77(false);

                    int numParams = X.getColumnDimension();

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
                    int[] method = {0, 3};
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
                    minimizer.optif9_f77(numParams, guess, this, typsiz, fscale, method, iexp, msg, ndigit, itnlim, iagflg, iahflg, dlt, gradtl, stepmx, steptl, xpls, fpls, gpls, itrmcd, a, udiag);

                    Jama.Matrix betaUncmin = new Jama.Matrix(X.getColumnDimension(), 1);
                    for (int i = 0; i < guess.length - 1; i++) {
                        betaUncmin.set(i, 0, xpls[i + 1]);
                    }

                    // produce the essentially same answer 
//            pmUtility.prettyPrintVector(betaOLS);
//            pmUtility.prettyPrintVector(betaUncmin);
//            System.exit(0);
                    beta = betaUncmin.copy();
                    double sse = 0;
                    Jama.Matrix fit = X.times(beta);
                    for (int i = 0; i < Y.getRowDimension(); i++) {
                        sse += Math.pow(Y.get(i, 0) - fit.get(i, 0), 2);
                    }
                    goodnessOfFit = sse;

//                    System.out.print("\t\tFound uncmin beta: ");
//                    pmUtility.prettyPrintVector(betaUncmin);
//                    Jama.Matrix betaOLS = pmUtility.OLSsvd(X, Y, false);
//                    System.out.print("\t\tCompared to betaOLS: ");
//                    pmUtility.prettyPrintVector(betaOLS);
                } else {
                    Jama.Matrix betaOLS = pmUtility.OLSsvd(X, Y, false);
                    beta = betaOLS.copy();
                    double sse = 0;
                    Jama.Matrix fit = X.times(beta);
                    for (int i = 0; i < Y.getRowDimension(); i++) {
                        sse += Math.pow(Y.get(i, 0) - fit.get(i, 0), 2);
                    }
                    goodnessOfFit = sse;
                    // objectiveFunctionValue = getMoment(beta);
                }

                if (debugVerbose) {
                    System.out.format("ContainerLinear.computeBetaAndErrors SSE: %g ", +goodnessOfFit);
                    pmUtility.prettyPrintVector(beta);
                }
            } catch (Exception e) {
                e.printStackTrace();
                if (debugVerbose) {
                    System.out.println("Matrix not invertible");
                }
                beta = null;
                goodnessOfFit = Double.POSITIVE_INFINITY;
            }
        }
    }

    public Jama.Matrix getGi(Jama.Matrix beta, int i) {
        Jama.Matrix X = lens.getX();
        Jama.Matrix Y = lens.getY();
        int numMoments = X.getColumnDimension();
        Jama.Matrix fittedY = X.times(beta);
        Jama.Matrix e = fittedY.minus(Y);

        Jama.Matrix gi = new Jama.Matrix(numMoments, 1);
        for (int k = 0; k < numMoments; k++) {
            gi.set(k, 0, e.get(i, 0) * X.get(i, k));
        }

        return gi;
    }

    private double getMoment(Jama.Matrix beta, boolean debugMoment) {
        /**
         * Let's implement the moment-based version of OLS here (need this for a
         * variety of reasons, also will extend nicely to other models more
         * directly this way)
         */
        Jama.Matrix X = lens.getX();
        Jama.Matrix Y = lens.getY();
        // Jama.Matrix runningTotal = new Jama.Matrix(X.getRowDimension(), X.getColumnDimension());
        int numMoments = X.getColumnDimension();
        Jama.Matrix g = new Jama.Matrix(numMoments, 1); // x'e, one row for each x
        Jama.Matrix fittedY = X.times(beta);
        Jama.Matrix e = fittedY.minus(Y);
        Jama.Matrix omega = new Jama.Matrix(numMoments, numMoments);

        // turns out using gi method here is crazy slow!        
        for (int i = 0; i < X.getRowDimension(); i++) {
            Jama.Matrix gi = new Jama.Matrix(numMoments, 1);
            for (int k = 0; k < numMoments; k++) {
                gi.set(k, 0, e.get(i, 0) * X.get(i, k));                
                // runningTotal.set(i, k, g.get(k, 0));
            }
//            pmUtility.prettyPrintVector(gi);
            g.plusEquals(gi);
            omega.plusEquals(gi.times(gi.transpose()));
        }
        omega.timesEquals(1.0 / Y.getRowDimension());
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
            System.out.println("omega inverse:");
            pmUtility.prettyPrint(omega.inverse());
            System.out.println("q: " + q + " norm2: " + e.norm2());

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
    public Matrix getVariance() {
        return variance;
    }

    @Override
    public double getMomentFunctionImposingHomogeneity(int k, double value
    ) {
        Jama.Matrix betaHomogeneous = beta.copy();
        betaHomogeneous.set(k, 0, value);
        return getMoment(betaHomogeneous, false);
    }

    @Override
    public double f_to_minimize(double[] x
    ) {
        Jama.Matrix b = new Jama.Matrix(x.length - 1, 1);
        for (int i = 0; i < x.length - 1; i++) {
            b.set(i, 0, x[i + 1]);
        }
        return getMoment(b, false);
    }

    @Override
    public void gradient(double[] x, double[] g
    ) {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    public void hessian(double[] x, double[][] h
    ) {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    public double getMomentFunctionValue(Jama.Matrix b
    ) {
        return getMoment(b, false);
    }

}
