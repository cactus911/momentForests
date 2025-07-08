/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package experimental;

import core.DataLens;
import optimization.Uncmin_f77;
import optimization.Uncmin_methods;

/**
 *
 * @author Stephen P. Ryan <stephen.p.ryan@wustl.edu>
 */
public class LinearMoment implements Uncmin_methods {

    Jama.Matrix X;
    Jama.Matrix Y;
    Jama.Matrix Z;
    boolean useCUE;

    public LinearMoment(DataLens lens, boolean useCUE) {
        this.useCUE = useCUE;
        this.X = lens.getX();
        this.Y = lens.getY();
        this.Z = lens.getZ();
    }

    public double[] solveBeta() {
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

//        Jama.Matrix betaUncmin = new Jama.Matrix(X.getColumnDimension(), 1);
//        for (int i = 0; i < guess.length - 1; i++) {
//            betaUncmin.set(i, 0, xpls[i + 1]);
//        }
        return xpls;
    }

    @Override
    public double f_to_minimize(double[] x) {
        Jama.Matrix g = new Jama.Matrix(X.getColumnDimension(), 1);
        Jama.Matrix omega = new Jama.Matrix(X.getColumnDimension(), X.getColumnDimension());
        Jama.Matrix beta = new Jama.Matrix(x.length - 1, 1);
        for (int i = 0; i < x.length - 1; i++) {
            beta.set(i, 0, x[i + 1]);
        }
        Jama.Matrix fit = X.times(beta);
        Jama.Matrix e = fit.minus(Y);
        for (int i = 0; i < X.getRowDimension(); i++) {
            Jama.Matrix gi = new Jama.Matrix(X.getColumnDimension(), 1);
            for (int k = 0; k < X.getColumnDimension(); k++) {
                gi.set(k, 0, X.get(i, k) * e.get(i, 0));
                g.set(k, 0, g.get(k, 0) + gi.get(k, 0));
            }
            omega.plusEquals(gi.times(gi.transpose()));
        }
        g.timesEquals(1.0 / X.getRowDimension());
        omega.timesEquals(1.0 / X.getRowDimension());

        // adjusting by number of observations seems important here (in the classical GMM framework, doesn't matter since N is constant)
//        double q = X.getRowDimension() * 2.0 * ((g.transpose()).times(g)).get(0, 0); // go with identity to start
//        if (useCUE) {
//            q = X.getRowDimension() * 2.0 * ((g.transpose()).times(omega.inverse()).times(g)).get(0, 0); // go with identity to start
//        }
        
        double q = ((g.transpose()).times(g)).get(0, 0); // go with identity to start
        if (useCUE) {
            q = ((g.transpose()).times(omega.inverse()).times(g)).get(0, 0); // go with identity to start
        }
        return q;
    }

    @Override
    public void gradient(double[] x, double[] g) {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    public void hessian(double[] x, double[][] h) {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }
}
