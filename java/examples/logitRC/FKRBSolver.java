/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/Classes/Class.java to edit this template
 */
package examples.logitRC;

import Jama.Matrix;
import com.mathworks.engine.MatlabEngine;
import java.io.StringWriter;
import mcmc.gibbsLTEGeneralized;
import optimization.Uncmin_f77;
import optimization.Uncmin_methods;
import utility.pmUtility;

/**
 *
 * @author stephen.p.ryan
 */
public class FKRBSolver implements Uncmin_methods, mcmc.mcmcFunction {

    Jama.Matrix Y;
    Jama.Matrix X;
    double penalty = 0;
    private double violationUpperBound;
    private double violationLowerBound;
    private double violationSummingConstraint;

    public FKRBSolver(Matrix Y, Matrix X) {
        this.Y = Y;
        this.X = X;
    }

    public Jama.Matrix solve() {
        boolean useMATLAB = true;
        if (useMATLAB) {
            /**
             * Convert this to use Matlab!
             */
            try {
                MatlabEngine eng;
                String[] engines = MatlabEngine.findMatlab();
                if (engines.length > 0) {
                    eng = MatlabEngine.connectMatlab(engines[0]);
                } else {
                    eng = MatlabEngine.startMatlab();
                }
                eng.putVariable("X", X.getArray());

                int numParams = X.getColumnDimension();

                // this is so cool!
                eng.putVariable("Y", Y.getArray());
                eng.eval("Aeq=ones(1," + numParams + ");");
                eng.eval("beq=1"); // parameters sum to one
                eng.eval("lb=zeros(" + numParams + ",1);");
                eng.eval("ub=ones(" + numParams + ",1);");

                // StringWriter writer = new StringWriter();
                // eng.eval("size(Aeq)", writer, null);
                // System.out.println(writer.toString());
                // System.out.println("X: "+X.getRowDimension()+" by "+X.getColumnDimension());
                // System.out.println("Y: "+X.getRowDimension()+" by "+Y.getColumnDimension());
                eng.eval("beta=lsqlin(X,Y,[],[],Aeq,beq,lb,ub)");

                double[] bhat = eng.getVariable("beta");
                pmUtility.prettyPrint(new Jama.Matrix(bhat, 1));

                eng.close();
                return new Jama.Matrix(bhat, 1).transpose();
            } catch (Exception e) {
                e.printStackTrace();
                System.exit(0);
            }
        }

        int numParams = X.getColumnDimension() - 1;

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
        int[] ndigit = {0, 8};
        int[] itnlim = {0, 150};
        int[] iagflg = {0, 0};
        int[] iahflg = {0, 0};
        double[] dlt = {0, 1};
        double[] gradtl = {0, 1E-8};
        double[] stepmx = {0, 1.0}; // size of maximum step (default is 1E8! Making this MUCH smaller to prevent this thing from blowing up into outer space)
        double[] steptl = {0, 1E-8};

        boolean useUncmin = false;
        if (useUncmin) {
            minimizer.optif9_f77(numParams, guess, this, typsiz, fscale, method, iexp, msg, ndigit, itnlim, iagflg, iahflg, dlt, gradtl, stepmx, steptl, xpls, fpls, gpls, itrmcd, a, udiag);
        }

        boolean useLTE = true;
        if (useLTE) {
            for (penalty = 1E7; penalty <= 1E7; penalty *= 10) {
                mcmc.gibbsLTEGeneralized lte = new gibbsLTEGeneralized(this, 200, 200, guess, false);
                xpls = lte.getLowestPoint();
                guess = xpls;
                System.out.print("penalty: " + penalty + " sum: " + violationSummingConstraint + " lower: " + violationLowerBound + " upper: " + violationUpperBound + " x: ");
                Jama.Matrix beta = new Jama.Matrix(xpls.length, 1);
                double sum = 0;
                for (int i = 1; i < beta.getRowDimension(); i++) {
                    beta.set(i, 0, xpls[i]);
                    sum += xpls[i];
                }
                beta.set(0, 0, 1.0 - sum); // impose summing up constrain
                pmUtility.prettyPrintVector(beta);
            }
        }

        f_to_minimize(xpls);

        double sum = 0;
        Jama.Matrix beta = new Jama.Matrix(xpls.length, 1);
        for (int i = 1; i < beta.getRowDimension(); i++) {
            beta.set(i, 0, xpls[i]);
            sum += xpls[i];
        }
        beta.set(0, 0, 1.0 - sum); // impose summing up constraint

        return beta;
    }

    @Override
    public double f_to_minimize(double[] x) {
        // try to impose summing up constraint here programmatically?
        // try that; here are bounds on variables
        for (int i = 1; i < x.length; i++) {
            if (x[i] > 1.0) {
                x[i] = 1.0;
            }
            if (x[i] < 0.0) {
                x[i] = 0.0;
            }
        }

        double sum = 0;
        Jama.Matrix beta = new Jama.Matrix(x.length, 1);
        for (int i = 1; i < beta.getRowDimension(); i++) {
            beta.set(i, 0, x[i]);
            sum += x[i];
        }
        beta.set(0, 0, 1.0 - sum); // impose summing up constraint
        // pmUtility.prettyPrintVector(beta);

        // so there is still a problem here when the first element gets set to something negative to make everything work
        // i'm going to leave this for later when i can actually use a proper optimizer for a constrained optimization
        Jama.Matrix fit = X.times(beta);
        Jama.Matrix error = fit.minus(Y);

        double sumWeights = 0;
        violationUpperBound = 0;
        violationLowerBound = 0;
        for (int i = 0; i < beta.getRowDimension(); i++) {
            double parameter = beta.get(i, 0);
            sumWeights += parameter;
            if (parameter > 1) {
                // System.out.println("Greater than one "+x[i+1]);
                violationUpperBound += Math.pow(parameter - 1.0, 2);
                // System.out.println("upper bound: "+violationUpperBound);
            }
            if (parameter < 0) {
                // System.out.println("Less than one "+x[i+1]);
                violationLowerBound += Math.pow(-parameter, 2);
                // System.out.println("lower bound: "+violationLowerBound);
            }
        }
        violationSummingConstraint = Math.pow(1.0 - sumWeights, 2);

        return error.norm2() + penalty * (violationSummingConstraint + violationLowerBound + violationUpperBound);
    }

    @Override
    public void gradient(double[] x, double[] g) {
        throw new UnsupportedOperationException("Not supported yet."); // Generated from nbfs://nbhost/SystemFileSystem/Templates/Classes/Code/GeneratedMethodBody
    }

    @Override
    public void hessian(double[] x, double[][] h) {
        throw new UnsupportedOperationException("Not supported yet."); // Generated from nbfs://nbhost/SystemFileSystem/Templates/Classes/Code/GeneratedMethodBody
    }

    @Override
    public double objectiveFunction(double[] x) {
        return -f_to_minimize(x);
    }

    @Override
    public double pi(double[] x) {
        // System.out.println("Inside Pi");
//        for (int i = 1; i < x.length; i++) {
//            if (x[i] > 1.0) {
//                return 0.0;
//            }
//            if (x[i] < 0.0) {
//                return 0.0;
//            }
//        }
//        double sum = 0;
//        for (int i = 1; i < x.length; i++) {
//            sum += x[i];
//        }
//        if (sum > 1.0) {
//            return 0.0;
//        }
        return 1.0;
    }

}
