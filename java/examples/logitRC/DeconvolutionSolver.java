/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/Classes/Class.java to edit this template
 */
package examples.logitRC;

import core.MomentSpecification;
import optimization.Uncmin_f77;
import optimization.Uncmin_methods;
import utility.pmUtility;

/**
 *
 * @author stephen.p.ryan
 */
public class DeconvolutionSolver implements Uncmin_methods {

    int numModels;
    MomentSpecification mySpecification;
    Jama.Matrix testX;
    Jama.Matrix testY;

    private double[][] betaList;
    private double[] betaWeights;

    public DeconvolutionSolver(Jama.Matrix testY, Jama.Matrix testX, MomentSpecification mySpecification) {
        this.testX = testX;
        this.testY = testY;

        this.mySpecification = mySpecification;

        boolean simpleTypes = false;
        if (simpleTypes) {
            double[][] b2 = {{-1, -1},
            {-1, 0},
            {-1, 1},
            {-1, 2},
            {-1, 3}};
            betaList = b2;
            numModels = betaList.length;
        } else {
            numModels = 10;
            betaList = new double[numModels][2];
            double lowerBeta = -5;
            double upperBeta = 5;
            double increment = (upperBeta - lowerBeta) / (numModels - 1);
            for (int i = 0; i < numModels; i++) {
                betaList[i][1] = lowerBeta + increment * i;
            }
        }
        betaWeights = new double[betaList.length];
    }

    public double solve() {
        int numParams = 1;

        Uncmin_f77 minimizer = new Uncmin_f77(true);

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

        boolean go = true;
        if (numParams == 1 && go) {
            double left = -5.0;
            double right = 5.0;
            return goldenRatioSearch(left, right)[1];
        } else {

//         guess[1] = -1;
//        System.out.println("fmin at truth: "+f_to_minimize(guess));
//        System.exit(0);
            minimizer.optif9_f77(numParams, guess, this, typsiz, fscale, method, iexp, msg, ndigit, itnlim, iagflg, iahflg, dlt, gradtl, stepmx, steptl, xpls, fpls, gpls, itrmcd, a, udiag);
        }

        f_to_minimize(xpls);

        return xpls[1];
    }

    private double[] goldenRatioSearch(double xLower, double xUpper) {
        double tol = 1E-5;
        double gold = ((Math.sqrt(5.0) - 1) / 2.0);
        double d = gold * (xUpper - xLower);
        double x1 = xLower + d;
        double x2 = xUpper - d;
        double[] fx1 = {0, x1};
        double[] fx2 = {0, x2};
        double f1 = f_to_minimize(fx1);
        double f2 = f_to_minimize(fx2);
        boolean go = true;
        while (go) {
            if (f1 < f2) {
                xLower = x2;
                x2 = x1;
                f2 = f1;
                x1 = xLower + gold * (xUpper - xLower);
                fx1[1] = x1;
                f1 = f_to_minimize(fx1);
            } else {
                xUpper = x1;
                x1 = x2;
                f1 = f2;
                x2 = xUpper - gold * (xUpper - xLower);
                fx2[1] = x2;
                f2 = f_to_minimize(fx2);
            }
            // System.out.println(xLower + " " + x1 + " " + x2 + " " + xUpper + " f(x1): " + f1 + " f(x2): " + f2 + " interval length: " + (xUpper - xLower));
            if (xUpper - xLower < tol) {
                go = false;
            }
        }
        double[] v = {xLower, xUpper}; // interval containing minimum
        return v;
    }

    @Override
    public double f_to_minimize(double[] x) {
        Jama.Matrix FKRB_X = new Jama.Matrix(testY.getRowDimension(), numModels);
        for (int i = 0; i < testY.getRowDimension(); i++) {
            Jama.Matrix xi = testX.getMatrix(i, i, 0, testX.getColumnDimension() - 1);

            for (int j = 0; j < numModels; j++) {
                Jama.Matrix betaJ = new Jama.Matrix(2, 1);

                betaJ.set(0, 0, x[1]);
                betaJ.set(1, 0, betaList[j][1]);
                FKRB_X.set(i, j, mySpecification.getPredictedY(xi, betaJ, null));
            }
        }

        Jama.Matrix weights;
        boolean useSolver = !true;
        if (useSolver) {
            FKRBSolver fsol = new FKRBSolver(testY, FKRB_X);
            weights = fsol.solve();
        } else {
            weights = pmUtility.OLSsvd(FKRB_X, testY, false);
        }

        for (int i = 0; i < weights.getRowDimension(); i++) {
            betaWeights[i] = weights.get(i, 0);
        }

        // pmUtility.prettyPrintVector(weights);
        // System.out.println("Compared to unconstrained OLS:");
        // pmUtility.prettyPrintVector(weightsOLS);
        // System.exit(0);
        Jama.Matrix fit = FKRB_X.times(weights);
        Jama.Matrix error = fit.minus(testY);
        return error.norm2();
    }

    @Override
    public void gradient(double[] x, double[] g) {
        throw new UnsupportedOperationException("Not supported yet."); // Generated from nbfs://nbhost/SystemFileSystem/Templates/Classes/Code/GeneratedMethodBody
    }

    @Override
    public void hessian(double[] x, double[][] h) {
        throw new UnsupportedOperationException("Not supported yet."); // Generated from nbfs://nbhost/SystemFileSystem/Templates/Classes/Code/GeneratedMethodBody
    }

    /**
     * @return the betaList
     */
    public double[][] getBetaList() {
        return betaList;
    }

    /**
     * @return the betaWeights
     */
    public double[] getBetaWeights() {
        return betaWeights;
    }

}
