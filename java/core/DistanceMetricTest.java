/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package core;

import Jama.Matrix;
import mcmc.gibbsLTEGeneralized;
import optimization.Uncmin_f77;
import optimization.Uncmin_methods;
import utility.pmUtility;

/**
 *
 * @author Stephen P. Ryan <stephen.p.ryan@wustl.edu>
 */
public class DistanceMetricTest implements Uncmin_methods, mcmc.mcmcFunction {

    Jama.Matrix leftX;
    Jama.Matrix rightX;
    Jama.Matrix leftY;
    Jama.Matrix rightY;
    int indexConstrainedParameter = -1;

    Jama.Matrix omega; // weighting matrix in GMM
    boolean useCUE = false; // utilize continuously-updated weighting matrix

    public DistanceMetricTest(Matrix leftX, Matrix rightX, Matrix leftY, Matrix rightY) {
        this.leftX = leftX;
        this.rightX = rightX;
        this.leftY = leftY;
        this.rightY = rightY;
    }

    public double computeStatistic(int indexConstrainedParameter) {

        System.out.println("------------------------ Testing parameter k = " + indexConstrainedParameter + " ------------------------");

        double dm = 0;

        // let's test this unconstrained first to make sure that we get the same parameters back
        int numParamsEachSplit = rightX.getColumnDimension();

        double[] unconstrainedX = computeParameters(numParamsEachSplit * 2);
        double fminUnconstrained = f_to_minimize(unconstrainedX);

        Jama.Matrix betaLeft = new Jama.Matrix(numParamsEachSplit, 1);
        Jama.Matrix betaRight = new Jama.Matrix(numParamsEachSplit, 1);
        for (int i = 0; i < numParamsEachSplit; i++) {
            betaLeft.set(i, 0, unconstrainedX[i + 1]);
            betaRight.set(i, 0, unconstrainedX[i + 1 + numParamsEachSplit]);
        }

        System.out.println("Unconstrained Estimates");
        System.out.print(" \tbeta left: ");
        pmUtility.prettyPrintVector(betaLeft);
        System.out.print("\tbeta right: ");
        pmUtility.prettyPrintVector(betaRight); // this works perfectly
        // then impose constraint on indexParameterConstrain, report twice the difference
        this.indexConstrainedParameter = indexConstrainedParameter;
        double[] constrainedX = computeParameters(numParamsEachSplit * 2 - 1);
        double fminConstrained = f_to_minimize(constrainedX);

        Jama.Matrix betaLeftC = new Jama.Matrix(numParamsEachSplit, 1);
        Jama.Matrix betaRightC = new Jama.Matrix(numParamsEachSplit, 1);
        for (int i = 0; i < numParamsEachSplit; i++) {
            betaLeftC.set(i, 0, constrainedX[i + 1]);
        }
        int counter = 0;
        for (int i = 0; i < numParamsEachSplit; i++) {
            if (i != indexConstrainedParameter) {
                betaRightC.set(i, 0, constrainedX[counter + 1 + numParamsEachSplit]);
                counter++;
            } else {
                betaRightC.set(i, 0, betaLeftC.get(i, 0));
            }
        }

        System.out.println("Constrained Estimates");
        System.out.print(" \tbeta left: ");
        pmUtility.prettyPrintVector(betaLeftC);
        System.out.print("\tbeta right: ");
        pmUtility.prettyPrintVector(betaRightC);

        dm = 2.0 * (leftY.getRowDimension() + rightY.getRowDimension()) * (fminConstrained - fminUnconstrained);
        System.out.println("k = " + indexConstrainedParameter + " unconstrained: " + fminUnconstrained + " constrained: " + fminConstrained + " dm: " + dm);

        return dm;
    }

    private double[] computeParameters(int numParams) {
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
        double[] stepmx = {0, 1E8};
        double[] steptl = {0, 1E-8};

        // with identity weighting matrix (Step 1)
        int numMoments = leftX.getColumnDimension() * 2;
        omega = Jama.Matrix.identity(numMoments, numMoments);
        minimizer.optif9_f77(numParams, guess, this, typsiz, fscale, method, iexp, msg, ndigit, itnlim, iagflg, iahflg, dlt, gradtl, stepmx, steptl, xpls, fpls, gpls, itrmcd, a, udiag);
        for (int i = 0; i < guess.length; i++) {
            guess[i] = xpls[i];
        }
        System.out.print("after first step: ");
        pmUtility.prettyPrint(new Jama.Matrix(xpls, 1));

        boolean useCUEInSecondStep = true;
        if (useCUEInSecondStep) {
            useCUE = true;
            minimizer = new Uncmin_f77(false);
            minimizer.optif9_f77(numParams, guess, this, typsiz, fscale, method, iexp, msg, ndigit, itnlim, iagflg, iahflg, dlt, gradtl, stepmx, steptl, xpls, fpls, gpls, itrmcd, a, udiag);
            System.out.print("after second step with CUE: ");
            pmUtility.prettyPrint(new Jama.Matrix(xpls, 1));
        }

        boolean useOptimalTwoStepWeightingMatrix = false;
        if (useOptimalTwoStepWeightingMatrix) {
            // with optimal weighting matrix (Step 2)
            computeOptimalOmega(xpls);
            minimizer = new Uncmin_f77(false);
            minimizer.optif9_f77(numParams, guess, this, typsiz, fscale, method, iexp, msg, ndigit, itnlim, iagflg, iahflg, dlt, gradtl, stepmx, steptl, xpls, fpls, gpls, itrmcd, a, udiag);
            System.out.print("after second step with fixed Omega: ");
            pmUtility.prettyPrint(new Jama.Matrix(xpls, 1));
        }

        double start = f_to_minimize(xpls);
        mcmc.gibbsLTEGeneralized lte = new gibbsLTEGeneralized(this, 1000, 1000, xpls, false);
        guess = lte.getLowestPoint();
        double end = f_to_minimize(guess);
        System.out.print("After LTE: ");
        pmUtility.prettyPrint(new Jama.Matrix(guess, 1));
        if (end < start) {
            // hmm, this is sometimes finding a lower point, which is really important here since we relying on these f_min values to compute test statistics
            System.out.println("LTE made an improvement start: " + start + " end: " + end);
            // System.exit(0);
            for (int i = 0; i < guess.length; i++) {
                xpls[i] = guess[i];
            }
        }

        return xpls;
    }

    /**
     * This is for the testing case where we need to impose a common parameter
     * Return the objective function for stacked moments (with CUE weighting)
     *
     * @param x
     * @return
     */
    @Override
    public double f_to_minimize(double[] x) {
        int numParamsEachSplit = leftX.getColumnDimension();
        Jama.Matrix leftBeta = new Jama.Matrix(numParamsEachSplit, 1);
        Jama.Matrix rightBeta = new Jama.Matrix(numParamsEachSplit, 1);

        if (indexConstrainedParameter < 0) {
            // this is the unconstrained case
            for (int i = 0; i < numParamsEachSplit; i++) {
                leftBeta.set(i, 0, x[i + 1]);
                rightBeta.set(i, 0, x[i + 1 + numParamsEachSplit]);
            }
        } else {
            // with a constraint
            for (int i = 0; i < numParamsEachSplit; i++) {
                leftBeta.set(i, 0, x[i + 1]);
            }
            int counter = 0;
            for (int i = 0; i < numParamsEachSplit; i++) {
                if (i != indexConstrainedParameter) {
                    rightBeta.set(i, 0, x[counter + 1 + numParamsEachSplit]);
                    counter++;
                } else {
                    rightBeta.set(i, 0, leftBeta.get(i, 0));
                }
            }
//            pmUtility.prettyPrint(new Jama.Matrix(x, 1));
//            pmUtility.prettyPrintVector(leftBeta);
//            pmUtility.prettyPrintVector(rightBeta);
        }

        int numMoments = 2 * leftX.getColumnDimension();
        Jama.Matrix g = new Jama.Matrix(numMoments, 1); // x'e, one row for each x, multiply by two for each split (left/right)
        Jama.Matrix leftFittedY = leftX.times(leftBeta);
        Jama.Matrix leftError = leftFittedY.minus(leftY);

        Jama.Matrix rightFittedY = rightX.times(rightBeta);
        Jama.Matrix rightError = rightFittedY.minus(rightY);

        int K = leftX.getColumnDimension();

        /**
         * Left split
         */
        for (int i = 0; i < leftX.getRowDimension(); i++) {
            Jama.Matrix gi = new Jama.Matrix(numMoments, 1);
            for (int k = 0; k < K; k++) {
                gi.set(k, 0, leftError.get(i, 0) * leftX.get(i, k));
            }
            g.plusEquals(gi);
        }

        /**
         * Right split
         */
        for (int i = 0; i < rightX.getRowDimension(); i++) {
            Jama.Matrix gi = new Jama.Matrix(numMoments, 1);
            for (int k = 0; k < K; k++) {
                gi.set(k + K, 0, rightError.get(i, 0) * rightX.get(i, k));
            }
            g.plusEquals(gi);
        }

        g.timesEquals(1.0 / (leftY.getRowDimension() + rightY.getRowDimension()));

        // omega may be messing this up
        // try just using the identity weighting matrix to see if it fixes things up
        // can turn off CUE using the boolean here
        if (useCUE) {
            computeOptimalOmega(x);
        }

        double q = 0.5 * (((g.transpose()).times(omega.inverse())).times(g)).get(0, 0);

        return q;
    }

    private void computeOptimalOmega(double[] x) {
        int numParamsEachSplit = leftX.getColumnDimension();
        Jama.Matrix leftBeta = new Jama.Matrix(numParamsEachSplit, 1);
        Jama.Matrix rightBeta = new Jama.Matrix(numParamsEachSplit, 1);

        if (indexConstrainedParameter < 0) {
            // this is the unconstrained case
            for (int i = 0; i < numParamsEachSplit; i++) {
                leftBeta.set(i, 0, x[i + 1]);
                rightBeta.set(i, 0, x[i + 1 + numParamsEachSplit]);
            }
        } else {
            // with a constraint
            for (int i = 0; i < numParamsEachSplit; i++) {
                leftBeta.set(i, 0, x[i + 1]);
            }
            int counter = 0;
            for (int i = 0; i < numParamsEachSplit; i++) {
                if (i != indexConstrainedParameter) {
                    rightBeta.set(i, 0, x[counter + 1 + numParamsEachSplit]);
                    counter++;
                } else {
                    rightBeta.set(i, 0, leftBeta.get(i, 0));
                }
            }
//            pmUtility.prettyPrint(new Jama.Matrix(x, 1));
//            pmUtility.prettyPrintVector(leftBeta);
//            pmUtility.prettyPrintVector(rightBeta);
        }

        int numMoments = 2 * leftX.getColumnDimension();

        Jama.Matrix leftFittedY = leftX.times(leftBeta);
        Jama.Matrix leftError = leftFittedY.minus(leftY);

        Jama.Matrix rightFittedY = rightX.times(rightBeta);
        Jama.Matrix rightError = rightFittedY.minus(rightY);

        omega = new Jama.Matrix(numMoments, numMoments);

        int K = leftX.getColumnDimension();

        /**
         * Left split
         */
        for (int i = 0; i < leftX.getRowDimension(); i++) {
            Jama.Matrix gi = new Jama.Matrix(numMoments, 1);
            for (int k = 0; k < K; k++) {
                gi.set(k, 0, leftError.get(i, 0) * leftX.get(i, k));
                // g.set(k, 0, g.get(k,0) + gi.get(k,0));
            }
            omega.plusEquals(gi.times(gi.transpose()));
        }

        /**
         * Right split
         */
        for (int i = 0; i < rightX.getRowDimension(); i++) {
            Jama.Matrix gi = new Jama.Matrix(numMoments, 1);
            for (int k = 0; k < K; k++) {
                gi.set(k + K, 0, rightError.get(i, 0) * rightX.get(i, k));
                // g.set(k, 0, g.get(k,0) + gi.get(k,0));
            }
            omega.plusEquals(gi.times(gi.transpose()));
        }

        omega.timesEquals(1.0 / (leftY.getRowDimension() + rightY.getRowDimension()));
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
    public double objectiveFunction(double[] x) {
        return -f_to_minimize(x);
    }

    @Override
    public double pi(double[] x) {
        return 1.0;
    }

}
