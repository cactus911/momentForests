/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package core;

import examples.linear.ContainerLinear;
import java.util.ArrayList;
import mcmc.gibbsLTEGeneralized;
import optimization.Uncmin_f77;
import optimization.Uncmin_methods;
import utility.pmUtility;

/**
 *
 * @author Stephen P. Ryan <stephen.p.ryan@wustl.edu>
 */
public class DistanceMetricTestWholeTree implements Uncmin_methods, mcmc.mcmcFunction {

    ArrayList<DataLens> v;
    int indexConstrainedParameter = -1;
    private double valueConstrainedParameter;

    Jama.Matrix omega; // weighting matrix in GMM
    boolean useCUE = false; // utilize continuously-updated weighting matrix
    boolean useSumOfSquaredErrors = false; // use SSE to get a starting value for GMM, which can be sensitive in small samples

    public DistanceMetricTestWholeTree(ArrayList<DataLens> v) {
        this.v = v;
    }

    public double computeStatistic(int indexConstrainedParameter) {

        System.out.println("------------------------ Testing parameter k = " + indexConstrainedParameter + " ------------------------");

        double dm = 0;

        // let's test this unconstrained first to make sure that we get the same parameters back
        int numParamsEachSplit = v.get(0).getX().getColumnDimension();

        double[] unconstrainedX = computeParameters(numParamsEachSplit * v.size());
        double fminUnconstrained = f_to_minimize(unconstrainedX);

        System.out.println("Unconstrained Estimates");
        pmUtility.prettyPrintVector(new Jama.Matrix(unconstrainedX, 1));
        System.exit(0);
        
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
                valueConstrainedParameter = betaLeftC.get(i, 0);
            }
        }

        System.out.println("Constrained Estimates");
        System.out.print(" \tbeta left: ");
        pmUtility.prettyPrintVector(betaLeftC);
        System.out.print("\tbeta right: ");
        pmUtility.prettyPrintVector(betaRightC);

        System.out.println("Need to get right number of observations here");
        dm = -666; // 2.0 * (leftY.getRowDimension() + rightY.getRowDimension()) * (fminConstrained - fminUnconstrained);
        System.exit(0);
        System.out.println("k = " + indexConstrainedParameter + " unconstrained: " + fminUnconstrained + " constrained: " + fminConstrained + " dm: " + dm);

        return dm;
    }

    private double[] computeParameters(int numParams) {
        
        System.out.println("Number of parameters: "+numParams);
        
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
        double[] stepmx = {0, 1.0}; // default is 1E8 (!!!), declining this to help ensure it doesn't shoot off into outer space
        double[] steptl = {0, 1E-8};

        // with identity weighting matrix (Step 1)
        int numMoments = v.get(0).getX().getColumnDimension() * v.size();
        omega = Jama.Matrix.identity(numMoments, numMoments);

        boolean minimizeSSEFirst = true;
        if (minimizeSSEFirst) {
            useSumOfSquaredErrors = true;
            minimizer.optif9_f77(numParams, guess, this, typsiz, fscale, method, iexp, msg, ndigit, itnlim, iagflg, iahflg, dlt, gradtl, stepmx, steptl, xpls, fpls, gpls, itrmcd, a, udiag);
            for (int i = 0; i < guess.length; i++) {
                guess[i] = xpls[i];
            }
            System.out.print("After using SSE+Uncmin: ");
            pmUtility.prettyPrint(new Jama.Matrix(xpls, 1));
            useSumOfSquaredErrors = false;
        }

        boolean useUncminFirst = true;
        if (useUncminFirst) {
            minimizer.optif9_f77(numParams, guess, this, typsiz, fscale, method, iexp, msg, ndigit, itnlim, iagflg, iahflg, dlt, gradtl, stepmx, steptl, xpls, fpls, gpls, itrmcd, a, udiag);
            for (int i = 0; i < guess.length; i++) {
                guess[i] = xpls[i];
            }
            System.out.print("After using GMM+Uncmin: ");
            pmUtility.prettyPrint(new Jama.Matrix(xpls, 1));
        }

        boolean useLTE = false;
        if (useLTE) {
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
        }

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
        ArrayList<Jama.Matrix> cellBetaList = new ArrayList<>();
        int numParamsEachLeaf = v.get(0).getX().getColumnDimension();

        /**
         * Fill up a list with betas that will match the moments in each cell
         */
        int counter = 1;
        for (int leaf = 0; leaf < v.size(); leaf++) {
            Jama.Matrix beta = new Jama.Matrix(numParamsEachLeaf, 1);

            if (indexConstrainedParameter < 0) {
                // this is the unconstrained case
                for (int i = 0; i < numParamsEachLeaf; i++) {
                    // System.out.println(i+" "+leaf+" "+v.size());
                    beta.set(i, 0, x[(leaf * numParamsEachLeaf) + i + 1]);
                }
            } else {
                // with a constraint
                for (int i = 0; i < numParamsEachLeaf; i++) {
                    if (i == indexConstrainedParameter) {
                        beta.set(i, 0, x[1]); // we are going to use the convention that the first parameter is the constrained one

                    } else {
                        beta.set(i, 0, x[counter + 1]);
                        counter++;
                    }
                }
//            pmUtility.prettyPrint(new Jama.Matrix(x, 1));
//            pmUtility.prettyPrintVector(leftBeta);
//            pmUtility.prettyPrintVector(rightBeta);
            }
            cellBetaList.add(beta);
        }

        /**
         * This is the number of moments in each cell (this is hardwired right
         * now for the linear case, need to come back to this and figure this
         * out in general when we have broader cases)
         */
        int K = v.get(0).getX().getColumnDimension();

        /**
         * The total size of the stacked moment vector is the dimensionality of
         * moments in each cell times number of cells
         */
        int numMoments = v.size() * K;
        Jama.Matrix g = new Jama.Matrix(numMoments, 1); // x'e, one row for each x, multiply by two for each split (left/right)

        /**
         * Go leaf by leaf and stack moments
         */
        int numObs = 0;
        for (int leaf = 0; leaf < v.size(); leaf++) {
            DataLens lens = v.get(leaf);
            numObs += lens.getNumObs();
            ContainerLinear c = new ContainerLinear(lens);
            Jama.Matrix leafG = c.getMomentG(cellBetaList.get(leaf));
            for (int j = 0; j < leafG.getRowDimension(); j++) {
                g.set(leaf * K + j, 0, leafG.get(j, 0));
            }
        }
        g.timesEquals(1.0 / numObs);

        // omega may be messing this up
        // try just using the identity weighting matrix to see if it fixes things up
        // can turn off CUE using the boolean here
        if (useCUE) {
            computeOptimalOmega(x);
        }

        double q = 0.5 * (((g.transpose()).times(omega.inverse())).times(g)).get(0, 0);

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

    @Override
    public double objectiveFunction(double[] x) {
        return -f_to_minimize(x);
    }

    @Override
    public double pi(double[] x) {
        return 1.0;
    }

    /**
     * @return the valueConstrainedParameter
     */
    public double getValueConstrainedParameter() {
        return valueConstrainedParameter;
    }

    private void computeOptimalOmega(double[] x) {
        throw new UnsupportedOperationException("Not supported yet."); // Generated from nbfs://nbhost/SystemFileSystem/Templates/Classes/Code/GeneratedMethodBody
    }

}
