/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package core;

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
    private final MomentSpecification spec;

    public DistanceMetricTestWholeTree(ArrayList<DataLens> v, MomentSpecification spec) {
        this.v = v;
        this.spec = spec;
    }

    public double computeStatistic(int indexConstrainedParameter) {

        System.out.println("------------------------ Testing parameter k = " + indexConstrainedParameter + " ------------------------");

        double dm = 0;

        // let's test this unconstrained first to make sure that we get the same parameters back
        int numParamsEachSplit = spec.getNumParams();

        double[] unconstrainedX = computeParameters(numParamsEachSplit * v.size());
        double fminUnconstrained = f_to_minimize(unconstrainedX);

        System.out.print("Unconstrained Estimates: ");
        pmUtility.prettyPrint(new Jama.Matrix(unconstrainedX, 1));
        for (Jama.Matrix beta : convertToBetaList(unconstrainedX)) {
            pmUtility.prettyPrintVector(beta);
        }
        System.out.println("SSE (unconstrained): " + computeSSE(unconstrainedX));

        // then impose constraint on indexParameterConstrain, report twice the difference
        this.indexConstrainedParameter = indexConstrainedParameter;
        double[] constrainedX = computeParameters((numParamsEachSplit - 1) * v.size() + 1);
        double fminConstrained = f_to_minimize(constrainedX);

        int numObs = 0;
        for (DataLens h : v) {
            numObs += h.getNumObs();
        }

        System.out.print("Constrained Estimates: ");
        pmUtility.prettyPrint(new Jama.Matrix(constrainedX, 1));

        valueConstrainedParameter = constrainedX[1];

        for (Jama.Matrix beta : convertToBetaList(constrainedX)) {
            pmUtility.prettyPrintVector(beta);
        }
        System.out.println("SSE (constrained): " + computeSSE(constrainedX));

        dm = 2.0 * numObs * (fminConstrained - fminUnconstrained);

        System.out.println("k = " + indexConstrainedParameter + " unconstrained: " + fminUnconstrained + " constrained: " + fminConstrained + " dm: " + dm);

        return dm;
    }

    private double[] computeParameters(int numParams) {
        System.out.println("Number of parameters: " + numParams);

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

        int xCounter = 0;
        boolean first = true;
        for (int leaf = 0; leaf < v.size(); leaf++) {
            DataLens lens = v.get(leaf);
            ContainerMoment cm = spec.getContainerMoment(lens);
            cm.computeBetaAndErrors();
            // pmUtility.prettyPrintVector(cm.getBeta());
            if (indexConstrainedParameter < 0) {
                // unconstrained case, just put in double[] x in order
                for (int i = 0; i < cm.getBeta().getRowDimension(); i++) {
                    guess[xCounter + 1] = cm.getBeta().get(i, 0);
                    xCounter++;
                }
            } else {
                // have to do something more complex here, figure that once i get the above working
                for (int i = 0; i < cm.getBeta().getRowDimension(); i++) {
                    if (indexConstrainedParameter == i) {
                        if (first) {
                            guess[1] = cm.getBeta().get(i, 0);
                            first = false;
                        }
                    } else {
                        guess[xCounter + 2] = cm.getBeta().get(i, 0);
                        xCounter++;
                    }
                }
            }
        }
        System.out.print("Seeded starting values: ");
        pmUtility.prettyPrint(new Jama.Matrix(guess, 1));
        for (int i = 0; i < xpls.length; i++) {
            xpls[i] = guess[i];
        }

        // with identity weighting matrix (Step 1)
        int numMoments = v.get(0).getX().getColumnDimension() * v.size();
        omega = Jama.Matrix.identity(numMoments, numMoments);

        boolean minimizeSSEFirst = false;
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

        boolean useUncminFirst = false;
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

        boolean useCUEInSecondStep = false;
        if (useCUEInSecondStep) {
            System.out.println("Recomputing with CUE...");
            useCUE = true;
            minimizer = new Uncmin_f77(true);
            minimizer.optif9_f77(numParams, guess, this, typsiz, fscale, method, iexp, msg, ndigit, itnlim, iagflg, iahflg, dlt, gradtl, stepmx, steptl, xpls, fpls, gpls, itrmcd, a, udiag);
            System.out.print("after second step with CUE: ");
            pmUtility.prettyPrint(new Jama.Matrix(xpls, 1));
        }

        boolean useOptimalTwoStepWeightingMatrix = true;
        if (useOptimalTwoStepWeightingMatrix) {
            for (int i = 0; i < 2; i++) {
                // with optimal weighting matrix (Step 2)
                System.out.print("Computing optimal weighting matrix...");
                omega = computeOptimalOmega(xpls);
                System.out.println("done.");

                System.out.println("With optimal weighting matrix:");
                pmUtility.prettyPrint(omega);

                minimizer = new Uncmin_f77(false);
                minimizer.optif9_f77(numParams, guess, this, typsiz, fscale, method, iexp, msg, ndigit, itnlim, iagflg, iahflg, dlt, gradtl, stepmx, steptl, xpls, fpls, gpls, itrmcd, a, udiag);
                System.out.print("after second step with fixed Omega: ");
                pmUtility.prettyPrint(new Jama.Matrix(xpls, 1));
                System.arraycopy(xpls, 0, guess, 0, guess.length);
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
        boolean oracle = false;
        if (oracle) {
            if (indexConstrainedParameter == 1) {
                x[1] = 1.0;
            }
        }

        ArrayList<Jama.Matrix> cellBetaList = convertToBetaList(x);

        /**
         * This is the number of moments in each cell (this is hardwired right
         * now for the linear case, need to come back to this and figure this
         * out in general when we have broader cases)
         */
        int K = spec.getNumParams();

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
            ContainerMoment c = spec.getContainerMoment(lens);
            Jama.Matrix leafG = c.getMomentGWithoutDivision(cellBetaList.get(leaf));
            for (int j = 0; j < leafG.getRowDimension(); j++) {
                g.set(leaf * K + j, 0, leafG.get(j, 0));
            }
        }
        g.timesEquals(1.0 / numObs);

        // omega may be messing this up
        // try just using the identity weighting matrix to see if it fixes things up
        // can turn off CUE using the boolean here
        if (useCUE) {
            omega = computeOptimalOmega(x);
        }

        double q = 0.5 * (((g.transpose()).times(omega.inverse())).times(g)).get(0, 0);

        if (useCUE) {
            System.out.println("-----------");
            System.out.print("f: " + q + " ");
            pmUtility.prettyPrint(new Jama.Matrix(x, 1));
            pmUtility.prettyPrint(omega);
        }

        return q;
    }

    public double computeSSE(double[] x) {
        ArrayList<Jama.Matrix> cellBetaList = convertToBetaList(x);

        /**
         * This is the number of moments in each cell (this is hardwired right
         * now for the linear case, need to come back to this and figure this
         * out in general when we have broader cases)
         */
        if(v.get(0).getX().getColumnDimension()!=spec.getNumParams()) {
            System.out.println("WARNING! Check your Container that you implemented SSE correctly!\nDistanceMetricTestWholeTree.java:296 Use of specification getNumParams gives different answer than previous hardwired columnDimension!");
            // System.exit(0);
        }
        int K = spec.getNumParams(); // v.get(0).getX().getColumnDimension();

        double SSE = 0;

        /**
         * Go leaf by leaf and stack moments
         */
        int numObs = 0;
        Jama.Matrix g = new Jama.Matrix(0, 1);
        for (int leaf = 0; leaf < v.size(); leaf++) {
            DataLens lens = v.get(leaf);
            numObs += lens.getNumObs();
            ContainerMoment c = spec.getContainerMoment(lens);
            g = pmUtility.stackMatrix(g, c.getMomentGWithoutDivision(cellBetaList.get(leaf)));
            SSE += c.computeMeasureOfFit(cellBetaList.get(leaf));
        }
        g.timesEquals(1.0 / numObs);
        System.out.println("Moment Vector:");
        pmUtility.prettyPrint(g);
        System.out.println("g'g: " + ((g.transpose()).times(g)).get(0, 0));

        return SSE;
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

    private ArrayList<Jama.Matrix> convertToBetaList(double[] x) {
        ArrayList<Jama.Matrix> cellBetaList = new ArrayList<>();
        int numParamsEachLeaf = spec.getNumParams();

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
            }
            cellBetaList.add(beta);
        }
        return cellBetaList;
    }

    private Jama.Matrix computeOptimalOmega(double[] x) {
        ArrayList<Jama.Matrix> cellBetaList = convertToBetaList(x);

        /**
         * This is the number of moments in each cell (this is hardwired right
         * now for the linear case, need to come back to this and figure this
         * out in general when we have broader cases)
         */
        int K = spec.getNumMoments();

        /**
         * The total size of the stacked moment vector is the dimensionality of
         * moments in each cell times number of cells
         */
        int numMoments = v.size() * K;
        Jama.Matrix W = new Jama.Matrix(numMoments, numMoments);

        /**
         * Go leaf by leaf and stack moments
         */
        int numObs = 0;
        for (int leaf = 0; leaf < v.size(); leaf++) {
            // need to get gi, which is a little weird here

            DataLens lens = v.get(leaf);
            numObs += lens.getNumObs();
            ContainerMoment c = spec.getContainerMoment(lens);
            for (int i = 0; i < lens.getNumObs(); i++) {
                Jama.Matrix gi = new Jama.Matrix(numMoments, 1);
                Jama.Matrix leafGi = c.getGi(cellBetaList.get(leaf), i);
                for (int j = 0; j < leafGi.getRowDimension(); j++) {
                    gi.set(leaf * K + j, 0, leafGi.get(j, 0));
                }
                W.plusEquals(gi.times(gi.transpose()));
            }

        }
        W.timesEquals(1.0 / numObs);

        return W;
    }

}
