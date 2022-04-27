/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package core;

import java.awt.BorderLayout;
import java.util.ArrayList;
import java.util.Random;
import javax.swing.JFrame;
import optimization.Uncmin_f77;
import optimization.Uncmin_methods;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;
import utility.pmUtility;

/**
 *
 * @author Stephen P. Ryan <stephen.p.ryan@wustl.edu>
 */
public class HomogeneousSearchContainer implements Uncmin_methods, mcmc.mcmcFunction {

    MomentSpecification mySpecification;
    int numberTreesInForest;
    boolean verbose;
    double minImprovement;
    int minObservationsPerLeaf;
    int maxTreeDepth;
    int numParams;
    long rngSeedBaseMomentForest;
    long rngSeedBaseOutOfSample;
    ArrayList<Integer> homogeneousParameterIndex;

    private Jama.Matrix estimatedHomogeneousParameters;
    boolean allParametersHomogeneous = true;

    long t1 = System.currentTimeMillis();
    DataLens homogenizedForestLens;

    public HomogeneousSearchContainer(MomentSpecification mySpecification, int numberTreesInForest, boolean verbose, double minImprovement, int minObservationsPerLeaf, int maxTreeDepth, ArrayList<Integer> homogeneousParameterIndex, long rngSeedBaseMomentForest, long rngSeedBaseOutOfSample) {
        this.mySpecification = mySpecification;
        this.numberTreesInForest = numberTreesInForest;
        this.verbose = verbose;
        this.minImprovement = minImprovement;
        this.minObservationsPerLeaf = minObservationsPerLeaf;
        this.maxTreeDepth = maxTreeDepth;
        this.numParams = homogeneousParameterIndex.size();
        this.rngSeedBaseMomentForest = rngSeedBaseMomentForest;
        this.rngSeedBaseOutOfSample = rngSeedBaseOutOfSample;
        this.homogeneousParameterIndex = homogeneousParameterIndex;

        homogenizedForestLens = new DataLens(mySpecification.getX(), mySpecification.getY(), mySpecification.getZ(), null);

        for (boolean b : mySpecification.getHomogeneousIndex()) {
            if (!b) {
                allParametersHomogeneous = false;
            }
        }
    }

    public void executeSearch() {
        // System.out.println("Inside executeSearch");

        Uncmin_f77 minimizer = new Uncmin_f77(true);

        System.out.println("Number of parameters: " + numParams);
        // System.exit(0);

        double[] guess = new double[numParams + 1];

        // use the starting values from the test
        for (int k = 0; k < homogeneousParameterIndex.size(); k++) {
            guess[k + 1] = mySpecification.getHomogeneousParameter(homogeneousParameterIndex.get(k));
        }
        System.out.print("Starting values taken from test: ");
        pmUtility.prettyPrint(new Jama.Matrix(guess, 1));

        System.out.println("F_min(x): " + f_to_minimize(guess));

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

        if (numParams == 1) {
            // just use the number passed from the legit optimizer that already did all the work in detecting the parameter
            xpls[1] = guess[1];
        }

        if (numParams == 1) {
            if (1 == 2) {
                t1 = System.currentTimeMillis();
                int numEvals = 9; // number evaluations within each bracket (MIN: 3)
                int R = 5; // number of times to bracket grid search
                double left = guess[1] - 1.0;
                double right = guess[1] + 1.0;
                double increment = (right - left) / numEvals;
                // for (int r = 0; r < R; r++) {
                boolean converged = false;
                int r = 0;
                boolean first = true;

                // plotFunction(-4, 0, 25);
                // see how well this works
                guess[1] = goldenRatioSearch(left, right)[1];
                // this appears to work incredibly well

                boolean useGridSearch = false;

                if (useGridSearch) {
                    while (!converged) {
                        double[] v = gridSearch(left, right, increment);
                        System.out.println("r = " + r + " Grid search produced best x = " + v[0] + " f(x) = " + v[1] + " [" + (2.0 * increment) + "]");
                        if (Math.abs(guess[1] - v[0]) < 1E-4) { // this criterion really just checks the length of the search interval, which is probably fine for now
                            if (first) {
                                first = false;
                            } else {
                                converged = true;
                            }
                        }
                        guess[1] = v[0];
                        left = guess[1] - increment;
                        right = guess[1] + increment;
                        increment = (right - left) / numEvals;
                        r++;
                    }
                    long t2 = System.currentTimeMillis();
                    System.out.println("Grid search with numEvals: " + numEvals + " R: " + R + " time: " + (t2 - t1) + " f: " + f_to_minimize(guess));
                }

                // plotFunction(left, right, numEvals);
                // System.exit(0);
                xpls[1] = guess[1];

                // try one newton step after this?
                // doesn't look like it is necessary
                // itnlim[1] = 2;
                // minimizer.optif9_f77(numParams, guess, this, typsiz, fscale, method, iexp, msg, ndigit, itnlim, iagflg, iahflg, dlt, gradtl, stepmx, steptl, xpls, fpls, gpls, itrmcd, a, udiag);
            }
        } else {

            // good for robustness, way too slow to actually use
//            long t1 = System.currentTimeMillis();
//            mcmc.gibbsLTEGeneralized lte = new gibbsLTEGeneralized(this, 10, 0, guess, false);
//            guess = lte.getLowestPoint();
//            long t2 = System.currentTimeMillis();
//            System.out.println("===================================== "+(t2-t1)+" ms to starting value: "+pmUtility.stringPrettyPrint(new Jama.Matrix(guess,1)));
// it is getting stuck internally (not even calling f_to_minimize!)
            minimizer.optif9_f77(numParams, guess, this, typsiz, fscale, method, iexp, msg, ndigit, itnlim, iagflg, iahflg, dlt, gradtl, stepmx, steptl, xpls, fpls, gpls, itrmcd, a, udiag);

            // System.out.print("Post-uncmin f: " + fpls[1] + " x = ");
            // pmUtility.prettyPrint(new Jama.Matrix(xpls, 1));
        }

        Jama.Matrix compactHomogeneousParameterVector = new Jama.Matrix(numParams, 1);
        for (int i = 0; i < numParams; i++) {
            compactHomogeneousParameterVector.set(i, 0, xpls[i + 1]);
        }
        setCompactEstimatedHomogeneousParameters(compactHomogeneousParameterVector); // need to account for the fact that this is potentially a shorter vector

    }

    private double[] gridSearch(double leftEnd, double rightEnd, double interval) {
        double bestF = 0;
        double bestX = leftEnd;

        boolean first = true;
        double h = 1E-5;
        // look in a bracket around starting value that came from DM test
        for (double xp = leftEnd; xp <= rightEnd + h; xp += interval) {
            double[] tfx = {0, xp};
            double f = f_to_minimize(tfx);
            if (first) {
                first = false;
                bestX = xp;
                bestF = f;
            }
            if (f < bestF) {
                bestF = f;
                bestX = xp;
            }
            // System.out.println("xp: " + xp + " f(xp): " + f);
        }
        double[] v = {bestX, bestF};
        return v;
    }

    @Override
    public double f_to_minimize(double[] x) {
        for (int i = 0; i < numParams; i++) {
            mySpecification.setHomogeneousParameter(homogeneousParameterIndex.get(i), x[i + 1]);
        }
        // System.out.print("Inside homogeneousSearchContainer f_to_min: ");
        // pmUtility.prettyPrint(new Jama.Matrix(x, 1));

        // i should seed this in some way with the estimates from the DistanceMetricTest results
        double v = computeObjectiveFunction();
//        long t2 = System.currentTimeMillis();
//        if(t2-t1>60000) { // if longer than a minute
//            System.out.print("f: "+v+" x: ");
//            pmUtility.prettyPrint(new Jama.Matrix(x,1));
//        }
        return v;
    }

    public double computeObjectiveFunction() {
        // verbose = true;

        /**
         * Two things: 1. need to be using in-sample Y,X,Z here, not
         * out-of-sample 2. do not need to grow any trees when all the
         * parameters are homogeneous
         */
        if (allParametersHomogeneous) {
            Jama.Matrix beta = new Jama.Matrix(mySpecification.getHomogeneousIndex().length, 1);
            for (int i = 0; i < mySpecification.getHomogeneousIndex().length; i++) {
                beta.set(i, 0, mySpecification.getHomogeneousParameter(i));
            }
            double f = 0;
            Jama.Matrix Y = mySpecification.getY();
            Jama.Matrix X = mySpecification.getX();
            for (int i = 0; i < Y.getRowDimension(); i++) {
                f += mySpecification.getGoodnessOfFit(Y.get(i, 0), X.getMatrix(i, i, 0, X.getColumnDimension() - 1), beta);
            }
            return f;
        }

        boolean testParameterHomogeneity = false;
        // System.out.println("Initializing forest");
        MomentForest myForest = new MomentForest(mySpecification, numberTreesInForest, rngSeedBaseMomentForest, homogenizedForestLens, verbose, new TreeOptions());
        TreeOptions cvOptions = new TreeOptions(0.01, minObservationsPerLeaf, minImprovement, maxTreeDepth, testParameterHomogeneity); // k = 1
        // System.out.println("Setting options");
        myForest.setTreeOptions(cvOptions);
        /**
         * Grow the moment forest
         */
        // System.out.println("Growing forest");
        myForest.growForest();

        // System.out.println("Computing out of sample fits");
        // Random rng = new Random(888);
        double outOfSampleFit = 0;
        for (int i = 0; i < mySpecification.getZ().getRowDimension(); i++) {
            Jama.Matrix fullXi = mySpecification.getX().getMatrix(i, i, 0, mySpecification.getX().getColumnDimension() - 1);
            Jama.Matrix zi = mySpecification.getZ().getMatrix(i, i, 0, mySpecification.getZ().getColumnDimension() - 1);

            // old code from the linear case (doesn't work in general)
//            double fitY = mySpecification.getPredictedY(fullXi, myForest.getEstimatedParameterForest(zi), rng);
//            boolean outputFits = false;
//            if (outputFits) {
//                System.out.print("z: " + pmUtility.stringPrettyPrint(zi) + " beta: " + pmUtility.stringPrettyPrintVector(myForest.getEstimatedParameterForest(zi)));
//                System.out.print(" x: " + pmUtility.stringPrettyPrint(fullXi));
//                double error = fitY - (testY.get(i, 0));
//                System.out.println(" fitY: " + fitY + " Y: " + testY.get(i, 0) + " sqErr: " + error * error);
//            }
//            double error = fitY - (testY.get(i, 0));
//
//            outOfSampleFit += error * error;
            // in principle, we could actually use something like GMM in here
            // we are going to!
            outOfSampleFit += mySpecification.getGoodnessOfFit(mySpecification.getY().get(i, 0), fullXi, myForest.getEstimatedParameterForest(zi));
        }
        return outOfSampleFit / mySpecification.getZ().getRowDimension();
    }

    @Override
    public void gradient(double[] x, double[] g) {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    public void hessian(double[] x, double[][] h) {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    /**
     * @return the estimatedHomogeneousParameters
     */
    public Jama.Matrix getEstimatedHomogeneousParameters() {
        return estimatedHomogeneousParameters;
    }

    /**
     * @param estimatedHomogeneousParameters the estimatedHomogeneousParameters
     * to set
     */
    public void setCompactEstimatedHomogeneousParameters(Jama.Matrix estimatedHomogeneousParameters) {
        this.estimatedHomogeneousParameters = estimatedHomogeneousParameters;
    }

    @Override
    public double objectiveFunction(double[] x) {
        return -f_to_minimize(x);
    }

    @Override
    public double pi(double[] x) {
        return 1.0;
    }

    private void plotFunction(double left, double right, int numEvals) {
        JFrame f = new JFrame("Function Plot");
        f.setBounds(100, 100, 500, 500);
        f.getContentPane().setLayout(new BorderLayout());
        XYSeries xy = new XYSeries("f");
        ChartPanel panel = new ChartPanel(ChartFactory.createXYLineChart("F", "x", "F", new XYSeriesCollection(xy)));
        f.getContentPane().add(panel);
        f.setVisible(true);

        for (double x = left; x <= right; x += (right - left) / numEvals) {
            double[] fx = {0, x};
            xy.add(x, f_to_minimize(fx));
        }
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

}
