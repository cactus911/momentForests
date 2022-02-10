/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package core;

import java.awt.BorderLayout;
import java.util.ArrayList;
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

    DataLens oosDataLens;
    Jama.Matrix testZ;
    Jama.Matrix residualizedX;
    Jama.Matrix testX;
    Jama.Matrix testY;
    long t1 = System.currentTimeMillis();

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

        oosDataLens = mySpecification.getOutOfSampleXYZ(2000, rngSeedBaseOutOfSample); // this should eventually be modified to come out of the data itself (or generalized in some way)
        testZ = oosDataLens.getZ();
        testX = oosDataLens.getX();
        residualizedX = mySpecification.residualizeX(testX);
        testY = oosDataLens.getY();
    }

    public void executeSearch() {
        // System.out.println("Inside executeSearch");

        Uncmin_f77 minimizer = new Uncmin_f77(true);

        double[] guess = new double[numParams + 1];

        // cheat this for speeding up testing purposes
//            double[] truth = {-1,1};
//            for(int k=0;k<homogeneousParameterIndex.size();k++) {
//                guess[k] = truth[homogeneousParameterIndex.get(k)];
//            }
        // use the starting values from the DM test
        for (int k = 0; k < homogeneousParameterIndex.size(); k++) {
            guess[k + 1] = mySpecification.getHomogeneousParameter(homogeneousParameterIndex.get(k));
        }
        System.out.print("Starting values taken from DM test: ");
        pmUtility.prettyPrint(new Jama.Matrix(guess, 1));

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
            // plot this function (is there something goofy going on here?)
            // there are kink points and nonconvexities in that function, that may be why things go off the rails very rarely
            // plotFunction(-2, 0, 50);
            // try using secant method
            // guess = secant(guess);

            // do a grid search (successively bracketing smaller intervals)
            // could make this adaptive so that if it picks the endpoints we move the endpoint and start over
            // good starting values should hopefully fix that issue
            // put a flag in here though if that's the case (maybe a system.exit) since that's no good
            // this appears to work
            // however, in more general cases this could be an issue with more parameters that we are searching over (in the mixed case)
            // cannot do a grid search in so many dimensions
            // the objective function looks like it has kinks and nonconvex parts, so that makes optimization a challenge (although min is nicely behaved)
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
            
            // see how well this works
            goldenRatioSearch(left, right);
            
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

            // two additional ideas: iterate on this until we get some stability? did the first
            // second: use the golden ratio approach
            

            // plotFunction(left, right, numEvals);
            // System.exit(0);
            xpls[1] = guess[1];

            // try one newton step after this?
            // doesn't look like it is necessary
            // itnlim[1] = 2;
            // minimizer.optif9_f77(numParams, guess, this, typsiz, fscale, method, iexp, msg, ndigit, itnlim, iagflg, iahflg, dlt, gradtl, stepmx, steptl, xpls, fpls, gpls, itrmcd, a, udiag);
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
        Jama.Matrix homogeneousParameterVector = new Jama.Matrix(numParams, 1);
        for (int i = 0; i < numParams; i++) {
            homogeneousParameterVector.set(i, 0, x[i + 1]);
            mySpecification.setHomogeneousParameter(homogeneousParameterIndex.get(i), x[i + 1]);
        }

        // i should seed this in some way with the estimates from the DistanceMetricTest results
        double v = computeOutOfSampleMSE();
//        long t2 = System.currentTimeMillis();
//        if(t2-t1>60000) { // if longer than a minute
//            System.out.print("f: "+v+" x: ");
//            pmUtility.prettyPrint(new Jama.Matrix(x,1));
//        }
        return v;
    }

    public double computeOutOfSampleMSE() {
        // long t1 = System.currentTimeMillis();
        /**
         * Need to regenerate a new DataLens with the homogeneity restrictions
         * imposed; also, getY depends on guess of homogeneous parameters, so
         * need to regenerate for each guess of X here
         */
        // System.out.println("Initializing lens");

        // if whole model is homogeneous this doens't work since it pulls nulls and grows a null tree
        // just skip directly to evaluating fit using completely specified homogeneous model (COME BACK TO THIS)
        // some mystery here as to why this isn't essentially instant when we have a fully saturated model (at least in the OLS case)
        MomentForest myForest = null;
        boolean allParametersHomogeneous = true;
        for (boolean b : mySpecification.getHomogeneousIndex()) {
            if (!b) {
                allParametersHomogeneous = false;
            }
        }
        if (allParametersHomogeneous) {
            // don't need to grow tree in this circumstance, obviously
        } else {
            DataLens homogenizedForestLens = new DataLens(mySpecification.getX(), mySpecification.getY(true), mySpecification.getZ(), null);

            boolean testParameterHomogeneity = false;
            // System.out.println("Initializing forest");
            myForest = new MomentForest(mySpecification, numberTreesInForest, rngSeedBaseMomentForest, homogenizedForestLens, verbose, new TreeOptions());
            TreeOptions cvOptions = new TreeOptions(0.01, minObservationsPerLeaf, minImprovement, maxTreeDepth, testParameterHomogeneity); // k = 1
            // System.out.println("Setting options");
            myForest.setTreeOptions(cvOptions);
            /**
             * Grow the moment forest
             */
            // System.out.println("Growing forest");
            myForest.growForest();
        }

        /**
         * Test vectors for assessment
         */
//        DataLens oosDataLens = mySpecification.getOutOfSampleXYZ(2000, rngSeedBaseOutOfSample); // this should eventually be modified to come out of the data itself (or generalized in some way)
//        Jama.Matrix testZ = oosDataLens.getZ();
//        Jama.Matrix residualizedX = mySpecification.residualizeX(oosDataLens.getX());
//        // pmUtility.prettyPrint(pmUtility.concatMatrix(pmUtility.concatMatrix(oosDataLens.getY(), oosDataLens.getX()), testZ));
        /**
         * Compute out-of-sample fit at current homogeneous parameter vector
         */
        // System.out.println("Computing out of sample fits");
        double outOfSampleFit = 0;
        for (int i = 0; i < testZ.getRowDimension(); i++) {
            Jama.Matrix fullXi = testX.getMatrix(i, i, 0, testX.getColumnDimension() - 1);

            double fitY = mySpecification.getHomogeneousComponent(fullXi);
            if (!allParametersHomogeneous) {
                Jama.Matrix zi = testZ.getMatrix(i, i, 0, testZ.getColumnDimension() - 1);
                Jama.Matrix b = myForest.getEstimatedParameterForest(zi);
                Jama.Matrix residualizedXi = residualizedX.getMatrix(i, i, 0, residualizedX.getColumnDimension() - 1);
                fitY += residualizedXi.times(b).get(0, 0);

                boolean outputFits = false;
                if (outputFits) {
                    System.out.print("z: " + pmUtility.stringPrettyPrint(zi) + " beta: " + pmUtility.stringPrettyPrintVector(b));
                    System.out.print(" x: " + pmUtility.stringPrettyPrint(fullXi));
                    System.out.print(" residualizedXb: " + residualizedXi.times(b).get(0, 0) + " hc: " + mySpecification.getHomogeneousComponent(fullXi));
                    double error = fitY - (testY.get(i, 0));
                    System.out.println(" fitY: " + fitY + " Y: " + testY.get(i, 0) + " sqErr: " + error * error);
                }
            }
            double error = fitY - (testY.get(i, 0));

            outOfSampleFit += error * error;

            // in principle, we could actually use something like GMM in here
        }
        double MSE = outOfSampleFit / testZ.getRowDimension();
        // System.out.println("Out-of-sample MSE: " + MSE);
        // long t2 = System.currentTimeMillis();
        // System.out.println("All parameter homogeneous: " + allParametersHomogeneous + " time: " + (t2 - t1));
        // is this slow because it keeps regenerating the data to fit each time? i bet it would be a LOT faster
        // to do that once and just pull that
        // OK something super weird happening here. The run time is increasing higher and higher as the program runs
        // something is leaking in here

        return MSE;
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
            System.out.println(xLower+" "+x1+" "+x2+" "+xUpper+" f1: "+f1+" f2: "+f2+" interval length: "+(xUpper-xLower));
            if (xUpper - xLower < tol) {
                go = false;
            }
        }
        double[] v = {xLower, xUpper}; // interval containing minimum
        return v;
    }

}
