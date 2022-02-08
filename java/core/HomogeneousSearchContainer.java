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

        Uncmin_f77 minimizer = new Uncmin_f77(false);

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

        // good for robustness, way too slow to actually use
//            long t1 = System.currentTimeMillis();
//            mcmc.gibbsLTEGeneralized lte = new gibbsLTEGeneralized(this, 10, 0, guess, false);
//            guess = lte.getLowestPoint();
//            long t2 = System.currentTimeMillis();
//            System.out.println("===================================== "+(t2-t1)+" ms to starting value: "+pmUtility.stringPrettyPrint(new Jama.Matrix(guess,1)));
        long t1 = System.currentTimeMillis();
        minimizer.optif9_f77(numParams, guess, this, typsiz, fscale, method, iexp, msg, ndigit, itnlim, iagflg, iahflg, dlt, gradtl, stepmx, steptl, xpls, fpls, gpls, itrmcd, a, udiag);
        long t2 = System.currentTimeMillis();

        boolean allParametersHomogeneous = true;
        for (boolean b : mySpecification.getHomogeneousIndex()) {
            if (!b) {
                allParametersHomogeneous = false;
            }
        }
        System.out.println("All parameters homogeneous: " + allParametersHomogeneous + " time elapsed: " + (t2 - t1));

        Jama.Matrix compactHomogeneousParameterVector = new Jama.Matrix(numParams, 1);
        for (int i = 0; i < numParams; i++) {
            compactHomogeneousParameterVector.set(i, 0, xpls[i + 1]);
        }
        setCompactEstimatedHomogeneousParameters(compactHomogeneousParameterVector); // need to account for the fact that this is potentially a shorter vector

    }

    @Override
    public double f_to_minimize(double[] x) {
        Jama.Matrix homogeneousParameterVector = new Jama.Matrix(numParams, 1);
        for (int i = 0; i < numParams; i++) {
            homogeneousParameterVector.set(i, 0, x[i + 1]);
            mySpecification.setHomogeneousParameter(homogeneousParameterIndex.get(i), x[i + 1]);
        }

        // i should seed this in some way with the estimates from the DistanceMetricTest results
        double v= computeOutOfSampleMSE();
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

}
