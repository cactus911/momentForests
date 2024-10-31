/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package core;

import java.util.ArrayList;
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
    int numHomogeneousParametersToSearchOver;
    long rngSeedBaseMomentForest;
    long rngSeedBaseOutOfSample;
    ArrayList<Integer> homogeneousParameterList;

    private Jama.Matrix estimatedHomogeneousParameters;
    boolean allParametersHomogeneous = true;

    boolean debug = false;

    long t1 = System.currentTimeMillis();
    DataLens homogenizedForestLens;

    public HomogeneousSearchContainer(MomentSpecification mySpecification, int numberTreesInForest, boolean verbose, double minImprovement, int minObservationsPerLeaf, int maxTreeDepth,
            ArrayList<Integer> homogeneousParameterList, long rngSeedBaseMomentForest, long rngSeedBaseOutOfSample) {
        this.mySpecification = mySpecification;
        this.numberTreesInForest = numberTreesInForest;
        this.verbose = verbose;
        this.minImprovement = minImprovement;
        this.minObservationsPerLeaf = minObservationsPerLeaf;
        this.maxTreeDepth = maxTreeDepth;
        numHomogeneousParametersToSearchOver = homogeneousParameterList.size();
        this.rngSeedBaseMomentForest = rngSeedBaseMomentForest;
        this.rngSeedBaseOutOfSample = rngSeedBaseOutOfSample;
        this.homogeneousParameterList = homogeneousParameterList;

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

        if (debug) {
            System.out.println("Number of parameters: " + numHomogeneousParametersToSearchOver);
        }
        // System.exit(0);

        double[] guess = new double[numHomogeneousParametersToSearchOver + 1];

        // use the starting values from the test
        for (int k = 0; k < homogeneousParameterList.size(); k++) {
            if (!Double.isNaN(mySpecification.getHomogeneousParameter(homogeneousParameterList.get(k)))) {
                guess[k + 1] = mySpecification.getHomogeneousParameter(homogeneousParameterList.get(k));
            }
        }
        if (debug) {
            System.out.print("Starting values taken from test: ");
            pmUtility.prettyPrint(new Jama.Matrix(guess, 1));

            System.out.println("F_min(x): " + f_to_minimize(guess));
        }

        double[] xpls = new double[numHomogeneousParametersToSearchOver + 1];
        double[] fpls = new double[2];
        double[] gpls = new double[numHomogeneousParametersToSearchOver + 1];
        int[] itrmcd = new int[2];
        double[][] a = new double[numHomogeneousParametersToSearchOver + 1][numHomogeneousParametersToSearchOver + 1];
        double[] udiag = new double[numHomogeneousParametersToSearchOver + 1];
        double[] typsiz = new double[numHomogeneousParametersToSearchOver + 1];
        for (int i = 1; i <= numHomogeneousParametersToSearchOver; i++) {
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

        if (numHomogeneousParametersToSearchOver == 1) {
            // just use the number passed from the legit optimizer that already did all the work in detecting the parameter
            xpls[1] = guess[1];
        } else {
            minimizer.optif9_f77(numHomogeneousParametersToSearchOver, guess, this, typsiz, fscale, method, iexp, msg, ndigit, itnlim, iagflg, iahflg, dlt, gradtl, stepmx, steptl, xpls, fpls, gpls, itrmcd, a, udiag);
        }

        if (itrmcd[1] == 4 || itrmcd[1] == 5) {
            System.out.println("***** CAUTION: Uncmin reported failing to find a minimum *****");
        }

        // System.out.print("<<<<<<<<< After search, parameter vector is: ");
        // pmUtility.prettyPrint(new Jama.Matrix(xpls,1));
        for (int i = 0; i < numHomogeneousParametersToSearchOver; i++) {
            mySpecification.setHomogeneousParameter(homogeneousParameterList.get(i), xpls[i + 1]);
        }

        Jama.Matrix compactHomogeneousParameterVector = new Jama.Matrix(numHomogeneousParametersToSearchOver, 1);
        for (int i = 0; i < numHomogeneousParametersToSearchOver; i++) {
            compactHomogeneousParameterVector.set(i, 0, xpls[i + 1]);
        }
        setCompactEstimatedHomogeneousParameters(compactHomogeneousParameterVector); // need to account for the fact that this is potentially a shorter vector

    }

    @Override
    public double f_to_minimize(double[] x) {
        for (int i = 0; i < numHomogeneousParametersToSearchOver; i++) {
            mySpecification.setHomogeneousParameter(homogeneousParameterList.get(i), x[i + 1]);
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

        /**
         * Use the moment functions stacked across all trees for the objective
         * function here This is generalizable and leverages the moments
         * (necessary in cases outside of OLS)
         */
        boolean useMoments = true;
        if (useMoments) {
            // average across all the trees

            // we grew a tree under the parameter restriction, should be able to find the sum of (unaveraged??? unweighted?) GMM functions across leaves
            double avgGMMObjectiveFunctionValue = 0;
            for (TreeMoment tree : myForest.forest) {
                avgGMMObjectiveFunctionValue += tree.getTreeMomentObjectiveFunctionAtComputedParameters(verbose);
            }
            avgGMMObjectiveFunctionValue /= myForest.getForestSize();
            if (verbose) {

            }
            return avgGMMObjectiveFunctionValue;
        } else {
            System.out.println("Deprecated and marked for deletion Oct 22 2024");
            System.exit(0);
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
