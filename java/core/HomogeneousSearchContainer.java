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
public class HomogeneousSearchContainer implements Uncmin_methods, mcmc.mcmcFunction  {

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
    }

    public void executeSearch() {
        // System.out.println("Inside executeSearch");

        boolean useGridSearch = false;

        // is this blowing up because i'm calling this even when there are zero parameters to search over?!?
        if (useGridSearch && numParams == 1) {
            // System.out.println("Starting grid search");
//            XYSeries xy = new XYSeries("Parameter 1");
//            XYSeriesCollection xyc = new XYSeriesCollection(xy);
//            boolean showGUI = false;
//            if (showGUI) {
//                JFrame fr = new JFrame("Grid Search");
//                fr.setBounds(200, 200, 500, 500);
//                fr.getContentPane().setLayout(new BorderLayout());
//                JFreeChart chart = ChartFactory.createXYLineChart("Grid Search Homogeneous Parameter", "Theta", "MSE", xyc);
//                fr.getContentPane().add(new ChartPanel(chart));
//                fr.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
//                fr.setVisible(true);
//            }

            boolean first = true;
            double bestF = 0;
            double bestX = 0;
            for (double x = -2.0; x <= -0.0; x += 0.1) {
                double[] xp = new double[2];
                xp[1] = x;
                // System.out.println("Calling f_to_minimize");
                double f = f_to_minimize(xp);
                System.out.println("x = " + x + " f_to_minimize: " + f);
//                xy.add(x, f);
                if (f < bestF || first) {
                    bestF = f;
                    bestX = x;
                    first = false;
                }
            }
            System.out.println("Optimal x = " + bestX);

            setCompactEstimatedHomogeneousParameters(new Jama.Matrix(1, 1, bestX));

        } else {
            Uncmin_f77 minimizer = new Uncmin_f77(true);

            double[] guess = new double[numParams + 1];
            
            // cheat this for speeding up testing purposes
            double[] truth = {-1,1};
            for(int k=0;k<homogeneousParameterIndex.size();k++) {
                guess[k] = truth[homogeneousParameterIndex.get(k)];
            }

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
            
            minimizer.optif9_f77(numParams, guess, this, typsiz, fscale, method, iexp, msg, ndigit, itnlim, iagflg, iahflg, dlt, gradtl, stepmx, steptl, xpls, fpls, gpls, itrmcd, a, udiag);

            Jama.Matrix compactHomogeneousParameterVector = new Jama.Matrix(numParams, 1);
            for (int i = 0; i < numParams; i++) {
                compactHomogeneousParameterVector.set(i, 0, xpls[i + 1]);
            }
            setCompactEstimatedHomogeneousParameters(compactHomogeneousParameterVector); // need to account for the fact that this is potentially a shorter vector
        }
    }

    @Override
    public double f_to_minimize(double[] x) {
        Jama.Matrix homogeneousParameterVector = new Jama.Matrix(numParams, 1);
        for (int i = 0; i < numParams; i++) {
            homogeneousParameterVector.set(i, 0, x[i + 1]);
            mySpecification.setHomogeneousParameter(homogeneousParameterIndex.get(i), x[i + 1]);
        }
        return computeOutOfSampleMSE();
    }

    public double computeOutOfSampleMSE() {
        /**
         * Need to regenerate a new DataLens with the homogeneity restrictions
         * imposed; also, getY depends on guess of homogeneous parameters, so
         * need to regenerate for each guess of X here
         */
        // System.out.println("Initializing lens");

        // if whole model is homogeneous this doens't work since it pulls nulls and grows a null tree
        // just skip directly to evaluating fit using completely specified homogeneous model (COME BACK TO THIS)
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
        DataLens oosDataLens = mySpecification.getOutOfSampleXYZ(2000, rngSeedBaseOutOfSample); // this should eventually be modified to come out of the data itself (or generalized in some way)
        Jama.Matrix testZ = oosDataLens.getZ();
        Jama.Matrix residualizedX = mySpecification.residualizeX(oosDataLens.getX());
        // pmUtility.prettyPrint(pmUtility.concatMatrix(pmUtility.concatMatrix(oosDataLens.getY(), oosDataLens.getX()), testZ));

        /**
         * Compute out-of-sample fit at current homogeneous parameter vector
         */
        // System.out.println("Computing out of sample fits");
        double outOfSampleFit = 0;
        for (int i = 0; i < testZ.getRowDimension(); i++) {
            Jama.Matrix fullXi = oosDataLens.getX().getMatrix(i, i, 0, oosDataLens.getX().getColumnDimension() - 1);

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
                    double error = fitY - (oosDataLens.getY().get(i, 0));
                    System.out.println(" fitY: " + fitY + " Y: " + oosDataLens.getY().get(i, 0) + " sqErr: " + error * error);
                }
            }
            double error = fitY - (oosDataLens.getY().get(i, 0));

            outOfSampleFit += error * error;
        }
        double MSE = outOfSampleFit / testZ.getRowDimension();
        // System.out.println("Out-of-sample MSE: " + MSE);
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