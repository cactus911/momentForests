/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package core;

import java.awt.BorderLayout;
import javax.swing.JFrame;
import optimization.Uncmin_f77;
import optimization.Uncmin_methods;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;
import utility.pmUtility;

/**
 *
 * @author Stephen P. Ryan <stephen.p.ryan@wustl.edu>
 */
public class HomogeneousSearchContainer implements Uncmin_methods {

    MomentSpecification mySpecification;
    int numberTreesInForest;
    boolean verbose;
    double minImprovement;
    int minObservationsPerLeaf;
    int maxTreeDepth;
    int numParams;
    long rngSeedBaseMomentForest;
    long rngSeedBaseOutOfSample;

    private Jama.Matrix estimatedHomogeneousParameters;

    public HomogeneousSearchContainer(MomentSpecification mySpecification, int numberTreesInForest, boolean verbose, double minImprovement, int minObservationsPerLeaf, int maxTreeDepth, int numParams, long rngSeedBaseMomentForest, long rngSeedBaseOutOfSample) {
        this.mySpecification = mySpecification;
        this.numberTreesInForest = numberTreesInForest;
        this.verbose = verbose;
        this.minImprovement = minImprovement;
        this.minObservationsPerLeaf = minObservationsPerLeaf;
        this.maxTreeDepth = maxTreeDepth;
        this.numParams = numParams;
        this.rngSeedBaseMomentForest = rngSeedBaseMomentForest;
        this.rngSeedBaseOutOfSample = rngSeedBaseOutOfSample;
    }

    public void executeSearch() {

        boolean useGridSearch = true;

        if (useGridSearch && numParams == 1) {
            XYSeries xy = new XYSeries("Parameter 1");
            XYSeriesCollection xyc = new XYSeriesCollection(xy);
            JFrame fr = new JFrame("Grid Search");
            fr.setBounds(200, 200, 500, 500);
            fr.getContentPane().setLayout(new BorderLayout());
            JFreeChart chart = ChartFactory.createXYLineChart("Grid Search Homogeneous Parameter", "Theta", "MSE", xyc);
            fr.getContentPane().add(new ChartPanel(chart));
            fr.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
            fr.setVisible(true);

            boolean first = true;
            double bestF = 0;
            double bestX = 0;
            for (double x = -2.0; x <= -0.0; x += 0.1) {
                double[] xp = new double[2];
                xp[1] = x;
                double f = f_to_minimize(xp);
                System.out.println("x = " + x + " f_to_minimize: " + f);
                xy.add(x, f);
                if (f < bestF || first) {
                    bestF = f;
                    bestX = x;
                    first = false;
                }
            }
            System.out.println("Optimal x = " + bestX);

            setEstimatedHomogeneousParameters(new Jama.Matrix(1, 1, bestX));

        } else {
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
            minimizer.optif9_f77(numParams, guess, this, typsiz, fscale, method, iexp, msg, ndigit, itnlim, iagflg, iahflg, dlt, gradtl, stepmx, steptl, xpls, fpls, gpls, itrmcd, a, udiag);

            Jama.Matrix homogeneousParameterVector = new Jama.Matrix(numParams, 1);
            for (int i = 0; i < numParams; i++) {
                homogeneousParameterVector.set(i, 0, xpls[i + 1]);
            }
            setEstimatedHomogeneousParameters(homogeneousParameterVector);
        }
    }

    @Override
    public double f_to_minimize(double[] x) {
        Jama.Matrix homogeneousParameterVector = new Jama.Matrix(numParams, 1);
        for (int i = 0; i < numParams; i++) {
            homogeneousParameterVector.set(i, 0, x[i + 1]);
        }
        mySpecification.setHomogeneousParameters(homogeneousParameterVector);

        /**
         * Need to regenerate a new datalens with the homogeneity restrictions
         * imposed; also, getY depends on guess of homogeneous parameters, so
         * need to regenerate for each guess of X here
         */
        DataLens homogenizedForestLens = new DataLens(mySpecification.getX(), mySpecification.getY(true), mySpecification.getZ(), null);

        MomentForest myForest = new MomentForest(mySpecification, numberTreesInForest, rngSeedBaseMomentForest, homogenizedForestLens, verbose, new TreeOptions());
        TreeOptions cvOptions = new TreeOptions(0.01, minObservationsPerLeaf, minImprovement, maxTreeDepth, false); // k = 1
        myForest.setTreeOptions(cvOptions);
        /**
         * Grow the moment forest
         */
        myForest.growForest();

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
        double outOfSampleFit = 0;
        for (int i = 0; i < testZ.getRowDimension(); i++) {
            Jama.Matrix zi = testZ.getMatrix(i, i, 0, testZ.getColumnDimension() - 1);
            Jama.Matrix b = myForest.getEstimatedParameterForest(zi);

            Jama.Matrix residualizedXi = residualizedX.getMatrix(i, i, 0, residualizedX.getColumnDimension() - 1);
            Jama.Matrix fullXi = oosDataLens.getX().getMatrix(i, i, 0, oosDataLens.getX().getColumnDimension() - 1);

            // note, xi is part of x that is not homogeneous
            // for the homogeneous component below, need to find the other parts of x not contained in xi
            // not only that, the way that that method is specified is that it takes the ENTIRE xi to get the homogeneous part, where I have pulled out the subset of X already in xi
            // so how to get the whole row?
            // i am going to do two things here: one, remove residualization from getOutOfSampleXYZ
            // added a new method to residualize an X matrix
            double fitY = residualizedXi.times(b).get(0, 0) + mySpecification.getHomogeneousComponent(fullXi);
            double error = fitY - (oosDataLens.getY().get(i, 0));

            outOfSampleFit += error * error;

            boolean outputFits = false;
            if (outputFits) {
                System.out.print("z: " + pmUtility.stringPrettyPrint(zi) + " beta: " + pmUtility.stringPrettyPrintVector(b));
                System.out.print(" x: " + pmUtility.stringPrettyPrint(fullXi));
                System.out.print(" residualizedXb: " + residualizedXi.times(b).get(0, 0) + " hc: " + mySpecification.getHomogeneousComponent(fullXi));
                System.out.println(" fitY: " + fitY + " Y: " + oosDataLens.getY().get(i, 0) + " SE: " + error * error);
            }
        }
        double MSE = outOfSampleFit / testZ.getRowDimension();
        System.out.println("Out-of-sample MSE: " + MSE);
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
    public void setEstimatedHomogeneousParameters(Jama.Matrix estimatedHomogeneousParameters) {
        this.estimatedHomogeneousParameters = estimatedHomogeneousParameters;
    }

}
