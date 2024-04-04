/*
 * The MIT License
 *
 * Copyright 2020 Stephen P. Ryan <stephen.p.ryan@wustl.edu>.
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */
package examples.logitRC;

// import JSci.maths.statistics.NormalDistribution;
import JSci.maths.statistics.NormalDistribution;
import Jama.Matrix;
import core.ChartGenerator;
import core.ContainerMoment;
import core.DataLens;
import core.IntegerPartition;
import core.MomentContinuousSplitObj;
import core.MomentForest;
import core.MomentPartitionObj;
import core.MomentSpecification;
import core.OutOfSampleStatisticsContainer;
import core.TreeOptions;
import java.io.FileReader;
import java.io.BufferedReader;
import java.util.ArrayList;
import java.util.Random;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;
import utility.pmUtility;

/**
 *
 * @author Stephen P. Ryan <stephen.p.ryan@wustl.edu>
 */
public class LogitRCMomentSpecification implements MomentSpecification {

    Jama.Matrix X;
    Jama.Matrix Y;
    Jama.Matrix Z;
    Jama.Matrix balancingVector; // is treatment status in the RCT setting
    int numObs;
    int numtrees;
    int[] variableSearchIndex; // this should be restricted to only Z
    Boolean[] DiscreteVariables; // also this should be restricted to only Z
    String filename;
    boolean MONTE_CARLO = true;
    NormalDistribution normal = new NormalDistribution();

    /**
     * We are going to control homogeneous parameters through these variables
     */
    private boolean[] homogeneityIndex; // = new boolean[X.getColumnDimension()];
    private Jama.Matrix homogeneousParameterVector; // how to keep everything straight--easier way is to just to make it the size of the parameter vector and go from there (and never read elements that aren't labeled as homogeneous)

    public LogitRCMomentSpecification(int numObs) {
        this.numObs = numObs;
        // all these indices are hard-coded; want to change that down the road!
        homogeneityIndex = new boolean[2];
        homogeneousParameterVector = new Jama.Matrix(2, 1);
        resetHomogeneityIndex();

        /**
         * Here we are going to search over X (I think this works like this!)
         */
        // int[] vsi = {1}; //Search over x2 only (x1 is a constant; have to think about what that means)
        int[] vsi = {0, 1}; // Going to change this to make X1 also a continuous variable to search over

        Boolean[] wvd = {false, false}; // x1, x2 continuous (this thing has to be filled out for all the Z's, apparently; the length of this array is used in momentTree!
        variableSearchIndex = vsi;
        DiscreteVariables = wvd;
    }

    public LogitRCMomentSpecification(String filename) {
        this.filename = filename;
    }

    @Override
    public boolean[] getHomogeneousIndex() {
        return homogeneityIndex;
    }

    @Override
    public double getHomogeneousParameter(int parameterIndex) {
        return homogeneousParameterVector.get(parameterIndex, 0);
    }

    @Override
    public void setHomogeneousParameter(int parameterIndex, double value) {
        homogeneousParameterVector.set(parameterIndex, 0, value);
    }

    @Override
    public Matrix getBalancingVector() {
        return balancingVector;
    }

    public LogitRCMomentSpecification(Jama.Matrix X, Jama.Matrix Y, Jama.Matrix Z, int numtrees, int[] variableSearchIndex, Boolean[] DiscreteVariables) {
        this.X = X;
        this.Y = Y;
        this.Z = Z;
        this.numtrees = numtrees;
        this.variableSearchIndex = variableSearchIndex;
        this.DiscreteVariables = DiscreteVariables;
        balancingVector = pmUtility.getColumn(X, 0); // treatment status is the first column of X
    }

    @Override
    public Double getPredictedY(Matrix xi, Jama.Matrix beta, Random rng) {
        return LogitRCDataGenerator.getLogitShare(xi, beta);
    }

    @Override
    public int[] getVariableIndicesToSearchOver() {
        return variableSearchIndex;
    }

    @Override
    public MomentContinuousSplitObj getFminObjective(DataLens lens, int indexSplitVariable, double minProportionEachPartition, int minCountEachPartition) {
        return new MomentContinuousSplitObjLogitRC(indexSplitVariable, lens, minProportionEachPartition, minCountEachPartition, this);
    }

    @Override
    public MomentPartitionObj getMomentPartitionObj(DataLens lens, int indexSplitVariable, IntegerPartition partition) {
        return new MomentPartitionObjLogitRC(partition, indexSplitVariable, lens, this);
    }

    @Override
    public ContainerMoment computeOptimalBeta(DataLens lens, boolean allParametersHomogeneous) {
        ContainerLogitRC l = new ContainerLogitRC(lens, homogeneityIndex, homogeneousParameterVector, allParametersHomogeneous);
        l.computeBetaAndErrors();
        return l;
    }

    @Override
    public Matrix getY() {
        return Y;
    }

    @Override
    public Matrix getX() {
        return X;
    }

    @Override
    public Matrix getZ() {
        // do I simply return the X's here?!?!?
        return X;
        // return Z;
    }

    @Override
    public int numberoftrees() {
        return numtrees;
    }

    @Override
    public Boolean[] getDiscreteVector() {
        return DiscreteVariables;
    }

    @Override
    public double getGoodnessOfFit(double yi, Matrix xi, Matrix beta) {
        // here, let's return the LLH of this observation
        return ContainerLogitRC.computeSSEi(yi, xi, beta);
    }

    //Return the true parameter vector for a given observation
    @Override
    public Matrix getBetaTruth(Matrix zi, Random rng) {
        Jama.Matrix beta = new Jama.Matrix(2, 1); // Beta is a scalar
        beta.set(0, 0, -1.0 + 0.0 * normal.inverse(rng.nextDouble()));
        beta.set(1, 0, 1.0 + 1.0 * normal.inverse(rng.nextDouble()));

//        double draw = rng.nextDouble();
//        if (draw < 0.7) {
//            beta.set(1, 0, beta.get(1, 0) + normal.inverse(rng.nextDouble()));
//        } else {
//            beta.set(1, 0, -2.3 + normal.inverse(rng.nextDouble()));
//        }
        return beta;
    }

    @Override
    public OutOfSampleStatisticsContainer computeOutOfSampleStatistics(int numberTreesInForest, long rngBaseSeedMomentForest, boolean verbose,
            int minObservationsPerLeaf, double minImprovement, int maxTreeDepth, long rngBaseSeedOutOfSample) {

        double outOfSampleResultsY = 0;
        double outOfSampleResultsBeta = 0;

        MomentForest myForest;

        DataLens homogenizedForestLens = new DataLens(getX(), getY(), getZ(), null);

        myForest = new MomentForest(this, numberTreesInForest, rngBaseSeedMomentForest, homogenizedForestLens, verbose, new TreeOptions());
        TreeOptions cvOptions = new TreeOptions(0.01, minObservationsPerLeaf, minImprovement, maxTreeDepth, false); // k = 1
        myForest.setTreeOptions(cvOptions);
        /**
         * Grow the moment forest
         */
        myForest.growForest();

        // myForest.getTree(0).printTree();
        // debugOutputArea.append(myForest.getTree(0).toString());
        /**
         * Test vectors for assessment
         */
        DataLens oosDataLens = getOutOfSampleXYZ(25000, rngBaseSeedOutOfSample); // this should eventually be modified to come out of the data itself (or generalized in some way)
        Jama.Matrix testZ = oosDataLens.getZ();
        Jama.Matrix testX = oosDataLens.getX();
        Jama.Matrix testY = oosDataLens.getY();

        Random rng = new Random(rngBaseSeedOutOfSample - 3);

        boolean computeBetaAndYFits = false;
        if (computeBetaAndYFits) {
            for (int i = 0; i < testZ.getRowDimension(); i++) {
                Jama.Matrix zi = testZ.getMatrix(i, i, 0, testZ.getColumnDimension() - 1);
                Jama.Matrix xi = testX.getMatrix(i, i, 0, testX.getColumnDimension() - 1);

                // going to compare directly to the true parameter vector in this method instead of using fit of Y
                Jama.Matrix bTruth = getBetaTruth(zi, rng);

                Jama.Matrix compositeEstimatedBeta = myForest.getEstimatedParameterForest(zi);

                outOfSampleResultsY += getGoodnessOfFit(testY.get(i, 0), xi, compositeEstimatedBeta);

                if (i < 10) {
                    String hString = "[ ";
                    for (int k = 0; k < bTruth.getRowDimension(); k++) {
                        if (getHomogeneousIndex()[k]) {
                            hString = hString + "X ";
                        } else {
                            hString = hString + "O ";
                        }
                    }
                    hString = hString + "]";
                    System.out.print("Composite estimated beta: " + pmUtility.stringPrettyPrintVector(compositeEstimatedBeta) + " " + hString + " " + testY.get(i, 0) + " " + getPredictedY(xi, compositeEstimatedBeta, rng) + "\n");
                }
                //pmUtility.prettyPrintVector(compositeEstimatedBeta);

                outOfSampleResultsBeta += pmUtility.sumSquaredElements((compositeEstimatedBeta.minus(bTruth)));
            }

            outOfSampleResultsBeta /= testZ.getRowDimension();
            outOfSampleResultsY /= testZ.getRowDimension();
        }

        /**
         * Awesome, this works with discrete types where we seed it with the
         * right set of possible types Need to do several things:
         *
         * 1. I am imposing the truth on the homogeneous parameter; need to have
         * a search over that with FKRB as a nested inside loop (DONE: WORKS)
         *
         * 2. Want to impose probability constraints on OLS using Gurobi (after
         * 2 hours, that isn't going to happen; the Java code for this is
         * insanely unwieldy) I am going to try some basic barrier constraint
         * stuff where I successively increase the penalty function on the
         * constraints and see if that works Looks like the quick and dirty
         * works; let's do point 4 next and see if it still works The short
         * answer is that it sort of works. Need a proper constrained optimizer
         * here in the future.
         *
         * 3. Want to extend to using normal distributions instead of point
         * masses to get smooth CDFs
         *
         * 4. Want to add a bunch more points; in a single dimension that
         * doesn't help (without a proper solver anyways) probably due to
         * collinearity of the predicted models as you add more and more points
         *
         * 5. Need to connect the X's that the tree splits on to automate the RC
         * testing / estimation procedure
         */
        // to point 5 above, query the moment forest to see which variables
        // are ever split on
        boolean rcDetected = false;
        
        ArrayList<Integer> rcVarIndices = new ArrayList<>();
        double[] countSplitVariables = myForest.getCountSplitVariables();
        System.out.print("Detected random coefficients on following X indices: ");
        for (int i = 0; i < countSplitVariables.length; i++) {
            if (countSplitVariables[i] > 0) {
                rcDetected = true;
                rcVarIndices.add(i);
                System.out.print(i + " " + countSplitVariables[i]);
            }
        }
        System.out.println("");
        

        if (rcDetected) {
            DeconvolutionSolver desolve = new DeconvolutionSolver(testY, testX, this, rcVarIndices);
            desolve.solve();
            ArrayList<double[]> betaList = desolve.getBetaList();
            double[] betaWeights = desolve.getBetaWeights();

            boolean plot = true;

            if (plot) {
                // let's plot the fitted distribution of F(\beta)
                XYSeries xy = new XYSeries("Estimated");
                XYSeries xytruth = new XYSeries("True");

                boolean first = true;
                double minBetaSupport = 0;
                double maxBetaSupport = 0;

                for (int i = 0; i < betaList.size(); i++) {
                    xy.add(betaList.get(i)[1], betaWeights[i]);
                    if (betaList.get(i)[1] < minBetaSupport || first) {
                        minBetaSupport = betaList.get(i)[1];
                        first = false;
                    }
                    if (betaList.get(i)[1] > maxBetaSupport || first) {
                        maxBetaSupport = betaList.get(i)[1];
                        first = false;
                    }
                }
                double numPointsPlot = 100;
                double plotIncrement = (maxBetaSupport - minBetaSupport) / (numPointsPlot - 1);
                double numDraws = 100000;

                for (double beta = minBetaSupport; beta < maxBetaSupport; beta += plotIncrement) {
                    double mixtureF = 0;

                    for (int k = 0; k < numDraws; k++) {
                        mixtureF += normal.probability(beta - getBetaTruth(null, rng).get(1, 0));
                    }
                    mixtureF /= numDraws;
                    xytruth.add(beta, mixtureF);
                }

                XYSeriesCollection xyc = new XYSeriesCollection(xy);
                xyc.addSeries(xytruth);
                ChartGenerator.makeXYLine(xyc, "Fitted f(\\beta)", "\\beta", "f(\\beta)");
            }

            for (int r = 0; r < betaList.size(); r++) {
                System.out.format("beta1: %.3f beta2: %.3f weight: %.3f %n", betaList.get(r)[0], betaList.get(r)[1], betaWeights[r]);
            }
            // pmUtility.prettyPrintVector(weights);
            for (int i = 0; i < 10; i++) {
                Jama.Matrix xi = testX.getMatrix(i, i, 0, testX.getColumnDimension() - 1);
                double fittedShare = 0;
                for (int r = 0; r < betaList.size(); r++) {
                    Jama.Matrix beta = new Jama.Matrix(testX.getColumnDimension(), 1);
                    for (int j = 0; j < testX.getColumnDimension(); j++) {
                        beta.set(j, 0, betaList.get(r)[j]);
                    }
                    fittedShare += betaWeights[r] * LogitRCDataGenerator.getLogitShare(xi, beta);
                }
                System.out.format("Y: %.4f Fitted: %.4f %n", testY.get(i, 0), fittedShare);
            }
        }

        // jt.append("betaMSE: " + (outOfSampleFit / testZ.getRowDimension()) + " \t [" + rngSeed + "]\n");
        return new OutOfSampleStatisticsContainer(outOfSampleResultsBeta, outOfSampleResultsY);
    }

    @Override
    public DataLens getOutOfSampleXYZ(int numObsOOS, long rngSeed) {
        LogitRCDataGenerator xyz = new LogitRCDataGenerator(numObsOOS, this, rngSeed);
        // again, the Z's are X's here
        return new DataLens(xyz.getX(), xyz.getY(), xyz.getX(), null);
    }

    @Override
    public void loadData(long rngSeed) {
        if (!MONTE_CARLO) {
            int numObsFile = 0;
            try {
                BufferedReader in = new BufferedReader(new FileReader(filename)); // Inputting data. What is the cd here?
                String line = in.readLine(); // headers
                while (line != null) { // Each line is an observation
                    line = in.readLine(); // Read in data line by line
                    if (line != null) {
                        numObsFile++;
                    }
                }
                in.close();
                numObsFile /= 1;
                numObs = numObsFile;
            } catch (Exception e) {
                e.printStackTrace();
            }

            System.out.format("Number of observations = %,d %n", numObsFile);
            Jama.Matrix dX = new Jama.Matrix(numObsFile, 3); // Used previous loop to create arrays of the correct size, memory saver?
            Jama.Matrix dY = new Jama.Matrix(numObsFile, 1);

            try {
                BufferedReader in = new BufferedReader(new FileReader(filename));
                String line = in.readLine(); // headers
                int i = 0;
                while (line != null) {
                    line = in.readLine();
                    if (line != null) {
                        int a = 0;
                        int b = line.indexOf(",", a); //Returns the index within this string of the first occurrence of "," starting at 0; this is a comma delimited file
                        // System.out.println(line);
                        dY.set(i, 0, Double.valueOf(line.substring(a, b))); // outcome, assuming lines are comma delimitted and y value begins at index 0 and ends before the comma

                        a = b + 1;
                        b = line.indexOf(",", a); // treatment status is the next outcome in the comma delimited line
                        dX.set(i, 0, Double.valueOf(line.substring(a, b))); // treatment status

                        a = b + 1;
                        b = line.indexOf(",", a); // X1 is the next outcome in the comma delimited line
                        dX.set(i, 1, Double.valueOf(line.substring(a, b))); //X1

                        a = b + 1;
                        b = line.indexOf(",", a); // X1 is the next outcome in the comma delimited line
                        dX.set(i, 2, Double.valueOf(line.substring(a))); //X2

                        i++;
                    }
                }
                X = dX;
                Y = dY;

                in.close();
            } catch (Exception e) {
                e.printStackTrace();
            }
        } else {

            LogitRCDataGenerator xyz = new LogitRCDataGenerator(numObs, this, rngSeed);
            X = xyz.getX();
            Y = xyz.getY();
            Z = xyz.getZ();

        }
//        pmUtility.prettyPrint(pmUtility.concatMatrix(Y,pmUtility.ooncatMatrix(X,Z)));
//        System.exit(0);
    }

    @Override
    public String getVariableName(int variableIndex) {
        if (variableIndex == 0) {
            return "X1";
        }
        if (variableIndex == 1) {
            return "X2";
        }
        return "Unknown";
    }

    @Override
    public String getFixedEffectName(int variableIndex, int fixedEffectIndex) {
        return "Group " + fixedEffectIndex;
    }

    @Override
    public String formatTreeLeafOutput(Matrix beta, Matrix variance) {
        if (beta == null) {
            return "null (shouldn't be here!)";
        }
        // double b = beta.get(0, 0);
//        double se = Math.sqrt(variance.get(0, 0));
//        String stars = "";
//        NormalDistribution normal = new NormalDistribution(0, 1);
//        if (Math.abs(b / se) > Math.abs(normal.inverse(0.90))) {
//            stars = "*";
//        }
//        if (Math.abs(b / se) > Math.abs(normal.inverse(0.95))) {
//            stars = "**";
//        }
//        if (Math.abs(b / se) > Math.abs(normal.inverse(0.99))) {
//            stars = "***";
//        }
//        return String.format("%.2f (%.2f) %s", b, se, stars);
        return pmUtility.stringPrettyPrintVector(beta);
    }

    /**
     * @return the homogeneityIndex
     */
    public boolean[] getHomogeneityIndex() {
        return homogeneityIndex;
    }

    @Override
    public void resetHomogeneityIndex() {
        for (int i = 0; i < homogeneityIndex.length; i++) {
            homogeneityIndex[i] = false;
        }
    }

    @Override
    public void setHomogeneousIndex(Integer i) {
        homogeneityIndex[i] = true;
    }

    /**
     * @return the homogeneousParameterVector
     */
    @Override
    public Jama.Matrix getHomogeneousParameterVector() {
        return homogeneousParameterVector;
    }

    @Override
    public ContainerMoment getContainerMoment(DataLens lens) {
        /**
         * Should this boolean ever be true???
         *
         * I think the answer is no since the only place that calls this is in
         * HomogeneousSearchContainer, and the reason that I added the boolean
         * (allParametersHomogeneous) is to estimate the stump parameters when
         * the parameter are already determined to be all homogeneous
         */
        return new ContainerLogitRC(lens, homogeneityIndex, homogeneousParameterVector, false);
    }

    @Override
    public int getNumMoments() {
        return 2;
    }

    @Override
    public int getNumParams() {
        throw new UnsupportedOperationException("Not supported yet."); // Generated from nbfs://nbhost/SystemFileSystem/Templates/Classes/Code/GeneratedMethodBody
    }

}
