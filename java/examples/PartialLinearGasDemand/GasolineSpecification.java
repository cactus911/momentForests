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
package examples.PartialLinearGasDemand;

// import JSci.maths.statistics.NormalDistribution;
import JSci.maths.statistics.NormalDistribution;
import Jama.Matrix;
import core.ContainerMoment;
import core.DataLens;
import core.IntegerPartition;
import core.MomentContinuousSplitObj;
import core.MomentPartitionObj;
import core.MomentSpecification;
import java.io.FileReader;
import java.io.BufferedReader;
import java.util.Random;
import utility.pmUtility;

/**
 *
 * @author Stephen P. Ryan <stephen.p.ryan@wustl.edu>
 */
public class GasolineSpecification implements MomentSpecification {

    Jama.Matrix X;
    Jama.Matrix Y;
    Jama.Matrix Z;
    Jama.Matrix balancingVector; // is treatment status in the RCT setting
    int numObs;
    int numtrees;
    int[] variableSearchIndex; // this should be restricted to only Z
    Boolean[] DiscreteVariables; // also this should be restricted to only Z
    String filename;
    boolean MONTE_CARLO = false;

    DataLens outSampleLens;

    String[] varNames = {"constant", "logHHSize", "logNumDrivers", "logAge", "urban", "hhfaminc", "census_d", "lif_cyc", "logCost"};

    /**
     * We are going to control homogeneous parameters through these variables
     */
    final private boolean[] homogeneityIndex; // = new boolean[X.getColumnDimension()];
    final private Jama.Matrix homogeneousParameterVector; // = new Jama.Matrix(X.getColumnDimension(), 1); // this is a compact vector (only consists of the parameters we are imposing for homogeneity)

    public GasolineSpecification(String filename) {
        this.filename = filename;

        // what specification am i going to use here?
        // y = x'b + h(t) + e
        // the baseline models is going to be y = \alpha + \beta' X + \epsilon
        // we are going to let Z = X
        // if we don't split, then we have OLS for the whole sample
        // if we have NP component on some X, then we keep splitting on X over and over again
        // will we have interaction of X with \beta'X though? is that weird?
        // we'll grow the unrestricted tree, detect the homogeneous parameters, then let the tree
        // go to town holding those fixed across leaves; should work
        // observation: fixed effects don't need to be in the OLS baseline model
        // as long as an intercept is there and we allow the tree to split on those categorical variables as Z's
        // we can build the fixed effects within the estimator
        loadData(787);

        homogeneityIndex = new boolean[X.getColumnDimension()];
        homogeneousParameterVector = new Jama.Matrix(X.getColumnDimension(), 1);
        resetHomogeneityIndex();

        /**
         * These all refer to Z (not X)!!!
         *
         * Z =
         *
         * 0.constant
         *
         * 1. logHHSize
         *
         * 2. logNumDrivers
         *
         * 3. logAge (continuous due to averaging within family)
         *
         * 4. urban
         *
         * 5. hhfaminc
         *
         * 6. census_d
         *
         * 7. lif_cyc
         *
         * 8. logCost (this is the only truly continuous variable in the data
         * set)
         */
        int[] vsi = {1, 2, 3, 4, 5, 6, 7};
        // Boolean[] wvd = {true, false, false, false, true, true, true, true, false};
        // Boolean[] wvd = {true, false, false, false, true, false, true, true, false};
        // Boolean[] wvd = {true, true, true, false, true, true, true, true, true};
        // since we can split along here, why not have it treat everything continuously, since it does a grid search?
        Boolean[] wvd = {false, false, false, false, false, false, false, false, false};
        variableSearchIndex = vsi;
        DiscreteVariables = wvd;
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
    public int getNumMoments() {
        return X.getColumnDimension();
    }

    @Override
    public void setHomogeneousParameter(int parameterIndex, double value) {
        homogeneousParameterVector.set(parameterIndex, 0, value);
    }

    @Override
    public Matrix getBalancingVector() {
        return balancingVector;
    }

    @Override
    public double getGoodnessOfFit(double yi, Matrix xi, Matrix beta) {
        double fit = (xi.times(beta)).get(0, 0);
        double error = yi - fit;
        // System.out.println(fit+" "+yi);
        return error * error;
    }

    @Override
    public Double getPredictedY(Matrix xi, Jama.Matrix beta, Random rng) {
        /**
         * This may have to be adjusted when we impose homogeneity, depending on
         * what this is used for
         */
        // pmUtility.prettyPrint(xi);
        if (beta != null) {
            double yhat = (xi.times(beta)).get(0, 0);
            return yhat;
        }
        return null;
    }

    @Override
    public int[] getVariableIndicesToSearchOver() {
        return variableSearchIndex;
    }

    @Override
    public MomentContinuousSplitObj getFminObjective(DataLens lens, int indexSplitVariable, double minProportionEachPartition, int minCountEachPartition) {
        return new MomentContinuousSplitObjLinear(indexSplitVariable, lens, minProportionEachPartition, minCountEachPartition, this);
    }

    @Override
    public MomentPartitionObj getMomentPartitionObj(DataLens lens, int indexSplitVariable, IntegerPartition partition) {
        return new MomentPartitionObjLinear(partition, indexSplitVariable, lens, this);
    }

    @Override
    public ContainerMoment computeOptimalBeta(DataLens lens, boolean allParametersHomogeneous) {
        ContainerLinear l = new ContainerLinear(lens, homogeneityIndex, homogeneousParameterVector, allParametersHomogeneous);
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
        return Z;
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
    public Matrix getBetaTruth(Matrix zi, Random rng) {
        // we don't know, this shouldn't be called in a real application
        throw new UnsupportedOperationException("Not supported yet."); // Generated from nbfs://nbhost/SystemFileSystem/Templates/Classes/Code/GeneratedMethodBody
    }

    @Override
    public DataLens getOutOfSampleXYZ(int numObsOutOfSample, long rngSeed) {
        return outSampleLens;
    }

    @Override
    public void loadData(long rngSeed) {

        int numObsFile = 0;
        try {
            BufferedReader in = new BufferedReader(new FileReader(filename));
            String line = in.readLine(); // headers
            while (line != null) { // Each line is an observation
                line = in.readLine(); // Read in data line by line
                if (line != null) {
                    numObsFile++;
                }
            }
            in.close();

            boolean subsample = true;
            if (subsample) {
                numObsFile = Math.floorDiv(numObsFile, 10);
            }

            numObs = numObsFile;
        } catch (Exception e) {
            e.printStackTrace();
        }

        System.out.format("Number of observations = %,d %n", numObsFile);

        // in the gasoline data, we have 9 X variables and one Y
        // in order, they are: logGallons (this is the Y) constant	logHHSize	logNumDrivers	logAge	urban	hhfaminc	census_d	lif_cyc	logCost
        Jama.Matrix dX = new Jama.Matrix(numObsFile, 9);
        Jama.Matrix dY = new Jama.Matrix(numObsFile, 1);

        try {
            BufferedReader in = new BufferedReader(new FileReader(filename));
            String line = in.readLine(); // headers
            int i = 0;
            for (int obs = 0; obs < numObs; obs++) {
                line = in.readLine();
                // System.out.println(line);
                if (line != null) {
                    int a = 0;
                    int b = line.indexOf(",", a); //Returns the index within this string of the first occurrence of "," starting at 0; this is a comma delimited file
                    // System.out.println(line);
                    dY.set(i, 0, Double.valueOf(line.substring(a, b))); // log gallons used

                    a = b + 1;
                    b = line.indexOf(",", a);
                    dX.set(i, 0, Double.valueOf(line.substring(a, b))); // constant

                    a = b + 1;
                    b = line.indexOf(",", a);
                    dX.set(i, 1, Double.valueOf(line.substring(a, b))); // log household size

                    a = b + 1;
                    b = line.indexOf(",", a);
                    dX.set(i, 2, Double.valueOf(line.substring(a, b))); // log number drivers

                    a = b + 1;
                    b = line.indexOf(",", a);
                    dX.set(i, 3, Double.valueOf(line.substring(a, b))); // log age

                    a = b + 1;
                    b = line.indexOf(",", a);
                    dX.set(i, 4, Double.valueOf(line.substring(a, b))); // urban (categorical)

                    a = b + 1;
                    b = line.indexOf(",", a);
                    dX.set(i, 5, Double.valueOf(line.substring(a, b))); // household family income (categorical)

                    a = b + 1;
                    b = line.indexOf(",", a);
                    dX.set(i, 6, Double.valueOf(line.substring(a, b))); // census district (categorical)

                    a = b + 1;
                    b = line.indexOf(",", a);
                    dX.set(i, 7, Double.valueOf(line.substring(a, b))); // life cycle (categorical)

                    a = b + 1;
                    b = line.indexOf(",", a);
                    dX.set(i, 8, Double.valueOf(line.substring(a))); // log cost per gallon

                    i++;
                }
            }

            // going to pull a subset of loaded X's to use in base linear model
            X = pmUtility.getColumn(dX, 0); // constant

            // when this is just a constant, we have the standard regression tree
            // the S&S specification has all of these variables in it
            X = pmUtility.concatMatrix(X, pmUtility.getColumn(dX, 1)); // log household size
            X = pmUtility.concatMatrix(X, pmUtility.getColumn(dX, 2)); // log number drivers
            // X = pmUtility.concatMatrix(X, pmUtility.getColumn(dX, 3)); // log age
            // X = pmUtility.concatMatrix(X, pmUtility.getColumn(dX, 4)); // catogorical: urban 
            // X = pmUtility.concatMatrix(X, pmUtility.getColumn(dX, 5)); // categorical: family income
            // X = pmUtility.concatMatrix(X, pmUtility.getColumn(dX, 6)); // categorical: census district
            // X = pmUtility.concatMatrix(X, pmUtility.getColumn(dX, 7)); // categorical: life cycle
            /**
             * MAJOR POINT: ContainerLinear has no idea how to deal with
             * categorical variables right now since they are stacked and not
             * expanded
             *
             * TODO: implement categorical function (expansion?)
             *
             * There is a strange observation here that you can actually treat
             * them as continuous variables; if they actually matter, the forest
             * should split on their discrete values and estimate separate
             * subtrees for each of those splits!
             *
             * Maybe try to just treat them all as continuous to begin
             */
            // X = pmUtility.concatMatrix(X, pmUtility.getColumn(dX, 1)); // cost per gallon
            Y = dY;

            boolean runMonteCarlo = false;
            if (runMonteCarlo) {
                /**
                 * Easy Monte Carlo here for testing purposes
                 */
                NormalDistribution normal = new NormalDistribution();
                Random rng = new Random(rngSeed);
                for (int w = 0; w < Y.getRowDimension(); w++) {
                    if (dX.get(w, 5) > 5) {
                        Y.set(w, 0, 4.5 + dX.get(w, 2) + 0.1 * normal.inverse(rng.nextDouble()));
                    } else {
                        Y.set(w, 0, 4.5 - dX.get(w, 2) + 0.1 * normal.inverse(rng.nextDouble()));
                    }
                }
            }
            System.out.println("Mean of Y: " + pmUtility.mean(Y, 0));

            Z = dX.copy();

            // can split these into in-sample and out-of-sample here
            int cutoff = (int) Math.round(X.getRowDimension() * 0.9);

            DataLens lens = new DataLens(X, Y, Z, null);
            DataLens inSample = lens.getSubsetData(0, cutoff);
            outSampleLens = lens.getSubsetData(cutoff, X.getRowDimension() - 1);
            X = inSample.getX();
            Y = inSample.getY();
            Z = inSample.getZ();

//            pmUtility.prettyPrint(dX, 10);
//            System.out.println("----");
//            pmUtility.prettyPrint(Z, 10);
//            System.exit(0);
            in.close();
        } catch (Exception e) {
            e.printStackTrace();
        }

        // have to decide what I'm going to split on (Z vector)
        // System.exit(0);
//        pmUtility.prettyPrint(pmUtility.concatMatrix(Y,pmUtility.concatMatrix(X,Z)));
//        System.exit(0);
    }

    @Override
    public String getVariableName(int variableIndex) {
        return varNames[variableIndex];
    }

    @Override
    public String getFixedEffectName(int variableIndex, int fixedEffectIndex) {
        // return "Group " + fixedEffectIndex;
        return varNames[variableIndex] + " " + fixedEffectIndex;
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
    public Jama.Matrix getHomogeneousParameterVector() {
        return homogeneousParameterVector;
    }

    @Override
    public ContainerMoment getContainerMoment(DataLens lens) {
        return new ContainerLinear(lens, homogeneityIndex, homogeneousParameterVector, false);
    }

}
