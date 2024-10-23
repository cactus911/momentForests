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
package examples.CardIV;

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
public class MomentSpecificationCardIV implements MomentSpecification {

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
    // NEED TO UPDATE
    String[] varNames = {"ed76", "constant", "exp76", "exp762", "black", "south76", "city76", "region_1966", "smsa66r", "daded", "momed", "nodaded", "nomomed", "famed", "momdad14", "sinmom14", "near4"};

    /**
     * We are going to control homogeneous parameters through these variables
     */
    final private boolean[] homogeneityIndex; // = new boolean[X.getColumnDimension()];
    final private Jama.Matrix homogeneousParameterVector; // = new Jama.Matrix(X.getColumnDimension(), 1); // this is a compact vector (only consists of the parameters we are imposing for homogeneity)

    public MomentSpecificationCardIV(String filename) {
        this.filename = filename;

        // what specification am i going to use here?
        // y = alpha(Z)X + epsilon
        loadData(787);

        homogeneityIndex = new boolean[getNumParams()];
        homogeneousParameterVector = new Jama.Matrix(getNumParams(), 1);
        resetHomogeneityIndex();

        // NEED TO UPDATE
        /**
         * These all refer to Z (not X)!!!
         *
         * Z =
         *
         * 1. constant 0. education (flipped these below) 2. experience 3.
         * experience^2 4. dummy = 1 if black 5. dummy = 1 if living in the
         * South in 1976 6. dummy = 1 if living in SMSA in 1976 7. region in
         * 1966, categorical 8. dummy = 1 if living in SMSA in 1966 9. father's
         * years of education 10. mother's years of education 11. dummy = 1 if
         * father's education is missing 12. dummy = 1 if mother's education is
         * missing 13. interactions of family education, categorical 14. dummy =
         * 1 if household contains both parents 15. dummy = 1 if household is a
         * single mother
         *
         */
        int[] vsi = {4, 5, 6, 8, 9, 11, 12, 13, 14, 15};
        Boolean[] whichVariablesAreDiscrete = {false, // 0
            false,
            false,
            false,
            true, // black
            true,
            true,
            true,
            true,
            false,
            false,
            true,
            true,
            false, // ordered (interaction of education and family)
            true,
            true
        };
        variableSearchIndex = vsi;
        DiscreteVariables = whichVariablesAreDiscrete;
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
        return X.getColumnDimension() - 1; // one endo X (first column), one instrument (last column of "X")
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
        double error = yi - ((xi.getMatrix(0, 0, 0, xi.getColumnDimension() - 2)).times(beta)).get(0, 0);
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
            double yhat = ((xi.getMatrix(0, 0, 0, xi.getColumnDimension() - 2)).times(beta)).get(0, 0);
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
        return new MomentContinuousSplitObjCardIV(indexSplitVariable, lens, minProportionEachPartition, minCountEachPartition, this);
    }

    @Override
    public MomentPartitionObj getMomentPartitionObj(DataLens lens, int indexSplitVariable, IntegerPartition partition) {
        return new MomentPartitionObjCardIV(partition, indexSplitVariable, lens, this);
    }

    @Override
    public ContainerMoment computeOptimalBeta(DataLens lens, boolean allParametersHomogeneous) {
        // System.out.println("Are all parameters homogeneous? " + allParametersHomogeneous);
        ContainerCardIV l = new ContainerCardIV(lens, homogeneityIndex, homogeneousParameterVector, allParametersHomogeneous, this);
        l.computeBetaAndErrors();
        return l;
    }

    @Override
    public Matrix getY() {
        return Y;
    }

    @Override
    public Matrix getZ() {
        return Z;
    }

    @Override
    public Matrix getX() {
        return X;
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

            boolean subsample = false;
            if (subsample) {
                numObsFile = Math.floorDiv(numObsFile, 5);
            }

            numObs = numObsFile;
        } catch (Exception e) {
            e.printStackTrace();
        }

        System.out.format("Number of observations = %,d %n", numObsFile);

        // NEED TO UPDATE
        Jama.Matrix dX = new Jama.Matrix(numObsFile, 17);
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
                    dY.set(i, 0, Double.valueOf(line.substring(a, b)));

                    a = b + 1;
                    b = line.indexOf(",", a);
                    dX.set(i, 0, Double.valueOf(line.substring(a, b))); 	// constant

                    a = b + 1;
                    b = line.indexOf(",", a);
                    dX.set(i, 1, Double.valueOf(line.substring(a, b))); 	// education (POTENTIALLY ENDOGENOUS)

                    a = b + 1;
                    b = line.indexOf(",", a);
                    dX.set(i, 2, Double.valueOf(line.substring(a, b))); 	// experience

                    a = b + 1;
                    b = line.indexOf(",", a);
                    dX.set(i, 3, Double.valueOf(line.substring(a, b))); 	// experience^2

                    a = b + 1;
                    b = line.indexOf(",", a);
                    dX.set(i, 4, Double.valueOf(line.substring(a, b))); 	// dummy = 1 if black

                    a = b + 1;
                    b = line.indexOf(",", a);
                    dX.set(i, 5, Double.valueOf(line.substring(a, b))); 	// dummy = 1 if living in the South in 1976  

                    a = b + 1;
                    b = line.indexOf(",", a);
                    dX.set(i, 6, Double.valueOf(line.substring(a, b))); 	// dummy = 1 if living in SMSA in 1976

                    a = b + 1;
                    b = line.indexOf(",", a);
                    dX.set(i, 7, Double.valueOf(line.substring(a, b))); 	// region in 1966, categorical  

                    a = b + 1;
                    b = line.indexOf(",", a);
                    dX.set(i, 8, Double.valueOf(line.substring(a, b))); 	// dummy = 1 if living in SMSA in 1966 

                    a = b + 1;
                    b = line.indexOf(",", a);
                    dX.set(i, 9, Double.valueOf(line.substring(a, b))); 	// father's years of education 

                    a = b + 1;
                    b = line.indexOf(",", a);
                    dX.set(i, 10, Double.valueOf(line.substring(a, b))); 	// mother's years of education  

                    a = b + 1;
                    b = line.indexOf(",", a);
                    dX.set(i, 11, Double.valueOf(line.substring(a, b))); 	// dummy = 1 if father's education is missing 

                    a = b + 1;
                    b = line.indexOf(",", a);
                    dX.set(i, 12, Double.valueOf(line.substring(a, b))); 	// dummy = 1 if mother's education is missing 

                    a = b + 1;
                    b = line.indexOf(",", a);
                    dX.set(i, 13, Double.valueOf(line.substring(a, b))); 	// interactions of family education, categorical 

                    a = b + 1;
                    b = line.indexOf(",", a);
                    dX.set(i, 14, Double.valueOf(line.substring(a, b))); 	// dummy = 1 if household contains both parents

                    a = b + 1;
                    b = line.indexOf(",", a);
                    dX.set(i, 15, Double.valueOf(line.substring(a, b))); 	// dummy = 1 if single mom

                    a = b + 1;
                    b = line.indexOf(",", a);
                    dX.set(i, 16, Double.valueOf(line.substring(a))); 		// dummy = 1 if near a 4-year college (instrument)

                    i++;
                }
            }

            /**
             * NATHAN: going to put the endogenous X here in the first column to
             * make keeping track of it easier
             */
            X = pmUtility.getColumn(dX, 1); // the endogenous variable
            X = pmUtility.concatMatrix(X, pmUtility.getColumn(dX, 0));
            X = pmUtility.concatMatrix(X, pmUtility.getColumn(dX, 2));
            //X = pmUtility.concatMatrix(X, pmUtility.getColumn(dX, 3)); 
            //X = pmUtility.concatMatrix(X, pmUtility.getColumn(dX, 4)); 
            //X = pmUtility.concatMatrix(X, pmUtility.getColumn(dX, 5)); 
            //X = pmUtility.concatMatrix(X, pmUtility.getColumn(dX, 6)); 
            //X = pmUtility.concatMatrix(X, pmUtility.getColumn(dX, 7)); 
            //X = pmUtility.concatMatrix(X, pmUtility.getColumn(dX, 8));
            //X = pmUtility.concatMatrix(X, pmUtility.getColumn(dX, 9)); 
            //X = pmUtility.concatMatrix(X, pmUtility.getColumn(dX, 10)); 
            //X = pmUtility.concatMatrix(X, pmUtility.getColumn(dX, 11)); 
            //X = pmUtility.concatMatrix(X, pmUtility.getColumn(dX, 12)); 
            //X = pmUtility.concatMatrix(X, pmUtility.getColumn(dX, 13)); 
            //X = pmUtility.concatMatrix(X, pmUtility.getColumn(dX, 14)); 
            //X = pmUtility.concatMatrix(X, pmUtility.getColumn(dX, 15));
            X = pmUtility.concatMatrix(X, pmUtility.getColumn(dX, 16));  	// the instrument

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
             * Maybe try to just treat them all as continuous to begin; little
             * experiment shows the exact same fit doing both approaches
             *
             * ALSO: categorical variables don't have to enter the OLS model
             * (X). They only have to enter on Z! The FE are absorbed into the
             * constant in each cell (conditioning on splitting on it!).
             */
            Y = dY;

            // SPLITTING VARIABLES FOR HETEROGENEITY
            //Z = dX.copy();  
            Z = pmUtility.getColumn(dX, 0);
            Z = pmUtility.concatMatrix(Z, pmUtility.getColumn(dX, 1));
            Z = pmUtility.concatMatrix(Z, pmUtility.getColumn(dX, 2));
            Z = pmUtility.concatMatrix(Z, pmUtility.getColumn(dX, 3));
            Z = pmUtility.concatMatrix(Z, pmUtility.getColumn(dX, 4));
            Z = pmUtility.concatMatrix(Z, pmUtility.getColumn(dX, 5));
            Z = pmUtility.concatMatrix(Z, pmUtility.getColumn(dX, 6));
            Z = pmUtility.concatMatrix(Z, pmUtility.getColumn(dX, 7));
            Z = pmUtility.concatMatrix(Z, pmUtility.getColumn(dX, 8));
            Z = pmUtility.concatMatrix(Z, pmUtility.getColumn(dX, 9));
            Z = pmUtility.concatMatrix(Z, pmUtility.getColumn(dX, 10));
            Z = pmUtility.concatMatrix(Z, pmUtility.getColumn(dX, 11));
            Z = pmUtility.concatMatrix(Z, pmUtility.getColumn(dX, 12));
            Z = pmUtility.concatMatrix(Z, pmUtility.getColumn(dX, 13));
            Z = pmUtility.concatMatrix(Z, pmUtility.getColumn(dX, 14));
            Z = pmUtility.concatMatrix(Z, pmUtility.getColumn(dX, 15));

            boolean FAKE_DATA = true;
            if (FAKE_DATA) {
                Random rng = new Random(rngSeed);
                NormalDistribution normal = new NormalDistribution();
                for (i = 0; i < numObs; i++) {
                    // generate X that is correlated with Z
                    X.set(i, 3, 3 * rng.nextDouble()); // this is now the instrument
                    double error = normal.inverse(rng.nextDouble());
                    X.set(i, 0, X.get(i, 3) + 0.3 * error); // this is the endogenous variable

                    // generate Y with correlation between one X and the error
                    double beta1 = -5.0;
                    double beta2 = 3.0;
                    boolean observableHeterogeneity = true;
                    if (observableHeterogeneity) {
                        if (Z.get(i, 4) == 1.0) {
                            // System.out.println("black");
                            beta1 = 5.0;
                        } else {
                            // System.out.println("non-black");
                        }

                        if (Z.get(i, 5) == 1.0 && Z.get(i, 4) == 1.0) {
                            beta2 = 2.0;
                        }
                    }
                    Y.set(i, 0, X.get(i, 0) * 1 + X.get(i, 1) * beta1 + X.get(i, 2) * beta2 + error);
                }
            }
            System.out.println("Mean of Y: " + pmUtility.mean(Y, 0));
            System.out.print("OLS: ");
            pmUtility.prettyPrintVector(pmUtility.OLS(X.getMatrix(0, X.getRowDimension() - 1, 0, X.getColumnDimension() - 2), Y, false));

            in.close();
        } catch (Exception e) {
            e.printStackTrace();
        }

        // have to decide what instrumentVariables'm going to split on (Z vector)
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
        return "" + fixedEffectIndex;
    }

    @Override
    public String formatTreeLeafOutput(Matrix beta, Matrix variance) {
        if (beta == null) {
            System.out.println("Null beta, shouldn't be here!");
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
        // System.out.print("Returning homogeneous parameter vector: ");
        // pmUtility.prettyPrintVector(homogeneousParameterVector);
        return homogeneousParameterVector;
    }

    @Override
    public ContainerMoment getContainerMoment(DataLens lens) {
        return new ContainerCardIV(lens, homogeneityIndex, homogeneousParameterVector, false, this);
    }

    @Override
    public int getNumParams() {
        return X.getColumnDimension() - 1;
    }
}
