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
package examples.linear;

// import JSci.maths.statistics.NormalDistribution;
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
public class LinearMomentSpecification implements MomentSpecification {

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
    
    /**
     * We are going to control homogeneous parameters through these variables
     */
    private boolean[] homogeneityIndex = {false, false};
    Jama.Matrix homogeneousParameterVector = null; // this is a compact vector (only consists of the parameters we are imposing for homogeneity)

    public LinearMomentSpecification(int numObs) {
        this.numObs = numObs;
        int[] vsi = {0, 1, 2}; //Search over z1, z2, z3 
        Boolean[] wvd = {false, false, true}; // z1, z2 continuous, z3 discrete
        variableSearchIndex = vsi;
        DiscreteVariables = wvd;
    }

    public LinearMomentSpecification(String filename) {
        this.filename = filename;
    }

    @Override
    public void setHomogeneousParameters(Matrix h) {
        this.homogeneousParameterVector = h;
    }

    @Override
    public Matrix getBalancingVector() {
        return balancingVector;
    }

    public LinearMomentSpecification(Jama.Matrix X, Jama.Matrix Y, Jama.Matrix Z, int numtrees, int[] variableSearchIndex, Boolean[] DiscreteVariables) {
        this.X = X;
        this.Y = Y;
        this.Z = Z;
        this.numtrees = numtrees;
        this.variableSearchIndex = variableSearchIndex;
        this.DiscreteVariables = DiscreteVariables;
        balancingVector = pmUtility.getColumn(X, 0); // treatment status is the first column of X
    }

    @Override
    public Double getPredictedY(Matrix xi, Jama.Matrix beta) {
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
        return new MomentContinuousSplitObjLinear(indexSplitVariable, lens, minProportionEachPartition, minCountEachPartition);
    }

    @Override
    public MomentPartitionObj getMomentPartitionObj(DataLens lens, int indexSplitVariable, IntegerPartition partition) {
        return new MomentPartitionObjLinear(partition, indexSplitVariable, lens);
    }

    @Override
    public ContainerMoment computeOptimalBeta(DataLens lens) {
        ContainerLinear l = new ContainerLinear(lens);
        l.computeBetaAndErrors();
        return l;
    }

    @Override
    public Matrix getY(boolean residualizeY) {
        /**
         * Want to return the residualized Y (taking out the part explained by
         * the globally-imposed homogeneous subset of parameters)
         *
         * I am going to impose the convention that the parameter vector of
         * homogeneous parameters is only as long as the number of restrictions
         * (as opposed to having it be a full length of X with some zeros for
         * non-restrictions)
         */
        if (!residualizeY) {
            return Y;
        }
        Jama.Matrix residualizedY = Y.copy();

        int homogeneousParameterIndex = 0;
        for (int k = 0; k < homogeneityIndex.length; k++) {
            if (homogeneityIndex[k]) {
                for (int i = 0; i < Y.getRowDimension(); i++) {
                    residualizedY.set(i, 0, residualizedY.get(i, 0) - X.get(i, k) * homogeneousParameterVector.get(homogeneousParameterIndex, 0));
                }
                homogeneousParameterIndex++;
            }
        }
        return residualizedY;
    }

    @Override
    public double getHomogeneousComponent(Jama.Matrix xi) {
        double homogeneousComponent = 0;
        int homogeneousParameterIndex = 0;
        for (int k = 0; k < homogeneityIndex.length; k++) {
            if (homogeneityIndex[k]) {
                homogeneousComponent += xi.get(0, k) * homogeneousParameterVector.get(homogeneousParameterIndex, 0);

                homogeneousParameterIndex++;
            }
        }
        return homogeneousComponent;
    }

    @Override
    public Matrix getX() {
        /**
         * We want to return the part of X that is not restricted via a global
         * parameter homogeneity assumption
         */

        Jama.Matrix residualizedX = null;
        for (int k = 0; k < homogeneityIndex.length; k++) {
            if (!homogeneityIndex[k]) {
                if (residualizedX == null) {
                    residualizedX = pmUtility.getColumn(X, k);
                } else {
                    residualizedX = pmUtility.concatMatrix(residualizedX, pmUtility.getColumn(X, k));
                }
            }
        }
        return residualizedX;
    }
    
    public Jama.Matrix residualizeX(Jama.Matrix Xp) {
        Jama.Matrix residualizedX = null;
        for (int k = 0; k < homogeneityIndex.length; k++) {
            if (!homogeneityIndex[k]) {
                if (residualizedX == null) {
                    residualizedX = pmUtility.getColumn(Xp, k);
                } else {
                    residualizedX = pmUtility.concatMatrix(residualizedX, pmUtility.getColumn(Xp, k));
                }
            }
        }
        return residualizedX;
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

    //Return the true parameter vector for a given observation
    @Override
    public Matrix getBetaTruth(Matrix zi) {
        Jama.Matrix beta = new Jama.Matrix(2, 1); // Beta is a scalar
        beta.set(0, 0, -1);
        beta.set(1, 0, 1);

        boolean singleBeta = true;
        if (singleBeta) {
            return beta;
        }

        boolean simplest = true;
        if (simplest) {
            if (zi.get(0, 0) > 0) {
                beta.set(1, 0, -2);
            }
            return beta;
        }

        boolean imposeUniformBeta1 = true;
        if (zi.get(0, 0) > 0) { // if z1 > 0 \beta_0 = 1
            if (!imposeUniformBeta1) {
                beta.set(0, 0, 2);
            }
            if (zi.get(0, 1) < 0.5) { // if also z2<0.5, \beta_1 = -1;
                beta.set(1, 0, -3);
            }
        }
        if (zi.get(0, 2) == 1) {
            if (!imposeUniformBeta1) {
                beta.set(0, 0, 3);
            }
            beta.set(1, 0, beta.get(1, 0) * 1.5);
        }
        return beta;
    }

    @Override
    public DataLens getOutOfSampleXYZ(int numObsOOS, long rngSeed) {
        LinearDataGenerator xyz = new LinearDataGenerator(numObsOOS,this, rngSeed);        
        return new DataLens(xyz.getX(), xyz.getY(), xyz.getZ(), null);
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

            LinearDataGenerator xyz = new LinearDataGenerator(numObs, this, rngSeed);
            X = xyz.getX();
            Y = xyz.getY();
            Z = xyz.getZ();

        }
//        pmUtility.prettyPrint(pmUtility.concatMatrix(Y,pmUtility.concatMatrix(X,Z)));
//        System.exit(0);
    }

    @Override
    public String getVariableName(int variableIndex) {
        if (variableIndex == 0) {
            return "Z1";
        }
        if (variableIndex == 1) {
            return "Z2";
        }
        if (variableIndex == 2) {
            return "Z3";
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

}
