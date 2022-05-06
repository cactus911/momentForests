/*
 * The MIT License
 *
 * Copyright 2022 Stephen P. Ryan <stephen.p.ryan@wustl.edu>.
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
    
    int dimensionX = 5;

    /**
     * We are going to control homogeneous parameters through these variables
     */
    private boolean[] homogeneityIndex; // = new boolean[X.getColumnDimension()];
    private Jama.Matrix homogeneousParameterVector; // = new Jama.Matrix(X.getColumnDimension(), 1); // this is a compact vector (only consists of the parameters we are imposing for homogeneity)
    // i'm going to change that, since that sounds like a recipe for disaster as the set of parameters that gets the homogeneous label shifts around
    // how to keep everything straight--easier way is to just to make it the size of the parameter vector and go from there (and never read elements that aren't labeled as homogeneous)

    public LinearMomentSpecification(int numObs, int dimX) {
        this.numObs = numObs;
        this.dimensionX = dimX;
        // all these indices are hard-coded; want to change that down the road!
        homogeneityIndex = new boolean[dimensionX];
        homogeneousParameterVector = new Jama.Matrix(dimensionX,1);
        resetHomogeneityIndex();
        int[] vsi = {0, 1, 2}; //Search over z1, z2, z3 
        Boolean[] wvd = {false, false, true}; // z1, z2 continuous, z3 discrete
        variableSearchIndex = vsi;
        DiscreteVariables = wvd;
    }

    public LinearMomentSpecification(String filename) {
        this.filename = filename;
    }

    @Override
    public boolean[] getHomogeneousIndex() {
        return homogeneityIndex;
    }

    @Override
    public double getHomogeneousParameter(int parameterIndex) {
        return homogeneousParameterVector.get(parameterIndex,0);
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
    public double getGoodnessOfFit(double yi, Matrix xi, Matrix beta) {
        double fit = (xi.times(beta)).get(0,0);
        double error = yi-fit;
        // System.out.println(fit+" "+yi);
        return error*error;
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

    //Return the true parameter vector for a given observation
    @Override
    public Matrix getBetaTruth(Matrix zi, Random rng) {
        Jama.Matrix beta = new Jama.Matrix(dimensionX, 1, 1); // Beta is a scalar
        beta.set(0, 0, -1);
        beta.set(1, 0, 1.0);
        
        boolean partiallyLinearModel = true;
        if(partiallyLinearModel) {
            // want to get the model y = x\beta + g(z), or x\beta+1*\beta(Z) where the second function is complex (like a cosine function?)
            beta.set(0, 0, 2.5*Math.sin(zi.get(0,0)) + 0.25*Math.pow(zi.get(0,0),2));
            return beta;
        }        

        boolean singleBeta = false;
        if (singleBeta) {
            return beta;
        }

        boolean oneDimensionHeterogeneity = false;
        if (oneDimensionHeterogeneity) {
            if (zi.get(0, 0) > 0) {
                beta.set(1, 0, -1.05);
            }
            return beta; 
        }
        
        boolean twoDimensionHeterogeneity = false;
        if (twoDimensionHeterogeneity) {
            if (zi.get(0, 0) > 0) {
                beta.set(0, 0, 0.33);
                beta.set(1, 0, -1.05);
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
        LinearDataGenerator xyz = new LinearDataGenerator(numObsOOS, this, rngSeed, dimensionX);
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

            LinearDataGenerator xyz = new LinearDataGenerator(numObs, this, rngSeed, dimensionX);
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
