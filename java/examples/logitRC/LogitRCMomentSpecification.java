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
        beta.set(0, 0, -1.0);
        beta.set(1, 0, 2.3);

        // what does this even mean in a RC logit setting?
        // these the mean parameters of the utility, I guess
        
        return beta;
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

}
