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
package examples.SimpleRCT;

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
import utility.pmUtility;

/**
 *
 * NOTE: REQUIRES THAT FIRST COLUMN OF X BE THE INDICATOR FOR TREATMENT / CONTROL
 * 
 * @author Stephen P. Ryan <stephen.p.ryan@wustl.edu>
 */
public class SimpleRCTMomentSpecification implements MomentSpecification {

    Jama.Matrix X;
    Jama.Matrix Y;
    Jama.Matrix balancingVector; // is treatment status in the RCT setting
    int numObs;
    int numtrees;
    int[] variableSearchIndex;
    Boolean[] DiscreteVariables;

    public SimpleRCTMomentSpecification(int numObs) {
        this.numObs = numObs;
        int[] vsi = {1, 2}; //Search over X1, X2 
        Boolean[] wvd = {true, true, true}; //Treatment, X1, X2 all discrete
        variableSearchIndex = vsi;
        DiscreteVariables = wvd;
    }

    @Override
    public Matrix getBalancingVector() {
        return balancingVector;
    }

    public void SimpleRCTMomentSpecification() {
        // this.numObs = numObs;
    }

    public SimpleRCTMomentSpecification(Jama.Matrix X, Jama.Matrix Y, int numtrees, int[] variableSearchIndex, Boolean[] DiscreteVariables) {
        this.X = X;
        this.Y = Y;
        this.numtrees = numtrees;
        this.variableSearchIndex = variableSearchIndex;
        this.DiscreteVariables = DiscreteVariables;
        balancingVector = pmUtility.getColumn(X, 0); // treatment status is the first column of X
    }

    @Override
    public Double getPredictedY(Matrix xi, Jama.Matrix beta) {
        // pmUtility.prettyPrint(xi);
        if (beta != null) {
            double yhat = xi.get(0, 0) * beta.get(0, 0); // There is a single independent variable "treatment" in the simple RCT
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
        return new MomentContinuousSplitObjRCT(indexSplitVariable, lens, minProportionEachPartition, minCountEachPartition);
    }

    @Override
    public MomentPartitionObj getMomentPartitionObj(DataLens lens, int indexSplitVariable, IntegerPartition partition) {
        return new MomentPartitionObjRCT(partition, indexSplitVariable, lens);
    }

    @Override
    public ContainerMoment computeOptimalBeta(DataLens lens) {
        return new ContainerRCT(lens);
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
    public Matrix getXoriginal() {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    public Matrix cvparameters() {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    public int numberoftrees() {
        return numtrees;
    }

    @Override
    public Boolean[] getDiscreteVector() {
        return DiscreteVariables;
    }

    //Return the true treatment effect for a given observation
    @Override
    public Matrix getBetaTruth(Matrix xi) {
        Jama.Matrix beta = new Jama.Matrix(1, 1); // Beta is a scalar
        beta.set(0, 0, 0);
        beta.set(0, 0, xi.get(0, 1) + xi.get(0, 1) * (xi.get(0, 2) - 1)); //Treatment effect is x1+x1*(x2-1)
        return beta;
    }

    @Override
    public Matrix getOutOfSampleX() {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    public void loadData() {

        int numObsFile = 0;
        try {
            BufferedReader in = new BufferedReader(new FileReader("../Monte_Carlo/saturated_heterogeneity_25000.csv")); // Inputting data. What is the cd here?
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
            BufferedReader in = new BufferedReader(new FileReader("../Monte_Carlo/saturated_heterogeneity_25000.csv"));
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
    }

    @Override
    public String getVariableName(int variableIndex) {
        if (variableIndex == 0) {
            return "Treatment Indicator";
        }
        if (variableIndex == 1) {
            return "Outcome";
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
        double b = beta.get(0, 0);
        double se = Math.sqrt(variance.get(0, 0));
        String stars = "";
        NormalDistribution normal = new NormalDistribution(0, 1);
        if (Math.abs(b / se) > Math.abs(normal.inverse(0.90))) {
            stars = "*";
        }
        if (Math.abs(b / se) > Math.abs(normal.inverse(0.95))) {
            stars = "**";
        }
        if (Math.abs(b / se) > Math.abs(normal.inverse(0.99))) {
            stars = "***";
        }
        return String.format("%.2f (%.2f) %s", b, se, stars);
    }
}
