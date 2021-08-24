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
import Jama.Matrix;
import core.ContainerMoment;
import core.DataLens;
import core.IntegerPartition;
import core.MomentContinuousSplitObj;
import core.MomentPartitionObj;
import core.MomentSpecification;
import core.NaiveContainer;
import java.util.ArrayList;
import java.io.FileReader;
import java.io.BufferedReader;
import utility.pmUtility;

/**
 *
 * @author Stephen P. Ryan <stephen.p.ryan@wustl.edu>
 */
public class SimpleRCTMomentSpecification implements MomentSpecification {

    Jama.Matrix X;
    Jama.Matrix Y;
    int numObs;
    int numtrees;
    int[] variableSearchIndex;
    Boolean[] DiscreteVariables;
    
    /*Begin new*/
    public SimpleRCTMomentSpecification(int numObs) {
        this.numObs = numObs;
        int[] vsi = {1, 2}; //Search over X1, X2 
        Boolean[] wvd = {true, true, true}; //Treatment, X1, X2 all discrete
        variableSearchIndex = vsi;
        DiscreteVariables = wvd;
    }
   /*End new*/
    
    public void SimpleRCTMomentSpecification() {
        // this.numObs = numObs;
    }
    
    public SimpleRCTMomentSpecification(Jama.Matrix X, Jama.Matrix Y, int numtrees, int[] variableSearchIndex, Boolean[] DiscreteVariables) {      
        this.X = X;
        this.Y = Y;
        this.numtrees = numtrees;
        this.variableSearchIndex = variableSearchIndex;
        this.DiscreteVariables = DiscreteVariables;
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
    public NaiveContainer computeNaiveStatistics() {
        // this is to run separate OLS regressions for each treatment
        System.out.println("\nIndependent Outcome Estimates"); // These estimates are specific to the simple RCT context
        for (int k = 0; k < 10; k++) {
            ArrayList<Jama.Matrix> typeX = new ArrayList<>();
            ArrayList<Double> typeY = new ArrayList<>();
            for (int i = 0; i < X.getRowDimension(); i++) {
                if (X.get(i, 1) == k) {
                    typeX.add(X.getMatrix(i, i, 0, X.getColumnDimension() - 1));
                    typeY.add(Y.get(i, 0));
                }
            }
            Jama.Matrix subX = new Jama.Matrix(typeX.size(), 1);
            Jama.Matrix subY = new Jama.Matrix(typeX.size(), 1);
            int countControl = 0;
                int countTreatment = 0;
            for (int i = 0; i < typeX.size(); i++) {
                for (int j = 0; j < 1; j++) {
                    subX.set(i, j, typeX.get(i).get(0, j));
                }
                subY.set(i, 0, typeY.get(i));
                if (subX.get(i, 0) == 0) {
                        countControl++;
                    } else {
                        countTreatment++;
                    }
            }
            // pmUtility.prettyPrint(pmUtility.concatMatrix(subY, subX));
            // Jama.Matrix[] bootOLS = pmUtility.bootstrapOLS(subX, subY, false, 500, 787);

            Jama.Matrix olsBeta = pmUtility.OLSsvd(subX, subY, false);
            Jama.Matrix olsVar = pmUtility.getOLSVariances(subY, subX, false);
            String sig = "";
            if (Math.abs(olsBeta.get(0, 0) / Math.sqrt(olsVar.get(0, 0))) > 1.98) {
                sig = "*";
            }
            System.out.format("OLS Formula  Group %d: %g (%g) %s %n", k, olsBeta.get(0, 0), Math.sqrt(olsVar.get(0, 0)), sig);
            System.out.format("Bootstrapped Group %d: [%d, %d] %g (%g) %s %n", k, countControl, countTreatment, sig); // bootOLS[0].get(0, 0), bootOLS[1].get(0, 0),
        }
        
        boolean computeOracle = true;
        if (computeOracle) { // This oracle estimator is also specific to the simple RCT context
            // let's run the oracle estimator and see how that compares
            // System.out.println("\nOracle Estimator");
            for (int k = 0; k < 4; k++) {
                ArrayList<Jama.Matrix> typeX = new ArrayList<>();
                ArrayList<Double> typeY = new ArrayList<>();
                for (int i = 0; i < X.getRowDimension(); i++) {
                    if (k == 0 && X.get(i, 1) < 4) {
                        typeX.add(X.getMatrix(i, i, 0, X.getColumnDimension() - 1));
                        typeY.add(Y.get(i, 0));
                    }
                    if (k == 1 && X.get(i, 1) >= 4 && X.get(i, 1) < 8) {
                        typeX.add(X.getMatrix(i, i, 0, X.getColumnDimension() - 1));
                        typeY.add(Y.get(i, 0));
                    }
                    if (k == 2 && X.get(i, 1) == 8) {
                        typeX.add(X.getMatrix(i, i, 0, X.getColumnDimension() - 1));
                        typeY.add(Y.get(i, 0));
                    }
                    if (k == 3 && X.get(i, 1) == 9) {
                        typeX.add(X.getMatrix(i, i, 0, X.getColumnDimension() - 1));
                        typeY.add(Y.get(i, 0));
                    }
                }
                Jama.Matrix subX = new Jama.Matrix(typeX.size(), 1);
                Jama.Matrix subY = new Jama.Matrix(typeX.size(), 1);
                int countControl = 0;
                int countTreatment = 0;
                for (int i = 0; i < typeX.size(); i++) {
                    for (int j = 0; j < 1; j++) {
                        subX.set(i, j, typeX.get(i).get(0, j));
                    }
                    subY.set(i, 0, typeY.get(i));
                    if (subX.get(i, 0) == 0) {
                        countControl++;
                    } else {
                        countTreatment++;
                    }
                }
                // pmUtility.prettyPrint(pmUtility.concatMatrix(subY, subX));
                Jama.Matrix olsBeta = pmUtility.OLS(subX, subY, false);
                Jama.Matrix olsVar = pmUtility.getOLSVariances(subY, subX, false);
                // Jama.Matrix[] bootOLS = pmUtility.bootstrapOLS(subX, subY, false, 5000, 787);

                String sig = "";
                if (Math.abs(olsBeta.get(0, 0) / Math.sqrt(olsVar.get(0, 0))) > 1.98) {
                    sig = "*";
                }
                System.out.format("Group %d: [%d, %d] %g (%g) %s %n", k, countControl, countTreatment, olsBeta.get(0, 0), Math.sqrt(olsVar.get(0, 0)), sig);
                System.out.format("Bootstrapped Group %d: [%d, %d] %g (%g) %s %n", k, countControl, countTreatment, sig); // bootOLS[0].get(0, 0), bootOLS[1].get(0, 0),
            }
        }

        return null;
    }
    /*
    @Override
    public void loadData() {
	throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }
    */
    
    @Override
    public void loadData() {
        
        int numObsFile = 0;
        try {
            BufferedReader in = new BufferedReader(new FileReader("C:/Users/Spare/Dropbox/MomentForests/SaturatedHeterogeneityRCT/saturated_heterogeneity_25000.csv")); // Inputting data. What is the cd here?
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

        //System.out.format("Number of observations = %,d %n", numObsFile);
        Jama.Matrix dX = new Jama.Matrix(numObsFile, 3); // Used previous loop to create arrays of the correct size, memory saver?
        Jama.Matrix dY = new Jama.Matrix(numObsFile, 1);

        try {
            BufferedReader in = new BufferedReader(new FileReader("C:/Users/Spare/Dropbox/MomentForests/SaturatedHeterogeneityRCT/saturated_heterogeneity_25000.csv"));
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
                    dX.set(i, 1, Double.valueOf(line.substring(a,b))); //X1

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

    /*
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
    */
}
