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

import JSci.maths.statistics.NormalDistribution;
import Jama.Matrix;
import core.ContainerMoment;
import core.DataLens;
import core.IntegerPartition;
import core.MomentContinuousSplitObj;
import core.MomentPartitionObj;
import core.MomentSpecification;
import core.NaiveContainer;
import examples.RCT.ContainerRCT;
import examples.RCT.MomentContinuousSplitObjRCT;
import examples.RCT.MomentPartitionObjRCT;
import java.io.BufferedReader;
import java.io.FileReader;
import java.util.ArrayList;
import utility.pmUtility;

/**
 *
 * @author Stephen P. Ryan <stephen.p.ryan@wustl.edu>
 */
public class SimpleRCTMomentSpecification implements MomentSpecification {

    private final boolean USE_STAR_DATA = false;
    private final boolean USE_CHARLIE_TEST_1 = true;
    private final boolean SPLIT_OBSERVABLES = false;

    int G = 10; // number of groups
    Jama.Matrix X;
    Jama.Matrix Y;
    int numObs;

    int[] variableSearchIndex = {1};
    Boolean[] whichVariablesDiscrete = {true, true};

    public SimpleRCTMomentSpecification(int numObs) {
        this.numObs = numObs;

        if (USE_STAR_DATA) {
            int[] vsi = {6}; // 
            if (SPLIT_OBSERVABLES) {
                int[] vsi2 = {3, 4, 5, 6}; // 6 is the outcome
                vsi = vsi2;
            }
            Boolean[] wvd = {true, true, true, true, true, true, true};
            variableSearchIndex = vsi;
            whichVariablesDiscrete = wvd;
        }
        if (USE_CHARLIE_TEST_1) {
            int[] vsi = {1}; // 
            Boolean[] wvd = {true, false};
            variableSearchIndex = vsi;
            whichVariablesDiscrete = wvd;
        }
    }

    @Override
    public Double getPredictedY(Matrix xi, Matrix beta) {
        // pmUtility.prettyPrint(xi);
        if (xi.get(0, getIndexBalancingVector()) == 1) {
            return beta.get(0, 0);
        }
        return 0.0;
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
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    public Boolean[] getDiscreteVector() {
        return whichVariablesDiscrete;
    }

    @Override
    public Matrix getBetaTruth(Matrix xi) {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    public Matrix getOutOfSampleX() {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    public NaiveContainer computeNaiveStatistics() {
        // this is to run separate OLS regressions for each treatment
        boolean runIndependentOLS = false;
        if (runIndependentOLS) {
            System.out.println("\nIndependent Outcome Estimates");
            for (int k = 0; k < G; k++) {
                ArrayList<Jama.Matrix> typeX = new ArrayList<>();
                ArrayList<Double> typeY = new ArrayList<>();
                for (int i = 0; i < X.getRowDimension(); i++) {
                    int indexOfGroupIndicator = 1;
                    int adjustment = 0;
                    if(USE_STAR_DATA) {
                        indexOfGroupIndicator = 6;
                        adjustment = 1;
                    }
                    if (X.get(i, indexOfGroupIndicator) == k + adjustment) {
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
                // Jama.Matrix[] bootOLS = pmUtility.bootstrapOLS(subX, subY, true, 500, 787);

//                pmUtility.prettyPrint(pmUtility.concatMatrix(subY, subX));
                
                Jama.Matrix olsBeta = pmUtility.OLSsvd(subX, subY, true);
                Jama.Matrix olsVar = pmUtility.getOLSVariances(subY, subX, true);
                
//                System.out.println("beta:");
//                pmUtility.prettyPrint(olsBeta);
//                System.out.println("var:");
//                pmUtility.prettyPrint(olsVar);
                
                String sig = "";
                if (Math.abs(olsBeta.get(1, 0) / Math.sqrt(olsVar.get(1, 0))) > 1.98) {
                    sig = "*";
                }
                System.out.format("OLS Formula  Group %d: %g (%g) %s %n", k, olsBeta.get(1, 0), Math.sqrt(olsVar.get(1, 0)), sig);
                // System.out.format("Bootstrapped Group %d: [%d, %d] %g (%g) %s %n", k, countControl, countTreatment, bootOLS[0].get(1, 0), bootOLS[1].get(1, 0), sig);
            }
        }

        boolean computeOracle = true;
        if (computeOracle && !USE_STAR_DATA) {
            // let's run the oracle estimator and see how that compares
            System.out.println("\nOracle Estimator");
            for (int k = 0; k < 1; k++) {
                ArrayList<Jama.Matrix> typeX = new ArrayList<>();
                ArrayList<Double> typeY = new ArrayList<>();
                for (int i = 0; i < X.getRowDimension(); i++) {
                    // if (k == 0 && X.get(i, 1) < 4) {
                    if (k == 0) {
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
                boolean useIntercept = true;
                // Jama.Matrix[] bootOLS = pmUtility.bootstrapOLS(subX, subY, useIntercept, 500, 787);

                // System.out.format("Group %d: [%d, %d] %g (%g) %s %n", k, countControl, countTreatment, olsBeta.get(0, 0), Math.sqrt(olsVar.get(0, 0)), sig);
                // System.out.format("Bootstrapped Group %d: [%d, %d] %g (%g) %n", k, countControl, countTreatment, bootOLS[0].get(1, 0), bootOLS[1].get(1, 0));
            }
            System.out.println("Single treatment effect theoretical SE (if imposing intercept is zero): " + Math.sqrt(2.0 / (numObs - 1)));
        }

        return null;
    }

    public String getSpecName() {
        String s = "_monteCarlo";
        if (USE_STAR_DATA) {
            s = "_star_splitOutcomes";
            if (isSplittingOnObservables()) {
                s = "_star_splitAllCovariates";
            }
        }
        if(USE_CHARLIE_TEST_1) {
            s = "_charlie1";
        }
        return s;
    }

    @Override
    public void loadData() {
        if (USE_STAR_DATA) {
            int numObsFile = 0;
            try {
                BufferedReader in = new BufferedReader(new FileReader("../star_grade3_scaled.csv"));
                String line = in.readLine(); // headers
                while (line != null) {
                    line = in.readLine();
                    if (line != null) {
                        numObsFile++;
                    }
                }
                in.close();
                numObs = numObsFile;
            } catch (Exception e) {
                e.printStackTrace();
            }
            System.out.format("Number of observations = %,d %n", numObsFile);
            Jama.Matrix dX = new Jama.Matrix(numObsFile, 7);
            Jama.Matrix dY = new Jama.Matrix(numObsFile, 1);

            try {
                BufferedReader in = new BufferedReader(new FileReader("../star_grade3_scaled.csv"));
                String line = in.readLine(); // headers
                int i = 0;
                while (line != null) {
                    line = in.readLine();
                    if (line != null) {
                        int a = 0;
                        int b = line.indexOf(",", a);
                        // student id

                        a = b + 1;
                        b = line.indexOf(",", a);
                        // grade

                        a = b + 1;
                        b = line.indexOf(",", a);
                        dY.set(i, 0, Double.valueOf(line.substring(a, b))); // outcome

                        a = b + 1;
                        b = line.indexOf(",", a);
                        dX.set(i, 0, Double.valueOf(line.substring(a, b))); // treatment status

                        a = b + 1;
                        b = line.indexOf(",", a);
                        dX.set(i, 1, Double.valueOf(line.substring(a, b))); // year star (?)

                        a = b + 1;
                        b = line.indexOf(",", a);
                        dX.set(i, 2, Double.valueOf(line.substring(a, b))); // year small (?)

                        a = b + 1;
                        b = line.indexOf(",", a);
                        dX.set(i, 3, Double.valueOf(line.substring(a, b))); // is white

                        a = b + 1;
                        b = line.indexOf(",", a);
                        dX.set(i, 4, Double.valueOf(line.substring(a, b))); // is free lunch

                        a = b + 1;
                        b = line.indexOf(",", a);
                        dX.set(i, 5, Double.valueOf(line.substring(a, b))); // is female

                        a = b + 1;
                        b = line.indexOf(",", a);
                        // dX.set(i, 5, Double.valueOf(line.substring(a, b))); // school wave FE

                        a = b + 1;
                        dX.set(i, 6, Double.valueOf(line.substring(a))); // outcome indicator (runs from 1 to 10)

                        i++;
                    }
                }
                X = dX;
                Y = dY;
                in.close();
            } catch (Exception e) {
                e.printStackTrace();
            }
        } else if (USE_CHARLIE_TEST_1) {
            int numObsFile = 0;
            try {
                BufferedReader in = new BufferedReader(new FileReader("../1_1default_continuous_backup.csv"));
                String line = in.readLine(); // headers
                while (line != null) {
                    line = in.readLine();
                    if (line != null) {
                        numObsFile++;
                    }
                }
                in.close();
                numObs = numObsFile;
            } catch (Exception e) {
                e.printStackTrace();
            }
            
            System.out.format("Number of observations = %,d %n", numObsFile);
            Jama.Matrix dX = new Jama.Matrix(numObsFile, 2);
            Jama.Matrix dY = new Jama.Matrix(numObsFile, 1);

            try {
                BufferedReader in = new BufferedReader(new FileReader("../1_1default_continuous_backup.csv"));
                String line = in.readLine(); // headers
                int i = 0;
                while (line != null) {
                    line = in.readLine();
                    if (line != null) {
                        int a = 0;
                        int b = line.indexOf(",", a);
                        // System.out.println(line);
                        dY.set(i, 0, Double.valueOf(line.substring(a, b))); // outcome

                        a = b + 1;
                        b = line.indexOf(",", a);
                        dX.set(i, 0, Double.valueOf(line.substring(a, b))); // treatment status

                        a = b + 1;
                        b = line.indexOf(",", a);
                        dX.set(i, 1, Double.valueOf(line.substring(a))); // continuous X1

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
            NormalDistribution normal = new NormalDistribution();
            int n = numObs;
            
            Y = new Jama.Matrix(n, 1);
            X = new Jama.Matrix(n, 2);

            for (int i = 0; i < n; i++) {
                X.set(i, 0, Math.floor(2 * Math.random())); // treatment indicator
                X.set(i, 1, Math.floor(G * Math.random())); // group number
                if (X.get(i, 0) == 1) { // treatment
                    // if (X.get(i, 1) < 2) {
                    Y.set(i, 0, 0.0);
//                    } else if (X.get(i, 1) == 8) {
//                        Y.set(i, 0, 10.0);
//                    } else if (X.get(i, 1) == 9) {
//                        Y.set(i, 0, -4.0);
//                    } else {
//                        Y.set(i, 0, 0.0);
//                    }
                }
                Y.set(i, 0, Y.get(i, 0) + 1.0 * normal.inverse(Math.random()));
            }
        }
    }

    public int getIndexBalancingVector() {
        return 0;
    }

    public boolean usingSTARData() {
        return USE_STAR_DATA;
    }

    @Override
    public String getVariableName(int variableIndex) {
        if (USE_STAR_DATA) {
            if (variableIndex == 3) {
                return "White";
            }
            if (variableIndex == 4) {
                return "Free Lunch";
            }
            if (variableIndex == 5) {
                return "Female";
            }
            if (variableIndex == 6) {
                return "Outcome";
            }
        }
        if(USE_CHARLIE_TEST_1) {
            if(variableIndex==0) {
                return "x1";
            }
        }
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

    /**
     * @return the SPLIT_OBSERVABLES
     */
    public boolean isSplittingOnObservables() {
        return SPLIT_OBSERVABLES;
    }

}
