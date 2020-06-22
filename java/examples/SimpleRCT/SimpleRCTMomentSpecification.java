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
import core.IntegerPartition;
import core.MomentContinuousSplitObj;
import core.MomentPartitionObj;
import core.MomentSpecification;
import core.NaiveContainer;
import examples.RCT.ContainerRCT;
import examples.RCT.MomentContinuousSplitObjRCT;
import examples.RCT.MomentPartitionObjRCT;
import java.util.ArrayList;
import utility.pmUtility;

/**
 *
 * @author Stephen P. Ryan <stephen.p.ryan@wustl.edu>
 */
public class SimpleRCTMomentSpecification implements MomentSpecification {

    Jama.Matrix X;
    Jama.Matrix Y;

    int[] variableSearchIndex = {1};
    Boolean[] whichVariablesDiscrete = {true, true};

    public SimpleRCTMomentSpecification() {
    }

    @Override
    public Double getPredictedY(Matrix xi, Matrix beta) {
        // pmUtility.prettyPrint(xi);
        if (xi.get(0, 0) == 1) {
            return beta.get(0, 0);
        }
        return 0.0;
    }

    @Override
    public int[] getVariableIndicesToSearchOver() {
        return variableSearchIndex;
    }

    @Override
    public MomentContinuousSplitObj getFminObjective(Matrix nodeY, Matrix nodeX, int indexSplitVariable, double minProportionEachPartition, int minCountEachPartition) {
        return new MomentContinuousSplitObjRCT(indexSplitVariable, nodeX, nodeY, minProportionEachPartition, minCountEachPartition);
    }

    @Override
    public MomentPartitionObj getMomentPartitionObj(Matrix nodeX, Matrix nodeY, int indexSplitVariable, IntegerPartition partition) {
        return new MomentPartitionObjRCT(partition, indexSplitVariable, nodeX, nodeY);
    }

    @Override
    public ContainerMoment computeOptimalBeta(Matrix nodeY, Matrix nodeX) {
        return new ContainerRCT(nodeX, nodeY);
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
        System.out.println("\nIndependent Outcome Estimates");
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
            for (int i = 0; i < typeX.size(); i++) {
                for (int j = 0; j < 1; j++) {
                    subX.set(i, j, typeX.get(i).get(0, j));
                }
                subY.set(i, 0, typeY.get(i));
            }
            // pmUtility.prettyPrint(pmUtility.concatMatrix(subY, subX));
            Jama.Matrix olsBeta = pmUtility.OLS(subX, subY, false);
            Jama.Matrix olsVar = pmUtility.getOLSVariances(subY, subX, false);
            String sig = "";
            if (Math.abs(olsBeta.get(0, 0) / Math.sqrt(olsVar.get(0, 0))) > 1.98) {
                sig = "*";
            }
            System.out.format("Group %d: %g (%g) %s %n", k, olsBeta.get(0, 0), Math.sqrt(olsVar.get(0, 0)), sig);
        }

        // let's run the oracle estimator and see how that compares
        System.out.println("\nOracle Estimator");
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
                if(subX.get(i,0)==0) {
                    countControl++;
                } else {
                    countTreatment++;
                }
            }
            // pmUtility.prettyPrint(pmUtility.concatMatrix(subY, subX));
            Jama.Matrix olsBeta = pmUtility.OLS(subX, subY, false);
            Jama.Matrix olsVar = pmUtility.getOLSVariances(subY, subX, false);
            String sig = "";
            if (Math.abs(olsBeta.get(0, 0) / Math.sqrt(olsVar.get(0, 0))) > 1.98) {
                sig = "*";
            }
            System.out.format("Group %d: [%d, %d] %g (%g) %s %n", k, countControl, countTreatment, olsBeta.get(0, 0), Math.sqrt(olsVar.get(0, 0)), sig);
        }

        return null;
    }

    @Override
    public void loadData() {
        NormalDistribution normal = new NormalDistribution();
        int n = 2000;
        int G = 10; // number of groups
        Y = new Jama.Matrix(n, 1);
        X = new Jama.Matrix(n, 2);

        for (int i = 0; i < n; i++) {
            X.set(i, 0, Math.floor(2 * Math.random())); // treatment indicator
            X.set(i, 1, Math.floor(G * Math.random())); // group number
            if (X.get(i, 0) == 1) {
                if (X.get(i, 1) < 4) {
                    Y.set(i, 0, -1.0);
                } else if (X.get(i, 1) > 7) {
                    Y.set(i, 0, 1.0 + X.get(i, 1) * 0.1);
                }
            }
            Y.set(i, 0, Y.get(i, 0) + 0.1 * normal.inverse(Math.random()));
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
