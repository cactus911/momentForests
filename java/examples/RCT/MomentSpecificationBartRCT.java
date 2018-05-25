/*
 * The MIT License
 *
 * Copyright 2018 Stephen P. Ryan <stephen.p.ryan@wustl.edu>.
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
package examples.RCT;

import Jama.Matrix;
import java.util.Random;
import JSci.maths.statistics.NormalDistribution;
import core.ContainerMoment;
import core.IntegerPartition;
import core.MomentContinuousSplitObj;
import core.MomentPartitionObj;
import core.MomentSpecification;
import core.NaiveContainer;
import java.io.BufferedReader;
import java.io.FileReader;
import java.util.TreeSet;

/**
 *
 * @author Stephen P. Ryan <stephen.p.ryan@wustl.edu>
 */
public class MomentSpecificationBartRCT implements MomentSpecification {

    boolean monteCarlo = false;
    
    Jama.Matrix X;
    Jama.Matrix Y;
    static NormalDistribution normal = new NormalDistribution();
    int[] searchArray = {1, 2, 3, 4, 5, 6, 7, 8, 9};
    // disceteIndicator includes the first column, which is the treatment effect?
    // okay, the way this works is that the indices above select variables from this list next
    // so this list +is+ defined across everything in the data set, but searchArray above controls
    // which of those X's are actually used when splitting the sample
    Boolean[] discreteIndicator = {true, false, true, true, true, true, true, true, true, false};
    final private String[] variableNames = {"treatment", "age", "symptomatic", "homosexual", "IVDrugUser", "PriorAZT", "Male", "White", "FutureDropout", "Baseline"};

    /**
     * Try the simplest possible thing here. Y = \tau W + \epsilon, where W is
     * randomly assigned
     *
     * So X is a vector of zeros and ones FIRST COLUMN: This is like W---the
     * treatment indicator We will skip that when evaluating splits add one
     * column to X as a set of group indicators
     */
    public MomentSpecificationBartRCT() {
    }

    @Override
    public Double getPredictedY(Matrix xi, Jama.Matrix beta) {
        // constant + beta*1(W=1)
//        System.out.println("****");
//        pmUtility.prettyPrintVector(beta);
//        pmUtility.prettyPrint(xi);
        if (beta != null) {
            double yhat = xi.get(0, 0) * beta.get(0, 0);
            return yhat;
        }
        return null;
    }

    @Override
    public MomentContinuousSplitObj getFminObjective(Matrix nodeY, Matrix nodeX, int indexSplitVariable, double minProportionEachPartition, int minCountEachPartition) {
        return new MomentContinuousSplitObjBart(indexSplitVariable, nodeX, nodeY, minProportionEachPartition, minCountEachPartition);
    }

    @Override
    public int[] getVariableIndicesToSearchOver() {
        return searchArray;
    }

    @Override
    public MomentPartitionObj getMomentPartitionObj(Matrix nodeX, Matrix nodeY, int indexSplitVariable, IntegerPartition partition) {
        return new MomentPartitionObjBartRCT(partition, indexSplitVariable, nodeX, nodeY);
    }

    @Override
    public ContainerMoment computeOptimalBeta(Matrix nodeY, Matrix nodeX) {
        return new ContainerRCT(nodeX, nodeY);
    }

    @Override
    public void generateData(int numObs, Random rng, boolean addNoise) {
        X = new Jama.Matrix(numObs, discreteIndicator.length);
        Y = new Jama.Matrix(numObs, 1);

        double meanTreatment = 0;
        double meanControl = 0;
        int countTreatment = 0;
        int countControl = 0;

        for (int i = 0; i < numObs; i++) {
            double draw = rng.nextDouble();
            if (draw < 0.5) {
                X.set(i, 0, 1);
            }
            for (int j = 0; j < discreteIndicator.length - 1; j++) {
                if (discreteIndicator[j + 1]) {
                    X.set(i, j + 1, 1 + Math.floor(rng.nextDouble() * 2.0));
                } else {
                    X.set(i, j + 1, rng.nextDouble() * 2.0);
                }
            }

            if (addNoise) {
                Y.set(i, 0, Y.get(i, 0) + 0 * normal.inverse(rng.nextDouble()));
            }

            Jama.Matrix obsBeta = getBetaTruth(X.getMatrix(i, i, 0, X.getColumnDimension() - 1));

            if (X.get(i, 0) == 1) {
                // treatment
//                Y.set(i, 0, Y.get(i, 0) + 1);
//                if (X.get(i, 2) == 2) {
//                    Y.set(i, 0, Y.get(i, 0) + 1);
//                }
//                if (X.get(i, 3) == 1) {
//                    Y.set(i, 0, Y.get(i, 0) + 2);
//                }
//                if (X.get(i, 4) == 1) {
//                    Y.set(i, 0, Y.get(i, 0) + 3);
//                }
                if (X.get(i, 1) > 1.0) {
                    Y.set(i, 0, Y.get(i, 0) + 1);
                }
                if (X.get(i, 1) < 0.50) {
                    Y.set(i, 0, Y.get(i, 0) + 2);
                }
                if (X.get(i, 1) > 0.5 && X.get(i, 2) == 1) {
                    Y.set(i, 0, Y.get(i, 0) - 3);
                }
                meanTreatment += Y.get(i, 0);
                countTreatment++;
            } else {
                meanControl += Y.get(i, 0);
                countControl++;
            }
        }
        meanTreatment /= countTreatment;
        meanControl /= countControl;
        System.out.println("meanControl: " + meanControl + " meanTreatment: " + meanTreatment + " diff: " + (meanTreatment - meanControl));
        // System.exit(0);
//        pmUtility.prettyPrint(pmUtility.concatMatrix(Y, X));
//        System.exit(0);
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
    public Boolean[] getDiscreteVector() {
        return discreteIndicator;
    }

    @Override
    public Matrix getBetaTruth(Matrix xi) {
        Jama.Matrix beta = new Jama.Matrix(1, 1);
        beta.set(0, 0, 0);
        /**
         * Add in some observable heterogeneity. I think this should work.
         */
        // if (xi.get(0, 1) == 1 && xi.get(0,2)==4) {
//        if (xi.get(0, 1) == 1) {
//            beta.set(0, 0, 20);
//        }
//        if (xi.get(0, 1) == 4 && xi.get(0,2)==3) {
//            beta.set(0, 0, -10);
//        }
        beta.set(0, 0, xi.get(0, 1) + xi.get(0, 1) * (xi.get(0, 2) - 1));
        return beta;
    }

    @Override
    public void loadData() {

        try {
            BufferedReader inCounter = new BufferedReader(new FileReader("ACT175_week56_2arms.csv"));
            String line = inCounter.readLine(); // headers
            int counter = 0;
            TreeSet<Integer> exclusionTree = new TreeSet<>();
            while (line != null) {
                if (line != null) {
                    try {
                        line = inCounter.readLine();
                        int a = 0;
                        for (int j = 0; j < 10; j++) {
                            int b = line.indexOf(",", a);
                            Double.parseDouble(line.substring(a, b));
                            a = b + 1;
                        }
                        Double.parseDouble(line.substring(a));
                    } catch (Exception e) {
                        System.out.println(e);
                        System.out.println("Missing data for observation " + counter);
                        exclusionTree.add(counter);
                    }
                    counter++;
                }
            }
            inCounter.close();

            X = new Jama.Matrix(counter - exclusionTree.size(), 10);
            Y = new Jama.Matrix(counter - exclusionTree.size(), 1);

            BufferedReader in = new BufferedReader(new FileReader("ACT175_week56_2arms.csv"));
            line = in.readLine(); // headers
            int counter2 = 0;
            for (int i = 0; i < counter; i++) {
                line = in.readLine();
                // System.out.println(line);
                if (!exclusionTree.contains(i)) {
                    int a = 0;
                    int b = line.indexOf(",", a);
                    // System.out.println(line.substring(a, b));
                    Y.set(counter2, 0, Double.parseDouble(line.substring(a, b)));
                    a = b + 1;
                    b = line.indexOf(",", a);
                    // System.out.println(line.substring(a, b));
                    X.set(counter2, 0, Double.parseDouble(line.substring(a, b)));
                    for (int j = 0; j < 9; j++) {
                        a = b + 1;
                        b = line.indexOf(",", a);
                        // System.out.println(line.substring(a, b));
                        if (j < 8) {
                            X.set(counter2, j + 1, Double.parseDouble(line.substring(a, b)));
                        } else {
                            X.set(counter2, j + 1, Double.parseDouble(line.substring(a)));
                        }
                    }
                    counter2++;
                }
            }
            System.out.println("Loaded " + counter2 + " observations after dropping " + exclusionTree.size() + " observations for missing data.");
            in.close();
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(0);
        }

        if (monteCarlo) {
            // generateData(X.getRowDimension(), new Random(), true);
            generateData(15000, new Random(), true);
        }
    }

    @Override
    public Matrix getOutOfSampleX() {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    public NaiveContainer computeNaiveStatistics() {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    public String getVariableName(int variableIndex) {
        return variableNames[variableIndex];
    }

    @Override
    public String getFixedEffectName(int variableIndex, int fixedEffectIndex) {
        // final private String[] variableNames = {"treatment", "age", "symptomatic", "homosexual", "IVDrugUser", "PriorAZT", "Male", "White", "FutureDropout", "Baseline"};
        String[][] feNames = {
            {"Control", "Treatment"}, // x0
            {"Age (continuous)"}, // x1
            {"Not Symptomatic", "Symptomatic"}, // x2
            {"Homosexual", "Heterosexual"}, // x3
            {"Not IV Drug User", "IV Drug User"}, // x4
            {"No Prior AZT", "Prior AZT"}, // x5
            {"Female", "Male"}, // x6
            {"Not White", "White"}, // x7
            {"Not Dropout", "Future Dropout"}, // x8
            {"Baseline (not coded)"} // x9
        };
        // System.out.println(variableIndex+" "+fixedEffectIndex);
        if(monteCarlo) {
            return "x[" + variableIndex + "]=" + fixedEffectIndex;
        }
        return feNames[variableIndex][fixedEffectIndex];
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
