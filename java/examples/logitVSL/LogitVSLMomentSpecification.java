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
package examples.logitVSL;

import JSci.maths.statistics.NormalDistribution;
import Jama.Matrix;
import core.ContainerMoment;
import core.DataLens;
import core.MomentSpecification;
import java.io.FileReader;
import java.io.BufferedReader;
import java.util.Random;

/**
 *
 * @author Stephen P. Ryan <stephen.p.ryan@wustl.edu>
 */
public class LogitVSLMomentSpecification extends MomentSpecification {

    Jama.Matrix X;
    Jama.Matrix Y;
    Jama.Matrix Z;

    int numObs;
    int numtrees;
    int[] variableSearchIndex; // this should be restricted to only Z
    Boolean[] DiscreteVariables; // also this should be restricted to only Z
    String filename;
    boolean MONTE_CARLO = true;
    NormalDistribution normal = new NormalDistribution();

    String[] varName = {"dep_prob", "afqsc", "entryyr", "term", "mos", "female", "race", "entryage", "edu", "state", "urate"};

    public LogitVSLMomentSpecification(String filename) {
        this.filename = filename;

        initializeHomogeneity(3);

        int[] vsi = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10}; //Search over z1, z2, z3
        // Boolean[] wvd = {false, false, false, true, true, true, true, true, true, true, true}; // z1, z2 continuous, z3 discrete
        Boolean[] wvd = {false, false, false, false, false, true, true, false, false, false, false}; // z1, z2 continuous, z3 discrete
        variableSearchIndex = vsi;
        DiscreteVariables = wvd;
    }

    @Override
    public Double getPredictedY(Matrix xi, Jama.Matrix beta, Random rng) {
        return LogitVSLDataGenerator.getLogitDiscreteOutcome(xi, beta, rng);
    }

    @Override
    public int[] getVariableIndicesToSearchOver() {
        return variableSearchIndex;
    }

    @Override
    public ContainerMoment computeOptimalBeta(DataLens lens, boolean allParametersHomogeneous) {
        ContainerLogitVSL l = new ContainerLogitVSL(lens, getHomogeneousIndex(), getHomogeneousParameterVector(), allParametersHomogeneous);
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
    public double getGoodnessOfFit(double yi, Matrix xi, Matrix beta) {
        // here, let's return the LLH of this observation
        return ContainerLogitVSL.computeLLHi(yi, xi, beta);
    }

    @Override
    public int getNumParams() {
        return 3;
    }

    @Override
    public int getNumMoments() {
        return 3;
    }

    @Override
    public void loadData(long rngSeed) {

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
        X = new Jama.Matrix(numObsFile, 3); // Used previous loop to create arrays of the correct size, memory saver?
        Z = new Jama.Matrix(numObsFile, 11);
        Y = new Jama.Matrix(numObsFile, 1);

        try {
            BufferedReader in = new BufferedReader(new FileReader(filename));
            String line = in.readLine(); // headers
            int i = 0;
            while (line != null) {
                line = in.readLine();
                if (line != null) {
                    int a = 0;
                    int b = line.indexOf(",", a);
                    double reupind = Double.parseDouble(line.substring(a, b));
                    Y.set(i, 0, reupind);
                    a = b + 1;

                    // constant
                    X.set(i, 0, 1.0);

                    // srb
                    b = line.indexOf(",", a);
                    double srb = Double.parseDouble(line.substring(a, b));
                    X.set(i, 1, srb);

                    // hazard
                    a = b + 1;
                    b = line.indexOf(",", a);
                    double hazard = Double.parseDouble(line.substring(a, b));
                    X.set(i, 2, hazard);

//                dep_prob;
                    a = b + 1;
                    b = line.indexOf(",", a);
                    Z.set(i, 0, Double.parseDouble(line.substring(a, b)));

                    // System.out.println(X.get(i,0)+" "+X.get(i,1)+" "+X.get(i,2));
//                entry_afqsc;
                    a = b + 1;
                    b = line.indexOf(",", a);
                    Z.set(i, 1, Double.parseDouble(line.substring(a, b)));

//                entryyr;
                    a = b + 1;
                    b = line.indexOf(",", a);
                    Z.set(i, 2, Double.parseDouble(line.substring(a, b)));

//                term;
                    a = b + 1;
                    b = line.indexOf(",", a);
                    Z.set(i, 3, Double.parseDouble(line.substring(a, b)));

//                mos;
                    a = b + 1;
                    b = line.indexOf(",", a);
                    Z.set(i, 4, Double.parseDouble(line.substring(a, b)));

//                female;
                    a = b + 1;
                    b = line.indexOf(",", a);
                    Z.set(i, 5, Double.parseDouble(line.substring(a, b)));

//                race;
                    a = b + 1;
                    b = line.indexOf(",", a);
                    Z.set(i, 6, Double.parseDouble(line.substring(a, b)));

//                entryage;
                    a = b + 1;
                    b = line.indexOf(",", a);
                    Z.set(i, 7, Double.parseDouble(line.substring(a, b)));

//                edu;
                    a = b + 1;
                    b = line.indexOf(",", a);
                    Z.set(i, 8, Double.parseDouble(line.substring(a, b)));

//                state;
                    a = b + 1;
                    b = line.indexOf(",", a);
                    Z.set(i, 9, Double.parseDouble(line.substring(a, b)));

//                urate;
                    a = b + 1;
                    b = line.indexOf(",", a);
                    Z.set(i, 10, Double.parseDouble(line.substring(a, b)));
                }

                i++;
            }

            in.close();
        } catch (Exception e) {
            e.printStackTrace();
        }

        boolean TEST_MC = true;
        if (TEST_MC) {
            Random mcRNG = new Random();
            for (int i = 0; i < numObs; i++) {
                double error = -Math.log(-Math.log(mcRNG.nextDouble()));
                double error0 = -Math.log(-Math.log(mcRNG.nextDouble()));
                double alpha = -0.5;
                double betaBonus = 1.0; // 0.0147;
                double betaHazard = -1.0; // -0.0035;
                if(Z.get(i,5)==1.0) {
                    alpha = 0.5;
                    // betaHazard = -2.0;
                }

                double utility = alpha + betaBonus * X.get(i, 1) + betaHazard * X.get(i, 2) + error;
                if (utility > error0) {
                    Y.set(i, 0, 1.0);
                } else {
                    Y.set(i, 0, 0.0);
                }
            }
        }

    }

    @Override
    public String getVariableName(int variableIndex) {
        return varName[variableIndex];
    }

}
