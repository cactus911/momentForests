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
package examples.linearIV;

import JSci.maths.statistics.NormalDistribution;
import Jama.Matrix;
import core.ContainerMoment;
import core.DataLens;
import core.MomentSpecification;

import java.io.FileReader;
import java.io.BufferedReader;
import java.util.Random;
import utility.pmUtility;

/**
 * IV (Instrumental Variables) moment specification.
 *
 * Key differences from LinearMomentSpecification (OLS):
 * - X matrix structure: [endogenous_var, exogenous_vars..., instrument]
 * - First column is the endogenous regressor
 * - Last column is the instrument (excluded from parameter estimation)
 * - Moment conditions: E[Z'e] = 0 where Z = instruments (columns 1..end of X)
 * - getNumParams() = X.columns - 1 (instrument is not a parameter)
 * - getNumMoments() = X.columns - 1 (endogenous var excluded from instruments)
 *
 * @author Stephen P. Ryan <stephen.p.ryan@wustl.edu>
 */
public class LinearIVMomentSpecification extends MomentSpecification {

    Jama.Matrix X;
    Jama.Matrix Y;
    Jama.Matrix Z;
    int numObs;
    int[] variableSearchIndex;
    Boolean[] DiscreteVariables;
    String filename;
    boolean failedEstimator = false;

    String[] varNames = {"constant", "ed76", "exp76", "exp762", "black", "south76", "city76", "region_1966", "smsa66r", "daded", "momed", "nodaded", "nomomed", "famed", "momdad14", "sinmom14", "near4"};

    public LinearIVMomentSpecification(String filename) {
        this.filename = filename;

        loadData(787);

        initializeHomogeneity(getNumParams());

        int[] vsi = {1, 2, 4, 5, 6, 7, 8, 9, 10, 11, 12, 14, 15};
        Boolean[] whichVariablesAreDiscrete = {false,
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
    public int getNumMoments() {
        return X.getColumnDimension() - 1; // one endo X (first column), one instrument (last column of "X")
    }

    @Override
    public double getGoodnessOfFit(double yi, Matrix xi, Matrix beta) {
        double error = yi - ((xi.getMatrix(0, 0, 0, xi.getColumnDimension() - 2)).times(beta)).get(0, 0);
        return error * error;
    }

    @Override
    public Double getPredictedY(Matrix xi, Jama.Matrix beta, Random rng) {
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
    public ContainerMoment computeOptimalBeta(DataLens lens, boolean allParametersHomogeneous) {
        ContainerLinearIV l = new ContainerLinearIV(lens, getHomogeneousIndex(), getHomogeneousParameterVector(), allParametersHomogeneous, this);
        l.computeBetaAndErrors();
        failedEstimator = l.didEstimatorFail();
        return l;
    }

    @Override
    public boolean didEstimatorFail() {
        return failedEstimator;
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
    public Boolean[] getDiscreteVector() {
        return DiscreteVariables;
    }

    @Override
    public void loadData(long rngSeed) {

        int numObsFile = 0;
        try {
            BufferedReader in = new BufferedReader(new FileReader(filename));
            String line = in.readLine(); // headers
            while (line != null) {
                line = in.readLine();
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

        Jama.Matrix dX = new Jama.Matrix(numObsFile, 17);
        Jama.Matrix dY = new Jama.Matrix(numObsFile, 1);

        try {
            BufferedReader in = new BufferedReader(new FileReader(filename));
            String line = in.readLine(); // headers
            int i = 0;
            for (int obs = 0; obs < numObs; obs++) {
                line = in.readLine();
                if (line != null) {
                    int a = 0;
                    int b = line.indexOf(",", a);
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
             * Put the endogenous X in the first column to make keeping track of
             * it easier
             */
            X = pmUtility.getColumn(dX, 1); // the endogenous variable
            X = pmUtility.concatMatrix(X, pmUtility.getColumn(dX, 0));
            X = pmUtility.concatMatrix(X, pmUtility.getColumn(dX, 2));
            X = pmUtility.concatMatrix(X, pmUtility.getColumn(dX, 16));  	// the instrument

            Y = dY;

            // SPLITTING VARIABLES FOR HETEROGENEITY
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

            boolean FAKE_DATA = false;
            if (FAKE_DATA) {
                Random rng = new Random(rngSeed);
                NormalDistribution normal = new NormalDistribution();
                for (i = 0; i < numObs; i++) {
                    X.set(i, 3, 3 * rng.nextDouble()); // this is now the instrument
                    double error = normal.inverse(rng.nextDouble());
                    X.set(i, 0, X.get(i, 3) + 0.3 * error); // this is the endogenous variable

                    double beta1 = -5.0;
                    double beta2 = 3.0;
                    boolean observableHeterogeneity = true;
                    if (observableHeterogeneity) {
                        if (Z.get(i, 4) == 1.0) { // black
                            beta1 = 5.0;
                        }

                        if (Z.get(i, 5) == 1.0 && Z.get(i, 4) == 1.0) { // south and black
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
    }

    @Override
    public String getVariableName(int variableIndex) {
        return varNames[variableIndex];
    }

    @Override
    public String getFixedEffectName(int variableIndex, int fixedEffectIndex) {
        if (variableIndex == 7) {
            String[] regionNames = {"New England", "Mid Atlantic", "East North Central", "West North Central", "South Atlantic", "East South Central", "West South Central", "Mountain", "Pacific"};
            return regionNames[fixedEffectIndex - 1];
        }
        return " " + fixedEffectIndex;
    }

    @Override
    public int getNumParams() {
        return X.getColumnDimension() - 1;
    }
}
