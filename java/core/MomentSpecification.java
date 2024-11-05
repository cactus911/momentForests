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
package core;

import Jama.Matrix;
import java.util.Random;
import utility.pmUtility;

/**
 *
 * @author Stephen P. Ryan <stephen.p.ryan@wustl.edu>
 */
public interface MomentSpecification {

    public Double getPredictedY(Jama.Matrix xi, Jama.Matrix beta, Random rng);

    public int[] getVariableIndicesToSearchOver();

    public MomentContinuousSplitObj getFminObjective(DataLens lens, int k, double minProportionEachPartition, int minCountEachPartition);

    public MomentPartitionObj getMomentPartitionObj(DataLens lens, int k, IntegerPartition get);

    public ContainerMoment computeOptimalBeta(DataLens lens, boolean allParametersHomogeneous);

    //  public void generateData(int numObs, Random rng, boolean addNoise);
    // public Jama.Matrix getY(boolean residualizeY);
    public Jama.Matrix getY();

    public Jama.Matrix getX();

    public Jama.Matrix getZ();

    public Jama.Matrix getBalancingVector();

    public int numberoftrees();

    public Boolean[] getDiscreteVector();

    /**
     * Return the true \beta at a given vector z_i
     * 
     * @param zi Observable vector that determines the parameters.
     * @param rng Random number generator for drawing random coefficients
     * @return 
     */
    public Matrix getBetaTruth(Matrix zi, Random rng);

    public DataLens getOutOfSampleXYZ(int numObsOutOfSample, long rngSeed);

    public void loadData(long rngSeed);
    
    public String getVariableName(int variableIndex);

    public String getFixedEffectName(int variableIndex, int fixedEffectIndex);

    public String formatTreeLeafOutput(Jama.Matrix beta, Jama.Matrix variance);
    
    public void setHomogeneousParameter(int parameterIndex, double value);
    public double getHomogeneousParameter(int parameterIndex);

    // public double getHomogeneousComponent(Jama.Matrix xi); // this doesn't make sense in general

    public void resetHomogeneityIndex();

    public void setHomogeneousIndex(Integer i);
    public boolean[] getHomogeneousIndex();
    public Jama.Matrix getHomogeneousParameterVector();
    
    // public Jama.Matrix residualizeX(Jama.Matrix Xp);

    public ContainerMoment getContainerMoment(DataLens lens);
    public abstract int getNumMoments();
    
    public double getGoodnessOfFit(double yi, Jama.Matrix xi, Jama.Matrix beta);
    
    default OutOfSampleStatisticsContainer computeOutOfSampleStatistics(int numberTreesInForest, long rngBaseSeedMomentForest, boolean verbose,
            int minObservationsPerLeaf, double minImprovement, int maxTreeDepth, long rngBaseSeedOutOfSample) {

        double outOfSampleResultsY = 0;
        double outOfSampleResultsBeta = 0;

        MomentForest myForest;

        DataLens homogenizedForestLens = new DataLens(getX(), getY(), getZ(), null);

        myForest = new MomentForest(this, numberTreesInForest, rngBaseSeedMomentForest, homogenizedForestLens, verbose, new TreeOptions());
        TreeOptions cvOptions = new TreeOptions(0.01, minObservationsPerLeaf, minImprovement, maxTreeDepth, false); // k = 1
        myForest.setTreeOptions(cvOptions);
        /**
         * Grow the moment forest
         */
        myForest.growForest();

        // myForest.getTree(0).printTree();
        // debugOutputArea.append(myForest.getTree(0).toString());
        /**
         * Test vectors for assessment
         */
        DataLens oosDataLens = getOutOfSampleXYZ(2000, rngBaseSeedOutOfSample); // this should eventually be modified to come out of the data itself (or generalized in some way)
        Jama.Matrix testZ = oosDataLens.getZ();
        Jama.Matrix testX = oosDataLens.getX();
        Jama.Matrix testY = oosDataLens.getY();

        Random rng = new Random(rngBaseSeedOutOfSample - 3);

        for (int i = 0; i < testZ.getRowDimension(); i++) {
            Jama.Matrix zi = testZ.getMatrix(i, i, 0, testZ.getColumnDimension() - 1);
            Jama.Matrix xi = testX.getMatrix(i, i, 0, testX.getColumnDimension() - 1);

            // going to compare directly to the true parameter vector in this method instead of using fit of Y
            Jama.Matrix bTruth = getBetaTruth(zi, rng);

            Jama.Matrix compositeEstimatedBeta = myForest.getEstimatedParameterForest(zi);

            outOfSampleResultsY += getGoodnessOfFit(testY.get(i, 0), xi, compositeEstimatedBeta);

            if (i < 10) {
                String hString = "[ ";
                for (int k = 0; k < bTruth.getRowDimension(); k++) {
                    if (getHomogeneousIndex()[k]) {
                        hString = hString + "X ";
                    } else {
                        hString = hString + "O ";
                    }
                }
                hString = hString + "]";
                System.out.print("Composite estimated beta: " + pmUtility.stringPrettyPrintVector(compositeEstimatedBeta) + " " + hString + " "+testY.get(i,0)+" "+getPredictedY(xi, compositeEstimatedBeta, rng)+ "\n");
            }
            //pmUtility.prettyPrintVector(compositeEstimatedBeta);

            outOfSampleResultsBeta += pmUtility.sumSquaredElements(compositeEstimatedBeta.minus(bTruth));
        }

        outOfSampleResultsBeta /= testZ.getRowDimension();
        outOfSampleResultsY /= testZ.getRowDimension();

        boolean manualCheck = false;
        if (manualCheck) {
            Jama.Matrix zi = testZ.getMatrix(0, 0, 0, testZ.getColumnDimension() - 1);

            zi.set(0, 0, -1);
            System.out.print("Forest z1 < 0: ");
            pmUtility.prettyPrintVector(myForest.getEstimatedParameterForest(zi));

            System.out.print("Forest z1 > 0: ");
            zi.set(0, 0, 1.0);
            pmUtility.prettyPrintVector(myForest.getEstimatedParameterForest(zi));
        }

        // jt.append("betaMSE: " + (outOfSampleFit / testZ.getRowDimension()) + " \t [" + rngSeed + "]\n");
        return new OutOfSampleStatisticsContainer(outOfSampleResultsBeta, outOfSampleResultsY);
    }

    public int getNumParams();

    public boolean didEstimatorFail();

}
