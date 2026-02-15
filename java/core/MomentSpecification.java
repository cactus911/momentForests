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
 * Abstract base class for moment specifications. Subclasses must provide:
 * - Model definition (computeOptimalBeta, getNumParams, getNumMoments)
 * - Data access (getX, getY, getZ)
 * - Variable configuration (getDiscreteVector, getVariableIndicesToSearchOver)
 *
 * Homogeneity management, split objective creation, and display formatting
 * are handled by default implementations that can be overridden if needed.
 *
 * @author Stephen P. Ryan <stephen.p.ryan@wustl.edu>
 */
public abstract class MomentSpecification {

    // ==================== Homogeneity Management ====================
    // Previously copy-pasted across every implementation; now managed here.

    private boolean[] homogeneityIndex;
    private Jama.Matrix homogeneousParameterVector;

    /**
     * Initialize homogeneity tracking arrays. Must be called by subclass
     * constructors after the number of parameters is known.
     */
    protected void initializeHomogeneity(int numParams) {
        homogeneityIndex = new boolean[numParams];
        homogeneousParameterVector = new Jama.Matrix(numParams, 1);
        resetHomogeneityIndex();
    }

    public void resetHomogeneityIndex() {
        if (homogeneityIndex != null) {
            for (int i = 0; i < homogeneityIndex.length; i++) {
                homogeneityIndex[i] = false;
            }
        }
    }

    public void setHomogeneousIndex(Integer i) {
        homogeneityIndex[i] = true;
    }

    public boolean[] getHomogeneousIndex() {
        return homogeneityIndex;
    }

    public Jama.Matrix getHomogeneousParameterVector() {
        return homogeneousParameterVector;
    }

    public void setHomogeneousParameter(int parameterIndex, double value) {
        homogeneousParameterVector.set(parameterIndex, 0, value);
    }

    public double getHomogeneousParameter(int parameterIndex) {
        return homogeneousParameterVector.get(parameterIndex, 0);
    }

    // ==================== Model Definition (MUST override) ====================

    /**
     * Create a ContainerMoment for the given data, estimate parameters.
     * This is the primary factory method that the framework uses.
     */
    public abstract ContainerMoment computeOptimalBeta(DataLens lens, boolean allParametersHomogeneous);

    /**
     * Number of parameters in the model.
     */
    public abstract int getNumParams();

    /**
     * Number of moment conditions.
     */
    public abstract int getNumMoments();

    // ==================== Data Access (MUST override) ====================

    public abstract Jama.Matrix getY();
    public abstract Jama.Matrix getX();
    public abstract Jama.Matrix getZ();

    // ==================== Variable Configuration (MUST override) ====================

    /**
     * Which Z variables are discrete (true) vs continuous (false).
     */
    public abstract Boolean[] getDiscreteVector();

    /**
     * Indices of Z variables to consider for splitting.
     */
    public abstract int[] getVariableIndicesToSearchOver();

    // ==================== Split Objectives (default: generic) ====================
    // These use the generic implementations. Override only if your model
    // needs custom split logic (e.g., custom effective observation counts).

    /**
     * Create the objective function for continuous variable splitting.
     * Default uses GenericContinuousSplitObj which calls computeOptimalBeta().
     */
    public MomentContinuousSplitObj getFminObjective(DataLens lens, int k, double minProportionEachPartition, int minCountEachPartition) {
        return new GenericContinuousSplitObj(k, lens, minProportionEachPartition, minCountEachPartition, this);
    }

    /**
     * Create the objective function for discrete variable partitioning.
     * Default uses GenericPartitionObj which calls computeOptimalBeta().
     */
    public MomentPartitionObj getMomentPartitionObj(DataLens lens, int k, IntegerPartition partition, int minCount) {
        return new GenericPartitionObj(partition, k, lens, this, minCount);
    }

    // ==================== Container Creation ====================
    // getContainerMoment is used in WaldTestWholeTree and HomogeneousSearchContainer.
    // Default delegates to computeOptimalBeta(lens, false).

    /**
     * Create a ContainerMoment for inference (Wald tests, etc.).
     * Default delegates to computeOptimalBeta with allParametersHomogeneous=false.
     * Override if you need different behavior for inference vs tree-building.
     */
    public ContainerMoment getContainerMoment(DataLens lens) {
        return computeOptimalBeta(lens, false);
    }

    // ==================== Display/Formatting (sensible defaults) ====================

    /**
     * Name of the Z variable at the given index. Override for descriptive names.
     */
    public String getVariableName(int variableIndex) {
        return "Z" + (variableIndex + 1);
    }

    /**
     * Name of a fixed effect level. Override for descriptive names.
     */
    public String getFixedEffectName(int variableIndex, int fixedEffectIndex) {
        return "Group " + fixedEffectIndex;
    }

    /**
     * Format leaf output for tree printing. Override for model-specific formatting.
     */
    public String formatTreeLeafOutput(Jama.Matrix beta, Jama.Matrix variance) {
        if (beta == null) {
            return "null";
        }
        return pmUtility.stringPrettyPrintVector(beta);
    }

    /**
     * Prefix labels for beta parameters. Override for model-specific labels.
     */
    public String getBetaPrefixes() {
        return "";
    }

    // ==================== Goodness of Fit ====================

    /**
     * Compute goodness-of-fit for a single observation.
     * Used in out-of-sample evaluation.
     */
    public abstract double getGoodnessOfFit(double yi, Jama.Matrix xi, Jama.Matrix beta);

    /**
     * Predicted Y for a single observation.
     */
    public abstract Double getPredictedY(Jama.Matrix xi, Jama.Matrix beta, Random rng);

    // ==================== Estimation Failure ====================

    /**
     * Whether the last call to computeOptimalBeta failed.
     * Default returns false; override if your Container can fail.
     */
    public boolean didEstimatorFail() {
        return false;
    }

    // ==================== Data Loading / Simulation (optional) ====================
    // These have do-nothing defaults. Override for Monte Carlo simulations.

    /**
     * Load or generate data. Override in your specification.
     */
    public void loadData(long rngSeed) {
        // Default: no-op. Override to load from file or generate Monte Carlo data.
    }

    /**
     * Return the true beta for Monte Carlo evaluation.
     * Only needed for simulation studies.
     */
    public Matrix getBetaTruth(Matrix zi, Random rng) {
        throw new UnsupportedOperationException("getBetaTruth is only available in Monte Carlo specifications");
    }

    /**
     * Generate out-of-sample test data.
     * Only needed for simulation studies.
     */
    public DataLens getOutOfSampleXYZ(int numObsOutOfSample, long rngSeed) {
        throw new UnsupportedOperationException("getOutOfSampleXYZ is only available in Monte Carlo specifications");
    }

    // ==================== Configuration (sensible defaults) ====================

    /**
     * Number of trees in the forest. Default 100.
     * Now a convenience method; prefer setting this on MomentForest directly.
     */
    public int numberoftrees() {
        return 100;
    }

    /**
     * Proportion of observations used for tree structure (vs honest estimation).
     * Default 0.5 for a balanced split.
     */
    public double getProportionObservationsToEstimateTreeStructure() {
        return 0.5;
    }

    /**
     * Balancing vector for treatment-control balance in resampling.
     * Default null (no balancing).
     */
    public Jama.Matrix getBalancingVector() {
        return null;
    }

    /**
     * Column indices for stratified sampling.
     * Default null (no stratification).
     */
    public int[] getStratificationIndex() {
        return null;
    }

    // ==================== Out-of-Sample Evaluation ====================

    /**
     * Compute out-of-sample statistics by growing a forest and evaluating it.
     */
    public OutOfSampleStatisticsContainer computeOutOfSampleStatistics(int numberTreesInForest, long rngBaseSeedMomentForest, boolean verbose,
            int minObservationsPerLeaf, double minImprovement, int maxTreeDepth, long rngBaseSeedOutOfSample) {

        double outOfSampleResultsY = 0;
        double outOfSampleResultsBeta = 0;

        DataLens originalDataLens = new DataLens(getX(), getY(), getZ(), null);
        DataLens[] twoLenses = originalDataLens.randomlySplitSample(0.8, rngBaseSeedOutOfSample);

        DataLens homogenizedForestLens = twoLenses[0];

        MomentForest myForest = new MomentForest(this, numberTreesInForest, rngBaseSeedMomentForest, homogenizedForestLens, verbose, new TreeOptions());
        TreeOptions cvOptions = new TreeOptions(0.01, minObservationsPerLeaf, minImprovement, maxTreeDepth, false);
        myForest.setTreeOptions(cvOptions);
        myForest.growForest();

        DataLens oosDataLens = twoLenses[1];
        Jama.Matrix testZ = oosDataLens.getZ();
        Jama.Matrix testX = oosDataLens.getX();
        Jama.Matrix testY = oosDataLens.getY();

        Random rng = new Random(rngBaseSeedOutOfSample - 3);

        for (int i = 0; i < testZ.getRowDimension(); i++) {
            Jama.Matrix zi = testZ.getMatrix(i, i, 0, testZ.getColumnDimension() - 1);
            Jama.Matrix xi = testX.getMatrix(i, i, 0, testX.getColumnDimension() - 1);

            Jama.Matrix compositeEstimatedBeta = myForest.getEstimatedParameterForest(zi);

            outOfSampleResultsY += getGoodnessOfFit(testY.get(i, 0), xi, compositeEstimatedBeta);

            Jama.Matrix bTruth = getBetaTruth(zi, rng);
            outOfSampleResultsBeta += pmUtility.sumSquaredElements(compositeEstimatedBeta.minus(bTruth));

            if (i < 10) {
                String hString = "[ ";
                for (int k = 0; k < getNumParams(); k++) {
                    if (getHomogeneousIndex()[k]) {
                        hString = hString + "X ";
                    } else {
                        hString = hString + "O ";
                    }
                }
                hString = hString + "]";
                System.out.print("Composite estimated beta: " + pmUtility.stringPrettyPrintVector(compositeEstimatedBeta) + " " + hString + " " + testY.get(i, 0) + " " + getPredictedY(xi, compositeEstimatedBeta, rng) + "\n");
            }
        }

        outOfSampleResultsY /= testZ.getRowDimension();
        outOfSampleResultsBeta /= testZ.getRowDimension();

        return new OutOfSampleStatisticsContainer(outOfSampleResultsBeta, outOfSampleResultsY);
    }
}
