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
package core;

/**
 * Configuration options for tree and forest construction.
 * Centralizes previously hard-coded constants from TreeMoment.
 *
 * @author Stephen P. Ryan <stephen.p.ryan@wustl.edu>
 */
public class TreeOptions {

    // Existing options
    private double minProportion = 0.001;
    private int minCount = 5;
    private double minMSEImprovement = 0.01;
    private int maxDepth = 100;
    private boolean testParameterHomogeneity = false;

    // Previously hard-coded in TreeMoment (Issue 8)
    private int randomForestMaxVariables = 5;
    private int gridSearchSteps = 100;
    private double gridSearchEpsilon = 1E-30;
    private double subsamplingExponent = 0.7;
    private int numSubsamples = 5000;
    private boolean useRandomForest = false;
    private boolean useSubsampling = true;
    private SubsampleListener subsampleListener = null;

    /**
     * Adaptive subsampling configuration.
     *
     * When enabled, the number of subsample draws adapts based on how decisive
     * the test outcome is. Drawing stops early when the Monte Carlo standard
     * error of the estimated p-value is small relative to the distance from
     * the significance level (alpha = 0.05), following:
     *
     *   Andrews, D.W.K. and Buchinsky, M. (2000). "A Three-Step Method Using
     *     Inequality and Equality Constraints for Minimum Number of Bootstrap
     *     Replications." Econometrica, 68(6), 1557-1574.
     *
     *   Andrews, D.W.K. and Buchinsky, M. (2001). "Evaluation of a Three-Step
     *     Method for Choosing the Number of Bootstrap Repetitions." Journal of
     *     Econometrics, 103, 345-386.
     *
     *   Davidson, R. and MacKinnon, J.G. (2000). "Bootstrap Tests: How Many
     *     Bootstraps?" Econometric Reviews, 19(1), 55-68.
     *
     * The stopping rule uses asymmetric multipliers to reflect the different
     * costs of false rejection vs. false non-rejection:
     *
     *   If p_hat < alpha: stop when (alpha - p_hat) > c_reject * se(p_hat)
     *   If p_hat > alpha: stop when (p_hat - alpha) > c_nonreject * se(p_hat)
     *
     * where se(p_hat) = sqrt(p_hat * (1 - p_hat) / B). The rejection-side
     * multiplier c_reject (default 3.0) gives ~99.7% confidence. The
     * non-rejection-side multiplier c_nonreject (default 5.0) is higher
     * because falsely concluding homogeneity (missing real heterogeneity)
     * is a substantive error, while falsely rejecting homogeneity is only
     * an efficiency cost.
     *
     * numSubsamples serves as B_max (the ceiling). minSubsamples is B_min.
     */
    private boolean useAdaptiveSubsampling = false;
    private int minSubsamples = 500;
    private double adaptiveStoppingMultiplier = 3.0;
    private double adaptiveNonRejectionMultiplier = 5.0;

    public TreeOptions() {
    }

    public TreeOptions(double minProportion, int minCount, double minMSEImprovement, int maxDepth, boolean testParameterHomogeneity) {
        this.maxDepth = maxDepth;
        this.minMSEImprovement = minMSEImprovement;
        this.minCount = minCount;
        this.minProportion = minProportion;
        this.testParameterHomogeneity = testParameterHomogeneity;
    }

    public void setTestParameterHomogeneity(boolean testParameterHomogeneity) {
        this.testParameterHomogeneity = testParameterHomogeneity;
    }

    public boolean isTestParameterHomogeneity() {
        return testParameterHomogeneity;
    }

    public double getMinProportion() {
        return minProportion;
    }

    public void setMinProportion(double minProportion) {
        this.minProportion = minProportion;
    }

    public int getMinCount() {
        return minCount;
    }

    public void setMinCount(int minCount) {
        this.minCount = minCount;
    }

    public double getMinMSEImprovement() {
        return minMSEImprovement;
    }

    public void setMinMSEImprovement(double minMSEImprovement) {
        this.minMSEImprovement = minMSEImprovement;
    }

    public int getMaxDepth() {
        return maxDepth;
    }

    public void setMaxDepth(int maxDepth) {
        this.maxDepth = maxDepth;
    }

    public int getRandomForestMaxVariables() {
        return randomForestMaxVariables;
    }

    public void setRandomForestMaxVariables(int randomForestMaxVariables) {
        this.randomForestMaxVariables = randomForestMaxVariables;
    }

    public int getGridSearchSteps() {
        return gridSearchSteps;
    }

    public void setGridSearchSteps(int gridSearchSteps) {
        this.gridSearchSteps = gridSearchSteps;
    }

    public double getGridSearchEpsilon() {
        return gridSearchEpsilon;
    }

    public void setGridSearchEpsilon(double gridSearchEpsilon) {
        this.gridSearchEpsilon = gridSearchEpsilon;
    }

    public double getSubsamplingExponent() {
        return subsamplingExponent;
    }

    public void setSubsamplingExponent(double subsamplingExponent) {
        this.subsamplingExponent = subsamplingExponent;
    }

    public int getNumSubsamples() {
        return numSubsamples;
    }

    public void setNumSubsamples(int numSubsamples) {
        this.numSubsamples = numSubsamples;
    }

    public boolean isUseRandomForest() {
        return useRandomForest;
    }

    public void setUseRandomForest(boolean useRandomForest) {
        this.useRandomForest = useRandomForest;
    }

    public boolean isUseSubsampling() {
        return useSubsampling;
    }

    public void setUseSubsampling(boolean useSubsampling) {
        this.useSubsampling = useSubsampling;
    }

    public SubsampleListener getSubsampleListener() {
        return subsampleListener;
    }

    public void setSubsampleListener(SubsampleListener subsampleListener) {
        this.subsampleListener = subsampleListener;
    }

    public boolean isUseAdaptiveSubsampling() {
        return useAdaptiveSubsampling;
    }

    public void setUseAdaptiveSubsampling(boolean useAdaptiveSubsampling) {
        this.useAdaptiveSubsampling = useAdaptiveSubsampling;
    }

    public int getMinSubsamples() {
        return minSubsamples;
    }

    public void setMinSubsamples(int minSubsamples) {
        this.minSubsamples = minSubsamples;
    }

    public double getAdaptiveStoppingMultiplier() {
        return adaptiveStoppingMultiplier;
    }

    public void setAdaptiveStoppingMultiplier(double adaptiveStoppingMultiplier) {
        this.adaptiveStoppingMultiplier = adaptiveStoppingMultiplier;
    }

    public double getAdaptiveNonRejectionMultiplier() {
        return adaptiveNonRejectionMultiplier;
    }

    public void setAdaptiveNonRejectionMultiplier(double adaptiveNonRejectionMultiplier) {
        this.adaptiveNonRejectionMultiplier = adaptiveNonRejectionMultiplier;
    }
}
