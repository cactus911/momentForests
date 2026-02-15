# Moment Forests: Architecture and User's Guide

This document describes the architecture of the Moment Forests framework after
the recent refactoring. It is intended for economics PhD students (and other
researchers) who know basic Java and econometrics but are new to this codebase.

The framework implements the moment forest estimator from:

> Nekipelov, D., Novosad, P., and Ryan, S. P. (2020). "Moment Forests."

The key idea: **you provide the moment conditions for your econometric model,
and the framework handles everything else** -- tree construction, splitting,
honest estimation, homogeneity testing, and inference.

---

## Table of Contents

1. [Overview](#1-overview)
2. [Core Concepts](#2-core-concepts)
3. [How to Implement a New Model](#3-how-to-implement-a-new-model)
4. [Reference Implementations](#4-reference-implementations)
5. [Configuration: TreeOptions](#5-configuration-treeoptions)
6. [Migration Guide (Old to New Architecture)](#6-migration-guide)

---

## 1. Overview

### What the Framework Does

A moment forest is an ensemble of moment trees. Each tree partitions the
covariate space by splitting on observable characteristics (the Z variables) to
find regions where the model parameters differ. Within each leaf of each tree,
the model is estimated using the moment conditions you provide. The forest
aggregates estimates across trees for stable, honest inference.

The estimation pipeline has three phases:

1. **Cross-validation** -- Select tuning parameters (tree depth, minimum leaf
   size, improvement threshold) by out-of-sample prediction.
2. **Homogeneity testing** -- Grow a forest with homogeneity testing enabled.
   Each tree tests whether each parameter is constant across leaves. Trees
   vote, and parameters classified as homogeneous by a majority are constrained
   to be constant globally.
3. **Final estimation** -- Grow the production forest with homogeneous
   parameters locked in. Extract parameter estimates, predictions, and
   confidence sets.

### What You Provide

To plug in a new econometric model, you write three files:

| File | Purpose |
|------|---------|
| `MyModelMomentSpecification.java` | Extends `MomentSpecification`. Defines the model: data access, number of parameters, number of moments, and a factory method that creates containers. |
| `ContainerMyModel.java` | Extends `ContainerMoment`. Holds the estimation logic for a single leaf: computes beta, evaluates moment conditions, reports goodness of fit. |
| `MyModelMain.java` | The entry point. Orchestrates cross-validation, homogeneity testing, and final estimation. |

The framework handles tree splitting, honest sample splitting, resampling,
parallel tree growth, and the voting mechanism for homogeneity classification.

### What Changed in the Refactoring

| Aspect | Old Architecture | New Architecture |
|--------|-----------------|-----------------|
| `MomentSpecification` | Interface with 34 methods | Abstract class with sensible defaults |
| Files per model | 6 files (~1,300 lines) | 3 files (~700 lines) |
| Boilerplate per model | ~350 lines of identical code | Zero boilerplate |
| Homogeneity management | Copy-pasted in every implementation | Managed by the base class |
| Split objective classes | Two per model (continuous + discrete) | Generic implementations shared by all models |
| Voting loop | Copy-pasted in every Main class | One method call: `myForest.applyHomogeneityVotes(verbose)` |
| Configuration constants | Hard-coded in `TreeMoment` | Centralized in `TreeOptions` |
| Error handling | `System.exit()` calls | Proper exceptions |

---

## 2. Core Concepts

This section describes each major class in the `core` package and how they
relate to one another.

### 2.1 MomentSpecification (abstract class)

**File:** `java/core/MomentSpecification.java`

This is the central abstraction. It defines your econometric model and serves
as the bridge between your model and the tree/forest machinery.

**Responsibilities:**
- Provides data matrices (X, Y, Z) to the framework.
- Creates `ContainerMoment` instances via the factory method
  `computeOptimalBeta()`.
- Manages homogeneity state (which parameters are homogeneous, their values).
- Supplies variable metadata (names, discrete/continuous classification, which
  variables to split on).

**Key design decision:** `MomentSpecification` is an abstract class, not an
interface. This means it can hold state (the homogeneity arrays) and provide
default implementations. You extend it and override only what your model
requires.

The class is organized into sections:

```
Homogeneity Management      -- initializeHomogeneity(), get/set methods
Model Definition (REQUIRED) -- computeOptimalBeta(), getNumParams(), getNumMoments()
Data Access (REQUIRED)      -- getX(), getY(), getZ()
Variable Config (REQUIRED)  -- getDiscreteVector(), getVariableIndicesToSearchOver()
Split Objectives (defaults) -- getFminObjective(), getMomentPartitionObj()
Display/Formatting          -- getVariableName(), formatTreeLeafOutput(), etc.
Goodness of Fit (REQUIRED)  -- getGoodnessOfFit(), getPredictedY()
Configuration (defaults)    -- numberoftrees(), getProportionObservationsToEstimateTreeStructure()
Data Loading (optional)     -- loadData(), getBetaTruth(), getOutOfSampleXYZ()
```

### 2.2 ContainerMoment (abstract class)

**File:** `java/core/ContainerMoment.java`

A `ContainerMoment` represents the estimated model for a single leaf (a subset
of the data). The framework calls your specification's `computeOptimalBeta()`
factory method, which creates one of these and calls `computeBetaAndErrors()`
to run the estimation.

**Methods you must implement:**

| Method | Purpose |
|--------|---------|
| `computeBetaAndErrors()` | Run the estimator. Populate beta and goodness-of-fit. |
| `getBeta()` | Return the estimated parameter vector (Jama.Matrix, K x 1). |
| `getGoodnessOfFit()` | Return a scalar measure of fit (e.g., SSE). The framework minimizes this when choosing splits. |
| `getMomentGWithoutDivision(beta)` | Return the stacked moment vector g(beta) = sum_i g_i(beta) **without** dividing by n. Used in the Wald test. |
| `getGi(beta, i)` | Return the moment contribution of observation i. Used to construct the variance matrix. |
| `getJacobianNoDivision(beta)` | Return the Jacobian of the moment conditions (sum, not averaged). Used in the Wald test. |
| `computeMeasureOfFit(beta)` | Evaluate fit at an arbitrary beta (used internally). |
| `getMomentFunctionValue(beta)` | Evaluate the GMM objective at an arbitrary beta. |
| `getMomentFunctionImposingHomogeneity(k, value)` | Evaluate the GMM objective with parameter k fixed to `value`. Used in homogeneity testing. |
| `getVariance(beta)` | Return the variance-covariance matrix. |

### 2.3 DataLens

**File:** `java/core/DataLens.java`

A `DataLens` is a lightweight, index-based view into the full data matrices.
Instead of copying subsets of the data when splitting or resampling, the
framework maintains an integer index array that points into the original
matrices.

**Why this matters for you:** When you receive a `DataLens` in
`computeOptimalBeta()`, call `lens.getX()`, `lens.getY()`, and `lens.getZ()`
to get the matrices for that leaf. These methods return freshly constructed
sub-matrices containing only the rows in the current view. You can also use
`lens.getNumObs()` for the current observation count, or `lens.getX(i, j)` for
element-level access without constructing a full matrix.

**Key operations performed by the framework:**
- `getResampledDataLens()` -- Bootstrap resampling for each tree.
- `randomlySplitSample()` -- Honest splitting (growing sample vs. estimation
  sample).
- `splitOnRule()` -- Partition observations by a split rule during tree
  traversal.
- `getSubsampledDataLens()` -- Subsampling for the homogeneity test.

### 2.4 TreeOptions

**File:** `java/core/TreeOptions.java`

A plain configuration object that centralizes all tuning parameters. Previously
these were hard-coded constants scattered throughout `TreeMoment`.

| Parameter | Default | Meaning |
|-----------|---------|---------|
| `minProportion` | 0.001 | Minimum fraction of observations in each child after a split. |
| `minCount` | 5 | Minimum number of observations in each child after a split. |
| `minMSEImprovement` | 0.01 | Minimum improvement in the objective function to justify a split. |
| `maxDepth` | 100 | Maximum tree depth. |
| `testParameterHomogeneity` | false | Whether to run homogeneity tests in this forest. |
| `gridSearchSteps` | 100 | Number of grid points when searching for the optimal continuous split. |
| `gridSearchEpsilon` | 1E-30 | Small constant added to grid search increments. |
| `subsamplingExponent` | 0.7 | Controls subsample size b = n^d for the homogeneity test. |
| `numSubsamples` | 5000 | Number of subsamples drawn for the homogeneity test. |
| `useRandomForest` | false | If true, randomly restrict the set of splitting variables at each node. |
| `randomForestMaxVariables` | 5 | Maximum number of variables to consider at each split (when `useRandomForest` is true). |

You typically create `TreeOptions` with the five-argument constructor and leave
the rest at defaults:

```java
TreeOptions opts = new TreeOptions(
    0.01,   // minProportion
    50,     // minCount
    1.0,    // minMSEImprovement
    5,      // maxDepth
    false   // testParameterHomogeneity
);
```

### 2.5 MomentForest

**File:** `java/core/MomentForest.java`

The forest is a collection of `TreeMoment` instances. It manages the
tree-building lifecycle:

1. **Construction:** Takes a `MomentSpecification`, the number of trees, a
   random seed, a `DataLens`, verbosity flag, and `TreeOptions`.

2. **`growForest()`:** For each tree:
   - Resamples the data (with or without treatment balancing / stratification).
   - Splits the resample into a growing sample and an honest estimation sample.
   - Constructs a `TreeMoment` and grows it using the growing sample.
   - Re-estimates parameters in every leaf using the honest sample.
   - Prunes invalid trees (where estimation failed).
   - Runs in parallel via Java 8 parallel streams.

3. **`getEstimatedParameterForest(zi)`:** Given an observation's Z vector,
   looks up the estimated beta in each tree and returns the average.

4. **`applyHomogeneityVotes(verbose)`:** The one-call replacement for the old
   boilerplate voting loop. It:
   - Calls `testHomogeneity()` on every tree (in parallel).
   - Collects votes: each tree reports which parameters it classified as
     homogeneous.
   - Applies a majority rule (>50% of trees must agree).
   - Sets the homogeneous flags and starting values on the specification.
   - Returns the list of parameter indices classified as homogeneous.

### 2.6 TreeMoment

**File:** `java/core/TreeMoment.java`

A single moment tree. This is a recursive data structure: each node is a
`TreeMoment` that is either terminal (a leaf) or has left and right children.

**Growth algorithm (`determineSplit()`):**
1. Estimate the model at the current node to get a baseline objective value.
2. For each candidate splitting variable:
   - If continuous: use Brent's method (`Fmin`) plus a two-pass grid search to
     find the optimal split point.
   - If discrete: enumerate all disjoint partitions and pick the best.
3. Compare the best split's objective to the baseline. If improvement exceeds
   the threshold, split and recurse. Otherwise, declare the node terminal.

**Honest estimation (`estimateHonestTree()`):**
After the tree structure is fixed, the honest sample is distributed to leaves.
Each leaf re-estimates its parameters using only the honest data. If a leaf
ends up empty or estimation fails, the node is pruned back to its parent.

**Homogeneity testing (`testHomogeneity()`):**
For each parameter k:
1. Compute a Wald-type test statistic across all leaves.
2. Use subsampling to build the null distribution of the test statistic.
3. Compute a p-value. Apply the Holm-Bonferroni correction for multiple
   testing.
4. Parameters that fail to reject are classified as homogeneous.

### 2.7 GenericContinuousSplitObj and GenericPartitionObj

**Files:** `java/core/GenericContinuousSplitObj.java`,
`java/core/GenericPartitionObj.java`

These are the generic split objective classes that replaced the per-model split
objectives from the old architecture. They work by delegating to
`spec.computeOptimalBeta()` for each candidate partition, so they work with
any model automatically.

You do not need to touch these classes. The default implementations of
`getFminObjective()` and `getMomentPartitionObj()` in `MomentSpecification`
instantiate them for you. Override only if your model has a non-standard notion
of "effective number of observations" (e.g., in a treatment-control design
where the effective sample depends on treatment counts, not raw counts).

---

## 3. How to Implement a New Model

This section walks through building a new model from scratch. We will use OLS
as the running example, referencing the actual code in `examples/linear/`.

### Step 1: Create Your MomentSpecification

Create a class that extends `MomentSpecification`. This is where you define
what your model is.

```java
package examples.mymodel;

import Jama.Matrix;
import core.ContainerMoment;
import core.DataLens;
import core.MomentSpecification;
import java.util.Random;

public class MyModelMomentSpecification extends MomentSpecification {

    private Matrix X, Y, Z;
    private int numObs;

    // -- Constructor: configure splitting variables and homogeneity --
    public MyModelMomentSpecification(int numObs) {
        this.numObs = numObs;

        // REQUIRED: initialize homogeneity tracking for your parameter count.
        // Call this AFTER you know how many parameters the model has.
        initializeHomogeneity(getNumParams());
    }

    // ===== REQUIRED: Model Definition =====

    @Override
    public ContainerMoment computeOptimalBeta(DataLens lens, boolean allParametersHomogeneous) {
        // Create your container, estimate, and return it.
        ContainerMyModel c = new ContainerMyModel(
            lens,
            getHomogeneousIndex(),
            getHomogeneousParameterVector(),
            allParametersHomogeneous
        );
        c.computeBetaAndErrors();
        return c;
    }

    @Override
    public int getNumParams() {
        return 2; // e.g., intercept + slope
    }

    @Override
    public int getNumMoments() {
        return 2; // for exactly-identified models, numMoments == numParams
    }

    // ===== REQUIRED: Data Access =====

    @Override
    public Matrix getX() { return X; }

    @Override
    public Matrix getY() { return Y; }

    @Override
    public Matrix getZ() { return Z; }

    // ===== REQUIRED: Variable Configuration =====

    @Override
    public Boolean[] getDiscreteVector() {
        // One entry per column of Z.
        // true = discrete (partitioned), false = continuous (split point).
        return new Boolean[] { false, false, true };
    }

    @Override
    public int[] getVariableIndicesToSearchOver() {
        // Which columns of Z to consider when searching for splits.
        return new int[] { 0, 1, 2 };
    }

    // ===== REQUIRED: Goodness of Fit =====

    @Override
    public double getGoodnessOfFit(double yi, Matrix xi, Matrix beta) {
        double predicted = (xi.times(beta)).get(0, 0);
        double error = yi - predicted;
        return error * error;
    }

    @Override
    public Double getPredictedY(Matrix xi, Matrix beta, Random rng) {
        if (beta == null) return null;
        return (xi.times(beta)).get(0, 0);
    }

    // ===== OPTIONAL: Data Loading =====

    @Override
    public void loadData(long rngSeed) {
        // Load from file or generate Monte Carlo data.
        // Populate X, Y, Z.
    }
}
```

**Required methods summary:**

| Method | What to return |
|--------|---------------|
| `computeOptimalBeta(lens, allHomogeneous)` | A new `ContainerMoment` with estimates for this data subset. |
| `getNumParams()` | Number of parameters (K). |
| `getNumMoments()` | Number of moment conditions (M >= K). |
| `getX()`, `getY()`, `getZ()` | The full data matrices. |
| `getDiscreteVector()` | Boolean array: which Z columns are discrete. |
| `getVariableIndicesToSearchOver()` | Int array: which Z columns to consider for splitting. |
| `getGoodnessOfFit(yi, xi, beta)` | Scalar loss for one observation. |
| `getPredictedY(xi, beta, rng)` | Predicted outcome for one observation. |

**Methods with sensible defaults (override only if needed):**

| Method | Default | Override when... |
|--------|---------|-----------------|
| `getFminObjective(...)` | Uses `GenericContinuousSplitObj` | Your model has a custom "effective N" |
| `getMomentPartitionObj(...)` | Uses `GenericPartitionObj` | Same as above |
| `getContainerMoment(lens)` | Delegates to `computeOptimalBeta(lens, false)` | You need different logic for inference vs. tree-building |
| `getVariableName(i)` | Returns "Z1", "Z2", ... | You want descriptive names in output |
| `getFixedEffectName(i, j)` | Returns "Group j" | You have named categories |
| `formatTreeLeafOutput(beta, var)` | Pretty-prints beta | You want custom formatting |
| `numberoftrees()` | 100 | You want a different default |
| `getProportionObservationsToEstimateTreeStructure()` | 0.5 | You need more/fewer growing observations |
| `getBalancingVector()` | null | You have treatment-control data |
| `getStratificationIndex()` | null | You need stratified sampling |
| `didEstimatorFail()` | false | Your estimator can fail |
| `loadData(seed)` | No-op | You load/generate data |
| `getBetaTruth(zi, rng)` | Throws exception | Monte Carlo only |
| `getOutOfSampleXYZ(n, seed)` | Throws exception | Monte Carlo only |

**Critical step:** Call `initializeHomogeneity(numParams)` in your constructor.
This sets up the boolean array and parameter vector that the framework uses to
track which parameters are homogeneous. Forgetting this call will cause a
`NullPointerException` at runtime.

### Step 2: Create Your ContainerMoment

This class holds the estimation logic for a single leaf. It extends
`ContainerMoment` and typically implements `Uncmin_methods` if you use the
built-in optimizer.

```java
package examples.mymodel;

import Jama.Matrix;
import core.ContainerMoment;
import core.DataLens;
import optimization.Uncmin_f77;
import optimization.Uncmin_methods;

public class ContainerMyModel extends ContainerMoment implements Uncmin_methods {

    private Matrix beta;
    private double goodnessOfFit = Double.POSITIVE_INFINITY;
    private Matrix X, Y;
    private boolean[] homogeneityIndex;
    private Matrix homogeneityParameters;
    private boolean allParametersHomogeneous;

    public ContainerMyModel(DataLens lens, boolean[] homogeneityIndex,
                            Matrix homogeneityParameters,
                            boolean allParametersHomogeneous) {
        this.X = lens.getX();
        this.Y = lens.getY();
        this.homogeneityIndex = homogeneityIndex;
        this.homogeneityParameters = homogeneityParameters;
        this.allParametersHomogeneous = allParametersHomogeneous;
    }

    @Override
    public void computeBetaAndErrors() {
        // Guard against tiny samples
        if (Y.getRowDimension() < 30) {
            beta = null;
            goodnessOfFit = Double.POSITIVE_INFINITY;
            return;
        }

        try {
            // Count free parameters (subtract homogeneous ones)
            int numFree = homogeneityIndex.length;
            if (!allParametersHomogeneous) {
                for (boolean b : homogeneityIndex) {
                    if (b) numFree--;
                }
            }

            // ... set up and run Uncmin_f77 or your own optimizer ...
            // ... or use a closed-form solution like OLS ...

            // After optimization, reconstruct the full beta vector:
            beta = new Matrix(homogeneityIndex.length, 1);
            int counter = 0;
            for (int i = 0; i < homogeneityIndex.length; i++) {
                if (homogeneityIndex[i] && !allParametersHomogeneous) {
                    // This parameter is constrained to be homogeneous.
                    beta.set(i, 0, homogeneityParameters.get(i, 0));
                } else {
                    // This parameter was estimated freely.
                    beta.set(i, 0, optimizedValues[counter]);
                    counter++;
                }
            }

            // Compute goodness of fit
            goodnessOfFit = computeSSE(beta);

        } catch (Exception e) {
            beta = null;
            goodnessOfFit = Double.POSITIVE_INFINITY;
        }
    }

    @Override
    public Matrix getBeta() { return beta; }

    @Override
    public double getGoodnessOfFit() { return goodnessOfFit; }

    @Override
    public Matrix getMomentGWithoutDivision(Matrix beta) {
        // Return sum_i g_i(beta) -- do NOT divide by n.
        // For OLS: g_i = x_i * (x_i'beta - y_i), so G = X'(Xbeta - Y)
        Matrix g = new Matrix(X.getColumnDimension(), 1);
        Matrix fitted = X.times(beta);
        Matrix errors = fitted.minus(Y);
        for (int i = 0; i < X.getRowDimension(); i++) {
            for (int k = 0; k < X.getColumnDimension(); k++) {
                g.set(k, 0, g.get(k, 0) + errors.get(i, 0) * X.get(i, k));
            }
        }
        return g;
    }

    @Override
    public Matrix getGi(Matrix beta, int i) {
        // Moment contribution of observation i.
        double fitted = (X.getMatrix(i, i, 0, X.getColumnDimension()-1)
                          .times(beta)).get(0, 0);
        double error = fitted - Y.get(i, 0);
        Matrix gi = new Matrix(X.getColumnDimension(), 1);
        for (int k = 0; k < X.getColumnDimension(); k++) {
            gi.set(k, 0, error * X.get(i, k));
        }
        return gi;
    }

    @Override
    public Matrix getJacobianNoDivision(Matrix beta) {
        // For OLS: Jacobian = X'X (sum, not averaged)
        return (X.transpose()).times(X);
    }

    @Override
    public double getMomentFunctionValue(Matrix beta) {
        // GMM objective: g'g where g = (1/n) * G
        Matrix g = getMomentGWithoutDivision(beta);
        g.timesEquals(1.0 / Y.getRowDimension());
        return (g.transpose().times(g)).get(0, 0);
    }

    @Override
    public double getMomentFunctionImposingHomogeneity(int k, double value) {
        Matrix betaCopy = beta.copy();
        betaCopy.set(k, 0, value);
        return getMomentFunctionValue(betaCopy);
    }

    @Override
    public double computeMeasureOfFit(Matrix beta) {
        Matrix fitted = X.times(beta);
        Matrix errors = fitted.minus(Y);
        return (errors.transpose().times(errors)).get(0, 0);
    }

    @Override
    public Matrix getVariance(Matrix beta) {
        // Standard OLS variance: sigma^2 * (X'X)^{-1}
        double sse = computeMeasureOfFit(beta);
        double sigma2 = sse / Y.getRowDimension();
        return (X.transpose().times(X)).inverse().times(sigma2);
    }

    // --- Uncmin_methods interface ---

    @Override
    public double f_to_minimize(double[] x) {
        // Reconstruct beta from x (skipping homogeneous params), evaluate GMM.
        // ...
        return getMomentFunctionValue(b);
    }

    @Override
    public void gradient(double[] x, double[] g) {
        throw new UnsupportedOperationException();
    }

    @Override
    public void hessian(double[] x, double[][] h) {
        throw new UnsupportedOperationException();
    }
}
```

**Important notes on the Container:**

- **Homogeneity handling.** When `homogeneityIndex[k]` is true and
  `allParametersHomogeneous` is false, parameter k has been classified as
  homogeneous. Fix it to `homogeneityParameters.get(k, 0)` and do not include
  it in the optimization. When `allParametersHomogeneous` is true, ignore the
  index and estimate all parameters freely (this happens when the tree is a
  stump).

- **`getMomentGWithoutDivision` must NOT divide by n.** The framework handles
  the division in different places depending on context. Returning the raw sum
  is essential for the Wald test to work correctly.

- **`getJacobianNoDivision` must also NOT divide by n.** Same reason.

- **Handle small samples gracefully.** If the leaf has too few observations to
  estimate the model, set `beta = null` and `goodnessOfFit =
  Double.POSITIVE_INFINITY`. The framework will prune or skip such leaves.

### Step 3: Create Your Main Class

The Main class orchestrates the three-phase pipeline. Here is a simplified
template based on `LinearTestMain`:

```java
package examples.mymodel;

import core.DataLens;
import core.HomogeneousSearchContainer;
import core.MomentForest;
import core.TreeOptions;
import java.util.ArrayList;

public class MyModelMain {

    public static void main(String[] args) {

        // ===== Phase 0: Setup =====
        int numObs = 2000;
        int numberTrees = 40;
        boolean verbose = false;
        long seed = 12345L;
        java.util.Random rng = new java.util.Random(seed);

        MyModelMomentSpecification spec = new MyModelMomentSpecification(numObs);
        spec.loadData(rng.nextLong());

        // Cross-validated tuning parameters (or set directly)
        int minObsPerLeaf = 50;
        double minImprovement = 1.0;
        int maxDepth = 5;

        long forestSeed = rng.nextLong();

        // ===== Phase 1: Homogeneity Testing =====
        spec.resetHomogeneityIndex();

        DataLens forestLens = new DataLens(
            spec.getX(), spec.getY(), spec.getZ(), null
        );

        TreeOptions testOpts = new TreeOptions(
            0.01, minObsPerLeaf, minImprovement, maxDepth, false
        );

        MomentForest testForest = new MomentForest(
            spec, numberTrees, forestSeed, forestLens, verbose, new TreeOptions()
        );
        testForest.setTreeOptions(testOpts);
        testForest.growForest();

        // One call replaces the old 30-line voting loop:
        ArrayList<Integer> homogeneousParams =
            testForest.applyHomogeneityVotes(verbose);

        System.out.println("Homogeneous parameters: " + homogeneousParams);

        // ===== Phase 2: Estimate Homogeneous Parameter Values =====
        if (!homogeneousParams.isEmpty()) {
            HomogeneousSearchContainer search = new HomogeneousSearchContainer(
                spec, numberTrees, verbose,
                minImprovement, minObsPerLeaf, maxDepth,
                homogeneousParams,
                forestSeed, rng.nextLong()
            );
            search.executeSearch();
            // Values are now set inside spec.
        }

        // ===== Phase 3: Final Estimation =====
        MomentForest finalForest = new MomentForest(
            spec, numberTrees, forestSeed, forestLens, verbose, new TreeOptions()
        );
        TreeOptions finalOpts = new TreeOptions(
            0.01, minObsPerLeaf, minImprovement, maxDepth, false
        );
        finalForest.setTreeOptions(finalOpts);
        finalForest.growForest();

        // Query the forest for estimated parameters at a point z:
        // Jama.Matrix betaHat = finalForest.getEstimatedParameterForest(zi);
    }
}
```

**Key workflow points:**

- Always call `spec.resetHomogeneityIndex()` before the homogeneity testing
  phase to ensure no stale flags from a previous run.
- `applyHomogeneityVotes(verbose)` does three things internally: (1) runs
  `testHomogeneity()` on every tree in parallel, (2) tallies votes across
  trees, (3) sets the homogeneous flags and starting values on `spec`. It
  returns the list of homogeneous parameter indices.
- `HomogeneousSearchContainer.executeSearch()` uses numerical optimization to
  find the best global values for the homogeneous parameters. It grows forests
  internally at each evaluation of the objective function.
- After setting homogeneous parameters, grow the final forest. The trees will
  automatically detect which parameters are constrained and produce stump trees
  if all parameters are homogeneous.

---

## 4. Reference Implementations

The `examples/` directory contains three complete models. Study them in this
order.

### 4.1 Linear Regression (OLS)

**Directory:** `java/examples/linear/`

**Files:**
- `LinearMomentSpecification.java` -- Extends `MomentSpecification`. The
  simplest possible specification.
- `ContainerLinear.java` -- OLS estimation. Uses a closed-form solution
  (`pmUtility.OLS()`) when there are no homogeneous parameters, and falls back
  to numerical optimization (`Uncmin_f77`) when some parameters are
  constrained.
- `LinearTestMain.java` -- Full Monte Carlo harness with cross-validation,
  homogeneity testing, and out-of-sample evaluation.
- `LinearDataGenerator.java` -- Generates synthetic data with known
  heterogeneity structure.

**What to learn from this example:**

- How `initializeHomogeneity(dimensionX)` is called in the constructor.
- How `computeOptimalBeta()` creates a `ContainerLinear`, calls
  `computeBetaAndErrors()`, and returns it.
- How the container handles the `homogeneityIndex` / `homogeneityParameters`
  / `allParametersHomogeneous` arguments to correctly fix or free parameters.
- The moment conditions for OLS: g_i(beta) = x_i * (x_i'beta - y_i).
- How `getGoodnessOfFit()` returns squared error for a single observation.

**Relevant code from `LinearMomentSpecification`:**

```java
public LinearMomentSpecification(int numObs, int dimX) {
    this.numObs = numObs;
    this.dimensionX = dimX;
    initializeHomogeneity(dimensionX);         // <-- required call
    int[] vsi = {0, 1, 2};
    Boolean[] wvd = {false, false, true};      // Z1 continuous, Z2 continuous, Z3 discrete
    variableSearchIndex = vsi;
    DiscreteVariables = wvd;
}

@Override
public ContainerMoment computeOptimalBeta(DataLens lens, boolean allParametersHomogeneous) {
    ContainerLinear l = new ContainerLinear(
        lens,
        getHomogeneousIndex(),
        getHomogeneousParameterVector(),
        allParametersHomogeneous,
        this
    );
    l.computeBetaAndErrors();
    failedEstimatorIndicator = l.didEstimatorFail();
    return l;
}
```

### 4.2 Instrumental Variables (IV / 2SLS via GMM)

**Directory:** `java/examples/linearIV/`

**Files:**
- `LinearIVMomentSpecification.java` -- Shows how to handle endogenous
  regressors and excluded instruments.
- `ContainerLinearIV.java` -- GMM estimation with instruments.
- `LinearIVTestMain.java` -- Applied example using Card (1995) returns to
  education data.

**What is different from OLS:**

1. **X matrix layout.** The first column of X is the endogenous regressor
   (education). The last column is the excluded instrument (proximity to
   college). The middle columns are exogenous regressors.

2. **`getNumParams()` differs from `X.getColumnDimension()`.** Parameters are
   the coefficients on the regressors, excluding the instrument column:
   ```java
   public int getNumParams() {
       return X.getColumnDimension() - 1;
   }
   ```

3. **Moment conditions use instruments, not regressors.** In OLS, the moment
   is E[x_i * e_i] = 0. In IV, it is E[z_i * e_i] = 0 where z_i includes
   the exogenous regressors and the instrument but excludes the endogenous
   regressor. See `ContainerLinearIV.addGi()`:
   ```java
   public void addGi(Matrix g, double error_i, int i) {
       // Skip the first column (endogenous regressor)
       for (int k = 1; k < X.getColumnDimension(); k++) {
           g.set(k - 1, 0, g.get(k - 1, 0) + error_i * X.get(i, k));
       }
   }
   ```

4. **Fitted values exclude the instrument.** When computing Xbeta, use only
   the first `X.getColumnDimension() - 1` columns:
   ```java
   Matrix fit = X.getMatrix(0, X.getRowDimension()-1, 0, X.getColumnDimension()-2)
                 .times(beta);
   ```

5. **Data is loaded from a real CSV file**, not generated synthetically.

### 4.3 Logit (Discrete Choice)

**Directory:** `java/examples/logit/`

**Files:**
- `LogitMomentSpecification.java` -- Binary choice model.
- `ContainerLogit.java` -- Maximum likelihood / GMM estimation of logit.
- `LogitTestMain.java` -- Monte Carlo harness.
- `LogitDataGenerator.java` -- Generates binary outcomes from the logit model.

**What is different from OLS:**

1. **Goodness of fit is log-likelihood based**, not SSE:
   ```java
   public double getGoodnessOfFit(double yi, Matrix xi, Matrix beta) {
       return ContainerLogit.computeLLHi(yi, xi, beta);
   }
   ```

2. **`getPredictedY` returns a discrete draw**, not a continuous prediction:
   ```java
   public Double getPredictedY(Matrix xi, Matrix beta, Random rng) {
       return LogitDataGenerator.getLogitDiscreteOutcome(xi, beta, rng);
   }
   ```

3. **The number of moments and parameters are both 2** (intercept and slope
   coefficient on a single regressor), matching the exactly-identified case.

---

## 5. Configuration: TreeOptions

This section provides guidance on choosing `TreeOptions` parameters.

### Creating TreeOptions

```java
// Default options (useful for the initial MomentForest constructor)
TreeOptions defaults = new TreeOptions();

// Customized options (commonly used for cross-validation and final estimation)
TreeOptions custom = new TreeOptions(
    0.01,   // minProportion: at least 1% of data in each child
    50,     // minCount: at least 50 observations in each child
    1.0,    // minMSEImprovement: split must improve objective by at least 1.0
    5,      // maxDepth: tree can have at most 5 levels
    false   // testParameterHomogeneity: do not test in this forest
);
```

### Parameter Guidance

**`minCount` (minimum observations per leaf):**
Start with 50 for moderate sample sizes (n ~ 1,000-5,000). Increase for
larger datasets or models with many parameters. Too low leads to noisy
estimates; too high prevents the tree from finding heterogeneity.

**`minMSEImprovement` (improvement threshold):**
This is the most impactful tuning parameter. A value of 0 allows all
splits that improve fit. Values of 1-10 are reasonable for most
applications. Cross-validate to choose.

**`maxDepth`:**
Controls model complexity. Values of 1-5 are typical. A depth of 0 forces
a stump (one leaf, no splits). The framework automatically sets depth to 0
when all parameters are classified as homogeneous.

**`subsamplingExponent`:**
Controls the subsample size b = n^d for the homogeneity test. The default
of 0.7 works well in practice. Values closer to 1.0 use larger subsamples
(more power, slower). Values closer to 0.5 use smaller subsamples (faster,
less power).

**`numSubsamples`:**
Number of subsamples for constructing the null distribution of the test
statistic. The default of 5,000 provides good resolution for the 95th
percentile critical value. Reduce for faster computation during
development; increase for final results.

### Applying Options

The pattern used throughout the codebase is:

```java
MomentForest forest = new MomentForest(
    spec, numTrees, seed, lens, verbose, new TreeOptions()  // defaults at construction
);
forest.setTreeOptions(customOptions);  // override before growing
forest.growForest();
```

The reason for this two-step pattern is that `MomentForest` stores the default
options at construction but allows them to be overridden before growing. This
supports workflows where you construct the forest once and grow it multiple
times with different options during cross-validation.

---

## 6. Migration Guide

If you have an existing model written against the old interface-based
architecture, follow these steps to migrate.

### Step 1: Change the class declaration

**Before:**
```java
public class MyModelSpecification implements MomentSpecification {
```

**After:**
```java
public class MyModelSpecification extends MomentSpecification {
```

### Step 2: Replace homogeneity boilerplate with one method call

**Before (approximately 30 lines per model):**
```java
private boolean[] homogeneityIndex;
private Jama.Matrix homogeneousParameterVector;

// constructor:
homogeneityIndex = new boolean[numParams];
homogeneousParameterVector = new Jama.Matrix(numParams, 1);

public void resetHomogeneityIndex() { /* ... */ }
public void setHomogeneousIndex(Integer i) { /* ... */ }
public boolean[] getHomogeneousIndex() { /* ... */ }
public void setHomogeneousParameter(int i, double v) { /* ... */ }
public double getHomogeneousParameter(int i) { /* ... */ }
public Jama.Matrix getHomogeneousParameterVector() { /* ... */ }
```

**After (one line in your constructor):**
```java
// In your constructor, after you know the number of parameters:
initializeHomogeneity(numParams);
```

Delete all the homogeneity fields and methods listed above. They are now
inherited from `MomentSpecification`.

### Step 3: Delete your split objective classes

**Before:** You had two files:
- `MyModelContinuousSplitObj.java` (extending `MomentContinuousSplitObj`)
- `MyModelPartitionObj.java` (extending `MomentPartitionObj`)

**After:** Delete both files. The generic implementations
(`GenericContinuousSplitObj` and `GenericPartitionObj`) handle this for all
models. The default `getFminObjective()` and `getMomentPartitionObj()` methods
in `MomentSpecification` instantiate them automatically.

### Step 4: Remove the split objective overrides from your specification

**Before:**
```java
@Override
public MomentContinuousSplitObj getFminObjective(DataLens lens, int k,
        double minProp, int minCount) {
    return new MyModelContinuousSplitObj(k, lens, minProp, minCount, this);
}

@Override
public MomentPartitionObj getMomentPartitionObj(DataLens lens, int k,
        IntegerPartition partition) {
    return new MyModelPartitionObj(partition, k, lens, this);
}
```

**After:** Delete both methods. The defaults in `MomentSpecification` do the
right thing.

### Step 5: Replace the voting loop with one method call

**Before (in your Main class, approximately 40 lines):**
```java
// Grow forest for homogeneity testing
myForest.testHomogeneity(verbose);
boolean[] votes = myForest.getHomogeneityVotes(jt, verbose);
double[] startingValues = myForest.getHomogeneityStartingValues();

ArrayList<Integer> homogeneousParameterList = new ArrayList<>();
ArrayList<Double> homogeneousStartingValues = new ArrayList<>();
for (int i = 0; i < votes.length; i++) {
    if (votes[i]) {
        homogeneousParameterList.add(i);
        homogeneousStartingValues.add(startingValues[i]);
    }
}
// ... sorting ...
spec.resetHomogeneityIndex();
for (int i = 0; i < homogeneousParameterList.size(); i++) {
    spec.setHomogeneousIndex(homogeneousParameterList.get(i));
    spec.setHomogeneousParameter(
        homogeneousParameterList.get(i),
        homogeneousStartingValues.get(i)
    );
}
```

**After (one line):**
```java
ArrayList<Integer> homogeneousParams = myForest.applyHomogeneityVotes(verbose);
```

The method handles testing, voting, sorting, resetting the specification, and
setting all the flags and values.

### Step 6: Update computeOptimalBeta signature

**Before:**
```java
public ContainerMoment computeOptimalBeta(DataLens lens);
```

**After:**
```java
public ContainerMoment computeOptimalBeta(DataLens lens, boolean allParametersHomogeneous);
```

The `allParametersHomogeneous` flag tells the container whether to bypass the
homogeneity index and estimate all parameters freely (used when the tree is a
stump). Pass it through to your Container constructor.

### Step 7: Replace System.exit() with exceptions

If your old code had `System.exit()` calls for error handling, replace them
with appropriate exceptions:

```java
// Before:
System.out.println("Error: matrix not invertible");
System.exit(1);

// After:
throw new IllegalStateException("Error: matrix not invertible");
```

### Step 8: Use TreeOptions for configuration

If you had hard-coded constants in your TreeMoment interactions, move them to
`TreeOptions`:

```java
// Before (scattered in Main class):
int gridSteps = 100;
double subsamplingD = 0.7;

// After:
TreeOptions opts = new TreeOptions();
opts.setGridSearchSteps(100);
opts.setSubsamplingExponent(0.7);
```

### Summary Checklist

- [ ] `implements MomentSpecification` changed to `extends MomentSpecification`
- [ ] Homogeneity fields and methods deleted; replaced with
      `initializeHomogeneity(numParams)` in constructor
- [ ] Model-specific `*ContinuousSplitObj.java` file deleted
- [ ] Model-specific `*PartitionObj.java` file deleted
- [ ] `getFminObjective()` override removed from specification
- [ ] `getMomentPartitionObj()` override removed from specification
- [ ] Voting loop in Main class replaced with
      `myForest.applyHomogeneityVotes(verbose)`
- [ ] `computeOptimalBeta` signature updated to include
      `boolean allParametersHomogeneous`
- [ ] `System.exit()` calls replaced with proper exceptions
- [ ] Hard-coded constants moved to `TreeOptions`

After migration, your model should consist of exactly three files (Specification,
Container, Main) totaling roughly half the line count of the old implementation.

---

## Appendix: Package Structure

```
java/
  core/
    MomentSpecification.java        -- Abstract base class for all models
    ContainerMoment.java            -- Abstract base class for leaf estimation
    MomentForest.java               -- Forest of trees, voting, parallel growth
    TreeMoment.java                 -- Single tree: splitting, honest estimation, tests
    DataLens.java                   -- Index-based view into data matrices
    TreeOptions.java                -- Configuration constants
    GenericContinuousSplitObj.java  -- Generic continuous split (replaces per-model classes)
    GenericPartitionObj.java        -- Generic discrete partition (replaces per-model classes)
    HomogeneousSearchContainer.java -- Numerical search for homogeneous parameter values
    SplitRule.java                  -- Encodes a single split decision
    SplitContainer.java             -- Holds left/right DataLens after a split
    IntegerPartition.java           -- Represents a partition of discrete values
    DisjointSet.java                -- Enumerates all partitions of a discrete set
    WaldTestWholeTree.java          -- Wald test for parameter homogeneity
    HomogeneousParameterSorter.java -- Utility for sorting parameter lists
    OutOfSampleStatisticsContainer.java -- Holds out-of-sample evaluation results
  examples/
    linear/                         -- OLS reference implementation
    linearIV/                       -- IV/2SLS reference implementation
    logit/                          -- Logit reference implementation
  optimization/
    Fmin.java                       -- Brent's method for 1-D minimization
    Uncmin_f77.java                 -- Quasi-Newton optimizer (BFGS/Newton)
    Uncmin_methods.java             -- Interface for Uncmin objective functions
  utility/
    pmUtility.java                  -- Matrix utilities (OLS, printing, sorting, etc.)
```
