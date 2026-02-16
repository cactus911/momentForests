# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Build & Run

No build system (no gradle/maven). Compile and run manually with javac/java. Classpath uses semicolons on Windows bash (must be quoted).

**Compile everything:**
```bash
cd /c/git/momentForests
javac -d bin -cp "jars/guava-31.1-jre.jar;jars/Jama-1.0.3.jar;jars/jsci-core.jar;jars/jfreechart-1.0.19.jar;jars/jcommon-1.0.23.jar;jars/itext-1.3.jar" \
  java/core/*.java java/examples/linear/*.java java/examples/linearIV/*.java java/examples/logit/*.java \
  java/pmUtility/src/utility/*.java java/optimization/src/linear_algebra/*.java \
  java/optimization/src/optimization/*.java java/optimization/src/mcmc/*.java \
  java/experimental/SmartScrollTextArea.java java/experimental/TestResamplingMatrix.java \
  java/experimental/DifferentialEvolution.java java/experimental/LinearMoment.java
```

Exclude `java/experimental/CallMatlab.java` (requires MATLAB engine JAR not present).

**Run an example (e.g., linear Monte Carlo):**
```bash
java -cp "bin;jars/guava-31.1-jre.jar;jars/Jama-1.0.3.jar;jars/jsci-core.jar;jars/jfreechart-1.0.19.jar;jars/jcommon-1.0.23.jar;jars/itext-1.3.jar" examples.linear.LinearTestMain
```

Replace `examples.linear.LinearTestMain` with other main classes: `examples.logit.LogitTestMain`, `examples.linearIV.LinearIVMain`.

**Testing:** No JUnit. Each example has a main class that runs Monte Carlo simulations. Verify correctness by checking classification rates and MSE output.

## Architecture

**Paper:** Nekipelov, Novosad, and Ryan (2020). "Moment Forests." Trees are universal approximators for parameter heterogeneity θ(Z).

### Three-class contract for new models

Every econometric model requires exactly three things:

1. **`MomentSpecification`** (abstract class) — defines the model: data loading, variable configuration, DGP truth (for Monte Carlo), goodness-of-fit. Call `initializeHomogeneity(numParams)` in constructor. Homogeneity state management is handled by the base class.

2. **`ContainerMoment`** (abstract class) — single-leaf estimator: `computeBetaAndErrors()`, `getBeta()`, `getGoodnessOfFit()`, moment functions (`getMomentGWithoutDivision`, `getGi`, `getJacobianNoDivision`), and `getMomentFunctionImposingHomogeneity(k, value)`. Critical: moment functions return **sums without dividing by n**.

3. **Main class** — orchestrates: create spec → `loadData()` → configure `TreeOptions` → build `MomentForest` → `growForest()` → optionally `applyHomogeneityVotes()` → evaluate.

### Core framework classes (java/core/)

- **`MomentForest`** — manages ensemble of `TreeMoment`s. Handles resampling, honest splitting, parallel growth, and the `applyHomogeneityVotes()` one-call homogeneity pipeline.
- **`TreeMoment`** — recursive tree structure. `determineSplit()` grows the tree (Brent's method for continuous, partition enumeration for discrete). `estimateHonestTree()` re-estimates using held-out data. `testHomogeneity()` runs Wald tests with subsampling.
- **`DataLens`** — index-based view into data matrices (no copying). Supports resampling, splitting, and subsampling operations.
- **`TreeOptions`** — centralized configuration (min observations, max depth, subsampling parameters, adaptive stopping, etc.).
- **`WaldTestWholeTree`** — Wald-type test for parameter homogeneity across leaves, with subsampling null distribution.
- **`HomogeneousSearchContainer`** — numerical optimization (BFGS) to find best global values for parameters classified as homogeneous.
- **`GenericContinuousSplitObj` / `GenericPartitionObj`** — generic split objectives that delegate to `computeOptimalBeta()`, eliminating per-model split classes.

### Homogeneity testing pipeline

1. Grow initial forest with `testParameterHomogeneity = true`
2. Each tree's `testHomogeneity()` computes Wald statistics per parameter, uses subsampling for p-values
3. `applyHomogeneityVotes()` does majority voting across trees, then searches for optimal homogeneous values
4. Grow final forest with homogeneous parameters fixed

**Adaptive subsampling:** Stops early when MC standard error of p-hat is small relative to |p_hat - alpha|. Asymmetric boundaries: `adaptiveStoppingMultiplier` (rejection side, default 3.0) vs `adaptiveNonRejectionMultiplier` (non-rejection side, default 5.0).

### Key data flow

- **X** = regressors (what multiplies beta), **Y** = outcome, **Z** = heterogeneity variables (tree splits on these)
- Honest splitting: `getProportionObservationsToEstimateTreeStructure()` fraction for growing, remainder for estimation
- Each tree resamples independently; `DataLens` tracks indices without copying matrices

## Examples as reference implementations

- **`examples/linear/`** — OLS, simplest reference. `LinearTestMain` has full Monte Carlo with homogeneity classification, progress bar, visualization.
- **`examples/linearIV/`** — GMM with instruments. Shows how `getNumParams()` differs from X column count.
- **`examples/logit/`** — MLE-based. Shows log-likelihood goodness of fit instead of SSE.

## Dependencies

All in `jars/`: Jama (matrix algebra), JSci (statistics/distributions), JFreeChart (plotting), Guava, iText (PDF export).

## Key conventions

- Matrix library is JAMA (`Jama.Matrix`). Row vectors are 1×K matrices, column vectors are K×1.
- `pmUtility` has helpers: `prettyPrint()`, `OLS()`, `getColumn()`, `concatMatrix()`, `sortByColumn()`.
- Parallel execution uses Java 8 parallel streams (`parallelStream().forEach()`).
- GUI uses Swing with `JFrame`, `JTextArea` for output, optional `JFreeChart` visualization panels.
