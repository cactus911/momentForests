## Welcome

This website contains Java source code related to the Moment Forest estimator from Nekipelov, Novosad, and Ryan (2020). Our goal is to  provide three sets of related source code:

1. Java source code implementing the core moment forest functionality, including bootstrapped standard errors, that works with moment-based objective functions. The goal is that, with minimal effort, researchers can modify these codes to use their own data sets and moment functions.
2. Java and Stata source code linking the moment forest estimator to Stata. This allows for an easy bridge between existing Stata data sets and workflow and the moment forest.
3. Worked examples of several common estimators, including randomized control trials (RCTs), regression discontinuity designs (RDDs), multiple outcomes, linear regression, and logit models.

This website is a work in progress but will be updated regularly.

### The Paper

The paper can be found here \[[PDF](https://sites.wustl.edu/stephenpryan/files/2016/10/momentTrees.pdf)\].

### How Does It Work?

The fundamental building block of our estimator is a moment:
\\[
E[Y-m(X;\theta)] = 0,
\\]
where \\( m(X;\theta) \\) is the data-generating process which maps \\( X \\) into \\( Y \\) as a function of the parameter \\( \theta \\). Empirical models specify that relationship and then seek to find a vector of parameters, \\( \theta \\), in order to match the outcomes as closely as possible in some metric.

The primary innovation of the moment forest is to replace that model with the following:
\\[
E[Y-m(X;\theta(Z))] = 0,
\\]
where the parameter vector now depends on the matrix \\( Z \\), which may or may not be a subset of \\( X \\). Essentially, the moment forest produces an estimate of how the structure of \\( \theta \\) changes with \\( Z \\), allowing for arbitrary observable heterogeneity in parameters.

### Java Source Code

### Hello World Example

### Link to STATA

### Worked Examples
