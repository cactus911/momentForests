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
where the parameter vector now depends on the matrix \\( Z \\), which may or may not be a subset of \\( X \\). The moment forest  estimates how the structure of \\( \theta \\) changes with \\( Z \\), allowing for arbitrary observable heterogeneity in parameters. Conditioning on that partitioning of the parameter space, one estimates the parameter \\( \theta \\) as usual through the generalized method of moments.

The basic element of the moment forest is a *moment tree*. Moment trees provide the partioning of parameters with \\( Z \\) and are grown recursively from a base *root* node. At the root node, a universal parameter is fitted to minimize the moment function on all the data. A search is then conducted along each dimension of the matrix \\( Z \\) by splitting the data into two disjoint partitions. A separate moment is fitted within each of those partitions for all possible splits. The moment tree branches along the split point which generates the largest improvement in fit. This process is then repeated recursively until a stopping condition is met. The end result is an assignment of parameters to subsets of \\( Z \\). A moment forest is grown by resampling the data and growing many moment trees. The moment forest returns an estimate of the parameter \\( \theta(Z) \\) by averaging parameters returned by each of constituent trees in the forest.

In our theory, we prove several useful theoretical properties of moment forests. As is well-known in the literature, the trees in our setting are universal approximators, and can recover any \(( \theta(Z) \)), subject to some modest complexity bounds and weak regularity conditions on the moment function. Second, we show consistency and prove uniform convergence, including an associated rate. These results are key for two reasons; uniform convergence allows researchers to conduct statistical tests at different points of \\( Z \\), while the rate of convergence allows our results to be used as plug-in estimators in more complex models.

### Java Source Code

The java source code is contained in two separate folders on the Github repository.

### Link to STATA

We have also developed a link between the functionality contained in the Java code and STATA.

### Worked Examples

A simple "Hello, world" style example is coming soon, both in the underlying Java code and using the associated STATA link.
