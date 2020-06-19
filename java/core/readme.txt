Directory for core files that will be used with any application.

We need three things here:

1. Load the data (either using a .csv file or using Stata)
2. Specify which variables are splitting covariates, and, of those, which are discrete and continuous (either in the specification file or via Stata options)
3. The objective function

We need to simplify the code so that this approach is accessible to non-Java experts.

We have two separate objective functions for discrete splits and continuous ones. Almost all the code is the same; we should figure out how to merge those for simplicity.

Also, there should be a default implementation for MomentSpecification to get rid of a lot of those methods.