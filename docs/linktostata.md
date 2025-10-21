# How to use Moment Forests on Stata

<br>
## I. Components

The source code is organized into two main directories: [jars and java](https://github.com/cactus911/momentForests)

The `jars` folder contains supporting Java utilities, the compiled `momentforests.jar` file (built from the Java source code), and a Stata .ado file for integration with Stata.

The `java` folder includes all Java source files that implement the Moment Forests framework. Within this folder, the `core` subdirectory contains general-purpose components applicable to a wide range of applications, while the `examples` subdirectory holds application-specific files, including the example implementation of Card (1995) in Stata.

The most recent and actively updated files are available in the Refactor branch of the repository.

## II. How to run Moment Forests on Stata

### Step 1/3. Find Stata personal directory

Stata personal directory can be found by typing “adopath” in the Stata command window. The one starts with “(PERSONAL)” is the personal directory. For example, "c:\ado\personal/".

<img src="./adopath.png" width="700" >

{% comment %} 
![](./adopath.png)
{% endcomment %}

### Step 2/3. Download related files to Stata personal directory

Download the needed files to the Stata personal directory mentioned in Step 1/3.

#### a. .jar files

Download the .jar files from [“jars”](https://github.com/cactus911/momentForests/tree/master/jars) and save them in the Stata personal directory. 

#### b. STATA files

Download the STATA .do and .ado files from [“for_Stata”](https://github.com/cactus911/momentForests/tree/master/for_Stata)

### Step 3/3.

Close then reopen Stata. Now you are ready to use Moment Forests on Stata!

[back](./index.md)

## III. Stata command “momentforest”

{% comment %} 
### Title
momentforest - Moment forest
### Description
momentforest performs a moment forest estimation proposed by Nekipelov, Novosad, and Ryan (2020).
{% endcomment %}

### Variables
Y : outcome variable

X : explanatory variables that have a direct effect on the outcome variable 

Z : observables that are a source of heterogeneity in the effects of X on Y

### Syntax
<!-- `momentforest y [x] [if] [in] [, options]` -->
`momentforest y x [, options]`

| **Option** | **Description** |
|:--|:--|
| **Required** |  |
| `z(varlist)` | List of variables in Z. |
| `numtrees(integer)` | Number of trees in the moment forest. |
| `seed(integer)` | Random-number seed. |
| **Optional** |  |
| `discretevars(varlist)` | Variables in Z that are treated as discrete. |
| `strata(varlist)` | One or more stratification variables in Z for stratified random sampling. |
| `propstructure(string)` | Proportion of observations used to estimate the structure of the trees. Default is `0.35`. |
| `testhomogeneity(string)` | Flag for whether the moment forest should test for homogeneity. Default is `true`. |
| `cv(string)` | Flag for whether the moment forest should perform cross-validation. If `cv(true)` is used, you must also supply all three `cvgrid()` options below. |
| `cvgridminleaf(string)` | Space-separated values for the minimum leaf size grid (e.g., `5 10 20`). Only allowed when `cv(true)` is set. If `cv()` is false, default is `25`. |
| `cvgridminimp(string)` | Space-separated values for the minimum MSE improvement grid for a split (e.g., `0.1 1 10`). Only allowed when `cv(true)` is set. If `cv()` is false, default is `0.1`. |
| `cvgridmaxdepth(string)` | Space-separated values for the maximum tree depth grid (e.g., `3 4 5 6`). Only allowed when `cv(true)` is set. If `cv()` is false, default is `5`. |
| `gen(string)` | Stub for names of generated parameter-estimate variables. |

[back](./index.md)


{% comment %} 

- If a user wants to check how the algorithm is structured, one needs to go through these java codes. 
- If a user wants to modify the codes, one should compile the modified codes into a jar file named “momentforests.jar”. Then, by replacing this jar file in their own Stata personal directory, one can use own version of Moment Forests.
{% endcomment %}

## IV. Monte Carlo simulations

This section works through the Monte Carlo simulations in the paper. Consider the following linear model:  

$$
Y = X'\beta(Z) + \epsilon.
$$

Let $X = \{X_1, X_2\}$, where $X_1$ is a vector of ones, $X_2 \sim N(0,2)$, and $\epsilon$ is an idiosyncratic error term that is distributed as a standard normal. Let $Z = \{Z_1, Z_2, Z_3 \}$, with $Z_1 \sim N(0,1)$, $Z_2 \sim U[0,1]$, and $Z_3$ is a discrete variable taking on two values (1,2) with probability (0.3,0.7). There are four sample sizes, $n \in$ {500, 1000, 2000, 4000}, and three different $\beta(Z)$ corresponding to the cases of no parameter heterogeneity, heterogeneity in only one dimension of $X$, and heterogeneity in both dimensions of $X$. For each configuration, the code contrast the performance of the estimator against the unrestricted model where no parameters are identified as being universally homogeneous. Each Monte Carlo experiment is repeated 500 times to obtain distributions of the statistics of interest.

The estimator faces several simultaneous challenges: first, the moment forest has to correctly classify which components of $Z$ best improve the fit of the model without overfitting. Second, given those splits, it has to consistently estimate the parameters that best match the empirical moments. Third, it has to correctly classify and estimate parameters that are restricted to be homogeneous across the state space. Errors in one stage directly lead to errors in the other stages, so this is a good test of how each of the stages of the estimator work separately and in conjunction.

To test these three stages, we simulate three scenarios. First, we consider a fully heterogeneous parameter specification:  

$$
\beta(Z) =
\begin{cases}
    (-1.0, 1.0) & \text{ if } Z_1 > 0 \\
    (0.33, -1.0) & \text{ else.}
\end{cases}
$$

In the second scenario, we restrict one of the parameters to be homogeneous across the state space:  

$$
\beta(Z) =
\begin{cases}
    (-1.0, 1.0) & \text{ if } Z_1 > 0 \\
    (-1.0, -1.0) & \text{ else.}
\end{cases}
$$

And in the third case, we fix parameters to be homogeneous everywhere:  

$$
\beta(Z) = (-1.0, 1.0).
$$



[back](./index.md)

## V. Worked example

This section works through Card (1995).

[back](./index.md)
