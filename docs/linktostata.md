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
| `cvgridminimp(string)` | Space-separated values for the minimum MSE improvement grid for a split (e.g., `0.0 1e-4 1e-3`). Only allowed when `cv(true)` is set. If `cv()` is false, default is `0.1`. |
| `cvgridmaxdepth(string)` | Space-separated values for the maximum tree depth grid (e.g., `3 4 5 6`). Only allowed when `cv(true)` is set. If `cv()` is false, default is `5`. |
| `gen(string)` | Stub for names of generated parameter-estimate variables. |

[back](./index.md)


{% comment %} 

- If a user wants to check how the algorithm is structured, one needs to go through these java codes. 
- If a user wants to modify the codes, one should compile the modified codes into a jar file named “momentforests.jar”. Then, by replacing this jar file in their own Stata personal directory, one can use own version of Moment Forests.
{% endcomment %}

## IV. Monte Carlo simulations

This section works through the Monte Carlo simulations in the paper.

[back](./index.md)

## V. Worked example

This section works through Card (1995).

[back](./index.md)
