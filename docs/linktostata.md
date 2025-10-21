# How to use Moment Forests on Stata

<br>
## I. Components

The codes are grouped into two directories: [jars and java](https://github.com/cactus911/momentForests)

Folder “jars” contains relevant java utilities, the moment forest jar file compiled from the Moment Forests java files, and a Stata ado file. 
Folder “java” contains all the java source codes that implement moment forests. The common components that can be applied to any appliations are stored in the subdirectory "core", while specific application files including the Card (1995) example are in the subdirectory "examples".


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
`momentforest y [x] [if] [in] [, options]`

<br>

| Options | Description |
|:-----------|:------------|
| | Required |
| `z(varlist)` | list of variables in Z | 
| `numtrees(#)` | number of trees in the moment forest | 
| `cv(boolean)` | `false` = do not perform cross-validation for hyper-parameters, `true` = do perform cross-validation for hyper-parameters |
| `gen(stub)` | stub for variables names of parameter estimates | 
| |
| | Optional | 
| `cvgridminleaf(range)` | grid search for the minimum number of observations in each leaf | 
| `cvgridminimp(range)` | grid search for the minimum MSE improvement for a split | 
| `cvgridmaxdepth(range)` | grid search for the maximum depth of a tree | 
| `discretevars(varlist)` | list of variables in Z that are discrete |
| `strata(varname)` | stratification variable for stratified random sampling  |

<br>
[back](./index.md)



{% comment %} 

- If a user wants to check how the algorithm is structured, one needs to go through these java codes. 
- If a user wants to modify the codes, one should compile the modified codes into a jar file named “momentforests.jar”. Then, by replacing this jar file in their own Stata personal directory, one can use own version of Moment Forests.
{% endcomment %}

## IV. Monte Carlo simulations

This section works through the Monte Carlo simulations in the paper

## V. Worked example

This section works through Card (1995)
