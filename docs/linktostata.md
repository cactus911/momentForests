# How to use Moment Forests on Stata


## I. Components

The codes are grouped into two directories: [jars and java](https://github.com/cactus911/momentForests)

Folder “jars” contains relevant java utilities, the Moment Forests jar file compiled from the Moment Forests java files, and a Stata ado file. 
Folder “java” contains all the java source codes that implement Moment Forests. The common components that can be applied to any appliations are stored in the subdirectory "core", while specific application files including RCT are in the subdirectory "examples".


## II. How to run Moment Forests on Stata

### Step 1/3. Figure out where the Stata personal directories are

The Stata personal directory can be found by typing “adopath” in the Stata command window. The one starts with “(PERSONAL)” is the personal directory. For example, "c:\ado\personal/".

<img src="./adopath.png" width="700" >

{% comment %} 
![](./adopath.png)
{% endcomment %}

### Step 2/3. Download files from the GitHub website.

Download all the files from [“jars”](https://github.com/cactus911/momentForests/tree/master/jars) and save them in the Stata personal directory. Make sure all 15 files (13 jar files, 1 readme text file, 1 Stata ado file) are properly downloaded and stored.

<img src="./jars.png" width="700" >


### Step 3/3.

Close then reopen Stata. Now you are ready to use Moment Forests on Stata!




## III. Stata command “momentforests”

# Title
momentforests - Moment forests

# Description
momentforests performs Moment forests estimation proposed by Nekipelov, Novosad, and Ryan (2020).

# Quick start
For example, 
momentforests y w x1 x2, mink_size(1) bootstrap(300) num_tree(150)

# Syntax
momentforests depvar treatment indepvars [if] [in] [, options]


syntax varlist [if] [in] [, num_tree(numlist)] [mink(numlist)] [mink_lower(numlist)] [mink_size(numlist)] [mink_upper(numlist)] [msebar(numlist)] [msebar_lower(numlist)] [msebar_size(numlist)] [msebar_upper(numlist)] [bootstrap(numlist) ] [cv(numlist) ]



[back](./index.md)







## How to modify the original codes

If a user wants to check how the algorithm is structured or wants to modify it, one needs to look at the original java codes in folder [“java”](https://github.com/cactus911/momentForests/tree/master/java). After making modifications, users should compile them into a jar file “momentforests.jar”. Then, by replacing this jar file in their own Stata personal directory, one can use own version of Moment Forests on Stata.
