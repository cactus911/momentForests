# How to use Moment Forests on Stata


## Components

The codes are grouped into two folders: jars and java. 

Folder “jars” contains relevant java utilities, the Moment Forests jar file that is compiled from the Moment Forests java files, and the Stata ado file. 
Folder “java” contains all the original java codes that implement Moment Forests. These java codes are grouped into two subdirectories. 


##  How to run the Stata code for Moment Forests 

### 1. Direct Download without using NetBeans.

#### Step 1/3. Find out Stata personal directories

The Stata personal directory can be found by typing “adopath” in the Stata command window. The one starts with “(PERSONAL)” is the personal directory. For example, "c:\ado\personal/".

<img src="./adopath.PNG" width="800" >

{% comment %} 
![](./adopath.PNG)
{% endcomment %}

#### Step 2/3. Download files from the GitHub website.

Download all the files from folder [“jars”](https://github.com/cactus911/momentForests/tree/master/jars) and save them in Stata personal directories. Make sure all 15 files (13 jar files, 1 readme text file, 1 Stata ado file) are properly downloaded and stored.

<img src="./jars.png" width="800" >


#### Step 3/3.

Close and reopen Stata. Now you are ready use Moment Forests on Stata.



### 2. Download using NetBeans.





### How to modify the original codes

If a user wants to check how the algorithm is structured or wants to modify some parts of it, one needs to look at the original java codes in folder “java”. After making modifications to these java codes, users can compile them into a jar file with the name “momentforests.jar”. Then, by replacing this jar file in their own Stata personal directory, users can use their own version of the Moment Forests Stata code.


### Stata function “momentforests” 

Syntax
momentforests depvar treatment indepvars [if] [in] [, options]

For example, 
momentforests y w x1 x2, mink_size(1) bootstrap(300) num_tree(150)


[back](./index.md)
