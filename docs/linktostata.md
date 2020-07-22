## How to use Moment Forests on Stata


### Components

The codes are grouped into two folders: jars and java. 

Folder “jars” contains relevant java utilities, the Moment Forests jar file that is compiled from the Moment Forests java files, and the Stata ado file. 
Folder “java” contains all the original java codes that implement Moment Forests. These java codes are grouped into two subdirectories. 


###  How to run the Stata code for Moment Forests 

Users should download all the files from folder “jars” and save them in their own Stata personal directories in order to run the Stata code “momentforests”. Users do not need to download the files from folder “java”.

The Stata personal directory can be found by typing “adopath” in the Stata command window. The one starts with “(PERSONAL)” is the personal directory. For example, "c:\ado\personal/".

<img scr="./adopath.PNG" width="10" >

{% comment %} 
![](./adopath.PNG)
{% endcomment %}

### How to modify the original codes

If a user wants to check how the algorithm is structured or wants to modify some parts of it, one needs to look at the original java codes in folder “java”. After making modifications to these java codes, users can compile them into a jar file with the name “momentforests.jar”. Then, by replacing this jar file in their own Stata personal directory, users can use their own version of the Moment Forests Stata code.


### Stata function “momentforests” 

Syntax
momentforests depvar treatment indepvars [if] [in] [, options]

For example, 
momentforests y w x1 x2, mink_size(1) bootstrap(300) num_tree(150)


