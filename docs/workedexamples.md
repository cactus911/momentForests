# Worked Examples of Moment Forests


<br>
## I. Data Generation Process
In order to show the overall performance of the algorithm, we proceeded with three simulations based on different types of parameters.

\\[
Y = X \beta + \epsilon
\\]
where \\( \epsilon \\) follows N(0,1)
 1. discrete case: \\( \beta = x1 +10*(x2 -1) \\), where \\( x1, x2 = {1,..,10} \\)
 2. continuous case: \\( \beta = sin(x) \\), where \\( x = (0,2 \pi) \\)
 3. hybrid case: \\( \beta = sin(x1)*(x2 -5) \\), where \\( x1= (0,2 \pi) and x2 = {1,..,10} \\)



<br>
## II. Detailed procedures to run Moment Forests on Stata

1. Load the [dataset](https://github.com/cactus911/momentForests/tree/master/Monte_Carlo) to Stata.
2. Proceed with Moment Forests with default values.
3. Plot the simulation results.



<br>
## III. Result 1/3 Discrete Case




<br>
## III. Result 2/3 Continuous Case





<br>
## III. Result 3/3 Hybrid Case






<br>
[back](./index.md)
















