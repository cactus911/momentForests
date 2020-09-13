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
 3. hybrid case: \\( \beta = sin(x1)*(x2 -5) \\), where \\( x1= (0,2 \pi) \\) & \\( x2 = {1,..,10} \\)



<br>
## II. Detailed procedures to run Moment Forests on Stata

1. Load the [simulation datasets](https://github.com/cactus911/momentForests/tree/master/Monte_Carlo) to Stata.
2. Proceed with Moment Forests estimation. (see the detailed specification below)
3. Report and plot the simulation results.



<br>
## III. Results
The estimation results are reported and plotted here. The full estimation results are stored [here](https://github.com/cactus911/momentForests/tree/master/Monte_Carlo)


### 1/3 Discrete Case

<img src="./1_discrete_screen.PNG" width="700" >

<img src="./1_discrete_graph.png" width="500" >


<br>
### 2/3 Continuous Case

<img src="./2_continuous_screen.PNG" width="700" >

<img src="./2_continuous_graph.png" width="500" >



<br>
### 3/3 Hybrid Case

<img src="./3_hybrid_screen.PNG" width="700" >

<img src="./3_hybrid_graph.png" width="500" >




<br>
[back](./index.md)
















