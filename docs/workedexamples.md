# Worked Examples of Moment Forests


<br>
## I. Data Generation Process
In order to show the overall performance of the algorithm, we proceeded with three simulations based on different types of parameters.

\\[
Y = W \beta + \epsilon
\\]
where \\( W \\) is a treatment dummy variable and \\( \epsilon \\) follows N(0,1)
 1. discrete case: \\( \beta(x1,x2) = x1 +10*(x2 -1) \\), where \\( x1, x2 = {1,..,10} \\)
 2. continuous case: \\( \beta(x) = sin(x) \\), where \\( x = (0,2 \pi) \\)
 3. hybrid case: \\( \beta(x1,x2) = sin(x1)*(x2 -5) \\), where \\( x1= (0,2 \pi) \\) & \\( x2 = {1,..,10} \\)



<br>
## II. Monte Carlo Simulation

#### Step 1/4. Download the simulation datasets and Stata do file in [here](https://github.com/cactus911/momentForests/tree/master/Monte_Carlo) to your own working directory.

<img src="./simulationdatacode.png" width="500" >

<br>
#### Step 2/4. Open Stata and open the Stata do file downloaded (workedexample.do).

<img src="./workedexample_do.png" width="700" >

<br>
#### Step 3/4. Set your working directory in the do file. 

For example, `cd "C:\Users\Valued Customer"`

<br>
#### Step 4/4. Run the rest of the code.


<br>
## III. Results
The estimation results are reported and plotted below.


### 1/3. Discrete Case: \\( \beta(x1,x2) = x1 +10*(x2 -1) \\)

| <img src="./1_discrete_screen.PNG" width="700" > |
|:--:| 
| *Results printed on Stata window* |


| <img src="./1_discrete_graph.png" width="500" > |
|:--:| 
| *Scatter plot, true beta against estimated beta* |


<br>
### 2/3. Continuous Case: \\( \beta(x) = sin(x) \\)

| <img src="./2_continuous_screen.PNG" width="700" > |
|:--:| 
| *Results printed on Stata window* |

| <img src="./2_continuous_graph.png" width="500" > |
|:--:| 
| *Scatter plot, x against estimated beta* |



<br>
### 3/3. Hybrid Case: \\( \beta(x1,x2) = sin(x1)*(x2 -5) \\)

| <img src="./3_hybrid_screen.PNG" width="700" > |
|:--:| 
| *Results printed on Stata window* |



| <img src="./3_hybrid_graph.png" width="500"  alt="{{ Scatter plot, true beta against estimated beta }}"> | 
|:--:| 
| *Scatter plot, true beta against estimated beta* |



<br>
[back](./index.md)




