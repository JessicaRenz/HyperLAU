# HyperLAU
Hypercubic transition paths with Linear Algebra for Uncertainty: Inference of accumulation pathways in discrete time, under models allowing different degree of dependencies, using observations that might include missing or uncertain states.

![graphical_abstract](https://github.com/user-attachments/assets/cb4aa11a-fdaa-493d-81dc-983b1167f10a)



## Requirements
For running HyperLAU, you need the ability to compile and run C++ code, and having the library `armadillo` installed. 

For running the plotting tools as well as the tools to reconstruct a HyperLAU input dataset from the output of PHYLIP, you need R with the following libraries: 
`ggplot2` , `stringr`, `ggraph`, `ggpubr`, `igraph`,`stringr`,`phytools` and `phangorn`.

## Running HyperLAU
HyperLAU can be compiled in the command line by the following command:
```
g++ HyperLAU.cpp -o HyperLAU -larmadillo
```
After compiling, HyperLAU can be runned directly from the command line by specifying the following input parameters:
```
./HyperLAU [name of the input file] [string length] [name of the output file] [number of bootstrap resamples] [number for a random seed] [model] [denom]
```
- **name of the input file** Name of the file that contains the input data, including the extension `.txt`. HyperLAU expects as an input a textfile containing a list of ancestor and descendant states, for example `01? 011`. Both states are encoded by binary strings, but can contain one or more `?` to mark missing or uncertain data. Every line is considered as a sample independent of the others. For using cross-sectional data, just set all ancestor states to the zero-string.
- **string length** Number of features to consider, has to be an integer.
- **name of the output file** HyperLAU will output several text-files. With this input parameter, you can specify the basis of the names of all these outputs. For every run, you will get the two output files `best_likelihood_[name of the output file].txt` and `transitions_[name of the output file].txt`. If the number of bootstrap resamples is specified as $>0$, you will also get two additional files called `mean_[name of the output file].txt`and `sd_[name of the output file].txt`.
- **number of bootstrap resamples** If specified as 0, no bootstrapping will be done.
- **number of random seed** Here you can specify the random seed you want to use for your simulations. Has to be an integer. Has to be specified, no underlying default value.
- **model** Has to be one of the following integers: $-1,1,2,3,4$. With this input parameter you choose which model, i.e. what degree of allowed interaction between the features should be used. Model $-1$ corresponds to the model of arbitrary dependencies, where all combinations of features can influence each other. In model $1$, every feature occurs with a fixed rate, independent of other features already obtained. In model $2$ there are already pairwise interactions allowed, and in model $3$ and $4$, pairs or triples can influence the probability of features to occur next, respectively.
- **denom** This parameter specifies, how fast the temperature in the Simulated Annealing Process should be decreased. After every optimization loop, the current temperature is devided by this parameter: `temp = temp/denom`. This parameter should be a double $>1$.

## Output
HyperLAU does always output the two text files `best_likelihood_[name of the output file].txt` and `transitions_[name of the output file].txt`. If the number of bootstrap resamples is chosen to be bigger than zero, it outputs two additional files `mean_[name of the output file].txt`and `sd_[name of the output file].txt`.
- **best_likelihood_[..].txt** In this output you can find the intermediate results for the currently best log-likelihood. Every line corresponds to a certain loop in the Simulated Annealing optimization and contains the log-likelihood that is stored as the best log-likelihood until the end of this loop. This file helps to track the progression of the optimization and can give some insight if the optimiziation process was running long enough for the log-likelihood to converge.
- **transitions_[...].txt** This file contains the four coloumns `From`, `To`, `Probability` and `Flux`. In the coloumns `From` and `To` are the origin and destination node of a transition, encoded in the integer representation of the binary string. The coloumn `Probability` contains the transition probability that HyperLAU assigned to this edge. The last coloumn reports the calculated flux for the corresponding edge, which reflects the proportion of the trajectories that makes this transition.
- **mean_[...].txt** This file contains three coloumns. The firts two ones give the origin and destination of the considered transition in the same way as in the file `transitions_[...].txt`. The last coloumn is the mean of the transition probabilities for all bootstrap resamples that were runned.
- **sd_[...].txt** This file contains three coloumns. The firts two ones give the origin and destination of the considered transition in the same way as in the file `transitions_[...].txt`. The last coloumn is the standard deviation of the transition probabilities for all bootstrap resamples that were runned.
  
## R Scripts
In this repository, we also provide a bunge of R scripts for plotting the results of HyperLAU, or for preparing or manipulating the input data. 

### random_qm.R
This script can be used to insert a certain amount of randomly uniform distributed "?" into an input dataset, as we did for the different case studies in the article introducing HyperLAU. The only required R library for running this script is `stringr`. Running the first part of the script inserts "?" in the whole dataset, i.e. at all positions in the string, whereas the second part does this only for a certain specified feature. In both cases one has to specify the input data set in the `read.table()` function, as well as the length of the binary strings `L`, the probability with which a position in the binary strings should be replaced by a "?" (this is done via the parameter `threshold`) and the name under which the manipulated data set should be saved (in the `write.table()` function). In the second case, when you want the "?" to be only in one certain feature, you have to specify the number of the position of this feature in the binary string (parameter `feature`). 

### reconstruct.R
This script can be used to reconstruct the input format HyperLAU expects from the output provided by PHYLIP. The required R libraries for running this script are `phytools` and `phangorn`. There are three parameters one has to specify at the top of the script: `input.trees`(output of PHYLIP which contains the information about the phylogenetic tree), `input.data` (the original cross-sectional data set that was used by PHYLIP to create the phylogenetic tree) and `output.data` (name under which the reconstructed file in the input format for HyperLAU should get). 

### plotting_HyperLAU_results.R
This script contains all function that we used to produce the figures in the article introducing HyperLAU. The required R libraries for running this script are `ggplot2`, `ggraph`. `igraph`, `stringr` and `ggpubr`. 

The first three functions `BinToDec`, `DecToBin` and `DecToBinV` are just functions that are needed by the actual plotting functions, for converting binary strings in integer and vice versa. 

The function `plot_embedded_hypercube` is the function that plots the learned transitions pathways embedded in the full hypercube. 

## Data
In the `data` folder in this repository you can find the data sets we used in the article introducing HyperLAU. 

In the subfolder `second_toyexample` there are the full dataset `full.txt` (all positions in all datapoints are specified as $0$ or $1$) we used for our second toy example among the ones where we inserted with a probability of $0.4$ uncertainty markers ("?") into a certain feature. 

The subfolder `tb` contains all three datasets we used for the tuberculosis case studies. `tb_data_9.txt` contains the phylogeny reconstructed by considering only the datapoints of string length nine containing no missing states. Into this dataset, we inserted uncertainty markers ("?") with a probability of $0.5$ at all positions, which is the dataset `tb_data_9_qm50.txt`. The phylogenetic reconstruction (by PHYLIP) of the original data set of string length ten including also data points with missing states is stored as `phylo_tb_10.txt`. 
