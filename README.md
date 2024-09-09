# HyperLAU
Hypercubic transition paths with Linear Algebra for Uncertainty: Inference of accumulation pathways in discrete time, under models allowing different degree of dependencies, using observations that might include missing or uncertain states.

![graphical_abstract](https://github.com/user-attachments/assets/25a74577-28a2-460e-a875-f6122375f53b)


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
