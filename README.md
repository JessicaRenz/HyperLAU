# HyperLAU
Hypercubic transition paths with Linear Algebra for Uncertainty: Inference of accumulation pathways in discrete time, under models allowing different degree of dependencies, using observations that might include missing or uncertain states.

![graphical_abstract-igj](https://github.com/user-attachments/assets/8a34d662-9124-46f4-8f7f-ad0301e6f266)



## Requirements
For running HyperLAU, you need the ability to compile and run C++ code, and having the library `armadillo` installed. 

For running the plotting tools as well as the tools to manipulate the input datasets, you need R with the following libraries: 
`ggplot2` , `stringr`, `ggraph`, `ggpubr`, `igraph`, `ggrepel`, `readxl`, `tidyr` and `dplyr`.

## Running HyperLAU
HyperLAU can be compiled in the command line, for example by the following command:
```
g++ HyperLAU.cpp -o HyperLAU -larmadillo
```
After compiling, HyperLAU can be runned directly from the command line by specifying the following input parameters:
```
./HyperLAU [name of input file] [label for output files] 
```
The name of the input file and the label for the output files are the only two necessary input arguments. 

Additionally, you can costumize the simulation by setting some additionally flags, that are set to a specific value by default:
```
./HyperLAU [name of input file] [label for output files] --bootstrap [number of bootstrap resamples] --seed [random seed] --model [model structure] --rate [annealing rate]
```

The command below, for example, will run a dataset called `data.txt` and stores the result in an output file called `transitions_data.txt`. No bootstrap resamples will be executed, the random seed is set to one, and we will use the fully parameterized model `F`. In every optimization loop, the current temperature will be divided by 1.001.
```
./HyperLAU data.txt data 0 1 -1 1.001
```

- **name of input file** Name of the file that contains the input data, including possible extensions like `.txt`. HyperLAU expects as an input a textfile containing a list of ancestor and descendant states separated by a blank space, for example `01? 011`. Both states are encoded by binary strings, but can contain one or more `?` to mark missing or uncertain data. Every line is considered as a sample independent of the others. For using cross-sectional data, just set all ancestor states to the zero-string.
- **label for output files** HyperLAU will output several text-files. With this input parameter, you can specify the basis of the names of all these outputs. For every run, you will get the two output files `best_likelihood_[name of the output file].txt` and `transitions_[name of the output file].txt`. If the number of bootstrap resamples is specified as $>0$, you will also get two additional files called `mean_[name of the output file].txt` and `sd_[name of the output file].txt`.
- **number of bootstrap resamples** Number of resamples to simulate. If specified as 0, no bootstrapping will be done. Optional argument. Default value: 0.
- **number of random seed** Specifies the integer random seed to use for simulations. Optional argument. Default value: 1.
- **model structure** Has to be one of the following integers: $-1,1,2,3,4$. With this input parameter you choose which model, i.e. what degree of allowed interaction between the features should be used. Model $-1$ corresponds to the model of arbitrary dependencies, where all combinations of features can influence each other. In the article introducing HyperLAU it is labeled as $F$ for "full". In model $1$, every feature occurs with a fixed rate, independent of other features already obtained. In model $2$ there are already pairwise interactions allowed, and in model $3$ and $4$, pairs or triples can influence the probability of features to occur next, respectively. Optional argument. Default value: -1.
- **annealing rate** This parameter specifies how fast the temperature in the Simulated Annealing Process should be decreased. After every optimization loop, the current temperature is devided by this parameter: `temp = temp/rate`. This parameter should be a double $>1$. Optional argument. Default value: 1.001.

## Output
HyperLAU does always output the two text files `best_likelihood_[name of the output file].txt` and `transitions_[name of the output file].txt`. If the number of bootstrap resamples is chosen to be bigger than zero, it outputs two additional files `mean_[name of the output file].txt`and `sd_[name of the output file].txt`, as well as a file `bootstrap_[bootstrap number]_[name of the outpu file].txt` for every bootstrap resample.
- **best_likelihood_[..].txt** In this output you can find the intermediate results for the currently best log-likelihood. Every line corresponds to a certain loop in the Simulated Annealing optimization and contains the log-likelihood that is stored as the best log-likelihood until the end of this loop. This file helps to track the progression of the optimization and can give some insight if the optimiziation process was running long enough for the log-likelihood to converge.
- **transitions_[...].txt** This file contains the four columns `From`, `To`, `Probability` and `Flux`. In the columns `From` and `To` are the origin and destination node of a transition, encoded in the integer representation of the binary string. The column `Probability` contains the transition probability that HyperLAU assigned to this edge. The last column reports the calculated flux for the corresponding edge, which reflects the proportion of the trajectories that makes this transition.
- **mean_[...].txt** This file contains four columns. The firts two ones give the origin and destination of the considered transition in the same way as in the file `transitions_[...].txt`. The column `Probabiliy` contains the mean of the transition probabilities for all bootstrap resamples that were runned, the column `Flux` the mean of all fluxes accordingly.
- **sd_[...].txt** This file contains four columns. The firts two ones give the origin and destination of the considered transition in the same way as in the file `transitions_[...].txt`. The column `Probability` contains the standard deviation of the transition probabilities for all bootstrap resamples that were runned, the column `Flux` the standard deviation of all fluxes accordingly.
- **bootstrap_[...]** These files contain the four columns `From`, `To`, `Probability` and `Flux`. In the columns `From` and `To` are the origin and destination node of a transition, encoded in the integer representation of the binary string. The column `Probability` contains the transition probability that HyperLAU assigned to this edge, based on the considered Bootstrap resample. The last column reports the corresponding calculated flux for the edge, which reflects the proportion of the trajectories that makes this transition.
  
## Bash Script
You can also find the bash script `run_examples.sh` in the repository. When executing the script, you will run example calculations on all data sets provided. In order to do this the folder `data`has to be sotred in the same direction as the bash script itself. Then you can make the bash script executable (if not already done) and run it, for example with
```
chmod +x run_examples.sh
./run_examples.sh
```

## R Scripts
In this repository, we also provide a collection of R scripts for plotting the results of HyperLAU, and for preparing or manipulating the input data. 

### random_qm.R
This script can be used to insert a certain amount of randomly uniform distributed "?" into an input dataset, as we did for the different case studies in the article introducing HyperLAU. The only required R library for running this script is `stringr`. Running the first part of the script inserts "?" in the whole dataset, i.e. at all positions in the string, whereas the second part does this only for a certain specified feature. In both cases one has to specify the input data set in the `read.table()` function, as well as the length of the binary strings `L`, the probability with which a position in the binary strings should be replaced by a "?" (this is done via the parameter `threshold`) and the name under which the manipulated data set should be saved (in the `write.table()` function). In the second case, when you want the "?" to be only in one certain feature, you have to specify the number of the position of this feature in the binary string (parameter `feature`). 


### plot_HyperLAU_results.R
This script gives examples how to create the different types of plots we presented in the preprint. The labels for the simulations for which plots are required must be specified in this R file. The corresponding functions are contained in the file `plotting_functions.R`, so make sure that this file is saved in the working directory. Then the different commands for creating the plots can directly be executed here in this file. 

### plotting_function.R
Contains all the functions that are needed to create the different plots we presented in the article. For examples how to execute them, see `plot_HyperLAU_results.R`.

### prepare_and_plot_phylogeny.R
This script was used to create and draw the phylogeny of the whole tuberculosis dataset (including uncertainties) that was used in the paper. It requires some source code from HyperTraPS that can be found here https://github.com/StochasticBiology/hypertraps-ct?tab=readme-ov-file. Make sure that you saved the files `hypertraps.R`, `hypertraps.c` and `hypertraps-r.cpp` in the same directory as the `prepare_and_plot_phylogeny.R`script.

## Data
In the `data` folder in this repository you can find the data sets we used in the article introducing HyperLAU. 

The file `first_toyexample.txt` contains the very small example with string length three, we discussed first in the article. 

In the subfolder `second_toyexample` there is the full dataset `full.txt` (all positions in all datapoints are specified as $0$ or $1$) we used for our second toy example among the ones where we inserted with a probability of $0.4$ uncertainty markers ("?") into a certain feature. 

The subfolder `tb` contains both datasets we used for the tuberculosis case studies. `tb_data_10.txt` contains the phyloegeny reconstructed by considering only the datapoints of string length ten containing no missing states. Into this dataset, we inserted uncertainty markers ("?") with a probability of $0.5$ at all positions, which is the dataset `tb_data_10_qm50.txt`. 
