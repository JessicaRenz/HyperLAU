source("plotting_functions.R")


# 1. plots the learned fluxes over a specified threshold embedded in the full hypercube

#input that needs to be specified
label = "toy5_m-1"  #f.e. for the HyperLAU output with "transitions_tb_data_10.txt", the label will be "tb_data_10"
L = 3  #string length
thresh = 0.05 # specified threshold, f.e. when set to 0.05, all fluxes that are bigger than 0.05 will be plotted. 
sd_yn = "Y" #do you want to plot the standard deviation derived from bootstrap results? to be specified with "Y" or "N"

plot.1 = plot_embedded_hypercube(label,L,thresh,sd_yn)
print(plot.1)



#2. creates the node to node graph with nodes labeled by the antibiotics considered in the tuberculosis example

#input that needs to be specified
file = "transitions_toy5_m-1.txt"  #name of the file containing the transition probabilities
L = 3 #string length
labels = c("INH","RIF","PZA")  #vector containing the labels of the feature in the corresponding order"

plot.2 = create_plot_node_to_node(file,L,labels)
plot.2



#3. plots the trajectories of two different likelihood progressions

#input that needs to be specified
file.1 = "best_likelihood_tb_data_red.txt"  #name of the first file containing a tracked best log-likelihood
file.2 = "best_likelihood_tb_data_10_qm50.txt" #name of the second file containing a tracked best log-likelihood
legend.1 ="tb_reduced" # label for the trajectory of the first file
legend.2 = "tb_reduced_50 % ?" # label for the trajectory of the second file
label = "TB data, L = 10" # title of the plot

plot.3 = plot_likelihood(file.1,file.2,legend.1,legend.2,label)
plot.3



# 4. Visualization of the features in the input data set

# input that needs to be specified
label = "tb_10"
L = 10
features = c("INH","RIF","PZA","EMB","STR","AMI","CAP","MOX","OFL","PRO")

plot.4 = data_vis(label,L,features)
plot.4



# 5. Visualization of the fluxes in a fixed structure

# input that needs to be specified
#note: f.e. for the HyperLAU output with "transitions_tb_data_10.txt", the label will be "tb_data_10"
label1 = "toy5_m-1" #label dataset one
label2 = "toy5_m1" #label dataset two
sd_yn = "Y" #do you want to plot the standard deviation derived from bootstrap results? to be specified with "Y" or "N"

plot.5 = skeleton_plot(label1,label2,sd_yn)
plot.5