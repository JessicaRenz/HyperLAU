# L: number of features
# prop: proportion of questionmarks in %
# average_loglik_time: average time of all likelihood calculations during the SA process
# runtime: time of the whole SA process
library(ggplot2)
library(tidyverse)
library(zoo)
library(readr)
library(ggpubr)

loglikelihoods = read.csv2("loglikelihoods.csv")
# here the columns L3-L10 contain the final loglikelihood

df_long <- pivot_longer(
  loglikelihoods,
  cols = -prop,
  names_to = "L",
  values_to = "loglik"
)



df_long$loglik <- as.numeric(df_long$loglik)

p = ggplot(df_long,
       aes(x = prop, y = loglik, colour = L)) +
  geom_point(size = 2) +
  geom_point() +
  scale_y_continuous(n.breaks = 5) +
  theme_minimal() +
  labs(
    x = "Proportion of '?' in %" ,
    y = "Log-Likelihood",
    colour = "number of features"
  )

ggsave(
  "loglikelihood_proportions.png",
  p,
  width = 8,
  height = 6,
  dpi = 600
)



runtime_features = read.csv2("runtime_features.csv")
runtime_features$average_loglik_time = as.numeric(runtime_features$average_loglik_time)
p1 = ggplot(runtime_features, aes(x = L, y = average_loglik_time)) +
  geom_point(size = 2) +
  theme_minimal()+
  scale_y_log10() +
  labs(
    x = "Number of features" ,
    y = "Average time for one loglikelihood calculation"
  ) + theme(axis.title = element_text(size = 13), 
            axis.text  = element_text(size = 13))

ggsave(
  "average_loglik_features.png",
  p1,
  width = 8,
  height = 6,
  dpi = 600
)
figS1= ggarrange(
  p1, p2,
  labels = c("A", "B"),
  ncol = 2
)

ggsave(
  "figS1.png",
  figS1,
  width = 10,
  height = 5.5,
  units = "in"
)
plot(runtime_prop_L3$prop,runtime_prop_L3$average_loglik_time)

runtime_prop_L3 = read.csv2("runtime_prop_L3.csv")
runtime_prop_L4 = read.csv2("runtime_prop_L4.csv")
runtime_prop_L5 = read.csv2("runtime_prop_L5.csv")
runtime_prop_L6 = read.csv2("runtime_prop_L6.csv")
runtime_prop_L7 = read.csv2("runtime_prop_L7.csv")
runtime_prop_L8 = read.csv2("runtime_prop_L8.csv")
runtime_prop_L9 = read.csv2("runtime_prop_L9.csv")
runtime_prop_L10 = read.csv2("runtime_prop_L10.csv")

runtime_prop_L3$average_loglik_time = as.numeric(runtime_prop_L3$average_loglik_time)
runtime_prop_L4$average_loglik_time = as.numeric(runtime_prop_L4$average_loglik_time)
runtime_prop_L5$average_loglik_time = as.numeric(runtime_prop_L5$average_loglik_time)
runtime_prop_L6$average_loglik_time = as.numeric(runtime_prop_L6$average_loglik_time)
runtime_prop_L7$average_loglik_time = as.numeric(runtime_prop_L7$average_loglik_time)
runtime_prop_L8$average_loglik_time = as.numeric(runtime_prop_L8$average_loglik_time)
runtime_prop_L9$average_loglik_time = as.numeric(runtime_prop_L9$average_loglik_time)
runtime_prop_L10$average_loglik_time = as.numeric(runtime_prop_L10$average_loglik_time)

runtime_prop_L3$L <- "L3"
runtime_prop_L4$L <- "L4"
runtime_prop_L5$L <- "L5"
runtime_prop_L6$L <- "L6"
runtime_prop_L7$L <- "L7"
runtime_prop_L8$L <- "L8"
runtime_prop_L9$L <- "L9"
runtime_prop_L10$L <- "L10"

runtime_prop <- rbind(
  runtime_prop_L3,
  runtime_prop_L4,
  runtime_prop_L5,
  runtime_prop_L6,
  runtime_prop_L7,
  runtime_prop_L8,
  runtime_prop_L9,
  runtime_prop_L10
)

p2 = ggplot(runtime_prop,
       aes(x = prop,
           y = average_loglik_time,
           colour = L,
           group = L)) +
  geom_point() +
  scale_y_log10() +
  theme_classic() +
  labs(
    x = "Proportion of '?' in %" ,
    y = "Average time for one loglikelihood calculation in seconds",
    colour = "Number of \n features"
  )+ theme(axis.title = element_text(size = 11), 
           axis.text  = element_text(size = 13),
           legend.title = element_text(size = 13),
           legend.text  = element_text(size = 13))

ggsave(
  "average_loglik_proportions.png",
  p2,
  width = 8,
  height = 6,
  dpi = 600
)

test_lik <- read.csv2(
  "likelihoods_testsets.csv",
  dec = "."
)

testlik_plot = ggplot(test_lik, aes(x = L)) +
  geom_line(aes(y = HyperHMM, colour = "HyperHMM")) +
  geom_line(aes(y = HyperTraPS, colour = "HyperTraPS")) +
  geom_line(aes(y = HyperLAU, colour = "HyperLAU")) +
  geom_line(aes(y = uniform, colour = "uniform")) +
  theme_bw() + 
  labs(
    x = "number of features",
    y = "log-likelihood of the test-dataset",
    colour = "model"
  )

ggsave(
  "loglik_testsets.png",
  testlik_plot,
  width = 8,
  height = 6,
  dpi = 600
)

trainlik = read.csv2("loglikelihoods_trainingsets.csv", dec = ".")
trainlik_plot = ggplot(trainlik, aes(x=L))+
  geom_line(aes(y=HyperTraPS, colour = "HyperTraPS")) + 
  geom_line(aes(y=HyperHMM, colour = "HyperHMM")) + 
  geom_line(aes(y= HyperLAU, colour = "HyperLAU")) + 
  theme_bw() + 
  labs(
    x="number of features",
    y = "log-likelohood of the trainings-dataset",
    colour = "model"
  )

ggsave(
  "loglik_trainingsset.png",
  trainlik_plot,
  width= 8,
  height= 6,
  dpi = 600
)

testlik_qm30 = read.csv2("lik_test_qm30.csv", dec = ".")
testlik_qm30_plot = ggplot(testlik_qm30, aes(x=L)) +
  geom_line(aes(y=HyperHMM, colour="HyperHMM"),
            na.rm=TRUE) +
  geom_point(aes(y=HyperHMM, colour="HyperHMM")) +
  geom_line(aes(y=HyperTraPS, colour="HyperTraPS")) +
  geom_point(aes(y=HyperTraPS, colour="HyperTraPS")) +
  geom_line(aes(y=HyperLAU, colour="HyperLAU")) +
  geom_point(aes(y=HyperLAU, colour="HyperLAU")) +
  geom_line(aes(y=uniform, colour="uniform")) +
  geom_point(aes(y=uniform, colour="uniform")) +
  theme_bw() + 
  labs(
    x="number of features",
    y = "log-likelohood of the test-dataset (learned with 30% '?'s)",
    colour = "model"
  ) +
  theme(axis.title = element_text(size = 13), 
        axis.text  = element_text(size = 11),
        legend.title = element_text(size = 13),
        legend.text  = element_text(size = 12)) 

figS2 = ggarrange(
  testlik_plot,testlik_qm10_plot,testlik_qm20_plot, testlik_qm30_plot,
  ncol = 2,
  nrow = 2,
  labels = c("A", "B", "C", "D"),
  common.legend = TRUE,
  legend = "bottom"
)
ggsave(
  "figS2.png",
  figS2,
  width = 10,
  height = 12
)

ggsave(
  "loglik_test.png",
  testlik_qm30_plot,
  width= 8,
  height= 6,
  dpi = 600
)