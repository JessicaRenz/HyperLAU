library(ggplot2)

label = "model 2 , L=6"
df.full = read.csv("best_likelihood_full_model2.txt", header=FALSE)
df.1 = read.csv("best_likelihood_feature1_40_model2.txt", header=FALSE)
df.2 = read.csv("best_likelihood_feature2_40_model2.txt", header=FALSE)
df.3 = read.csv("best_likelihood_feature3_40_model2.txt", header=FALSE)
df.4 = read.csv("best_likelihood_feature4_40_model2.txt", header=FALSE)
df.5 = read.csv("best_likelihood_feature5_40_model2.txt", header=FALSE)
df.6 = read.csv("best_likelihood_feature6_40_model2.txt", header=FALSE)

df.full$x = 1:nrow(df.full)
df.1$x = 1:nrow(df.1)
df.2$x = 1:nrow(df.2)
df.3$x = 1:nrow(df.3)
df.4$x = 1:nrow(df.4)
df.5$x = 1:nrow(df.5)
df.6$x = 1:nrow(df.6)

res = ggplot() + geom_line(data=df.2, aes(x=x, y=V1,color = "feature 2"), show.legend = TRUE) + 
  ggtitle(label)+
  xlab("Iteration") + ylab("log-likelihood")+
  geom_line(data=df.1,  aes(x=x, y=V1, color = "feature 1"), show.legend = TRUE) +
  geom_line(data=df.3,  aes(x=x, y=V1, color = "feature 3"), show.legend = TRUE) +
  geom_line(data=df.4,  aes(x=x,y=V1, color = "feature 4" ), show.legend = TRUE)+
  geom_line(data=df.5,  aes(x=x,y=V1, color = "feature 5" ), show.legend = TRUE)+
  geom_line(data=df.6,  aes(x=x,y=V1, color = "feature 6" ), show.legend = TRUE)+
  geom_line(data=df.full,  aes(x=x,y=V1, color = "full " ), show.legend = TRUE)
  scale_color_manual(values = c("feature 1" = "red", "feature 2" = "blue", "feature 3" = "#33CCFF", "feature 4" = "#FF0099", "feature 5" = "#336666", "feature 6" = "green", "full" = "orange"))

plot(res)

png(file= paste(c("full_qm_model2.png"), collapse = ""),
    width=500, height=400)
plot(res)
dev.off()