#require(Rcpp)
require(ggplot2)
require(ggpubr)
require(ggraph)
require(igraph)
require(stringr)
require(dplyr)
#require(phangorn)
#require(phytools)
#require(ggtree)

file1 = "c4-curated"




DecToBin <- function(x, len) {
  s = c()
  for(j in (len-1):0)
  {
    if(x >= 2**j) { s=c(s,1); x = x-2**j } else { s=c(s,0)}
  }
  return(paste(s, collapse=""))
}

BinToDec <- function(state) {
  this.ref = 0
  for(j in 1:length(state)) {
    this.ref = this.ref + state[j]*(2**(length(state)-j))
  }
  return(this.ref)
}

create_plot <- function(file){
tblabels = c("1", "2", "3", "4", "5", "6", "7", "8", "9","10")
df = read.table(paste(c("transitions_", file, ".txt"), collapse = ""), header=TRUE)
df.u = read.table(paste(c("sd_", file, ".txt"), collapse = ""), header = TRUE)
names(df.u)[4] = "sd_value"
df <- df %>%
  left_join(df.u %>% select(From, To, sd_value), by = c("From", "To"))
df$sd_value[is.na(df$sd_value)] <- 0
df$frombin = sapply(df$From, DecToBin, len=10)
df$tobin = sapply(df$To, DecToBin, len=10)
df$change = 0
for(i in 1:nrow(df)) {
  src = strsplit(df$frombin[i], split="")[[1]]
  dest = strsplit(df$tobin[i], split="")[[1]]
  ref = which(src != dest)
  df$change[i] = ref
}
df = df[df$Flux > 0.015,]
gdf = data.frame(From=df$frombin, To=df$tobin, Flux=df$Flux, Change=df$change, sd.value = df$sd_value)
trans.g = graph_from_data_frame(gdf)
bs = V(trans.g)$name
V(trans.g)$binname = bs
layers = str_count(bs, "1")
this.plot =  ggraph(trans.g, layout="sugiyama", layers=layers) + 
  geom_edge_link(aes(edge_width=Flux, label=tblabels[Change], edge_color = sd.value/Flux), alpha=0.75,label_size=2.5,check_overlap = TRUE) +  
  scale_edge_width(limits=c(0,NA)) + # scale_edge_alpha(limits=c(0,NA)) +
  scale_edge_color_gradientn(name = "CV",colors = c("blue", "pink","red"),values = c(0,0.2,1), limits = c(0,5))+
  theme_graph(base_family="sans") #aes(label=bs)) + theme_graph() 
}

plotB = create_plot(file1)
print(plotB)


sf = 4
png("mro_0.01_c4_0.015.png", width = 800*sf, height = 400*sf, res = 75*sf)
ggarrange(plotB,plotC,labels = c("A","B"))
dev.off()
