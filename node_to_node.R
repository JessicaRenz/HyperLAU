require(Rcpp)
require(ggplot2)
require(ggpubr)
require(ggraph)
require(igraph)
require(stringr)
require(phangorn)
require(phytools)
require(ggtree)

file1 = "transitions_tb_data_9_model-1.txt"
file2 = "transitions_tb_data_9_50_model-1.txt"

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
tblabels = c("INH", "RIF", "PZA", "EMB", "STR", "AMI", "CAP", "MOX", "OFL")
df = read.table(file, header=TRUE)
df$frombin = sapply(df$From, DecToBin, len=9)
df$tobin = sapply(df$To, DecToBin, len=9)
df$change = 0
for(i in 1:nrow(df)) {
  src = strsplit(df$frombin[i], split="")[[1]]
  dest = strsplit(df$tobin[i], split="")[[1]]
  ref = which(src != dest)
  df$change[i] = ref
}
df = df[df$Flux > 0.05,]
gdf = data.frame(From=df$frombin, To=df$tobin, Flux=df$Flux, Change=df$change)
trans.g = graph_from_data_frame(gdf)
bs = V(trans.g)$name
V(trans.g)$binname = bs
layers = str_count(bs, "1")
this.plot =  ggraph(trans.g, layout="sugiyama", layers=layers) + 
  geom_edge_link(aes(edge_width=Flux, edge_alpha=Flux, label=tblabels[Change]),color="#AAAAFF") +  
  #geom_edge_link(aes(edge_width=Flux, edge_alpha=Flux),color="#AAAAFF") + 
  scale_edge_width(limits=c(0,NA)) + scale_edge_alpha(limits=c(0,NA)) +
  theme_graph(base_family="sans") #aes(label=bs)) + theme_graph() 
}

plot1 = create_plot(file1)

plot2 = create_plot(file2)

sf = 4
png("tb_full_qm.png", width = 500*sf, height = 300*sf, res = 72*sf)
print(ggarrange(plot1,plot2,labels=c("A","B"),nrow = 1))
dev.off()
