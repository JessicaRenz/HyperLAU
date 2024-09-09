require(ggplot2)
require(ggraph)
require(igraph)
require(stringr)
require(ggpubr)

label_toy5.1 = "toy5_model-1"
label_toy5.2 = "toy5_model1"
label_toy5.3 = "toy5_model2"
L.toy5 = 3

label.ad.m1.1 = "full_model-1"
label.ad.m1.2 = "feature1_40_model-1"

label.ad.2.1 = "full_model2"
label.ad.2.4 = "feature3_40_model2"

L.ad = 6

 # binary to decimal function
  BinToDec <- function(x) {
    sum(2^(which(rev(unlist(strsplit(as.character(x), "")) == 1))-1))
  }

  # decimal to binary function
  DecToBin <- function(x, len) {
    s = c()
    for(j in (len-1):0)
    {
      if(x >= 2**j) { s=c(s,1); x = x-2**j } else { s=c(s,0)}
    }
    return(paste(s, collapse=""))
  }
  
   # decimal to binary function, returning a numerical vector
  DecToBinV <- function(x, len) {
    s = c()
    for(j in (len-1):0)
    {
      if(x >= 2**j) { s=c(s,1); x = x-2**j } else { s=c(s,0)}
    }
    return(s)
  }
 
    
plot_model <-function(label,L){
  

  thresh = 0.05
  
  data = read.table(paste(c("transitions_", label, ".txt"), collapse = ""), header = TRUE)
  sd.data = read.table(paste(c("sd_", label, ".txt"), collapse = ""), header = FALSE)
  names(sd.data) = c("From", "To", "SD")
  
  rel.data = data.frame()
  for( i in 1:nrow(data)){
    if (data$Flux[i] >= thresh){
      vec = c(DecToBin(data$From[i],L),DecToBin(data$To[i],L),data$Flux[i])
      rel.data = rbind(rel.data,vec)
    }
  }
  names(rel.data) = c("From", "To", "Flux")
  
  rel.sd.data = data.frame()
  for( i in 1:nrow(sd.data)){
    vec = c(DecToBin(sd.data$From[i],L),DecToBin(sd.data$To[i],L),sd.data$SD[i])
    rel.sd.data = rbind(rel.sd.data,vec)
  }
  
  names(rel.sd.data) = c("From", "To", "SD")
  
  learned.edges <- rel.data[, 1:2]
  learned.edges <- as.matrix(learned.edges)
  flux.values = as.numeric(rel.data$Flux)
  
  
  # graphA stores the set of edges we want to highlight
  graphA = graph_from_edgelist(learned.edges)
  
  
  # construct full L-hypercube for comparison
  pow2 = 2**((L-1):0)
  am = matrix(ncol=2)
  # produce list of decimal edges
  for(i in 1:(2**L-1)) {
    anc = DecToBinV(i-1, len=L)
    to.1 = which(anc == 0)
    for(j in 1:length(to.1)) {
      desc = i-1+pow2[to.1[j]]
      am = rbind(am, c(i-1, desc))
    }
  }
  
  # convert to graph with binary labels
  ambin = apply(am[2:nrow(am),], c(1,2), DecToBin, len=L)
  graphO <- graph_from_edgelist(as.matrix(ambin),directed = T)
  # this union isn't necessary unless our learned edge set involves multi-step transitions
  graphO = graph.union(graphO, graphA)
  graphO.layers = sapply(V(graphO)$name, str_count, "1")
  
  # figure out which edges in the complete hypercube are those that we found in graph A
  rowsA <- apply(get.edgelist(graphA), 1, paste, collapse = ",")
  rowsO <- apply(get.edgelist(graphO), 1, paste, collapse = ",")
  common_rows_A <- intersect(rowsA, rowsO)
  indexes_A <- which(rowsO %in% common_rows_A)
  
  # set a variable for those graphO edges included in graphA
  E(graphO)$skeleton_A = 0
  E(graphO)[indexes_A]$skeleton_A = 1
  # get the graphO vertices involved in graphA
  skeleton_A_v = unique(as.vector(ends(graphO, indexes_A)))
  # set a name label for these, and a blank label for all others
  V(graphO)$labels = V(graphO)$name
  V(graphO)$labels[!(V(graphO)$name %in% skeleton_A_v)] = "" 
  
  #Add Flux values to the edges of graph0
  E(graphO)$flux = 0
  E(graphO)[indexes_A]$flux = flux.values
  
  # IGJ attempt
  E(graphO)$flux = 0
  E(graphO)$sd = NA
  edgenames = as_edgelist(graphO)
  for(i in 1:nrow(rel.data)) {
    ref = which(edgenames[,1] == rel.data$From[i] & edgenames[,2] == rel.data$To[i])
    E(graphO)$flux[ref] = as.numeric(rel.data$Flux[i])
    E(graphO)$sd[ref] = 0
  }
  for(i in 1:nrow(rel.sd.data)) {
    ref = which(edgenames[,1] == rel.sd.data$From[i] & edgenames[,2] == rel.sd.data$To[i])
    E(graphO)$sd[ref] = as.numeric(rel.sd.data$SD[i])
  }
  
  plotted_graph =ggraph(graphO) + 
    #    geom_edge_link(aes(edge_color=factor(skeleton_A), edge_alpha=factor(skeleton_A))) +
    geom_edge_link(aes(edge_color=sd/flux, edge_alpha=factor(1-skeleton_A)),show.legend = c(edge_color = TRUE, edge_alpha = FALSE, flux = TRUE)) +
    geom_edge_link(aes(edge_color=sd/flux, edge_alpha=factor(skeleton_A), edge_width = flux),show.legend = c(edge_color = FALSE, edge_alpha = FALSE, flux = TRUE)) +
    #     geom_edge_link(aes(edge_width = flux)) +
    geom_node_text(aes(label=labels), angle=45, hjust=0, size=3,check_overlap = TRUE) +
    scale_edge_alpha_manual(values=c("0"=0, "1"=1)) + 
    # scale_edge_color_manual(values=c("0"="grey", "1" = "red")) +
    scale_edge_color_gradient(low = "blue", high = "red", na.value = "lightgrey") +
    theme_graph() + labs(edge_color = "CV", edge_alpha = NULL, flux = "Prob flux") 
    #scale_y_continuous(expand = c(0.4,0)) #+ scale_x_continuous(expand = c(-0.00001,-0.00001))
  return(plotted_graph)
}


graph_toy5.1 = plot_model(label_toy5.1,L.toy5) 
graph_toy5.2 = plot_model(label_toy5.2,L.toy5) 
graph_toy5.3 = plot_model(label_toy5.3,L.toy5) 

graph_ad.m1.1 = plot_model(label.ad.m1.1,L.ad) 
graph_ad.m1.2 = plot_model(label.ad.m1.2,L.ad) 


graph_ad.2.1 = plot_model(label.ad.2.1,L.ad) 
graph_ad.2.4 = plot_model(label.ad.2.4,L.ad) 

print(ggarrange(ggarrange(graph_toy5.1,graph_toy5.2,graph_toy5.3, labels = c("A","",""), nrow = 3),
          ggarrange(ggarrange(graph_ad.m1.1,graph_ad.m1.2,labels = c("B","D"),nrow = 2),
                    ggarrange(graph_ad.2.1,graph_ad.2.4,labels = c("C","E"),nrow = 2),nrow = 1),widths = c(1,3), nrow = 1))

print(ggarrange(graph_ad.m1.1,graph_ad.2.1, labels = c("B","C"), nrow = 1))

sf = 4
png("figure1_a.png", width = 500*sf, height = 300*sf, res = 72*sf)
ggarrange(graph_toy5.1,graph_toy5.2,graph_toy5.3, labels = c("A","",""), nrow = 1)
dev.off()

sf = 4
png("figure1_b.png", width = 500*sf, height = 300*sf, res = 72*sf)
ggarrange(graph_ad.m1.1,graph_ad.2.1, labels = c("B","C"), nrow = 1)
dev.off()

sf = 4
png("figure1_c.png", width = 500*sf, height = 300*sf, res = 72*sf)
ggarrange(graph_ad.m1.2,graph_ad.2.4, labels = c("D","E"), nrow = 1)
dev.off()
