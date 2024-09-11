require(ggplot2)
require(ggraph)
require(igraph)
require(stringr)
require(ggpubr)



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


#plots the learned transition above a certain threshold embedded in the full hypercube
plot_embedded_hypercube <-function(label,L,thresh){
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


#creates the node to node graph with nodes labeled by the antibiotics considered in the tuberculosis example
create_plot_node_to_node <- function(file){
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

#plots the trajectories of six different likelihood progressions
plot_likelihood <- function(df.1,df.2,df.3,df.4,df.5,df.6){
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
    scale_color_manual(values = c("feature 1" = "red", "feature 2" = "blue", "feature 3" = "#33CCFF", "feature 4" = "#FF0099", "feature 5" = "#336666", "feature 6" = "green"))

    plot(res)
}