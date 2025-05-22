require(ggplot2)
require(ggraph)
require(igraph)
require(stringr)
require(ggpubr)
library(dplyr)
library(ggrepel)



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

plotter = function(g,node_positions,sd_yn){
  # make edits here to style plots as required. the essential feature is that we're reusing the same layout
  if(sd_yn=="Y"){
    graph = ggraph(g, x = node_positions$x, y=node_positions$y) +
      geom_edge_link(aes(edge_width=Flux, edge_alpha=Flux,edge_color=sd_value/Flux)) +
      scale_edge_width(range=c(0,max(E(g)$Flux)/0.2)) +
      scale_edge_alpha(range=c(0,max(E(g)$Flux)/0.2)) +
      scale_edge_color_gradient(name = "CV", low = "blue", high = "red", na.value = "lightgrey",limits = c(0,5)) +
      theme_void() 
  }
  
  if(sd_yn == "N"){
    graph = ggraph(g, x= node_positions$x, y=node_positions$y) +
      geom_edge_link(aes(edge_width=Flux,edge_alpha=Flux)) +
      scale_edge_width(range=c(0,max(E(g)$Flux)/0.2)) +
      scale_edge_alpha(range=c(0,max(E(g)$Flux)/0.2)) +
      theme_void()
  } 
  return(graph)
}

#plots the learned transition above a certain threshold embedded in the full hypercube
plot_embedded_hypercube <-function(label,L,thresh,sd_yn){
  data = read.table(paste(c("transitions_", label, ".txt"), collapse = ""), header = TRUE)
  if(sd_yn == "Y"){
    sd.data = read.table(paste(c("sd_", label, ".txt"), collapse = ""), header = TRUE)
    names(sd.data) = c("From", "To", "SD")
  }
  
  rel.data = data.frame()
  for( i in 1:nrow(data)){
    if (data$Flux[i] >= thresh){
      vec = c(DecToBin(data$From[i],L),DecToBin(data$To[i],L),data$Flux[i])
      rel.data = rbind(rel.data,vec)
    }
  }
  names(rel.data) = c("From", "To", "Flux")
  
  if(sd_yn=="Y"){
    rel.sd.data = data.frame()
    for( i in 1:nrow(sd.data)){
      vec = c(DecToBin(sd.data$From[i],L),DecToBin(sd.data$To[i],L),sd.data$SD[i])
      rel.sd.data = rbind(rel.sd.data,vec)
    }
  
    names(rel.sd.data) = c("From", "To", "SD")
  }
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
  if(sd_yn == "Y"){
    E(graphO)$sd = NA
  }
  edgenames = as_edgelist(graphO)
  for(i in 1:nrow(rel.data)) {
    ref = which(edgenames[,1] == rel.data$From[i] & edgenames[,2] == rel.data$To[i])
    E(graphO)$flux[ref] = as.numeric(rel.data$Flux[i])
    if(sd_yn == "Y"){
      E(graphO)$sd[ref] = 0
    }
  }
  
  if(sd_yn == "Y"){
    for(i in 1:nrow(rel.sd.data)) {
      ref = which(edgenames[,1] == rel.sd.data$From[i] & edgenames[,2] == rel.sd.data$To[i])
      E(graphO)$sd[ref] = as.numeric(rel.sd.data$SD[i])
    }
  }
  
  if(sd_yn=="N"){
    plotted_graph =ggraph(graphO) + 
      scale_edge_color_manual(values=c("0"="grey", "1" = "red")) +
      geom_edge_link(aes(edge_color=factor(skeleton_A), edge_width = flux)) +
      geom_node_text(aes(label=labels), angle=45, hjust=0, size=3,check_overlap = TRUE) +
      theme_graph() + labs( edge_color = "flux > thresh", edge_width = "Prob flux") 
  }
  
  if(sd_yn == "Y"){
    plotted_graph =ggraph(graphO) + 
      geom_edge_link(aes(edge_color=sd/flux, edge_alpha=factor(1-skeleton_A)),show.legend = c(edge_color = TRUE, edge_alpha = FALSE, flux = TRUE)) +
      geom_edge_link(aes(edge_color=sd/flux, edge_alpha=factor(skeleton_A), edge_width = flux),show.legend = c(edge_color = FALSE, edge_alpha = FALSE, flux = TRUE)) +
      geom_node_text(aes(label=labels), angle=45, hjust=0, size=3,check_overlap = TRUE) +
      scale_edge_alpha_manual(values=c("0"=0, "1"=1)) + 
      scale_edge_color_gradient(low = "blue", high = "red", na.value = "lightgrey", limits = c(0,5)) +
      theme_graph() + labs(edge_color = "CV", edge_alpha = NULL, edge_width = "Prob flux") 
  }
  return(plotted_graph)
}



#creates the node to node graph with nodes labeled by the antibiotics considered in the tuberculosis example
create_plot_node_to_node <- function(file, L,labels){
  df = read.table(file, header=TRUE)
  df$frombin = sapply(df$From, DecToBin, len=L)
  df$tobin = sapply(df$To, DecToBin, len=L)
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
    geom_edge_link(aes(edge_width=Flux, edge_alpha=Flux, label=labels[Change]),color="#AAAAFF") +  
    #geom_edge_link(aes(edge_width=Flux, edge_alpha=Flux),color="#AAAAFF") + 
    scale_edge_width(limits=c(0,NA)) + scale_edge_alpha(limits=c(0,NA)) +
    theme_graph(base_family="sans") #aes(label=bs)) + theme_graph() 
  return(this.plot)
}

#plots the trajectories of two different likelihood progressions
plot_likelihood <- function(file.1,file.2,legend.1,legend.2,label){
  
  df.1 = read.table(file.1, header = FALSE)
  df.2 = read.table(file.2, header = FALSE)
  
  df.1$x = 1:nrow(df.1)
  df.2$x = 1:nrow(df.2)

  res = ggplot() + geom_line(data=df.1, aes(x=x, y=V1,color = legend.1), show.legend = TRUE) + 
    ggtitle(label)+
    xlab("Iteration") + ylab("log-likelihood")+
    geom_line(data=df.2,  aes(x=x, y=V1, color = legend.2), show.legend = TRUE) +
    labs(color = NULL)
    
  return(res)
}

#data visualisation

data_vis = function(label, L, features){
  df = read.table(label, stringsAsFactors = FALSE)
  strings = df[2:nrow(df), 2]

  split_strings <- strsplit(strings, "")

  # Create a data frame with the string reference, position, and character
  df_long <- data.frame(
    string_reference = rep(seq_along(strings), sapply(split_strings, length)),
    position = sequence(sapply(split_strings, length)),
    character = unlist(split_strings)
  )
  num = c(1:L)
  tbl = ggplot(df_long, aes(x=position,y=string_reference, fill=character)) + geom_tile()+
      scale_fill_manual(values = c("?" = "red", "0" = "lightgrey","1" = "black"))+
    xlab("Feature") + ylab("Datapoint") + 
    scale_x_continuous(breaks=num,labels= features)
  
  return(tbl)
}


#Visualization of the fluxes in a fixed structure to compare two results

skeleton_plot = function(label1,label2,sd_yn){
  
  dfs = list()
  
  dfs[[1]] = read.csv(paste(c("transitions_",label1,".txt"),collapse=""), sep = " ")
  dfs[[2]] = read.csv(paste(c("transitions_",label2,".txt"),collapse=""), sep = " ")
  
  if(sd_yn == "Y"){
      sd1 = read.csv(paste(c("sd_", label1, ".txt"), collapse = ""),sep = " ")
      sd2 = read.csv(paste(c("sd_", label2, ".txt"), collapse = ""),sep = " ")
    
      names(sd1)[4] = "sd_value"
      names(sd2)[4] = "sd_value"
    
      sd1 <- as.data.frame(sd1)
      sd2 <- as.data.frame(sd2)
    
      dfs[[1]] <- dfs[[1]] %>%
        left_join(sd1 %>% select(From, To, sd_value), by = c("From", "To"))
      dfs[[2]] <- dfs[[2]] %>%
        left_join(sd2 %>% select(From, To, sd_value), by = c("From", "To"))
      
      dfs[[1]]$sd_value[is.na(dfs[[1]]$sd_value)] <- 0
      dfs[[2]]$sd_value[is.na(dfs[[2]]$sd_value)] <- 0
  }
  
  g = list()
  for(i in 1:length(dfs)) {
    g[[i]] = graph_from_data_frame(dfs[[i]])
  }
  
  # create a layout to use across subplots. the "kk" arrangement gives a reasonable view here (IMO)
  graph_layout <- create_layout(g[[1]], layout = "kk")
  node_positions <- graph_layout[, c("x", "y")]
  
  graph1 = plotter(g[[1]],node_positions,sd_yn)
  graph2 = plotter(g[[2]],node_positions,sd_yn)
  
  return(ggarrange(graph1,graph2))
}


pie_charts = function(data_file,L){
  
  data = read.table(data_file)
  trans_0_0 = c(rep(0,10))
  trans_0_1 = c(rep(0,10))
  trans_0_Q = c(rep(0,10))
  trans_1_1 = c(rep(0,10))
  trans_Q_0 = c(rep(0,10))
  trans_Q_Q = c(rep(0,10))
  trans_Q_1 = c(rep(0,10))

  for( i in 1:nrow(data)){
    before = data[i,1]
    after = data[i,2]
    for (j in 1:L){
      letter_1 = substr(before,j,j)
      letter_2 = substr(after,j,j)
      if ((letter_1 == "0")&&(letter_2 == "0")){
        trans_0_0[j] = trans_0_0[j] +1
      }
      if ((letter_1 == "0")&&(letter_2 == "1")){
        trans_0_1[j] = trans_0_1[j] +1
      }
      if ((letter_1 == "0")&&(letter_2 == "?")){
        trans_0_Q[j] = trans_0_Q[j] +1
      }
      if ((letter_1 == "1")&&(letter_2 == "1")){
        trans_1_1[j] = trans_1_1[j] +1
      }
      if ((letter_1 == "?")&&(letter_2 == "0")){
        trans_Q_0[j] = trans_Q_0[j] +1
      }
      if ((letter_1 == "?")&&(letter_2 == "?")){
        trans_Q_Q[j] = trans_Q_Q[j] +1
      }
      if ((letter_1 == "?")&&(letter_2 == "1")){
        trans_Q_1[j] = trans_Q_1[j] +1
      }
    }
  }

  features <- vector("list",L)
  
  for (i in 1:L){
    features[[i]] <- data.frame(category = c("0","1","?","1", "0", "1", "?"),
                value = c(trans_0_0[i], trans_0_1[i], trans_0_Q[i], trans_1_1[i],trans_Q_0[i],trans_Q_1[i],trans_Q_Q[i]),
                From = c("0","0","0","1","?","?","?")) 
    features[[i]] = features[[i]][features[[i]]$value != 0,]
  }
  
  plots = vector("list",L)
  
  for (i in 1:L){
    plots[[i]] = ggplot(features[[i]], aes(x = "", y = value, fill = From))+
                  geom_col(color = "white")+
                  geom_text(aes(x = 1.6,label = category), position = position_stack(vjust = 0.5))+
                  coord_polar(theta = "y")+
                  theme_void()+
                  ggtitle(paste("Feature ", i))+
                  scale_fill_manual(values = c("red","grey","black"))
    
  }


  plotB = do.call(ggarrange,plots)
  return(plotB)
}
