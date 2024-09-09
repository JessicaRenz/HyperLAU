library(phytools)
library(phangorn)

input.trees = "outtree_tb_9"
input.data = "tb_9"
output.data = "phylo_tb_9.txt"

# reconstruct an ancestral state based on two specified daughters
# need to generalise to more daughters
make.ancestor = function(s1, s2) {
  ss1 = strsplit(s1, split="")[[1]]
  ss2 = strsplit(s2, split="")[[1]]
  news = rep("0", length(ss1))
  # truth table
  # s1  000111???
  # s2  01?01?01?
  # new 00001?0??
  for(i in 1:length(ss1)) {
    if(ss1[i]=="1") {
      if(ss2[i]=="1") { news[i] = "1" }
      if(ss2[i]=="?") { news[i] = "?" }
    } else if(ss2[i] == "1" && ss1[i] == "?") { 
      news[i] = "?" 
    } else if(ss1[i]=="?" && ss2[i] == "?") { news[i] = "?" }
  }
  return(paste(news, collapse=""))
}

# read our set of trees
my.trees = read.tree(input.trees)
trans.df = data.frame()
if("multiPhylo" %in% class(my.trees)) {
  n.trees = length(my.trees)
} else {
  n.trees = 1
}


for(treeref in 1:length(my.trees)) {
  # pick out one
  if("multiPhylo" %in% class(my.trees)) {
    my.tree = my.trees[[treeref]]
  } else {
    my.tree = my.trees
  }
  # assign arbitrary labels to nodes for bookkeeping
  my.tree$node.label = as.character(1:my.tree$Nnode)
  tree.labels = c(my.tree$tip.label, my.tree$node.label)
  tree.states = rep("", length(tree.labels))
  
  # read in data
  my.data = readLines(input.data)
  these.lines = strsplit(my.data, split=" ")
  df = data.frame()
  # loop through lines and pull label and state from each line
  # our goal is to populate "tree.states", not just with these tips, but also the nodes
  for(i in 2:length(these.lines)) {
    this.line = these.lines[[i]]
    df = rbind(df, data.frame(label=this.line[1], state=this.line[length(this.line)]))
    ref = which(tree.labels == this.line[1])
    if(length(ref) > 0) {
      tree.states[ref] = this.line[length(this.line)]
    }  
  }
  # assume we've got work to do reconstructing ancestors
  change = TRUE
  # loop until we're done
  while(change == TRUE) {
    change = FALSE
    # which nodes don't yet have states
    todo = which(tree.states == "")
    if(length(todo) == 0) {
      change = FALSE
      break
    }
    # loop through empty nodes
    for(i in 1:length(todo)) {
      # get children of this node, and figure out if their states are specified
      ref = todo[i]
      kids = Children(my.tree, ref)
      kid.empty = which(tree.states[kids] == "")
      if(length(kids) == 2 && length(kid.empty) == 0) {
        # if the childrens' states are specified, reconstruct this ancestral state and add transitions to the dataframe
        print(paste("Tree", treeref, ": populating", ref, "based on kids", kids[1], "and", kids[2], collapse = " "))
        tree.states[ref] = make.ancestor(tree.states[kids[1]], tree.states[kids[2]])
        trans.df = rbind(trans.df, data.frame(tree=treeref, before=tree.states[ref], after=tree.states[kids[1]]))
        trans.df = rbind(trans.df, data.frame(tree=treeref, before=tree.states[ref], after=tree.states[kids[2]]))
        change = TRUE
      }
    }
  }
}

keep <- c("before","after")
trans.df = trans.df[keep]

names(trans.df) <- NULL
write.table(trans.df, output.data, sep = " ", row.names = FALSE, quote = FALSE)
