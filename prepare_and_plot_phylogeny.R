## example of evolutionary accumulation modelling (EvAM) for multidrug resistance (MDR)
# uses HyperTraPS-CT to explore tuberculosis MDR evolution

##########
##### Import libraries

library(readxl)       # to read Excel datafile
library(dplyr)        # for data wrangling
library(tidyr)

# get this from https://github.com/StochasticBiology/hypertraps-ct
source("hypertraps.R")

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

##########
##### Import and curate data

# get profiles of drug resistance/susceptibility
# Supplementary Table 4 of Casali et al., https://www.nature.com/articles/ng.2878.s3
# https://static-content.springer.com/esm/art%3A10.1038%2Fng.2878/MediaObjects/41588_2014_BFng2878_MOESM35_ESM.xls
o.df = read_excel("41588_2014_BFng2878_MOESM35_ESM.xls")

# extract isolate ID and resistance profiles to our ten drugs (discarding mutation info)
# remove any incomplete profiles and recast as binary strings
col.interest = c("Isolate", "INH", "RIF", "PZA", "EMB", "STR", "AMI", "CAP", "MOX", "OFL", "PRO")
df = o.df[,which(colnames(o.df) %in% col.interest)]
missing.rows = unique(which(df == ".", arr.ind=TRUE)[,1])
#final.df = df[-missing.rows,]
final.df = df %>%
  mutate(across(-Isolate, ~ case_when(
    . == "R" ~ "1",
    . == "S" ~ "0",
    . == "." ~ "?")
  )) %>%
  unite("concatenated", -1, sep = "", remove = TRUE) 
final.df = as.data.frame(final.df)
final.df$Isolate = as.character(final.df$Isolate)
final.df$concatenated = as.character(final.df$concatenated)
my.data = final.df 
colnames(my.data)[1] = "label"

# get phylogeny linking isolates
# Supplementary Data Set 1 of Casali et al., https://www.nature.com/articles/ng.2878.s3
# https://static-content.springer.com/esm/art%3A10.1038%2Fng.2878/MediaObjects/41588_2014_BFng2878_MOESM34_ESM.txt
my.tree = read.tree("41588_2014_BFng2878_MOESM34_ESM.txt")
my.rooted.tree = my.tree

# check for missing data in both barcodes and tips... from curate.tree
match.set = match(my.data$label, my.rooted.tree$tip.label)
if(any(is.na(match.set))) {
  message("Found observations that didn't correspond to tips of the tree!")
  my.data = my.data[-which(is.na(match.set)),]
  match.set = match(my.data$label, my.rooted.tree$tip.label)
}

if(any(duplicated(my.data$label))) {
  message("Duplicates in observation set!")
}

# prune tree to include only those tips in the barcode dataset
my.tree = drop.tip(my.rooted.tree,
                   my.rooted.tree$tip.label[-match.set])

final.df = my.data

# assign arbitrary labels to nodes for bookkeeping
my.tree$node.label = as.character(1:my.tree$Nnode)
tree.labels = c(my.tree$tip.label, my.tree$node.label)
tree.states = rep("", length(tree.labels))

##########
## assign states (with uncertainty) to tips and then nodes of tree

# populate tree states with states known from tip data
df = data.frame()
for(i in 1:nrow(final.df)) {
  this.line = final.df[i,]
  df = rbind(df, data.frame(label=this.line[1], state=this.line[length(this.line)]))
  ref = which(tree.labels == as.character(this.line[1]))
  if(length(ref) > 0) {
    tree.states[ref] = as.character(this.line[2])
  }  
}

# first, reconstruct full tree
trans.df = data.frame()
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
      print(paste("Tree", 1, ": populating", ref, "based on kids", kids[1], "and", kids[2], collapse = " "))
      tree.states[ref] = make.ancestor(tree.states[kids[1]], tree.states[kids[2]])
      trans.df = rbind(trans.df, data.frame(tree=1, b = ref, a = kids[1], before=tree.states[ref], after=tree.states[kids[1]]))
      trans.df = rbind(trans.df, data.frame(tree=1, b = ref, a = kids[2], before=tree.states[ref], after=tree.states[kids[2]]))
      change = TRUE
    }
  }
}

###########
### pull set of independent transitions 

# now, consider which transitions we need to drop from this set to make sure we're dealing with indepedent behaviours

# first, consider only those transitions to a tip state from its immediate ancestor...
prune.df = trans.df[trans.df$a < length(my.tree$tip.label),]
# ... and those that actually contain information (i.e. before and after states are different)
prune.df = prune.df[prune.df$before != prune.df$after,]

unsafe = prune.df
pairs = data.frame()
# go through every remaining transition
for(i in 1:nrow(unsafe)) {
  this.safe = TRUE
  if(nrow(pairs) > 0) {
    # if the ancestral *node* (not state) is already in our working set, ignore this transition
    # (avoid double counting P(0->a))
    if(unsafe$b[i] %in% pairs$b) {
      this.safe = FALSE
    }
    # go through transitions in our working set, comparing current transition with each
    for(j in 1:nrow(pairs)) {
      a1 = unsafe$b[i]
      a2 = pairs$b[j]
      # if the state of the ancestors' ancestor isn't compatible with 0^L, ignore this transition
      a12 = getMRCA(my.tree, c(a1, a2))
      if(grepl("1|\\?", tree.states[a12])) {
        this.safe = FALSE
      }
    }
  }
  # if we're not coupled to any existing transitions in the working set, add to working set
  if(this.safe == TRUE) {
    pairs = rbind(pairs, unsafe[i,])
  }
}

# confirm that we're only looking at ancestors of transition pairs that are compatible with 0^L
a.s = c()
ref = 2
for(i in 1:nrow(pairs)) {
  if(i != ref) {
    a.s = c(a.s, tree.states[getMRCA(my.tree, c(pairs$b[ref], pairs$b[i]))])
  }
}
a.s

# output to file
write.table(pairs[,c(4,5)], "independent-transitions_all0.txt", quote=FALSE, row.names = FALSE)

#### IGJ addition Apr 28 for plotting
require(ggtree)
require(ggtreeExtra)

plotHypercube.curated.unc.tree = function(tree.set, 
                                          scale.fn = geom_treescale(y=20, linesize=3, width =0.01),
                                          names = FALSE,
                                          font.size=4) {
  data.m = tree.set$data[1:length(tree.set$tree$tip.label), 2:ncol(tree.set$data)]
  data.m[which(data.m=="?", arr.ind = TRUE)] = 0.5
  data.m = apply(data.m, c(1,2), as.numeric)
  rownames(data.m) = tree.set$data[,1]
  
  g.core = ggtree(tree.set$tree) + scale.fn
  if(names == TRUE) {
    g.core = ggtree(tree.set$tree) + scale.fn + geom_tiplab(size=3, alpha=0.8, nudge_y=0.4, hjust=1)
  } else {
    ggtree(tree.set$tree) + scale.fn
  }
  this.plot = gheatmap(g.core, data.m, low="white", high="#000000",
                       colnames_angle=90, hjust=0, font.size=font.size) +
    theme(legend.position="none")
  
  return(this.plot)
}

tmp.df = strsplit(tree.states, split="")
new.df = data.frame()
for(i in 1:length(my.tree$tip.label)) {
  new.df = rbind(new.df, matrix(tmp.df[[i]], ncol=10, nrow=1))
}
new.df = cbind(my.tree$tip.label, new.df)
colnames(new.df) = col.interest
tree.plot = list()
tree.plot[["tree"]] = my.tree
tree.plot[["data"]] = new.df

plotA = plotHypercube.curated.unc.tree(tree.plot) + scale_fill_gradientn(colors = c("white", "red", "black"))

sf = 4
png("phylogeny.png", width = 800*sf, height = 400*sf, res = 75*sf)
ggarrange(print(plotA),labels = "A")
dev.off()