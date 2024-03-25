# kinda_date.R (thanks to Zach Griebenow for inspiring the name)
# Verison 2024-03-25
# by Marek Borowiec


# disable scientific notation
options(scipen=999)

# load needed libraries
library("ape")
library("phytools")

trees_dir <- file.path("./")
trees_files <- dir(path=trees_dir, pattern="*.tre")
tree_regex <- "(uce-[0-9]+).tre"

### AVERAGE BRANCH LENGTHS ###

# This takes a Newick tree with branch lengths
# and returns the number of tips for each tree,
# and calculates average branch length.

Br_length.trees <- function(file) {
  
  # read the phylogenetic tree
  tree <- read.tree(paste(trees_dir, file, sep=""))
  # gets number of tips
  no_tips <- length(tree$tip.label)
  # calculate avg branch length
  avg_br_length <- mean(tree$edge.length)
  # get UCE number from filename
  locus_no <- sub(tree_regex, "\\1", perl=TRUE, x=file)
  return(c(locus_no, avg_br_length))
  
}

### CLOCKLIKENESS ###

# define a function to calculate coefficient of variation
# of root to tip distances for the root minimizing the coefficient

Clocklikeness <- function(file) {
  
  # get UCE number from filename
  locus_no <- sub(tree_regex, "\\1", perl=TRUE, x=file)
  # read in tree
  tree <- read.tree(paste(trees_dir, file, sep=""))
  # record coefficient for all possible outgroups
  CV <- c()
  taxa <- tree$tip.label
  for (taxon in taxa) {
    # root tree
    rooted_tr <- root(phy=tree, outgroup=taxon, resolve.root=T)
    # get matrix diagonal of phylogenetic variance-covariance matrix
    # these are your distances from root
    root_dist <- diag(vcv.phylo(rooted_tr))
    std_dev_root_dist <- sd(root_dist)
    mean_root_dist <- mean(root_dist)
    CV <- c(CV, (std_dev_root_dist/mean_root_dist)*100)
  }
  # get lowest CV
  minCV <- min(CV)
  print(c(locus_no, minCV))
  return(c(locus_no, minCV))
}

### ROBINSON-FOULDS DISTANCES ###

# define a function to compute Robinson-Foulds distance (topology only)
# between gene tree and reference tree
# this assumes that reference tree has full set of taxa and gene tree a subset
# reference tree is pruned to match the gene tree prior to RF comparison

RF_distance <- function(file, reference) {
  # get UCE number from filename
  locus_no <- sub(tree_regex, "\\1", perl=TRUE, x=file)
  # read in tree
  tree <- read.tree(paste(trees_dir, file, sep=""))
  # unroot and resolve polytomies in gene tree
  mutree <- unroot(multi2di(tree))
  print(c(locus_no, is.rooted(mutree)))
  # find taxa missing in gene tree
  common_taxa <- intersect(mutree$tip.label, reference$tip.label)
  missing_taxa <- setdiff(union(mutree$tip.label, reference$tip.label), common_taxa)
  # prune missing taxa from reference tree for Robinson-Foulds comparison
  dreference <- drop.tip(reference, missing_taxa)
  # compute Robinson-Foulds distance using Penny and Hendy 1985 method
  # this does not consider branch lengths
  RF_dist <- dist.topo(mutree, dreference, method="PH85")
  #print(c(locus_no, RF_dist))
  return(c(locus_no, RF_dist))
}

### PLOTTING GENE TREES ### 

# This will save png files of png plots
# of all unrooted phylogenetic trees with tip labels
# and support values--useful for checking ranking results
# this is calibrated for a large phylogeny of ~300 taxa
# cex and width/height should be adjusted for other datasets

Plot_trees <- function(file) {
  
  # read the phylogenetic tree
  tree <- read.tree(paste(trees_dir, file, sep=""))
  # root at midpoint
  rtree <- midpoint.root(tree)
  # extract plot name (locus) from file name 
  plot_name <- sub(tree_regex, "\\1", perl=TRUE, x=file)
  # open png file
  png(file=paste(tree_plots_dir, plot_name, "-tree.png", sep=""), width=1800, height=3600)
  plot.phylo(rtree, show.node.label=T, cex=0.7)
  # give title as locus number
  title(main=plot_name)
  # close png file
  dev.off()
  
}


# note reference tree filename needs to be supplied
reference <- read.tree("ponerinae-792t-spruce-75p-iqtree-swscmerge-mfp_v2_v2.tre")
umreference <- unroot(multi2di(reference))
rf_distances <- lapply(trees_files, RF_distance, reference=umreference)
# loop over all files
br_lengths <- lapply(trees_files, Br_length.trees)
br_lengths <- data.frame(matrix(unlist(br_lengths), nrow=(length(br_lengths)), byrow=T))
colnames(br_lengths) <- c("Locus", "Average_branch_length")
cv_clocklikeness <- lapply(trees_files, Clocklikeness)
cv_clocklikeness <- data.frame(matrix(unlist(cv_clocklikeness), nrow=(length(cv_clocklikeness)), byrow=T))
colnames(cv_clocklikeness) <- c("Locus", "Clocklikeness")
rf_distances <- data.frame(matrix(unlist(rf_distances), nrow=(length(rf_distances)), byrow=T))
colnames(rf_distances) <- c("Locus", "RF_distance")
# putting together all the data
dtemp <- merge(cv_clocklikeness, br_lengths, by="Locus")
all_loci_stats <- merge(dtemp, rf_distances, by="Locus")
# convert characters to numbers
char_columns <- c("Clocklikeness", "Average_branch_length", "RF_distance")
all_loci_stats[, char_columns] <- apply(all_loci_stats[, char_columns], 2, as.numeric)
# add ranks
# default is smallest value = smallest (best) rank
# this needs to be reversed to smallest value = highest rank to properly weight
all_loci_stats$Clocklikeness_rank <- rank(all_loci_stats$Clocklikeness)
all_loci_stats$RF_rank <- rank(all_loci_stats$RF_distance)
# high value for average branch length = low (better) rank (because it indicates signal)
all_loci_stats$Br_length_rank <- rank(-all_loci_stats$Average_branch_length, ties.method="min")
all_loci_stats$Sum_Ranks <- rowSums(all_loci_stats[, c("Clocklikeness_rank", "Br_length_rank", "RF_rank")])
# define weights for each rank
weight_clock <- 0.5
weight_brlength <- 0.3
weight_rf <- 0.2
all_loci_stats$Weighted_Sum_Ranks <- NA  # Initialize the new column
for (i in 1:nrow(all_loci_stats)) {
  all_loci_stats$Weighted_Sum_Ranks[i] <- sum(all_loci_stats[i, "Clocklikeness_rank"] * weight_clock, 
  	all_loci_stats[i, "Br_length_rank"] * weight_brlength,
  	all_loci_stats[i, "RF_rank"] * weight_rf
  	)
}
write.csv(all_loci_stats, "tree_props_table.csv")
# create directory for tree plots
dir.create("./tree_plots")
tree_plots_dir <- file.path(getwd(), "tree_plots/")
#loop over all files to plot locus trees
lapply(trees_files, Plot_trees)
