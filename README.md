# KINDA DATE
This is a script for loci selection for divergence dating balancing criteria of clock-likeness and information content. It is suitable for unrooted locus trees. It requires `ape` and `phytools` (for plotting trees, optional). As input it requires a reference (species) tree, locus trees directory, and a regex pattern for locus tree names. See code comments for more details. The concept is similar to ["SortaDate"](https://github.com/FePhyFoFum/SortaDate), described [here](https://doi.org/10.1371/journal.pone.0197433), but this script has fewer dependencies, can be used with unrooted trees, and correctly weighs the three locus selection criteria. To run `kinda_date`, you will need 1) gene trees (best in `NEWICK` format) for each of your loci/genetic markers, 2) reference species tree, also in `NEWICK`.

## Installation
There is no installation required. You just need to copy the script using `git` or simply downloading it. You will need to have `R` installed with packages `ape` and `phytools`.

## Interface
The script has several hard-coded variables that need to be set before use:
```
trees_dir <- file.path("./trees/")
```
Change `./trees/` to the relative or absolute path to the directory with your gene trees.
```
trees_files <- dir(path=trees_dir, pattern="*.treefile")
``` 
Change `*.treefile` to whatever extension your gene tree files have. For example, if your gene trees are named `locus1.tre`, you want to change this to `*.tre`.
```
tree_regex <- "(uce-[0-9]+).treefile"
```
This line specifies the regular expression to extract the locus name from your tree file names. The way this one is set up, it will extract the number of locus `uce-12345` from tree file called `uce-12345.treefile`
```
reference <- read.tree("consensus.tre")
```
Here `consensus.tre` refers to the path to your reference tree file. This can be a species tree you will compare your gene trees against for similarity. Change it so that it refers to your file.
```
tree_plots_dir <- "./tree_plots/"
``` 
`kinda_date` will plot each of your gene trees so that you can examine them visually. `./tree_plots/` is the path of the directory in which these plots will be placed. Adjust accordingly.
```
weight_clock <- 0.5
weight_brlength <- 0.3
weight_rf <- 0.2
```
`kinda_date` will compute three gene tree characteristics: 

1) Clock-likeness, defined as minimum coefficient of variation of root to tip distances. This is a proxy of how "clock-like" each gene tree is. Weight is this criterion is set with `weight_clock`.

2) Average branch length. Higher is better, as we want to select loci with some information content. Weight is set with `weight_brlength`.

3) Robinson-Foulds distances to the reference tree. Lower is better, and this indicates how similar each gene tree is to your species tree. Weight is set with `weight_rf`.

Once you set all the variables correctly, you can run `kinda_date` from the command line with `Rscript kinda_date.R`. You will see values printed for each gene tree, first indicating whether the trees have been successfully unrooted (if you see all are `FALSE`, you can ignore the warning about some trees being rooted), then showing progress computing coefficients of variation of root to tip distances. Finally, the script will plot all your gene phylogenies in the directory you specified.

The results will be printed into a file called `tree_props_table.csv`. You can sort it by the `Weighted_Sum_Ranks` column to find your "best" genes. Lower numbers are better. Good luck and have fun!