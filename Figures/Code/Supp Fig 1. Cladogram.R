## -------------------------------------- ##
## ----    INITIALIZE ENVIRONMENT   ----- ##
## -------------------------------------- ##

# Clean environment
rm( list = ls() )

# load libraries
library(ggtree)
library(ape)
library(ggplot2)
library(phytools)


## -------------------------------------- ##
## ----    IMPORT METADATA FOR TREE ----- ##
## -------------------------------------- ##

# load metadata
setwd("~/Harvard/Research/GBS/Figures")
metadata_BCH <- read.csv("metadata_BCH_05292024.csv")


## -------------------------------------- ##
## ----    IMPORT TREE & DROP TIPS  ----- ##
## -------------------------------------- ##

# Import tree
setwd("~/Harvard/Research/GBS/For github/BCH-GBS-project/BCH_and_pediatric_CDC/Phylogeny")
tree_GBS_BCH_and_CDC <- read.tree("RAxML_bestTree.Core_phylogeny")
tree_GBS_BCH_and_CDC <- midpoint_root(tree_GBS_BCH_and_CDC)
tree_GBS_BCH <- keep.tip(tree_GBS_BCH_and_CDC, metadata_BCH$Sample)


## ------------------------------ ##
## ---- PREPARE & PLOT TREE ----- ##
## ------------------------------ ##

tree <- full_join(tree_GBS_BCH, 
                  metadata_BCH |> 
                    dplyr::select(Sample, CC, Serotype , ST) |> 
                    rename(label=Sample)  , by = 'label')


# plot tree (normal)
plot1 <- ggtree(tree, ladderize = T, lwd=0.5, layout="daylight", branch.length = 'none')

plot1 + 
  geom_tippoint(aes(colour= Serotype),size=2)+
  scale_colour_brewer(name="Serotype",palette = "Paired") 


# For adding CC and ST information
set.seed(1234)
CC_colors <- sample(RColorBrewer::brewer.pal(12, "Set3"))

plot1 + 
  geom_tippoint(aes(colour= as.factor(CC)),size=2)+
  scale_colour_manual(name="CC", values=CC_colors)



plot1 + 
  geom_tippoint(aes(colour= as.factor(ST)),size=2)



