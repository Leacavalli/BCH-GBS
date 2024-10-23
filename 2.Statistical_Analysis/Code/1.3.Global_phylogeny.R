## -------------------------------------- ##
## ----    INITIALIZE ENVIRONMENT   ----- ##
## -------------------------------------- ##

# Clean environment
rm( list = ls() )

# load libraries
library(ggtree)
library(ape)
library(dplyr)

# Load data
setwd("~/GitHub/BCH-GBS/2.Statistical_Analysis/Data")
Global_tree_df <- read.csv("Global_phylogeny_metadata.csv")

# Load tree
setwd("~/GitHub/BCH-GBS/1. Bioinformatic_Analysis/Outputs")
tree <- read.tree("Mashtree_Global.dnd")
tree_Final <- keep.tip(tree,  Global_tree_df$label)

# Join tree with data
tree_data <- full_join(tree_Final, Global_tree_df , by = 'label')


# plot tree (normal)
plot_tree <-ggtree(tree_data, ladderize = T, lwd=0.5) + 
  geom_tippoint(aes(colour= as.factor(clonal_complex)),size=2)+
  scale_colour_manual(name="clonal_complex", values=color("light")(9))
plot_tree

# plot tree (No Legend)
plot_tree_noLeg <-ggtree(tree_data, ladderize = T, lwd=0.5) 


## ------------------------------------- ##
## ---- PREPARE METADATA  HEATMAPS ----- ##
## ------------------------------------- ##


# SEROTYPE
Metadata_serotype = data.frame(Global_tree_df$Serotype)
row.names(Metadata_serotype)= Global_tree_df$label
colnames(Metadata_serotype)= c("1")


# Location
Metadata_Location = data.frame(Global_tree_df$Location)
row.names(Metadata_Location)= Global_tree_df$label
colnames(Metadata_Location)= c("3")



## ------------------------------ ##
## ---- PLOT TREE  HEATMAPS ----- ##
## ------------------------------ ##

p0 <- plot_tree + new_scale_fill()
p1 <- gheatmap(p0, Metadata_serotype,  width=0.1, offset=0.0002,
               font.size = 3, colnames_offset_y = -10, color=NA) +
  scale_fill_manual(values = c(color("vibrant")(7), color("medium contrast")(4)))+ 
  theme(legend.position="none")


# create extra gheatmap for collection date to extract legend 
p0x <- plot_tree_noLeg+ new_scale_fill()
p1x <- gheatmap(p0x, Metadata_serotype,  width=0.05,
                font.size = 3, colnames_offset_y = -5, color="white") +
  scale_fill_manual(values =c(color("vibrant")(7), color("medium contrast")(4)), name="1: Serotype") 


# add dates to the phylogenic tree as continuous heatmap
p2 <- p1 + new_scale_fill()
p3 <- gheatmap(p2, Metadata_Location,  width=0.1, offset=0.002, 
               font.size = 3, colnames_offset_y = -10, color=NA) +
  scale_fill_manual(values = color("light")(6))+ 
  theme(legend.position="none")


# create extra gheatmap for collection date to extract legend 
p2x <- plot_tree_noLeg + new_scale_fill()
p3x <- gheatmap(p2x, Metadata_Location,  width=0.05, 
                font.size = 3, colnames_offset_y = -5, color="white")  +
  scale_fill_manual(values = color("light")(6), name="2: Location of isolation")



## ------------------------------ ##
## ----     ADD LEGENDS     ----- ##
## ------------------------------ ##

# extract legends
pgp_leg   <- get_legend(plot_tree)
p1x_leg   <- get_legend(p1x)
p3x_leg   <- get_legend(p3x)

# plot object with legends
leg_grid0 = plot_grid(NULL,pgp_leg, NULL, ncol = 1, rel_heights= c(-0.15, 1, 0.5))
leg_grid1 = plot_grid(NULL,p1x_leg, NULL, ncol = 1, rel_heights= c(-0.15, 1, 0.5))
leg_grid2 = plot_grid(NULL,p3x_leg, NULL, ncol = 1, rel_heights= c(-0.15, 1, 2.5))

plot_grid1 <- plot_grid(leg_grid1,  NULL,
          leg_grid2,
          ncol = 1,  
          rel_heights  = c(0.5, -0.2, 0.2))

plot_grid2 <- plot_grid(NULL, leg_grid0, NULL,
          plot_grid1,  NULL,
          ncol = 5,  
          rel_widths = c(-0.15, 0.5, 0.1, 0.5, -0.15))

grod_plot_new <-  plot_grid(p3, NULL,
          plot_grid2,  NULL,
          ncol = 4,  
          rel_widths = c(0.5, 0.05, 0.2, 0.1))

grod_plot_new
