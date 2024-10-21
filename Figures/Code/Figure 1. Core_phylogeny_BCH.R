## -------------------------------------- ##
## ----    INITIALIZE ENVIRONMENT   ----- ##
## -------------------------------------- ##

# Clean environment
rm( list = ls() )

# load libraries
library(ggtree)
library(ape)
library(dplyr)
library(ggnewscale)
library(ggplot2)
library(phytools)
library(RColorBrewer)
library(stringr)
library(cowplot)
library(readr)
library(ggpubr)
library(adephylo)





## -------------------------------------- ##
## ----    IMPORT METADATA FOR TREE ----- ##
## -------------------------------------- ##

# load metadata
setwd("~/Harvard/Research/GBS/Figures")
metadata_BCH <- read.csv("metadata_BCH_05292024.csv") |> 
  mutate(Age_cat=factor(Age_cat, levels=c("EOD", "LOD","VLOD", "older children (1-17yr)", "adults (>18yr)")))

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

Metadata_CC <- metadata_BCH |> 
  dplyr::select(Sample, CC) |> 
  rename(label=Sample) 
# remove extra sample 

tree1 <- drop.tip(tree_GBS_BCH, tree_GBS_BCH$tip.label[which(! tree_GBS_BCH$tip.label %in% Metadata_CC$label )] )

Metadata_CC <- Metadata_CC[which(Metadata_CC$label %in% tree1$tip.label), ] 


# join metadata and tree 
tree <- full_join(tree1, Metadata_CC , by = 'label')

set.seed(1234)
CC_colors <- sample(RColorBrewer::brewer.pal(12, "Set3"))

# plot tree (normal)
plot_GBS_ST_legend <- ggtree(tree, ladderize = T, lwd=0.5) + 
  geom_tippoint(aes(colour= as.factor(CC)),size=2)+
  scale_colour_manual(name="CC", values=CC_colors)



# plot tree (no legend)
plot_GBS_ST <- ggtree(tree, ladderize = T, lwd=0.5) 



## -------------------------------------------- ##
## ---- PREPARE COLLECTION DATE METADATA  ----- ##
## -------------------------------------------- ##

# Make dataframe in the appropriate format to integrate it to the tree
Metadata_date = data.frame(metadata_BCH$Diagnosis_date)
row.names(Metadata_date)= metadata_BCH$Sample
colnames(Metadata_date)= c("2")


## ------------------------------------- ##
## ---- PREPARE SEROTYPE METADATA  ----- ##
## ------------------------------------- ##
# Make dataframe in the appropriate format to integrate it to the tree
Metadata_serotype = data.frame(metadata_BCH$Serotype)
row.names(Metadata_serotype)= metadata_BCH$Sample
colnames(Metadata_serotype)= c("3")

## ------------------------------------- ##
## ---- PREPARE CC METADATA  ----- ##
## ------------------------------------- ##
# Make dataframe in the appropriate format to integrate it to the tree
Metadata_Age_cat = data.frame(metadata_BCH$Age_cat)
row.names(Metadata_Age_cat)= metadata_BCH$Sample
colnames(Metadata_Age_cat)= c("1")


## ------------------------------- ##
## --- prepare surface protein --- ##
## ------------------------------- ##

setwd("~/Harvard/Research/GBS/For github/BCH-GBS-project/BCH_and_pediatric_CDC/Metadata")
metadata_w_sip <- read.csv("metadata_BCH_CDC_children_w_sip_gene.csv")
metadata_w_sip_BCH <- metadata_w_sip |>
  filter(Source=="BCH")
metadata_w_sip_BCH$Sample <- metadata_BCH$Sample

metadata_BCH$Sip <- metadata_w_sip_BCH$Sip_gene_presence_absence

metadata_VF <- metadata_BCH |>
  mutate(GGNN2= ifelse(Surface_protein.ALP1+Surface_protein.ALP23+Surface_protein.ALPHA+Surface_protein.RIB !=0, 1, 0) ,
         SRR= ifelse(Surface_protein.SRR1+Surface_protein.SRR2 !=0, 1, 0), 
         PI= ifelse(Surface_protein.PI1+Surface_protein.PI2A1+Surface_protein.PI2A2+Surface_protein.PI2B !=0, 1, 0))|>
  dplyr::select(Sample, GGNN2, PI,Sip, contains("C5a"), Virulence_factor.Lmb, Virulence_factor.FbsB,SRR) 



# Replace row names by the tree's tip names
row.names(metadata_VF) <- metadata_VF$Sample
# Transform each column in the table to a character vector
for(i in 1:ncol(metadata_VF)){
  metadata_VF[,i]= as.character(metadata_VF[,i])
}
metadata_VF <- metadata_VF[,-1]
# Create vector to store AB class names (for plot legend)
Surface_proteins <- colnames(metadata_VF)
# replace column names in the metadata by 1-5 (for heatmap names)
colnames(metadata_VF) <- LETTERS[1:ncol(metadata_VF)]


## -------------------------------------- ##
## ----  IMPORT AMR METADATA  ----- ##
## -------------------------------------- ##

metadata_AMR <- metadata_BCH |>
  mutate(Aminoglycosides= ifelse(Antibiotic_Gene.ANT6IA+Antibiotic_Gene.APH3III !=0, 1, 0), 
         Tetracyclines= ifelse(Antibiotic_Gene.TETM+Antibiotic_Gene.TETO !=0, 1, 0) ,
         M_type= ifelse(Antibiotic_Gene.MEFA+Antibiotic_Gene.MSRD !=0, 1, 0), 
         MLSB= ifelse(Antibiotic_Gene.ERMA+Antibiotic_Gene.ERMB+Antibiotic_Gene.ERMT!=0, 1, 0))|>
  dplyr::select(Sample, Aminoglycosides, Tetracyclines, MLSB, M_type) 


# Replace row names by the tree's tip names
row.names(metadata_AMR) <- metadata_AMR$Sample
# Transform each column in the table to a character vector
for(i in 1:ncol(metadata_AMR)){
  metadata_AMR[,i]= as.character(metadata_AMR[,i])
}
metadata_AMR <- metadata_AMR[,-1]
# Create vector to store AB class names (for plot legend)
AMRs <- colnames(metadata_AMR)
# replace column names in the metadata by 1-5 (for heatmap names)
colnames(metadata_AMR) <- LETTERS[(ncol(metadata_VF)+1):(ncol(metadata_VF)+ncol(metadata_AMR))]


## ------------------------------ ##
## ---- PREPARE & PLOT TREE ----- ##
## ------------------------------ ##

# add source to the phylogenic tree as heatmap

Source_colors <- RColorBrewer::brewer.pal(5, "Blues")

p0 <- plot_GBS_ST_legend + new_scale_fill()
p1 <- gheatmap(p0, Metadata_Age_cat,  width=0.15, 
               font.size = 5, colnames_offset_y = -10, color=NA)+
  scale_fill_manual(values = Source_colors) +
  theme(legend.position="none")+ 
  ggtree::vexpand(.03, -1)

# create extra gheatmap for collection date to extract legend 
p0x <- plot_GBS_ST+ new_scale_fill()
p1x <- gheatmap(p0x, Metadata_Age_cat,  width=0.15, 
                font.size = 5, colnames_offset_y = -10, color=NA)+
  scale_fill_manual(values = Source_colors, name="1: Age group")


# add dates to the phylogenic tree as continuous heatmap
p2 <- p1 + new_scale_fill()
p3 <- gheatmap(p2, Metadata_date,  width=0.15, offset=0.002, 
               font.size = 5, colnames_offset_y = -10, color=NA) +
  scale_fill_viridis_c(option="A", name="2: Collection\nDate") +
  theme(legend.position="none")


# create extra gheatmap for collection date to extract legend 
p2x <- plot_GBS_ST + new_scale_fill()
p3x <- gheatmap(p2x, Metadata_date,  width=0.15, 
                font.size = 5, colnames_offset_y = -10, color=NA) +
  scale_fill_viridis_c(option="A", name="2: Collection date")


# add serotype to the phylogenic tree as continuous heatmap
p4 <- p3 + new_scale_fill()
p5 <- gheatmap(p4, Metadata_serotype,  width=0.15, offset=0.004,
               font.size = 5, colnames_offset_y = -10, color=NA) +
  scale_fill_manual(values = color("vibrant")(6)) + 
  theme(legend.position="none")

# create extra gheatmap for collection date to extract legend 
p4x <- plot_GBS_ST + new_scale_fill()
p5x <- gheatmap(p4x, Metadata_serotype,  width=0.15, offset=0.004,
                font.size = 5, colnames_offset_y = -10, color=NA)+
  scale_fill_manual(values = color("vibrant")(6), name="3: Serotype") 

# add heatmap for Surface proteins
p8 <- p5 + new_scale_fill()
p9 <- gheatmap(p8, metadata_VF, offset=0.008, width=1, 
               font.size = 5, colnames_offset_y = -10, color=NA) +
  scale_fill_manual(values=c("1"= "chartreuse4", "0"="grey91")) +
  theme(legend.position="none")


# extract legend 
p8x <- plot_GBS_ST + new_scale_fill()
p9x <- gheatmap(p8x, metadata_VF, offset=0.01, width=1, 
                font.size = 5, colnames_offset_y = -10, color=NA) +
  scale_fill_manual(values=c("1"= "chartreuse4", "0"="grey91"), name='Potential Vaccine Targets')


# add heatmap for virulence factors
p10 <- p9 + new_scale_fill()
p11_new <- gheatmap(p10, metadata_AMR, offset=0.02, width=0.5, 
                font.size = 5, colnames_offset_y = -10, color=NA) +
  scale_fill_manual(values=c("1"= "orange", "0"="grey91")) +
  theme(legend.position="none")




# extract legend 
p10x <- plot_GBS_ST + new_scale_fill()
p11_newx <- gheatmap(p10x, metadata_AMR, offset=0.25, width=0.5, 
                 font.size = 5, colnames_offset_y = -10, color=NA) +
  scale_fill_manual(values=c("1"= "orange", "0"="grey91"), name='AMR Genotypes')


# extract legends
pgp_leg   <- get_legend(plot_GBS_ST_legend)
p1x_leg   <- get_legend(p1x)
p3x_leg   <- get_legend(p3x)
p5x_leg   <- get_legend(p5x)


# plot object with legends
leg_grid0 = plot_grid(NULL,pgp_leg, NULL, ncol = 1, rel_heights= c(-0.4, 1, 0.5))
leg_grid1 = plot_grid(NULL,p1x_leg, NULL, ncol = 1, rel_heights= c(-0.4, 1, 0.5))
leg_grid2 = plot_grid(NULL,p3x_leg, NULL, ncol = 1, rel_heights= c(-0.2, 1, 2.5))
leg_grid3 = plot_grid(NULL,p5x_leg, NULL, ncol = 1, rel_heights= c(-0.15, 1, 2.5))


grod_plot_new= plot_grid(p11_new, NULL,
                         leg_grid0, NULL,
                         leg_grid1,  NULL,
                         leg_grid2,NULL,
                         leg_grid3, NULL,
                         ncol = 10,  
                         rel_widths = c(1, -0.01, 0.2, 0.01, 0.2, 0.01, 0.2, 0.01,0.2, 0.02))






# add legend for virulence factor and surface protein
grod_plot_new+
  annotate(geom = "text", x = 0.7, y = 0.6, label = "A", hjust = "left", size=3)+ 
  annotate(geom = "text", x = 0.7, y = 0.57, label = "B", hjust = "left", size=3)+ 
  annotate(geom = "text", x = 0.7, y = 0.54, label = "C", hjust = "left", size=3)+ 
  annotate(geom = "text", x = 0.7, y = 0.51, label = "D", hjust = "left", size=3)+ 
  annotate(geom = "text", x = 0.7, y = 0.48, label = "E", hjust = "left", size=3)+ 
  annotate(geom = "text", x = 0.7, y = 0.45, label = "F", hjust = "left", size=3)+ 
  annotate(geom = "text", x = 0.7, y = 0.42, label = "G", hjust = "left", size=3)+ 
  annotate(geom = "text", x = 0.74, y = 0.6, label = "AlphaC/RIB/ALP1/ALP23", hjust = "left", size=3)+ 
  annotate(geom = "text", x = 0.74, y = 0.57, label = "PI1/PI2A/PI2B", hjust = "left", size=3)+ 
  annotate(geom = "text", x = 0.74, y = 0.54, label = "Sip", hjust = "left", size=3)+
  annotate(geom = "text", x = 0.74, y = 0.51, label = "C5a", hjust = "left", size=3)+ 
  annotate(geom = "text", x = 0.74, y = 0.48, label = "Lmb", hjust = "left", size=3)+ 
  annotate(geom = "text", x = 0.74, y = 0.45, label = "FbsB", hjust = "left",size =3)+ 
  annotate(geom = "text", x = 0.74, y = 0.42, label = "SRR1/SRR2", hjust = "left", size=3)+
  annotate(geom = "point", x = 0.9, y = 0.6, colour = "chartreuse4", shape=15, size=3)+ 
  annotate(geom = "point", x = 0.90, y = 0.57, colour = "chartreuse4", shape=15, size=3)+ 
  annotate(geom = "point", x = 0.90, y = 0.54, colour = "chartreuse4", shape=15, size=3)+ 
  annotate(geom = "point", x = 0.90, y = 0.51, colour = "chartreuse4", shape=15, size=3)+ 
  annotate(geom = "point", x = 0.90, y = 0.48, colour = "chartreuse4", shape=15, size=3)+ 
  annotate(geom = "point", x = 0.90, y = 0.45, colour = "chartreuse4", shape=15, size=3)+ 
  annotate(geom = "point", x = 0.90, y = 0.42, colour = "chartreuse4", shape=15, size=3)+ 
  annotate(geom = "text", x = 0.7, y = 0.26, label = "I", hjust = "left", size=3)+ 
  annotate(geom = "text", x = 0.7, y = 0.23, label = "J", hjust = "left", size=3)+ 
  annotate(geom = "text", x = 0.7, y = 0.20, label = "K", hjust = "left", size=3)+
  annotate(geom = "text", x = 0.7, y = 0.17, label = "L", hjust = "left", size=3)+
  annotate(geom = "text", x = 0.74, y = 0.26, label = "Aminoglycosides", hjust = "left", size=3)+ 
  annotate(geom = "text", x = 0.74, y = 0.23, label = "Tetracyclines", hjust = "left", size=3)+ 
  annotate(geom = "text", x = 0.74, y = 0.20, label = "MLSB", hjust = "left", size=3)+ 
  annotate(geom = "text", x = 0.74, y = 0.17, label = "M-type", hjust = "left", size=3)+
  annotate(geom = "point", x = 0.90, y = 0.26, colour = "orange", shape=15, size=3)+ 
  annotate(geom = "point", x = 0.90, y = 0.23, colour = "orange", shape=15, size=3)+ 
  annotate(geom = "point", x = 0.90, y = 0.20, colour = "orange", shape=15, size=3)+ 
  annotate(geom = "point", x = 0.90, y = 0.17, colour = "orange", shape=15, size=3)+ 
  annotate(geom = "text", x = 0.815, y = 0.64,  label = "Potential Vaccine Targets", hjust = "right")+  
  annotate(geom = "text", x = 0.815 , y =0.3,  label = "Resistance Phenotypes", hjust = "right")
