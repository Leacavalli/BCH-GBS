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

# load metadata
setwd("~/GitHub/BCH-GBS/2.Statistical_Analysis")
metadata_BCH <- read.xlsx("Clinical_Genomic_data.xlsx", sheetIndex = 1) |>
  filter(! ID %in% c("1_S1","32_S32","35_S35","81_S81"))|> 
  mutate(label = paste0("S", as.numeric(sub("_.*", "", ID))))

# Import tree
setwd("~/GitHub/BCH-GBS/1. Bioinformatic_Analysis/Outputs")
tree_GBS_BCH <- read.tree("RAxML_bipartitions.T3")
tree_GBS_BCH_clean <- keep.tip(tree_GBS_BCH, c(metadata_BCH$label, "Outgroup"))

## ---- CHECK BEST OUTGROUP ROOT ----- ##

# GAS 
tree_data <- full_join(tree_GBS_BCH_clean, metadata_BCH , by = 'label')
tree_data_GAS_root <- full_join(root(tree_GBS_BCH_clean, "Outgroup"), metadata_BCH , by = 'label')

ggtree(tree_data, ladderize = T, lwd=0.5) + 
  geom_tippoint(aes(colour= as.factor(CC)),size=2)

ggtree(tree_data_GAS_root, ladderize = T, lwd=0.5) + 
  geom_tippoint(aes(colour= as.factor(CC)),size=2)

# Spneumo 
tree_GBS_BCH_Spneumo <- read.tree("RAxML_bipartitions_Spneumo.T3")
tree_GBS_BCH_Spneumo_clean <- keep.tip(tree_GBS_BCH_Spneumo, c(metadata_BCH$label, "Spneumo_Outgroup"))
tree_data_Spneumo <- full_join(tree_GBS_BCH_Spneumo_clean, metadata_BCH , by = 'label')
tree_data_Spneumo_root <- full_join(root(tree_GBS_BCH_Spneumo_clean, "Spneumo_Outgroup"), metadata_BCH , by = 'label')

ggtree(tree_data_Spneumo, ladderize = T, lwd=0.5) + 
  geom_tippoint(aes(colour= as.factor(CC)),size=2)

ggtree(tree_data_Spneumo_root, ladderize = T, lwd=0.5) + 
  geom_tippoint(aes(colour= as.factor(CC)),size=2)



## ---- CHECK SAME PATIENT DISTRIBUTION ----- ##
tree_GBS_BCH_clean2 <- midpoint.root(drop.tip(tree_GBS_BCH,  c("1_S1","32_S32","35_S35","81_S81", "Reference", "Outgroup")))
# load data
setwd("~/GitHub/BCH-GBS/0.Raw_Data")
metadata <- read.xlsx("GBS_mastersheet_LC.xlsx", sheetIndex = 1) |>
  filter(! ID %in% c("1_S1","32_S32","35_S35","81_S81"))|> 
  mutate(label = paste0("S", as.numeric(sub("_.*", "", ID)))) |> 
  select(label, Same_patient)
#Join 
tree_Same_patient <- full_join(tree_GBS_BCH_clean2, metadata , by = 'label')
# Plot
ggtree(tree_Same_patient, ladderize = T, lwd=0.5) + 
  geom_tippoint(aes(colour= as.factor(Same_patient)),size=2)

## ------------------------------ ##
## ---- PREPARE & PLOT TREE ----- ##
## ------------------------------ ##
tree_GBS_BCH_Final <- (drop.tip(root(tree_GBS_BCH_clean, "Outgroup"),  c("Outgroup")))


Metadata_CC <- metadata_BCH |> select(label, CC)
tree_CC <- full_join(tree_GBS_BCH_Final, Metadata_CC , by = 'label')


# plot tree (normal)
plot_GBS_CC_legend <- ggtree(tree_CC, ladderize = T, lwd=0.5) + 
  geom_tippoint(aes(colour= as.factor(CC)),size=2)+
  scale_colour_manual(name="CC", values=color("light")(7))


# plot tree (no legend)
plot_GBS_CC <- ggtree(tree_CC, ladderize = T, lwd=0.5) 



## -------------------------------------------- ##
## ---- PREPARE METADATA  ----- ##
## -------------------------------------------- ##

# SEROTYPE
Metadata_serotype = data.frame(metadata_BCH$Serotype)
row.names(Metadata_serotype)= metadata_BCH$label
colnames(Metadata_serotype)= c("1")

# COLLECTION DATE
Metadata_date = data.frame(metadata_BCH$Date.of.cx)
row.names(Metadata_date)= metadata_BCH$label
colnames(Metadata_date)= c("2")

# Age cat
Metadata_Age_cat = data.frame(metadata_BCH$Age_cat)
row.names(Metadata_Age_cat)= metadata_BCH$label
colnames(Metadata_Age_cat)= c("3")

# Clinical Outcomes
# ICU
Metadata_ICU = data.frame(metadata_BCH$Admission_loc)
row.names(Metadata_ICU)= metadata_BCH$label
colnames(Metadata_ICU)= c("4")
# Meningitis
Metadata_Meningitis = data.frame(factor(metadata_BCH$Meningitis, levels= c(0,1)))
row.names(Metadata_Meningitis)= metadata_BCH$label
colnames(Metadata_Meningitis)= c("5")

# ALP Fam
Metadata_ALP_Fam = metadata_BCH |> 
  select(Surface_protein.ALP1, Surface_protein.ALP23, Surface_protein.ALPHA, Surface_protein.RIB)|>
  mutate(across(everything(), as.factor)) 
row.names(Metadata_ALP_Fam)= metadata_BCH$label
colnames(Metadata_ALP_Fam) <- LETTERS[1:ncol(Metadata_ALP_Fam)]

# Pilus islands
Metadata_Pilus = metadata_BCH |> 
  select(contains("Surface_protein.PI"))|>
  mutate(across(everything(), as.factor)) 
row.names(Metadata_Pilus)= metadata_BCH$label
colnames(Metadata_Pilus) <- LETTERS[(ncol(Metadata_ALP_Fam)+1):(ncol(Metadata_ALP_Fam)+ncol(Metadata_Pilus))]

# Other virulence protein
Metadata_Other_Vir = metadata_BCH |> 
  select(contains(c("Surface_protein.", "Virulence_Gene.")))|>
  select(-contains("Surface_protein.PI"), -Surface_protein.Sip.1a, -Surface_protein.Sip.3a, -Surface_protein.ALP1, -Surface_protein.ALP23, -Surface_protein.ALPHA, -Surface_protein.RIB)|>
  mutate(across(everything(), as.factor)) 
row.names(Metadata_Other_Vir)= metadata_BCH$label
colnames(Metadata_Other_Vir) <- LETTERS[(ncol(Metadata_ALP_Fam)+ncol(Metadata_Pilus)+1):(ncol(Metadata_ALP_Fam)+ncol(Metadata_Pilus)+ncol(Metadata_Other_Vir))]


## ------------------------------ ##
## ---- PLOT TREE  HEATMAPS ----- ##
## ------------------------------ ##

p0 <- plot_GBS_CC_legend + new_scale_fill()
p1 <- gheatmap(p0, Metadata_serotype,  width=0.05, offset=0.0002,
               font.size = 3, colnames_offset_y = -5, color="white") +
  scale_fill_manual(values = color("vibrant")(6))+ 
  theme(legend.position="none")


# create extra gheatmap for collection date to extract legend 
p0x <- plot_GBS_CC+ new_scale_fill()
p1x <- gheatmap(p0x, Metadata_serotype,  width=0.05,
                font.size = 3, colnames_offset_y = -5, color="white") +
  scale_fill_manual(values = color("vibrant")(6), name="1: Serotype") 


# add dates to the phylogenic tree as continuous heatmap
p2 <- p1 + new_scale_fill()
p3 <- gheatmap(p2, Metadata_date,  width=0.05, offset=0.0008, 
               font.size = 3, colnames_offset_y = -5, color="white") +
  scale_fill_viridis_c(option="A", name="2: Collection\nDate") +
  theme(legend.position="none")


# create extra gheatmap for collection date to extract legend 
p2x <- plot_GBS_CC + new_scale_fill()
p3x <- gheatmap(p2x, Metadata_date,  width=0.05, 
                font.size = 3, colnames_offset_y = -5, color="white") +
  scale_fill_viridis_c(option="A", name="2: Collection date")

# add source to the phylogenic tree as heatmap
Source_colors <- RColorBrewer::brewer.pal(5, "Blues")

p4 <- p3 + new_scale_fill()
p5 <- gheatmap(p4, Metadata_Age_cat,  width=0.05, offset=0.0014, 
               font.size = 3, colnames_offset_y = -5, color="white")+
  scale_fill_manual(values = Source_colors) +
  theme(legend.position="none")+ 
  ggtree::vexpand(.03, -1)

# create extra gheatmap for collection date to extract legend 
p4x <- plot_GBS_CC+ new_scale_fill()
p5x <- gheatmap(p0x, Metadata_Age_cat,  width=0.05, 
                font.size = 3, colnames_offset_y = -5, color="white")+
  scale_fill_manual(values = Source_colors, name="3: Age group")

# add source to the phylogenic tree as heatmap
p6 <- p5 + new_scale_fill()
p7 <- gheatmap(p6, Metadata_ICU,  width=0.05, offset=0.0022, 
               font.size = 3, colnames_offset_y = -5, color="white")+
  scale_fill_manual(values = c("red", "grey91")) +
  theme(legend.position="none")+ 
  ggtree::vexpand(.03, -1)

# create extra gheatmap for collection date to extract legend 
p6x <- plot_GBS_CC + new_scale_fill()
p7x <- gheatmap(p0x, Metadata_ICU,  width=0.05, 
                font.size = 3, colnames_offset_y = -5, color="white")+
  scale_fill_manual(values = c("red", "grey91"), name="4: ICU")

# add source to the phylogenic tree as heatmap
p8 <- p7 + new_scale_fill()
p9 <- gheatmap(p8, Metadata_Meningitis,  width=0.05, offset=0.0026, 
               font.size = 3, colnames_offset_y = -5, color="white")+
  scale_fill_manual(values = c("grey91", "orange")) +
  theme(legend.position="none")+ 
  ggtree::vexpand(.03, -1)

# create extra gheatmap for collection date to extract legend 
p8x <- plot_GBS_CC + new_scale_fill()
p9x <- gheatmap(p0x, Metadata_Meningitis,  width=0.05, 
                font.size = 3, colnames_offset_y = -5, color="white")+
  scale_fill_manual(values = c("grey91", "orange"), name="5: Meningitis")


# add heatmap for ALP proteins
p10 <- p9 + new_scale_fill()
p11 <- gheatmap(p10, Metadata_ALP_Fam,  width=0.20, offset=0.0034, 
                font.size = 3, colnames_offset_y = -5, color="white") +
  scale_fill_manual(values=c("1"= "blue", "0"="grey91")) +
  theme(legend.position="none")

# extract legend 
p10x <- plot_GBS_CC + new_scale_fill()
p11x <- gheatmap(p10x, Metadata_ALP_Fam, offset=0.01, width=1, 
                font.size = 3, colnames_offset_y = -5, color=NA) +
  scale_fill_manual(values=c("1"= "blue", "0"="grey91"), name='ALP proteins')


# add heatmap for Pilus Islands
p12 <- p11 + new_scale_fill()
p13 <- gheatmap(p12, Metadata_Pilus,  width=0.20, offset=0.0058, 
                font.size = 3, colnames_offset_y = -5, color="white") +
  scale_fill_manual(values=c("1"= "chartreuse4", "0"="grey91")) +
  theme(legend.position="none")

# extract legend 
p12x <- plot_GBS_CC + new_scale_fill()
p13x <- gheatmap(p12x, Metadata_Pilus, offset=0.01, width=1, 
                 font.size = 3, colnames_offset_y = -5, color=NA) +
  scale_fill_manual(values=c("1"= "chartreuse4", "0"="grey91"), name='Pilus Islands')


# add heatmap for Other virulence factors
p14 <- p13 + new_scale_fill()
p15 <- gheatmap(p14, Metadata_Other_Vir,  width=0.40, offset=0.0082, 
                font.size = 3, colnames_offset_y = -5, color="white") +
  scale_fill_manual(values=c("1"= "orange3", "0"="grey91")) +
  theme(legend.position="none")

# extract legend 
p14x <- plot_GBS_CC + new_scale_fill()
p15x <- gheatmap(p14x, Metadata_Other_Vir, offset=0.01, width=1, 
                 font.size = 3, colnames_offset_y = -5, color=NA) +
  scale_fill_manual(values=c("1"= "orange3", "0"="grey91"), name='Other Virulence Factors')


## ------------------------------ ##
## ----     ADD LEGENDS     ----- ##
## ------------------------------ ##


# extract legends
pgp_leg   <- get_legend(plot_GBS_CC_legend)
p1x_leg   <- get_legend(p1x)
p3x_leg   <- get_legend(p3x)
p5x_leg   <- get_legend(p5x)


# plot object with legends
leg_grid0 = plot_grid(NULL,pgp_leg, NULL, ncol = 1, rel_heights= c(-0.35, 1, 0.5))
leg_grid1 = plot_grid(NULL,p1x_leg, NULL, ncol = 1, rel_heights= c(-0.37, 1, 0.5))
leg_grid2 = plot_grid(NULL,p3x_leg, NULL, ncol = 1, rel_heights= c(-0.15, 1, 2.5))
leg_grid3 = plot_grid(NULL,p5x_leg, NULL, ncol = 1, rel_heights= c(-0.13, 1, 2.5))


grod_plot_new= plot_grid(p15, NULL,
                         leg_grid0, NULL,
                         leg_grid1,  NULL,
                         leg_grid2,NULL,
                         leg_grid3, NULL,
                         ncol = 10,  
                         rel_widths = c(0.5, -0.01, 0.1, 0.01, 0.05, 0.01, 0.15, 0.01,0.1, 0.02))

grod_plot_new




# add legend for virulence factor and surface protein
grod_plot_new+
  annotate(geom = "text", x = 0.6, y = 0.6, label = "A", hjust = "left", size=3)+ 
  annotate(geom = "text", x = 0.6, y = 0.57, label = "B", hjust = "left", size=3)+ 
  annotate(geom = "text", x = 0.6, y = 0.54, label = "C", hjust = "left", size=3)+ 
  annotate(geom = "text", x = 0.6, y = 0.51, label = "D", hjust = "left", size=3)+ 
  annotate(geom = "text", x = 0.6, y = 0.45, label = "E", hjust = "left", size=3)+ 
  annotate(geom = "text", x = 0.6, y = 0.42, label = "F", hjust = "left", size=3)+ 
  annotate(geom = "text", x = 0.6, y = 0.39, label = "G", hjust = "left", size=3)+ 
  annotate(geom = "text", x = 0.6, y = 0.36, label = "H", hjust = "left", size=3)+ 
  annotate(geom = "text", x = 0.82, y = 0.6, label = "I", hjust = "left", size=3)+ 
  annotate(geom = "text", x = 0.82, y = 0.57, label = "J", hjust = "left", size=3)+ 
  annotate(geom = "text", x = 0.82, y = 0.54, label = "K", hjust = "left", size=3)+
  annotate(geom = "text", x = 0.82, y = 0.51, label = "L", hjust = "left", size=3)+ 
  annotate(geom = "text", x = 0.82, y = 0.48, label = "M", hjust = "left", size=3)+ 
  annotate(geom = "text", x = 0.82, y = 0.45, label = "N", hjust = "left", size=3)+ 
  annotate(geom = "text", x = 0.82, y = 0.42, label = "O", hjust = "left", size=3)+ 
  annotate(geom = "text", x = 0.82, y = 0.39, label = "P", hjust = "left", size=3)+ 
  annotate(geom = "text", x = 0.64, y = 0.6, label = "ALP1", hjust = "left", size=3)+ 
  annotate(geom = "text", x = 0.64, y = 0.57, label = "ALP23", hjust = "left", size=3)+ 
  annotate(geom = "text", x = 0.64, y = 0.54, label = "Alpha", hjust = "left", size=3)+
  annotate(geom = "text", x = 0.64, y = 0.51, label = "RIB", hjust = "left", size=3)+
  annotate(geom = "text", x = 0.64, y = 0.45, label = "PI-1", hjust = "left", size=3)+ 
  annotate(geom = "text", x = 0.64, y = 0.42, label = "PI-2a1", hjust = "left", size=3)+ 
  annotate(geom = "text", x = 0.64, y = 0.39, label = "PI-2a2", hjust = "left",size =3)+ 
  annotate(geom = "text", x = 0.64, y = 0.36, label = "PI-2b", hjust = "left", size=3)+
  annotate(geom = "text", x = 0.86, y = 0.6, label = "HVGA", hjust = "left", size=3)+ 
  annotate(geom = "text", x = 0.86, y = 0.57, label = "SRR1", hjust = "left", size=3)+ 
  annotate(geom = "text", x = 0.86, y = 0.54, label = "SRR2", hjust = "left", size=3)+ 
  annotate(geom = "text", x = 0.86, y = 0.51, label = "Sip", hjust = "left", size=3)+ 
  annotate(geom = "text", x = 0.86, y = 0.48, label = "lmb", hjust = "left", size=3)+ 
  annotate(geom = "text", x = 0.86, y = 0.45, label = "scpB", hjust = "left", size=3)+ 
  annotate(geom = "text", x = 0.86, y = 0.42, label = "hylB", hjust = "left", size=3)+ 
  annotate(geom = "text", x = 0.86, y = 0.39, label = "fbsB", hjust = "left", size=3)+
  annotate(geom = "point", x = 0.7, y = 0.6, colour = "blue", shape=15, size=3)+ 
  annotate(geom = "point", x = 0.7, y = 0.57, colour = "blue", shape=15, size=3)+ 
  annotate(geom = "point", x = 0.7, y = 0.54, colour = "blue", shape=15, size=3)+ 
  annotate(geom = "point", x = 0.7, y = 0.51, colour = "blue", shape=15, size=3)+ 
  annotate(geom = "point", x = 0.7, y = 0.45, colour = "chartreuse4", shape=15, size=3)+ 
  annotate(geom = "point", x = 0.7, y = 0.42, colour = "chartreuse4", shape=15, size=3)+ 
  annotate(geom = "point", x = 0.7, y = 0.39, colour = "chartreuse4", shape=15, size=3)+
  annotate(geom = "point", x = 0.7, y = 0.36, colour = "chartreuse4", shape=15, size=3)+ 
  annotate(geom = "point", x = 0.92, y = 0.6, colour = "orange", shape=15, size=3)+ 
  annotate(geom = "point", x = 0.92, y = 0.57, colour = "orange", shape=15, size=3)+ 
  annotate(geom = "point", x = 0.92, y = 0.54, colour = "orange", shape=15, size=3)+ 
  annotate(geom = "point", x = 0.92, y = 0.51, colour = "orange", shape=15, size=3)+ 
  annotate(geom = "point", x = 0.92, y = 0.48, colour = "orange", shape=15, size=3)+ 
  annotate(geom = "point", x = 0.92, y = 0.45, colour = "orange", shape=15, size=3)+ 
  annotate(geom = "point", x = 0.92, y = 0.42, colour = "orange", shape=15, size=3)+ 
  annotate(geom = "point", x = 0.92, y = 0.39, colour = "orange", shape=15, size=3)+ 
  annotate(geom = "text", x = 0.54, y = 0.58,  label = "ALP\nProteins", hjust = "left")+  
  annotate(geom = "text", x = 0.54 , y =0.44,  label = "Pilus\nIslands", hjust = "left")+  
  annotate(geom = "text", x = 0.75 , y =0.57,  label = "Other\nVirulence\nFactors", hjust = "left")

