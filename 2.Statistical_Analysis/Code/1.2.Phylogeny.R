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
setwd("~/GitHub/BCH-GBS/2.Statistical_Analysis/Data")
metadata_BCH <- read.xlsx("Clinical_Genomic_data.xlsx", sheetIndex = 1) |>
  filter(! ID %in% c("1_S1","32_S32","35_S35","81_S81"))|> 
  mutate(label = paste0("S", as.numeric(sub("_.*", "", ID))))|>
  mutate(Age_cat= ifelse(grepl("Older",Age_cat), "Older Children", Age_cat))

# Import tree
setwd("~/GitHub/BCH-GBS/1. Bioinformatic_Analysis/Outputs")
tree_GBS_BCH <- read.tree("RAxML_bipartitions.T3")
tree_GBS_BCH_clean <- keep.tip(tree_GBS_BCH, c(metadata_BCH$label, "Outgroup"))

## ------------------------------ ##
## ---- PREPARE & PLOT TREE ----- ##
## ------------------------------ ##
tree_GBS_BCH_Final <- (drop.tip(root(tree_GBS_BCH_clean, "Outgroup"),  c("Outgroup")))


Metadata_CC <- metadata_BCH |> select(label, CC)
tree_CC <- full_join(tree_GBS_BCH_Final, Metadata_CC , by = 'label')


# plot tree (normal)
plot_GBS_CC_legend <- ggtree(tree_CC, ladderize = T, lwd=0.5) + 
  geom_tippoint(aes(colour= as.factor(CC)),size=2)+
  scale_colour_manual(name="Clonal\nComplex", values=color("light")(7))


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
Metadata_date = data.frame(year(metadata_BCH$Date.of.cx))
row.names(Metadata_date)= metadata_BCH$label
colnames(Metadata_date)= c("2")

# Age cat
Metadata_Age_cat = data.frame(metadata_BCH$Age_cat)
row.names(Metadata_Age_cat)= metadata_BCH$label
Metadata_Age_cat$metadata_BCH.Age_cat =  factor(Metadata_Age_cat$metadata_BCH.Age_cat, levels = c("EOD", "LOD", "VLOD", "Older Children", "Adults"))
colnames(Metadata_Age_cat)= c("3")

# Clinical Outcomes
# ICU
Metadata_ICU = data.frame(metadata_BCH$Admission_loc)
row.names(Metadata_ICU)= metadata_BCH$label
Metadata_ICU$metadata_BCH.Admission_loc = ifelse(Metadata_ICU$metadata_BCH.Admission_loc == "ICU", "Yes", "No")
colnames(Metadata_ICU)= c("4")

# Meningitis
Metadata_Meningitis = data.frame(factor(metadata_BCH$Meningitis, levels= c(0,1)))
row.names(Metadata_Meningitis)= metadata_BCH$label
Metadata_Meningitis$factor.metadata_BCH.Meningitis..levels...c.0..1.. =  ifelse(Metadata_Meningitis$factor.metadata_BCH.Meningitis..levels...c.0..1.. ==1, "Yes", "No")
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
p3x <- gheatmap(p2x,  Metadata_date,  width=0.05, 
                font.size = 3, colnames_offset_y = -5, color="white") +
  scale_fill_viridis_c(option="A", name="2: Collection date")

# add source to the phylogenic tree as heatmap
p4 <- p3 + new_scale_fill()
p5 <- gheatmap(p4, Metadata_Age_cat,  width=0.05, offset=0.0014, 
               font.size = 3, colnames_offset_y = -5, color="white")+
  scale_fill_manual(values = c("#A01813", "#125A56", "#52B2D9", "#FD9A44", "#AAD1EE")) +
  theme(legend.position="none")+ 
  ggtree::vexpand(.03, -1)

# create extra gheatmap for collection date to extract legend 
p4x <- plot_GBS_CC+ new_scale_fill()
p5x <- gheatmap(p4x, Metadata_Age_cat,  width=0.05, 
                font.size = 3, colnames_offset_y = -5, color="white")+
  scale_fill_manual(values =  c("#A01813", "#125A56", "#52B2D9", "#FD9A44", "#AAD1EE"), 
                    name="3: Age group")

# add source to the phylogenic tree as heatmap
p6 <- p5 + new_scale_fill()
p7 <- gheatmap(p6, Metadata_ICU,  width=0.05, offset=0.0022, 
               font.size = 3, colnames_offset_y = -5, color="white")+
  scale_fill_manual(values = c("grey91", "red")) +
  theme(legend.position="none")+ 
  ggtree::vexpand(.03, -1)

# create extra gheatmap for collection date to extract legend 
p6x <- plot_GBS_CC + new_scale_fill()
p7x <- gheatmap(p0x, Metadata_ICU,  width=0.05, 
                font.size = 3, colnames_offset_y = -5, color="white")+
  scale_fill_manual(values = c("grey91", "red"), name="4: ICU")

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
p7x_leg   <- get_legend(p7x)
p9x_leg   <- get_legend(p9x)


# plot object with legends
leg_grid0 = plot_grid(NULL,pgp_leg, NULL, ncol = 1, rel_heights= c(-0.35, 1, 0.5))
leg_grid1 = plot_grid(NULL,p1x_leg, NULL, ncol = 1, rel_heights= c(-0.37, 1, 0.5))
leg_grid2 = plot_grid(NULL,p3x_leg, NULL, ncol = 1, rel_heights= c(-0.15, 1, 2.5))
leg_grid3 = plot_grid(NULL,p5x_leg, NULL, ncol = 1, rel_heights= c(-0.13, 1, 2.5))
leg_grid4 = plot_grid(NULL,p7x_leg, NULL, ncol = 1, rel_heights= c(-0.13, 1, 2.5))
leg_grid5 = plot_grid(NULL,p9x_leg, NULL, ncol = 1, rel_heights= c(-0.13, 1, 2.5))

# Add legends 1-5 to phylogeny
legend_grid1_BCH <- plot_grid(leg_grid0, NULL,
          leg_grid1,  NULL,
          leg_grid2,NULL,
          ncol = 6,  
          rel_widths = c(0.1, 0.01, 0.05, 0.01, 0.15, 0.01))

legend_grid2_BCH <- plot_grid(leg_grid3,NULL,
                          leg_grid4,NULL,
                          leg_grid5, NULL,
                          ncol = 6,  
                          rel_widths = c(0.15, 0.01, 0.05, 0.01, 0.15, 0.01))

legend_grid_BCH <-  plot_grid(NULL, legend_grid1_BCH, NULL, legend_grid2_BCH, NULL,
                           ncol=1,
                           rel_heights = c(0.3,0.5, -0.01,0.5,0.8))

grod_plot_new_BCH = plot_grid(p15, NULL, legend_grid_BCH,
          ncol = 3,  
          rel_widths = c(0.5, -0.01, 0.5))


# add legend for virulence factor and surface protein (A-P)
BCH_phylo <- grod_plot_new_BCH +
  annotate(geom = "text", x = 0.6, y = 0.4, label = "A", hjust = "left", size = 3) + 
  annotate(geom = "text", x = 0.6, y = 0.37, label = "B", hjust = "left", size = 3) + 
  annotate(geom = "text", x = 0.6, y = 0.34, label = "C", hjust = "left", size = 3) + 
  annotate(geom = "text", x = 0.6, y = 0.31, label = "D", hjust = "left", size = 3) + 
  annotate(geom = "text", x = 0.6, y = 0.25, label = "E", hjust = "left", size = 3) + 
  annotate(geom = "text", x = 0.6, y = 0.22, label = "F", hjust = "left", size = 3) + 
  annotate(geom = "text", x = 0.6, y = 0.19, label = "G", hjust = "left", size = 3) + 
  annotate(geom = "text", x = 0.6, y = 0.16, label = "H", hjust = "left", size = 3) + 
  annotate(geom = "text", x = 0.82, y = 0.4, label = "I", hjust = "left", size = 3) + 
  annotate(geom = "text", x = 0.82, y = 0.37, label = "J", hjust = "left", size = 3) + 
  annotate(geom = "text", x = 0.82, y = 0.34, label = "K", hjust = "left", size = 3) +
  annotate(geom = "text", x = 0.82, y = 0.31, label = "L", hjust = "left", size = 3) + 
  annotate(geom = "text", x = 0.82, y = 0.28, label = "M", hjust = "left", size = 3) + 
  annotate(geom = "text", x = 0.82, y = 0.25, label = "N", hjust = "left", size = 3) + 
  annotate(geom = "text", x = 0.82, y = 0.22, label = "O", hjust = "left", size = 3) + 
  annotate(geom = "text", x = 0.82, y = 0.19, label = "P", hjust = "left", size = 3) + 
  annotate(geom = "text", x = 0.64, y = 0.4, label = "ALP1", hjust = "left", size = 3) + 
  annotate(geom = "text", x = 0.64, y = 0.37, label = "ALP23", hjust = "left", size = 3) + 
  annotate(geom = "text", x = 0.64, y = 0.34, label = "Alpha", hjust = "left", size = 3) +
  annotate(geom = "text", x = 0.64, y = 0.31, label = "RIB", hjust = "left", size = 3) +
  annotate(geom = "text", x = 0.64, y = 0.25, label = "PI-1", hjust = "left", size = 3) + 
  annotate(geom = "text", x = 0.64, y = 0.22, label = "PI-2a1", hjust = "left", size = 3) + 
  annotate(geom = "text", x = 0.64, y = 0.19, label = "PI-2a2", hjust = "left", size = 3) + 
  annotate(geom = "text", x = 0.64, y = 0.16, label = "PI-2b", hjust = "left", size = 3) +
  annotate(geom = "text", x = 0.86, y = 0.4, label = "HVGA", hjust = "left", size = 3) + 
  annotate(geom = "text", x = 0.86, y = 0.37, label = "SRR1", hjust = "left", size = 3) + 
  annotate(geom = "text", x = 0.86, y = 0.34, label = "SRR2", hjust = "left", size = 3) + 
  annotate(geom = "text", x = 0.86, y = 0.31, label = "Sip", hjust = "left", size = 3) + 
  annotate(geom = "text", x = 0.86, y = 0.28, label = "lmb", hjust = "left", size = 3) + 
  annotate(geom = "text", x = 0.86, y = 0.25, label = "scpB", hjust = "left", size = 3) + 
  annotate(geom = "text", x = 0.86, y = 0.22, label = "hylB", hjust = "left", size = 3) + 
  annotate(geom = "text", x = 0.86, y = 0.19, label = "fbsB", hjust = "left", size = 3) +
  annotate(geom = "point", x = 0.7, y = 0.4, colour = "blue", shape = 15, size = 3) + 
  annotate(geom = "point", x = 0.7, y = 0.37, colour = "blue", shape = 15, size = 3) + 
  annotate(geom = "point", x = 0.7, y = 0.34, colour = "blue", shape = 15, size = 3) + 
  annotate(geom = "point", x = 0.7, y = 0.31, colour = "blue", shape = 15, size = 3) + 
  annotate(geom = "point", x = 0.7, y = 0.25, colour = "chartreuse4", shape = 15, size = 3) + 
  annotate(geom = "point", x = 0.7, y = 0.22, colour = "chartreuse4", shape = 15, size = 3) + 
  annotate(geom = "point", x = 0.7, y = 0.19, colour = "chartreuse4", shape = 15, size = 3) +
  annotate(geom = "point", x = 0.7, y = 0.16, colour = "chartreuse4", shape = 15, size = 3) + 
  annotate(geom = "point", x = 0.92, y = 0.4, colour = "orange", shape = 15, size = 3) + 
  annotate(geom = "point", x = 0.92, y = 0.37, colour = "orange", shape = 15, size = 3) + 
  annotate(geom = "point", x = 0.92, y = 0.34, colour = "orange", shape = 15, size = 3) + 
  annotate(geom = "point", x = 0.92, y = 0.31, colour = "orange", shape = 15, size = 3) + 
  annotate(geom = "point", x = 0.92, y = 0.28, colour = "orange", shape = 15, size = 3) + 
  annotate(geom = "point", x = 0.92, y = 0.25, colour = "orange", shape = 15, size = 3) + 
  annotate(geom = "point", x = 0.92, y = 0.22, colour = "orange", shape = 15, size = 3) + 
  annotate(geom = "point", x = 0.92, y = 0.19, colour = "orange", shape = 15, size = 3) + 
  annotate(geom = "text", x = 0.54, y = 0.38, label = "ALP Proteins", hjust = "left") +  
  annotate(geom = "text", x = 0.54, y = 0.24, label = "Pilus Islands", hjust = "left") +  
  annotate(geom = "text", x = 0.75, y = 0.37, label = "Other Virulence Factors", hjust = "left")

# save landscape: 12 x 8

BCH_phylo





## ------------------------------ ##
## ----   Global Phylogeny  ----- ##
## ------------------------------ ##



# Load data
setwd("~/GitHub/BCH-GBS/2.Statistical_Analysis/Data")
Global_tree_df <- read.csv("Global_phylogeny_metadata.csv") |>
  mutate(label = str_replace(label, "BCH", "S"))|>
  mutate(clonal_complex = str_remove_all(clonal_complex, "cc"))|>
  filter(! is.na(clonal_complex))|>
  mutate(Serotype = ifelse(grepl("III", Serotype), "III", Serotype))|>
  mutate(Location = factor(Location, levels = c("BCH", "USA", "other")))



#  Global_tree_df <- Global_tree_df|>  filter(Country == "BCH") 

# Load tree
setwd("~/GitHub/BCH-GBS/1. Bioinformatic_Analysis/Outputs")
# tree <- read.tree("Mashtree_Global.dnd")
tree <- read.tree("T1.raxml.bestTree")
tree_Final <- keep.tip(tree,  Global_tree_df$label)

# Join tree with data
tree_data <- full_join(tree_Final, Global_tree_df , by = 'label')


# plot tree (normal)
plot_tree <-ggtree(tree_data, ladderize = T, lwd=0.5) + 
  geom_tippoint(aes(colour= as.factor(clonal_complex)),size=2)+
  scale_colour_manual(name="Clonal\nComplex", values=c("#77AADD", "#EE8866", "#EEDD88", "#FFAABB", "#BBCC33", "#99DDFF", "#DDDDDD", "#AAAA00", "#44BB99"))



# plot tree (No Legend)
plot_tree_noLeg <-ggtree(tree_data, ladderize = T, lwd=0.5) 


## ------------------------------------- ##
## ---- PREPARE METADATA  HEATMAPS ----- ##
## ------------------------------------- ##

# ST
Metadata_serotype = data.frame(Global_tree_df$Serotype)
row.names(Metadata_serotype)= Global_tree_df$label
colnames(Metadata_serotype)= c("1")


# SEROTYPE
Metadata_serotype = data.frame(Global_tree_df$Serotype)
row.names(Metadata_serotype)= Global_tree_df$label
colnames(Metadata_serotype)= c("1")


# Location
Metadata_Location = data.frame(Global_tree_df$Location)
row.names(Metadata_Location)= Global_tree_df$label
colnames(Metadata_Location)= c("6")



## ------------------------------ ##
## ---- PLOT TREE  HEATMAPS ----- ##
## ------------------------------ ##

p0 <- plot_tree + new_scale_fill()
p1 <- gheatmap(p0, Metadata_serotype,  width=0.1, offset=0.0002,
               font.size = 3, colnames_offset_y = -35, color=NA) +
  scale_fill_manual(values = c(color("vibrant")(7), color("medium contrast")(4)))+ 
  theme(legend.position="none")


# create extra gheatmap for collection date to extract legend 
p0x <- plot_tree_noLeg+ new_scale_fill()
p1x <- gheatmap(p0x, Metadata_serotype,  width=0.05,
                font.size = 3, colnames_offset_y = -35, color="white") +
  scale_fill_manual(values =c(color("vibrant")(7), color("medium contrast")(4)), name="1: Serotype") 


# add dates to the phylogenic tree as continuous heatmap
p2 <- p1 + new_scale_fill()
p3 <- gheatmap(p2, Metadata_Location,  width=0.1, offset=0.003, 
               font.size = 3, colnames_offset_y = -35, color=NA) +
  scale_fill_manual(values = color("light")(6))+ 
  theme(legend.position="none")


# create extra gheatmap for collection date to extract legend 
p2x <- plot_tree_noLeg + new_scale_fill()
p3x <- gheatmap(p2x, Metadata_Location,  width=0.1, 
                font.size = 3, colnames_offset_y = -35, color="white")  +
  scale_fill_manual(values = color("light")(6), name="6: Location\nof isolation")



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
                        rel_heights  = c(0.5, -0.3, 0.2))

plot_grid2 <- plot_grid(NULL, leg_grid0, NULL,
                        plot_grid1,  NULL,
                        ncol = 5,  
                        rel_widths = c(-0.15, 0.5, -0.4, 0.5, -0.15))

Global_phylo <-  plot_grid(p3, NULL,
                            plot_grid2,  NULL,
                            ncol = 4,  
                            rel_widths = c(0.5, 0.05, 0.2, 0.1))

Global_phylo


## ------------------------------ ##
## ----  JOIN PHYLOGENIES   ----- ##
## ------------------------------ ##

plots_join <- plot_grid(p15, NULL, p3 + ggtree::vexpand(.12, -1)+ ggtree::vexpand(.01, 1),
          ncol = 3,  
          rel_widths = c(0.5, -0.01, 0.2))

legend_plots_join <- plot_grid(legend_grid_BCH, NULL, plot_grid2,
          ncol = 3,  
          rel_widths = c(0.5, -0.1, 0.2))

plot_grid(plots_join, NULL, legend_plots_join,
          ncol = 1,  
          rel_heights = c(0.5, -0.01, 0.2))
