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
  mutate(label = paste0("S", as.numeric(sub("_.*", "", ID))))

# Import tree
setwd("~/GitHub/BCH-GBS/1. Bioinformatic_Analysis/Outputs")
# GAS
tree_GBS_BCH_GAS <- read.tree("RAxML_bipartitions.T3")
tree_GBS_BCH_GAS_clean <- keep.tip(tree_GBS_BCH_GAS, c(metadata_BCH$label, "Outgroup"))
# S.pneumo
tree_GBS_BCH_Spneumo <- read.tree("RAxML_bipartitions_Spneumo.T3")
tree_GBS_BCH_Spneumo_clean <- keep.tip(tree_GBS_BCH_Spneumo, c(metadata_BCH$label, "Spneumo_Outgroup"))


## -------------------------------------- ##
## ----   CHECK BEST OUTGROUP ROOT  ----- ##
## -------------------------------------- ##

# GAS 
tree_data_GAS_root <- full_join((drop.tip(root(tree_GBS_BCH_GAS_clean, "Outgroup"),  c("Outgroup"))), metadata_BCH , by = 'label')
plot_tree_GAS_root <- ggtree(tree_data_GAS_root, ladderize = T, lwd=0.5) + 
  geom_tippoint(aes(colour= as.factor(CC)),size=2)+
  scale_colour_manual(name="CC", values=color("light")(7))+
  ggtitle("GAS Outgroup")
plot_tree_GAS_root

# Spneumo 
tree_data_Spneumo_root <- full_join((drop.tip(root(tree_GBS_BCH_Spneumo_clean, "Spneumo_Outgroup"),  c("Spneumo_Outgroup"))), metadata_BCH , by = 'label')
plot_tree_Spneumo_root <- ggtree(tree_data_Spneumo_root, ladderize = T, lwd=0.5) + 
  geom_tippoint(aes(colour= as.factor(CC)),size=2)+
  scale_colour_manual(name="CC", values=color("light")(7))+
  ggtitle("S.pneumo Outgroup")
plot_tree_Spneumo_root

## -------------------------------------- ##
## -- CHECK SAME PATIENT DISTRIBUTION  -- ##
## -------------------------------------- ##

# load data
setwd("~/GitHub/BCH-GBS/0.Raw_Data")
metadata <- read.xlsx("GBS_mastersheet_LC.xlsx", sheetIndex = 1) |>
  filter(! ID %in% c("1_S1","32_S32","35_S35","81_S81"))|> 
  mutate(label = paste0("S", as.numeric(sub("_.*", "", ID)))) |> 
  select(label, Same_patient)
metadata$pair_id <- NA
metadata$pair_id[which(metadata$Same_patient =="Yes")] <- rep(seq(from = 1, to = 9), each = 2)

# prepare tree
tree_samepatient_GAS_root <- full_join((drop.tip(root(tree_GBS_BCH_GAS, "Outgroup"),  c("Outgroup"))), metadata_BCH , by = 'label')

#Join 
tree_Same_patient <- full_join(tree_samepatient_GAS_root, metadata , by = 'label')
# Plot
plot_tree_Same_patient <- ggtree(tree_Same_patient, ladderize = T, lwd=0.5) + 
  geom_tippoint(aes(colour= as.factor(pair_id)),size=3)+
  scale_colour_manual(name="Patient Pair", values=color("light")(9))+
  ggtitle("Blood and CSF isolates from the same patient")
plot_tree_Same_patient



## -------------------------------------- ##
## --     CHECK TWIN DISTRIBUTION      -- ##
## -------------------------------------- ##
# prepare data
metadata_twins <- metadata |> mutate(Twins=ifelse(label %in% c("S92", "S94"), "Twin ICU", ifelse(label =="S93", "Twin not-ICU", NA)))
# prepare tree
tree_twins_GAS_root <- full_join((drop.tip(root(tree_GBS_BCH_GAS, "Outgroup"),  c("Outgroup"))), metadata_BCH , by = 'label')
#Join 
tree_twins_patient <- full_join(tree_samepatient_GAS_root, metadata_twins , by = 'label')
# Plot
plot_tree_twins_patient <- ggtree(tree_twins_patient, ladderize = T, lwd=0.5) + 
  geom_tippoint(aes(colour= as.factor(Twins)),size=3)+
  scale_colour_manual(name="Twins", values=color("light")(9))+
  ggtitle("Twins")
plot_tree_twins_patient


# prepare tree
tree_twins_GAS_root2 <- full_join((keep.tip(root(tree_GBS_BCH_GAS, "Outgroup"),  c("S92", "S94","S93", "S77"))), metadata_BCH , by = 'label')
#Join 
tree_twins_patient2 <- full_join(tree_twins_GAS_root2, metadata_twins , by = 'label')
# Plot
plot_tree_twins_patient2 <- ggtree(tree_twins_patient2, ladderize = T, lwd=0.5) + 
  geom_tippoint(aes(colour= as.factor(Twins)),size=3)+
  scale_colour_manual(name="Twins", values=color("light")(9))+
  ggtitle("Twins")
plot_tree_twins_patient2

# Check genomic differences
metadata_BCH_twins <- metadata_BCH |> filter(ID %in% c("92_S92", "93_S93"))
diff_cols <- sapply(metadata_BCH_twins, function(col) col[1] != col[2])
# Keep only the columns where the values are different
metadata_BCH_twins_diff <- metadata_BCH_twins[, diff_cols]
# Check the result
print(metadata_BCH_twins_diff)
# There are no genomic differences.

