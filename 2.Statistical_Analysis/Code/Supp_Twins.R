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


# load data
setwd("~/GitHub/BCH-GBS/0.Raw_Data")
metadata <- read.xlsx("GBS_mastersheet_LC.xlsx", sheetIndex = 1) |>
  filter(! ID %in% c("1_S1","32_S32","35_S35","81_S81"))|>
  mutate(label = paste0("S", as.numeric(sub("_.*", "", ID)))) |>
  select(label, Same_patient, Source)
metadata$pair_id <- NA
metadata$pair_id[which(metadata$Same_patient =="Yes")] <- rep(seq(from = 1, to = 9), each = 2)

# Import tree
setwd("~/GitHub/BCH-GBS/1. Bioinformatic_Analysis/Outputs")
# GAS
tree_GBS_BCH_GAS <- read.tree("RAxML_bipartitions.T3")
tree_GBS_BCH_GAS_clean <- keep.tip(tree_GBS_BCH_GAS, c(metadata$label, "Outgroup"))


# prepare data
metadata_twins <- metadata |> 
  mutate(Patient=ifelse(label %in% c("S92", "S94"), "ICU Twin", ifelse(label =="S93", "In-patient Twin", "Other Patients")))
# prepare tree
tree_twins_GAS_root <- full_join((keep.tip(root(tree_GBS_BCH_GAS, "Outgroup"),  metadata$label)), metadata , by = 'label')
#Join 
tree_twins_patient <- full_join(tree_twins_GAS_root, metadata_twins , by = 'label')
# Plot
plot_tree_twins_patient <- ggtree(tree_twins_patient, ladderize = T, lwd=0.5) + 
  geom_tippoint(aes(colour= as.factor(Patient)),size=3)+
  scale_colour_manual(name="Patient", values=c(color("light")(2), "grey"))+
  theme(legend.position = "none")
plot_tree_twins_patient


# prepare tree
tree_twins_GAS_root2 <- full_join((keep.tip(root(tree_GBS_BCH_GAS, "Outgroup"),  c("S92", "S94","S93", "S77"))), 
                                  metadata_twins |> mutate(Patient= ifelse(label== "S77", "Other", Patient)), by = 'label')
#Join 
tree_twins_patient2 <- full_join(tree_twins_GAS_root2, metadata_twins , by = 'label')
# Plot
plot_tree_twins_patient2 <- ggtree(tree_twins_GAS_root2, ladderize = T, lwd=0.5) + 
  geom_tippoint(aes(colour= as.factor(Patient), shape= Source),size=5)+
  scale_colour_manual(name="Patient", values=c(color("light")(2), "grey"))+
  theme(legend.title =element_text(size=10), 
        legend.text = element_text(size=10), 
        legend.direction = "horizontal")

p2 <- plot_tree_twins_patient2 +
  theme(legend.position = "none")
legend <- get_legend(plot_tree_twins_patient2)


plot_tree_twins_patient+
  geom_rect(aes(xmin=0.006, xmax=0.0064, ymin=88, ymax= 92), alpha=0, colour = "red")+
  geom_segment(aes(x=0.0045, xend =0.006, y=47, yend = 88), lty="dashed", colour = "red")+
  geom_segment(aes(x=0.0097, xend =0.0064, y=47, yend = 88), lty="dashed", colour = "red")+
  geom_rect(aes(xmin=0.0045, xmax=0.0097, ymin=1, ymax= 47), alpha=0, colour = "red")+
  inset_element(p2, left = 0.55, bottom = 0.1, right = 0.9, top = 0.48)+ 
  inset_element(legend, left = 0.6, bottom = 0.05, right = 0.85, top = 0.1)

#save 10x9 

# Check genomic differences
metadata_BCH_twins <- metadata |> filter(ID %in% c("92_S92", "93_S93"))
diff_cols <- sapply(metadata_BCH_twins, function(col) col[1] != col[2])
# Keep only the columns where the values are different
metadata_BCH_twins_diff <- metadata_BCH_twins[, diff_cols]
# Check the result
print(metadata_BCH_twins_diff)
# There are no genomic differences.