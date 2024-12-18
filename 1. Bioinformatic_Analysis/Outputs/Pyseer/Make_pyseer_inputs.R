# Clean environment
rm( list = ls() )

# load libraries
library(xlsx)
library(tidyverse)
library(stringr)

# load data
setwd("~/GitHub/BCH-GBS/2.Statistical_Analysis/Data")
data <- read.xlsx("Clinical_Genomic_data.xlsx", sheetIndex = 1) 

# Remove contaminated samples
data_clean <- data|>
  filter(! ID %in% c("1_S1","32_S32","35_S35","81_S81"))|> 
  mutate(ICU = ifelse(Admission_loc=="ICU", 1,0))|> 
  mutate(ID = sub("^(\\d+)_S(\\d+)$", "S\\1", ID))

# Subset to infants alone 
data_clean_infants <- data_clean |> 
  filter(days_at_dx <= 365)

#####################################
#    Infants vs Older patients      #
#####################################

# gene_presence_absence
gene_presence_absence_Infants_vs_Older <- data_clean |> 
  select(ID, contains("Surface_protein.") , contains("Virulence_Gene."))|> 
  pivot_longer(
    cols = -ID,  # All columns except ID
    names_to = "Gene",  # Name of the new column for variable names
    values_to = "Value"  # Name of the column for the values
  ) |> 
  pivot_wider(
    names_from = ID,  # Pivot columns based on ID values
    values_from = Value  # Values to be filled in the new table
  )
write_delim(gene_presence_absence_Infants_vs_Older, "~/GitHub/BCH-GBS/2.Statistical_Analysis/Temp/Pyseer/Infants_vs_Older/gene_presence_absence.txt", delim = "\t")

# phenotypes
phenotypes_Infants_vs_Older <- data_clean |> 
  select(ID, infant) |>
  rename("samples" = "ID", 
         "binary" = "infant")
setwd("~/GitHub/BCH-GBS/2.Statistical_Analysis/Temp/Pyseer")
write_delim(phenotypes_Infants_vs_Older, "~/GitHub/BCH-GBS/2.Statistical_Analysis/Temp/Pyseer/Infants_vs_Older/phenotypes.txt", delim = "\t")

# Import tree
setwd("~/GitHub/BCH-GBS/1. Bioinformatic_Analysis/Outputs")
tree_GBS_BCH <- read.tree("RAxML_bipartitions.T3")
write.tree(tree_GBS_BCH, "~/GitHub/BCH-GBS/2.Statistical_Analysis/Temp/Pyseer/Infants_vs_Older/tree_GBS_BCH.tre")


#######################
#   LOD vs VLOD       #
#######################

# Subset data
data_LOD <- data_clean_infants |> 
  filter(Age_cat != "EOD") |> 
  mutate(LOD= ifelse(Age_cat == "LOD", 1, 0))

# gene_presence_absence
gene_presence_absence_LOD_vs_VLOD <- data_LOD |> 
  select(ID, contains("Surface_protein.") , contains("Virulence_Gene."))|> 
  pivot_longer(
    cols = -ID,  # All columns except ID
    names_to = "Gene",  # Name of the new column for variable names
    values_to = "Value"  # Name of the column for the values
  ) |> 
  pivot_wider(
    names_from = ID,  # Pivot columns based on ID values
    values_from = Value  # Values to be filled in the new table
  )
write_delim(gene_presence_absence_LOD_vs_VLOD, "~/GitHub/BCH-GBS/2.Statistical_Analysis/Temp/Pyseer/LOD_vs_VLOD/gene_presence_absence.txt", delim = "\t")

# phenotypes
phenotypes_LOD_vs_VLOD <- data_LOD |> 
  select(ID, LOD) |>
  rename("samples" = "ID", 
         "binary" = "LOD")
setwd("~/GitHub/BCH-GBS/2.Statistical_Analysis/Temp/Pyseer")
write_delim(phenotypes_LOD_vs_VLOD, "~/GitHub/BCH-GBS/2.Statistical_Analysis/Temp/Pyseer/LOD_vs_VLOD/phenotypes.txt", delim = "\t")

# Import tree
setwd("~/GitHub/BCH-GBS/1. Bioinformatic_Analysis/Outputs")
tree_GBS_BCH <- read.tree("RAxML_bipartitions.T3")
tree_GBS_BCH_LOD_vs_VLOD <- keep.tip(tree_GBS_BCH, data_LOD$ID)
write.tree(tree_GBS_BCH_LOD_vs_VLOD, "~/GitHub/BCH-GBS/2.Statistical_Analysis/Temp/Pyseer/LOD_vs_VLOD/tree_GBS_BCH.tre")




#######################
#    ICU admission    #
#######################

# gene_presence_absence
gene_presence_absence_infants <- data_clean_infants |> 
  select(ID, contains("Surface_protein.") , contains("Virulence_Gene."))|> 
  pivot_longer(
    cols = -ID,  # All columns except ID
    names_to = "Gene",  # Name of the new column for variable names
    values_to = "Value"  # Name of the column for the values
  ) |> 
  pivot_wider(
    names_from = ID,  # Pivot columns based on ID values
    values_from = Value  # Values to be filled in the new table
  )
write_delim(gene_presence_absence_infants, "~/GitHub/BCH-GBS/2.Statistical_Analysis/Temp/Pyseer/ICU_Infants/gene_presence_absence.txt", delim = "\t")

# phenotypes
phenotypes_ICU <- data_clean_infants |> 
  select(ID, ICU) |>
  rename("samples" = "ID", 
         "binary" = "ICU")
write_delim(phenotypes_ICU, "~/GitHub/BCH-GBS/2.Statistical_Analysis/Temp/Pyseer/ICU_Infants/phenotypes.txt", delim = "\t")

# Import tree
setwd("~/GitHub/BCH-GBS/1. Bioinformatic_Analysis/Outputs")
tree_GBS_BCH <- read.tree("RAxML_bipartitions.T3")
tree_GBS_BCH_infants <- keep.tip(tree_GBS_BCH, data_clean_infants$ID)
write.tree(tree_GBS_BCH_infants, "~/GitHub/BCH-GBS/2.Statistical_Analysis/Temp/Pyseer/ICU_Infants/tree_GBS_BCH.tre")




#######################
#    Meningitis       #
#######################

# gene_presence_absence
write_delim(gene_presence_absence_infants, "~/GitHub/BCH-GBS/2.Statistical_Analysis/Temp/Pyseer/Meningitis_Infants/gene_presence_absence.txt", delim = "\t")

# phenotypes
phenotypes_Meningitis <- data_clean_infants |> 
  select(ID, Meningitis) |>
  rename("samples" = "ID", 
         "binary" = "Meningitis")
write_delim(phenotypes_Meningitis, "~/GitHub/BCH-GBS/2.Statistical_Analysis/Temp/Pyseer/Meningitis_Infants/phenotypes.txt", delim = "\t")

# Import tree
write.tree(tree_GBS_BCH_infants, "~/GitHub/BCH-GBS/2.Statistical_Analysis/Temp/Pyseer/Meningitis_Infants/tree_GBS_BCH.tre")


#######################
#    Neutropenia       #
#######################

# gene_presence_absence
write_delim(gene_presence_absence_infants, "~/GitHub/BCH-GBS/2.Statistical_Analysis/Temp/Pyseer/Neutropenia_Infants/gene_presence_absence.txt", delim = "\t")

# phenotypes
phenotypes_Neutropenia <- data_clean_infants |> 
  select(ID, Neutropenia) |>
  rename("samples" = "ID", 
         "binary" = "Neutropenia")
write_delim(phenotypes_Neutropenia, "~/GitHub/BCH-GBS/2.Statistical_Analysis/Temp/Pyseer/Neutropenia_Infants/phenotypes.txt", delim = "\t")

# Import tree
write.tree(tree_GBS_BCH_infants, "~/GitHub/BCH-GBS/2.Statistical_Analysis/Temp/Pyseer/Neutropenia_Infants/tree_GBS_BCH.tre")

#######################
#    Leukopenia       #
#######################

# gene_presence_absence
write_delim(gene_presence_absence_infants, "~/GitHub/BCH-GBS/2.Statistical_Analysis/Temp/Pyseer/Leukopenia_Infants/gene_presence_absence.txt", delim = "\t")

# phenotypes
phenotypes_Leukopenia <- data_clean_infants |> 
  select(ID, Leukopenia) |>
  rename("samples" = "ID", 
         "binary" = "Leukopenia")
write_delim(phenotypes_Leukopenia, "~/GitHub/BCH-GBS/2.Statistical_Analysis/Temp/Pyseer/Leukopenia_Infants/phenotypes.txt", delim = "\t")

# Import tree
write.tree(tree_GBS_BCH_infants, "~/GitHub/BCH-GBS/2.Statistical_Analysis/Temp/Pyseer/Leukopenia_Infants/tree_GBS_BCH.tre")

#######################
#    Leukocytosis       #
#######################

# gene_presence_absence
write_delim(gene_presence_absence_infants, "~/GitHub/BCH-GBS/2.Statistical_Analysis/Temp/Pyseer/Leukocytosis_Infants/gene_presence_absence.txt", delim = "\t")

# phenotypes
phenotypes_Leukocytosis <- data_clean_infants |> 
  select(ID, Leukocytosis) |>
  rename("samples" = "ID", 
         "binary" = "Leukocytosis")
write_delim(phenotypes_Leukocytosis, "~/GitHub/BCH-GBS/2.Statistical_Analysis/Temp/Pyseer/Leukocytosis_Infants/phenotypes.txt", delim = "\t")

# Import tree
write.tree(tree_GBS_BCH_infants, "~/GitHub/BCH-GBS/2.Statistical_Analysis/Temp/Pyseer/Leukocytosis_Infants/tree_GBS_BCH.tre")
