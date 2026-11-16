# This script is for TaxSEA package testing taxon enrichment analysis with results from OTU table test
# 
# 
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("TaxSEA", version = "devel")


library(TaxSEA)
library(tidyverse)
library(dplyr)
# library(qiime2R)
data <- read.csv("H:\\OneDrive\\OneDrive - Indiana University\\Toshare\\Studies\\pet gut and cvd\\Enrichment\\07202026 result\\res_pair.csv")
taxa_info <- read.csv("H:\\OneDrive\\OneDrive - Indiana University\\Toshare\\Studies\\pet gut and cvd\\Enrichment\\07202026 result\\taxa_name.csv")
# implement those missing taxanomy information
colnames(taxa_info) <- c("Name",'kingdom', 'phylum', 'class', 
                         'order', 'family', 'genus', 'species','ID')
taxa_info <- taxa_info%>%
  mutate(kingdom = str_remove(kingdom,"k__"),
         phylum = str_remove(phylum,"p__"),
         class = str_remove(class,"c__"),
         order = str_remove(order,"o__"),
         family = str_remove(family,"f__"),
         genus = str_remove(genus,"g__"),
         species = str_remove(species,"s__"))
rownames(taxa_info) <- taxa_info$Name

test_taxa <- taxa_info%>%
  mutate(species=Name)
selected_taxon_level <- 'class'
my_custom_taxon_sets <- test_taxa %>%
  group_by(.data[[selected_taxon_level]]) %>% 
  summarise(species = list(species), .groups = "drop") %>%
  deframe()

test_col <- 3  # change here from 3 through 8 for each pair of comparison.
test <- data[,test_col]
colnames(data)[test_col]
names(test) <- data$taxon

custom_taxsea_results <- TaxSEA(taxon_ranks = test, custom_db = my_custom_taxon_sets)# , max_set_size = 500
custom_taxsea_results <- custom_taxsea_results$custom_sets
median(test <- data[,test_col])
