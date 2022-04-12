# Here we generate pheno and geno matrices for ML models
# The geno matrices could be created using the prewas package, but this needs debugging first
# Therefore, we simply upload the gene burden tests based on pyseer

# Library
library(tidyverse)
rm(list = ls())

setwd(paste0(here::here(), "/ALL/ML/no_variants_data/model2_lineages"))
getwd()

# Pheno and lineages file ---------------------------------------------------------------------------------------
f <- paste0(
  here::here(),
  "/ALL/Phenotypes/processed_data/df_mortality_all.Rda")
pheno <- readRDS(f) %>%
  drop_na(mortality, vancoetest, sab_duration, gender, age, hca, cci, pitt, pat, mecA, CC) %>%
  select(sample_id, mortality, CC) %>%
  mutate(CC = as.factor(CC))

# Geno file ----------------------------------------------------------------------------------------
# Raw data: pyseer gene burden output for VANANZ baseline mortality
# We have generated two sample x gene matrices:
# a matrix with rare mutations only 
# a matrix with all mutations
# 
# We coded mutated genes as 0 - non mutated, 1 - missense only, 2 - at least one truncation
f <- paste0(
  here::here(),
  "/VANANZ/baseline_strain_mortality/GWAS/processed_data/genotype_matrices/rare.genes.geno.Rda"
)
geno <- readRDS(f) 
colnames(geno) <- sort(colnames(geno))
# pheatmap::pheatmap(geno, show_rownames = F, show_colnames = F)
geno <- geno %>%
  as_tibble(rownames = "sample_id")

# Lineages
f <- paste0(
  here::here(),
  "/ALL/Lineages/processed_data/MDS/df_mash_mds_all_strains.Rda"
)
lineages <- readRDS(f)

# Merge into one dataframe -------------------------------------------------------------------------
pheno_geno <- pheno %>%
  # inner_join(geno) %>%
  inner_join(lineages) %>%
  select(-CC) %>%
  select(-sample_id) 
# quick check for missing values
map_int(pheno_geno, function(x) sum(is.na(x))) %>%
  head()
map_int(pheno_geno, function(x) sum(is.na(x))) %>%
  unique()

# Save processed data ------------------------------------------------------------------------------
pheno_geno %>%
  saveRDS("processed_data/pheno_geno_data/VANANZ.lineages.pheno.geno.Rda")

pheno_geno %>%
  saveRDS("~/Documents/Transfer_with_server/VANANZ.lineages.pheno.geno.Rda")

