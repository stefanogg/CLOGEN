# Here we generate 10 train-split replicates for the GWAS features selection model

# Library and working directory
library(tidyverse)
library(tidymodels)

setwd(paste0(here::here(), "/Cell_reports_revision/ML/model1_all_genes"))

# Prepare input data. We can't use the input data of the last run because they don't have isolates id
f <- list.files("input_data", pattern = "*pheno.geno.id.Rda", full.names = T)
pheno_geno <- readRDS(f)

for (i in 1:10){
  set.seed(i)
  print(glue::glue("Generating replicate {i}"))
  name <- str_c("replicate_", formatC(i, width = 2, format = "d", flag = "0"))
  
  # Data splitting
  splits <- rsample::initial_split(pheno_geno, prop = .8, strata = mortality)
  df_train <- training(splits)
  df_test <- testing(splits)
  
  # Save output
  saveRDS(splits, str_c("split_replicates/", name, ".Rda"))
  df_train %>%
    arrange(sample_id) %>%
    transmute(sample_id, mortality = as.integer(mortality == "Died")) %>%
    write_tsv(str_c("~/Documents/Transfer_with_server/", name, ".mortality_phenotype.tab"))
}
# split <- readRDS("split_replicates/replicate_01.Rda")
