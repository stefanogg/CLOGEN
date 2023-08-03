# Repeat pure genomics model 5 times (100 times on the server)

# Library and working directory
library(tidyverse)
library(tidymodels)
library(themis)
library(vip)
library(patchwork)

# Prepare input data
f <- list.files("input_data", pattern = "*pheno.geno.Rda", full.names = T)
pheno_geno <- readRDS(f)

for (i in 1:100){
  
  set.seed(i)
  print(glue::glue("Starting replicate with seed {i}"))
  dir <- str_c("seed_", i)
  dir.create(dir)
  
  setwd(dir)
  
  # Data splitting
  splits <- rsample::initial_split(pheno_geno, prop = .8, strata = mortality)
  df_train <- training(splits)
  df_test <- testing(splits)
  
  # Construct workflow
  rf_rec <- recipes::recipe(mortality ~ ., data = df_train) %>%
    themis::step_downsample(mortality)
  rf_mod <- parsnip::rand_forest(
    mtry = tune(),
    min_n = tune(),
    trees = 500
  ) %>%
    set_mode("classification") %>%
    set_engine("ranger")
  rf_wf <- workflows::workflow() %>%
    add_recipe(rf_rec) %>%
    add_model(rf_mod) 
  
  # Resampling 
  set.seed(i)
  folds <- rsample::vfold_cv(df_train, v = 10, strata = mortality)
  folds
  
  # Train model
  doParallel::registerDoParallel()
  
  rf_grid <- grid_regular(
    mtry(range = c(1,5)),
    min_n(range = c(30,40)),
    levels = 5
  )
  rf_fit_rs <- 
    rf_wf %>% 
    tune_grid(resamples = folds,
              grid = rf_grid,
              control = control_grid(save_pred = TRUE),
              metrics = metric_set(roc_auc))
  
  # Collect metrics
  df_metrics <- tune::collect_metrics(rf_fit_rs)
  metrics_train <- df_metrics %>%
    slice_max(order_by = mean, n = 1) %>%
    transmute(dataset = "train",
              .metric,
              .estimate = mean)
  
  p1 <- rf_fit_rs %>%
    collect_metrics() %>%
    filter(.metric == "roc_auc") %>%
    mutate(min_n = factor(min_n)) %>%
    ggplot(aes(mtry, mean, color = min_n)) +
    geom_line(alpha = 0.5, size = 1.5) +
    geom_point() +
    labs(y = "AUC") +
    theme_bw()
  
  # Collect predictions
  df_roc <- rf_fit_rs %>%
    tune::collect_predictions() 
  
  p2 <- df_roc %>%
    group_by(.config) %>%
    yardstick::roc_curve(truth = mortality, .pred_Died) %>%
    autoplot()
  
  # Chose the best model and fit it 
  best <- rf_fit_rs %>%
    tune::select_best()
  final_wf <- rf_wf %>%
    tune::finalize_workflow(best) 
  final_rf <- final_wf %>%
    extract_spec_parsnip() %>%
    set_engine(engine = "ranger", importance = "impurity")
  final_wf <- final_wf %>%
    workflows::update_model(final_rf)
  final_fit <- final_wf %>%
    tune::last_fit(splits)

  # Final assessment on the test dataset
  
  ## Metrics
  df_metrics <- final_fit %>%
    collect_metrics() %>%
    filter(.metric == "roc_auc") %>%
    mutate(dataset = "test") %>%
    select(-c(.config, .estimator)) %>%
    bind_rows(metrics_train)
    
  ## Confusion matrix
  
  rf_testing_pred <- final_fit %>%
    collect_predictions()
  
  rf_testing_pred <- rf_testing_pred %>%
    mutate(dataset = "test")
  
  p3 <- rf_testing_pred %>% 
    conf_mat(mortality, .pred_class) %>% 
    autoplot(type = "heatmap")
  
  ## ROC curve
 
  df_roc <- df_roc %>%
    mutate(dataset = "train") %>%
    semi_join(best) %>%
    bind_rows(rf_testing_pred)
  
  p4 <- df_roc %>%
    group_by(dataset) %>%
    yardstick::roc_curve(truth = mortality, .pred_Died) %>%
    autoplot() +
    scale_colour_manual(values = c("red", "blue"))
  
  # Here we construct the precision recall curve
  
  p5 <- df_roc %>%
    group_by(dataset) %>%
    yardstick::pr_curve(truth = mortality, .pred_Died) %>%
    autoplot() +
    scale_colour_manual(values = c("red", "blue"))
  
  df_metrics <- df_roc %>%
    group_by(dataset) %>%
    yardstick::pr_auc(truth = mortality, .pred_Died) %>%
    select(-.estimator) %>%
    bind_rows(df_metrics) %>%
    arrange(desc(.metric), desc(dataset))
  
  
  ## Feature importance
  
  p6 <- final_fit %>%
    purrr::pluck(".workflow", 1) %>%
    workflowsets::extract_fit_parsnip() %>%
    vip::vip() +
    theme_bw()
  
  # Output 
  
  ## Models
  
  name <- str_c("ML_mortality.with_lineages.replicate_", i)
  
  dir <- str_c("models/")
  dir.create(dir)
  
  rf_fit_rs %>%
    saveRDS(str_c(dir, name, ".train_resample.Rda"))
  
  final_fit %>%
    saveRDS(str_c(dir, name, ".last_fit.Rda"))
  
  ## Metrics
  dir <- "processed_data/"
  dir.create(dir)
  
  df_metrics %>%
    write_tsv(str_c(dir, name, ".metrics.tab"))

  df_roc %>%
    write_tsv(str_c(dir, name, ".predictions.tab"))
  
  ## Importance
  df_importance <- final_fit %>%
    extract_workflow() %>%
    extract_fit_parsnip() %>%
    vi()
  
  df_importance %>%
    write_tsv(str_c(dir, name, ".importance.tab"))
  
  
  ## Figures
  dir <- "figures/"
  dir.create(dir)
  
  p7 <- wrap_plots(p1, p2, p3, p4, p5, p6, nrow = 3)
  ggsave(str_c(dir, name, ".plots.pdf"), plot = p7, width = 11.9, height = 7.18, units = "in")

    
  # Go back to the working directory
  setwd("..")
}



