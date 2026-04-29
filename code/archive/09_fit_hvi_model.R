# ================================================================================================
# HVI 2.0 | 07_fit_hvi_model_cv.R
# Fit endpoint-specific HVI models with grouped k-fold CV by community
#
# Inputs:
#   - hvi_model_matrix_2019_2022.csv
#   - hvi_endpoint_metadata.csv
#
# Outputs:
#   - cv_metrics_by_endpoint.csv
#   - cv_predictions_all.csv
#   - hvi_scores_ca_year_endpoint.csv
#   - hvi_scores_ca_year_composite.csv
#   - hvi_model_coefficients.csv
#
# Modeling strategy:
#   outcome ~ heat_dose + z_vuln vars + heat_dose:z_vuln vars + seasonality + DOW + year
#
# Validation:
#   grouped k-fold CV by community
# ================================================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(mgcv)
  library(lubridate)
  library(janitor)
  library(readr)
  library(stringr)
})

source(file.path(if (dir.exists("code")) "code" else ".", "00_config.R"))

# -----------------------------
# CONFIG
# -----------------------------
project_dir <- HVI_PATHS$private

model_matrix_path <- file.path(project_dir, "hvi_model_matrix_2019_2022.csv")
endpoint_meta_path <- file.path(project_dir, "hvi_endpoint_metadata.csv")

legacy_out_dir <- HVI_PATHS$private_outputs$model_outputs
hvi_dir_create(legacy_out_dir)
out_cv_metrics      <- file.path(legacy_out_dir, "cv_metrics_by_endpoint.csv")
out_cv_predictions  <- file.path(legacy_out_dir, "cv_predictions_all.csv")
out_hvi_endpoint    <- file.path(legacy_out_dir, "hvi_scores_ca_year_endpoint.csv")
out_hvi_composite   <- file.path(legacy_out_dir, "hvi_scores_ca_year_composite.csv")
out_coef            <- file.path(legacy_out_dir, "hvi_model_coefficients.csv")

# Cross-validation
k_folds <- 5
seed    <- 20260401

# optional count family
# choices: "quasipoisson", "nb"
family_choice <- "quasipoisson"

# Minimum mean events/day required to fit endpoint
min_total_events <- 50

# Use offset(log(pop_offset)) only if pop_offset exists and is positive
use_offset_if_available <- TRUE

# Restrict vulnerability variables used in interactions if desired
# If NULL, all z_ vars in model matrix will be used
interaction_vuln_vars <- NULL

# -----------------------------
# Helper functions
# -----------------------------
make_group_folds <- function(groups, k = 5, seed = 1) {
  set.seed(seed)
  ug <- unique(groups)
  ug <- sample(ug, length(ug))
  fold_id <- rep(seq_len(k), length.out = length(ug))
  tibble(group = ug, fold = fold_id)
}

poisson_deviance <- function(obs, pred, eps = 1e-10) {
  pred <- pmax(pred, eps)
  obs  <- pmax(obs, 0)
  term <- ifelse(obs == 0, 0, obs * log(obs / pred))
  2 * sum(term - (obs - pred), na.rm = TRUE)
}

rescale_0_100 <- function(x) {
  rng <- range(x, na.rm = TRUE)
  if (anyNA(rng) || diff(rng) == 0) return(rep(50, length(x)))
  100 * (x - rng[1]) / diff(rng)
}

safe_cor <- function(x, y, method = "spearman") {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 3) return(NA_real_)
  suppressWarnings(cor(x[ok], y[ok], method = method))
}

get_interaction_betas <- function(model, heat_var, z_vars) {
  cf <- coef(model)
  nm <- names(cf)
  
  out <- tibble(
    z_var = z_vars,
    beta_interaction = map_dbl(z_vars, function(zv) {
      cand1 <- paste0(heat_var, ":", zv)
      cand2 <- paste0(zv, ":", heat_var)
      hit <- nm[nm %in% c(cand1, cand2)]
      if (length(hit) == 0) return(0)
      unname(cf[hit[1]])
    })
  )
  out
}

build_formula <- function(outcome_var, heat_var, z_vars, use_offset = FALSE) {
  main_z   <- paste(z_vars, collapse = " + ")
  inter_z  <- paste0(heat_var, ":", z_vars, collapse = " + ")
  
  rhs <- paste(
    c(
      heat_var,
      main_z,
      inter_z,
      "s(doy, bs = 'cc', k = 10)",
      "factor(dow)",
      "factor(year)"
    ),
    collapse = " + "
  )
  
  if (use_offset) {
    as.formula(paste0(outcome_var, " ~ ", rhs, " + offset(log(pop_offset))"))
  } else {
    as.formula(paste0(outcome_var, " ~ ", rhs))
  }
}

fit_one_model <- function(dat, form, family_choice = "quasipoisson") {
  if (family_choice == "nb") {
    mgcv::bam(
      formula = form,
      data = dat,
      family = nb(),
      method = "fREML",
      discrete = TRUE
    )
  } else {
    mgcv::bam(
      formula = form,
      data = dat,
      family = quasipoisson(link = "log"),
      method = "fREML",
      discrete = TRUE
    )
  }
}

# -----------------------------
# 1. Read inputs
# -----------------------------
dat <- read_csv(model_matrix_path, show_col_types = FALSE) %>%
  clean_names() %>%
  mutate(
    date = as.Date(date),
    community = as.character(community),
    year = as.integer(year),
    doy = as.integer(doy),
    dow = factor(dow, levels = c("Sun", "Mon", "Tue", "Wed", "Thu", "Fri", "Sat"))
  )

endpoint_meta <- read_csv(endpoint_meta_path, show_col_types = FALSE) %>%
  clean_names()

z_vars_all <- names(dat)[str_detect(names(dat), "^z_")]

if (is.null(interaction_vuln_vars)) {
  z_vars_use <- z_vars_all
} else {
  z_vars_use <- interaction_vuln_vars[interaction_vuln_vars %in% z_vars_all]
}

if (length(z_vars_use) == 0) {
  stop("No z_ vulnerability variables found for interactions.")
}

message("Using z variables in model:")
print(z_vars_use)

# -----------------------------
# 2. Create grouped folds by community
# -----------------------------
fold_tbl <- make_group_folds(dat$community, k = k_folds, seed = seed)

dat <- dat %>%
  left_join(fold_tbl, by = c("community" = "group"))

# -----------------------------
# 3. Loop over endpoints
# -----------------------------
cv_metrics_list <- list()
cv_pred_list    <- list()
coef_list       <- list()
hvi_list        <- list()

for (i in seq_len(nrow(endpoint_meta))) {
  
  endpoint_key <- endpoint_meta$endpoint_key[i]
  outcome_var  <- endpoint_meta$panel_outcome_col[i]
  heat_var     <- paste0("heat_dose__", endpoint_key)
  
  if (!outcome_var %in% names(dat)) next
  if (!heat_var %in% names(dat)) next
  
  message("\n==========================")
  message("Endpoint: ", endpoint_key)
  message("Outcome : ", outcome_var)
  message("Heat var: ", heat_var)
  
  dat_ep <- dat %>%
    select(
      community, date, year, doy, dow, fold,
      pop_offset,
      all_of(outcome_var),
      all_of(heat_var),
      all_of(z_vars_use)
    ) %>%
    rename(outcome = all_of(outcome_var), heat_dose = all_of(heat_var)) %>%
    filter(!is.na(outcome), !is.na(heat_dose))
  
  # Require at least some information
  if (sum(dat_ep$outcome, na.rm = TRUE) < min_total_events) {
    message("Skipping: too few total events.")
    next
  }
  
  use_offset <- FALSE
  if (use_offset_if_available &&
      "pop_offset" %in% names(dat_ep) &&
      sum(is.finite(dat_ep$pop_offset) & dat_ep$pop_offset > 0, na.rm = TRUE) > 0) {
    dat_ep <- dat_ep %>%
      mutate(pop_offset = ifelse(is.na(pop_offset) | pop_offset <= 0, NA_real_, pop_offset)) %>%
      drop_na(pop_offset)
    use_offset <- nrow(dat_ep) > 0
  }
  
  if (nrow(dat_ep) == 0) next
  
  # ---- Cross-validation ----
  endpoint_cv_metrics <- list()
  endpoint_cv_preds   <- list()
  
  for (f in seq_len(k_folds)) {
    train_dat <- dat_ep %>% filter(fold != f)
    test_dat  <- dat_ep %>% filter(fold == f)
    
    if (nrow(train_dat) == 0 || nrow(test_dat) == 0) next
    if (sum(train_dat$outcome, na.rm = TRUE) < min_total_events) next
    
    form <- build_formula(
      outcome_var = "outcome",
      heat_var    = "heat_dose",
      z_vars      = z_vars_use,
      use_offset  = use_offset
    )
    
    fit <- fit_one_model(train_dat, form, family_choice = family_choice)
    
    test_dat <- test_dat %>%
      mutate(pred = as.numeric(predict(fit, newdata = test_dat, type = "response")))
    
    # fold metrics
    fold_metrics <- tibble(
      endpoint_key = endpoint_key,
      fold = f,
      n_test = nrow(test_dat),
      obs_sum = sum(test_dat$outcome, na.rm = TRUE),
      pred_sum = sum(test_dat$pred, na.rm = TRUE),
      mae = mean(abs(test_dat$outcome - test_dat$pred), na.rm = TRUE),
      rmse = sqrt(mean((test_dat$outcome - test_dat$pred)^2, na.rm = TRUE)),
      poisson_deviance = poisson_deviance(test_dat$outcome, test_dat$pred),
      cor_daily_spearman = safe_cor(test_dat$outcome, test_dat$pred, method = "spearman")
    )
    
    # community-level ranking performance inside fold
    by_comm <- test_dat %>%
      group_by(community) %>%
      summarise(
        obs = sum(outcome, na.rm = TRUE),
        pred = sum(pred, na.rm = TRUE),
        .groups = "drop"
      )
    
    fold_metrics <- fold_metrics %>%
      mutate(
        cor_comm_spearman = safe_cor(by_comm$obs, by_comm$pred, method = "spearman"),
        cor_comm_pearson  = safe_cor(by_comm$obs, by_comm$pred, method = "pearson")
      )
    
    endpoint_cv_metrics[[f]] <- fold_metrics
    endpoint_cv_preds[[f]] <- test_dat %>%
      mutate(endpoint_key = endpoint_key, fold = f)
  }
  
  cv_metrics_list[[endpoint_key]] <- bind_rows(endpoint_cv_metrics)
  cv_pred_list[[endpoint_key]]    <- bind_rows(endpoint_cv_preds)
  
  # ---- Refit on full endpoint dataset ----
  form_full <- build_formula(
    outcome_var = "outcome",
    heat_var    = "heat_dose",
    z_vars      = z_vars_use,
    use_offset  = use_offset
  )
  
  fit_full <- fit_one_model(dat_ep, form_full, family_choice = family_choice)
  
  # store coefficients
  coef_tbl <- tibble(
    endpoint_key = endpoint_key,
    term = names(coef(fit_full)),
    estimate = as.numeric(coef(fit_full))
  )
  coef_list[[endpoint_key]] <- coef_tbl
  
  # Extract interaction coefficients and compute endpoint-specific HVI
  int_betas <- get_interaction_betas(
    model = fit_full,
    heat_var = "heat_dose",
    z_vars = z_vars_use
  ) %>%
    mutate(endpoint_key = endpoint_key)
  
  # CA-year HVI
  base_ep <- dat %>%
    select(community, year, all_of(z_vars_use)) %>%
    distinct()
  
  beta_vec <- int_betas$beta_interaction
  names(beta_vec) <- int_betas$z_var
  
  base_ep <- base_ep %>%
    rowwise() %>%
    mutate(
      hvi_raw = sum(c_across(all_of(z_vars_use)) * beta_vec[z_vars_use], na.rm = TRUE)
    ) %>%
    ungroup() %>%
    mutate(
      hvi_0_100 = rescale_0_100(hvi_raw),
      endpoint_key = endpoint_key
    )
  
  hvi_list[[endpoint_key]] <- base_ep
}

# -----------------------------
# 4. Bind outputs
# -----------------------------
cv_metrics_all <- bind_rows(cv_metrics_list)
cv_preds_all   <- bind_rows(cv_pred_list)
coef_all       <- bind_rows(coef_list)
hvi_endpoint   <- bind_rows(hvi_list)

# composite HVI across endpoints
hvi_composite <- hvi_endpoint %>%
  group_by(community, year) %>%
  summarise(
    hvi_raw_mean = mean(hvi_raw, na.rm = TRUE),
    hvi_0_100_mean = mean(hvi_0_100, na.rm = TRUE),
    n_endpoints = sum(!is.na(hvi_raw)),
    .groups = "drop"
  ) %>%
  mutate(
    hvi_composite_0_100 = rescale_0_100(hvi_raw_mean)
  )

# summarize CV performance
cv_summary <- cv_metrics_all %>%
  group_by(endpoint_key) %>%
  summarise(
    folds = n(),
    mean_mae = mean(mae, na.rm = TRUE),
    mean_rmse = mean(rmse, na.rm = TRUE),
    mean_poisson_deviance = mean(poisson_deviance, na.rm = TRUE),
    mean_cor_daily_spearman = mean(cor_daily_spearman, na.rm = TRUE),
    mean_cor_comm_spearman = mean(cor_comm_spearman, na.rm = TRUE),
    mean_cor_comm_pearson = mean(cor_comm_pearson, na.rm = TRUE),
    .groups = "drop"
  )

# -----------------------------
# 5. Write outputs
# -----------------------------
write_csv(cv_summary, out_cv_metrics)
write_csv(cv_preds_all, out_cv_predictions)
write_csv(hvi_endpoint, out_hvi_endpoint)
write_csv(hvi_composite, out_hvi_composite)
write_csv(coef_all, out_coef)

message("\nDone.")
message("CV metrics: ", out_cv_metrics)
message("CV predictions: ", out_cv_predictions)
message("Endpoint HVI scores: ", out_hvi_endpoint)
message("Composite HVI: ", out_hvi_composite)
message("Model coefficients: ", out_coef)

message("\nCV summary:")
print(cv_summary)
