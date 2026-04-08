# ================================================================================================
# HVI 2.0 | 09a_fit_endpoint_models.R
# Fit endpoint-specific heat-vulnerability models and export performance, coefficients,
# model objects, and endpoint-specific variable maps.
#
# ASSUMES THESE OBJECTS ARE ALREADY IN MEMORY:
#   - hvi_model_matrix
#   - hvi_endpoint_metadata
#
# Preferred selected-variable sources:
#   - selected_variables_final
#   - final_selected_vulnerability_vars
#   - saved artifacts from 08_variable_selection_hvi.R in:
#       {project_dir}/variable_selection/
# ================================================================================================

# -----------------------------
# CONFIG
# -----------------------------
project_dir <- "C:/Users/Peter Graffy/Box/HVI2.0"
code_dir <- file.path(project_dir, "code")


out_dir <- ensure_output_dir(project_dir, "09_model_outputs")

model_matrix_obj  <- "hvi_model_matrix"
endpoint_meta_obj <- "hvi_endpoint_metadata"

family_choice <- "quasipoisson"   # "nb" or "quasipoisson"
k_folds_spatial <- 5
seed <- 20260402
min_total_events <- 50
use_offset_if_available <- TRUE

doy_k <- 10
perform_temporal_cv <- TRUE
perform_spatial_cv  <- TRUE
save_model_rds      <- TRUE

# -----------------------------
# LOAD OBJECTS
# -----------------------------
if (!exists(model_matrix_obj, envir = .GlobalEnv)) {
  stop("Object not found: ", model_matrix_obj)
}
if (!exists(endpoint_meta_obj, envir = .GlobalEnv)) {
  stop("Object not found: ", endpoint_meta_obj)
}

hvi_model_matrix <- get(model_matrix_obj, envir = .GlobalEnv) %>%
  clean_names() %>%
  mutate(
    date = as.Date(date),
    community = as.character(community),
    year = factor(as.character(year)),
    doy = as.integer(doy),
    dow = factor(as.character(dow), levels = c("Sun", "Mon", "Tue", "Wed", "Thu", "Fri", "Sat"))
  ) %>%
  as.data.frame()

hvi_endpoint_metadata <- get(endpoint_meta_obj, envir = .GlobalEnv) %>%
  clean_names()

if (!"community" %in% names(hvi_model_matrix)) {
  stop("hvi_model_matrix is missing community column.")
}

global_selected_vars <- infer_global_selected_vars(
  hvi_model_matrix,
  project_dir = project_dir
)

if (length(global_selected_vars) == 0) {
  stop("No selected vulnerability variables were identified.")
}

fold_tbl <- make_group_folds(
  hvi_model_matrix$community,
  k = k_folds_spatial,
  seed = seed
)

hvi_model_matrix <- hvi_model_matrix %>%
  left_join(fold_tbl, by = "community")

# -----------------------------
# STORAGE OBJECTS
# -----------------------------
spatial_metrics_list <- list()
temporal_metrics_list <- list()
spatial_preds_list <- list()
temporal_preds_list <- list()
coef_list <- list()
interaction_beta_list <- list()
main_beta_list <- list()
model_list <- list()
endpoint_var_map <- list()
fit_log <- list()

# -----------------------------
# FIT LOOP
# -----------------------------
for (i in seq_len(nrow(hvi_endpoint_metadata))) {
  ep_key      <- hvi_endpoint_metadata$endpoint_key[i]
  outcome_var <- hvi_endpoint_metadata$panel_outcome_col[i]
  heat_var    <- get_heat_var(ep_key, hvi_model_matrix)
  
  if (!outcome_var %in% names(hvi_model_matrix)) {
    fit_log[[length(fit_log) + 1]] <- tibble(
      endpoint_key = ep_key,
      status = "skipped",
      reason = "missing outcome column",
      n_rows = NA_integer_,
      n_events = NA_real_,
      n_vars = NA_integer_,
      vars_used = NA_character_
    )
    next
  }
  
  if (length(heat_var) == 0 || is.na(heat_var)) {
    fit_log[[length(fit_log) + 1]] <- tibble(
      endpoint_key = ep_key,
      status = "skipped",
      reason = "missing heat_dose column",
      n_rows = NA_integer_,
      n_events = NA_real_,
      n_vars = NA_integer_,
      vars_used = NA_character_
    )
    next
  }
  
  z_vars_use <- infer_endpoint_selected_vars(
    endpoint_key = ep_key,
    dat = hvi_model_matrix,
    fallback_vars = global_selected_vars,
    project_dir = project_dir
  )
  
  z_vars_use <- z_vars_use[z_vars_use %in% names(hvi_model_matrix)]
  
  if (length(z_vars_use) == 0) {
    fit_log[[length(fit_log) + 1]] <- tibble(
      endpoint_key = ep_key,
      status = "skipped",
      reason = "no selected z vars",
      n_rows = NA_integer_,
      n_events = NA_real_,
      n_vars = 0L,
      vars_used = NA_character_
    )
    next
  }
  
  endpoint_var_map[[ep_key]] <- tibble(
    endpoint_key = ep_key,
    variable = z_vars_use
  )
  
  dat_ep <- hvi_model_matrix %>%
    select(
      community, date, year, doy, dow, fold, pop_offset,
      all_of(outcome_var), all_of(heat_var), all_of(z_vars_use)
    ) %>%
    rename(
      outcome = all_of(outcome_var),
      heat_dose = all_of(heat_var)
    ) %>%
    filter(!is.na(outcome), !is.na(heat_dose))
  
  if (sum(dat_ep$outcome, na.rm = TRUE) < min_total_events) {
    fit_log[[length(fit_log) + 1]] <- tibble(
      endpoint_key = ep_key,
      status = "skipped",
      reason = "too few total events",
      n_rows = nrow(dat_ep),
      n_events = sum(dat_ep$outcome, na.rm = TRUE),
      n_vars = length(z_vars_use),
      vars_used = paste(z_vars_use, collapse = "; ")
    )
    next
  }
  
  use_offset <- FALSE
  if (
    use_offset_if_available &&
    "pop_offset" %in% names(dat_ep) &&
    any(is.finite(dat_ep$pop_offset) & dat_ep$pop_offset > 0, na.rm = TRUE)
  ) {
    dat_ep <- dat_ep %>%
      mutate(
        pop_offset = ifelse(is.na(pop_offset) | pop_offset <= 0, NA_real_, pop_offset)
      ) %>%
      drop_na(pop_offset)
    
    use_offset <- nrow(dat_ep) > 0
  }
  
  dat_ep <- dat_ep %>%
    drop_na(all_of(z_vars_use), doy, dow, year) %>%
    mutate(
      community = as.character(community),
      date = as.Date(date),
      doy = as.integer(doy),
      year = factor(as.character(year)),
      dow = factor(as.character(dow), levels = c("Sun", "Mon", "Tue", "Wed", "Thu", "Fri", "Sat"))
    ) %>%
    as.data.frame()
  
  if (nrow(dat_ep) == 0) {
    fit_log[[length(fit_log) + 1]] <- tibble(
      endpoint_key = ep_key,
      status = "skipped",
      reason = "no complete cases after filtering",
      n_rows = 0L,
      n_events = NA_real_,
      n_vars = length(z_vars_use),
      vars_used = paste(z_vars_use, collapse = "; ")
    )
    next
  }
  
  form_spatial_full <- build_formula(
    outcome_var  = "outcome",
    heat_var     = "heat_dose",
    z_vars       = z_vars_use,
    use_offset   = use_offset,
    include_year = TRUE,
    doy_k        = doy_k
  )
  
  form_temporal <- build_formula(
    outcome_var  = "outcome",
    heat_var     = "heat_dose",
    z_vars       = z_vars_use,
    use_offset   = use_offset,
    include_year = FALSE,
    doy_k        = doy_k
  )
  
  # ---- Spatial CV ----
  if (perform_spatial_cv) {
    for (f in seq_len(k_folds_spatial)) {
      train_dat <- dat_ep %>% filter(fold != f)
      test_dat  <- dat_ep %>% filter(fold == f)
      
      if (nrow(train_dat) == 0 || nrow(test_dat) == 0) next
      if (sum(train_dat$outcome, na.rm = TRUE) < min_total_events) next
      
      train_dat <- train_dat %>%
        mutate(
          year = factor(as.character(year), levels = levels(dat_ep$year)),
          dow  = factor(as.character(dow), levels = levels(dat_ep$dow))
        ) %>%
        as.data.frame()
      
      test_dat <- test_dat %>%
        mutate(
          year = factor(as.character(year), levels = levels(dat_ep$year)),
          dow  = factor(as.character(dow), levels = levels(dat_ep$dow))
        ) %>%
        as.data.frame()
      
      fit_cv <- tryCatch(
        fit_endpoint_model(
          dat = train_dat,
          form = form_spatial_full,
          family_choice = family_choice
        ),
        error = function(e) NULL
      )
      
      if (is.null(fit_cv)) next
      
      test_dat$pred <- as.numeric(
        predict(fit_cv, newdata = test_dat, type = "response")
      )
      
      by_comm <- test_dat %>%
        group_by(community) %>%
        summarise(
          obs = sum(outcome, na.rm = TRUE),
          pred = sum(pred, na.rm = TRUE),
          .groups = "drop"
        )
      
      spatial_metrics_list[[length(spatial_metrics_list) + 1]] <- tibble(
        endpoint_key = ep_key,
        cv_type = "spatial",
        fold = f,
        n_test = nrow(test_dat),
        obs_sum = sum(test_dat$outcome, na.rm = TRUE),
        pred_sum = sum(test_dat$pred, na.rm = TRUE),
        mae = mean(abs(test_dat$outcome - test_dat$pred), na.rm = TRUE),
        rmse = sqrt(mean((test_dat$outcome - test_dat$pred)^2, na.rm = TRUE)),
        poisson_deviance = poisson_deviance(test_dat$outcome, test_dat$pred),
        cor_daily_spearman = safe_cor(test_dat$outcome, test_dat$pred, method = "spearman"),
        cor_comm_spearman = safe_cor(by_comm$obs, by_comm$pred, method = "spearman"),
        cor_comm_pearson = safe_cor(by_comm$obs, by_comm$pred, method = "pearson")
      )
      
      spatial_preds_list[[length(spatial_preds_list) + 1]] <- test_dat %>%
        mutate(
          endpoint_key = ep_key,
          cv_type = "spatial",
          split_id = paste0("fold_", f)
        )
    }
  }
  
  # ---- Temporal CV ----
  if (perform_temporal_cv) {
    years_holdout <- sort(unique(as.character(dat_ep$year)))
    
    for (yy in years_holdout) {
      train_dat <- dat_ep %>% filter(as.character(year) != yy)
      test_dat  <- dat_ep %>% filter(as.character(year) == yy)
      
      if (nrow(train_dat) == 0 || nrow(test_dat) == 0) next
      if (sum(train_dat$outcome, na.rm = TRUE) < min_total_events) next
      
      train_dat <- train_dat %>%
        mutate(
          dow = factor(as.character(dow), levels = levels(dat_ep$dow))
        ) %>%
        as.data.frame()
      
      test_dat <- test_dat %>%
        mutate(
          dow = factor(as.character(dow), levels = levels(dat_ep$dow))
        ) %>%
        as.data.frame()
      
      fit_cv <- tryCatch(
        fit_endpoint_model(
          dat = train_dat,
          form = form_temporal,
          family_choice = family_choice
        ),
        error = function(e) NULL
      )
      
      if (is.null(fit_cv)) next
      
      test_dat$pred <- as.numeric(
        predict(fit_cv, newdata = test_dat, type = "response")
      )
      
      by_comm <- test_dat %>%
        group_by(community) %>%
        summarise(
          obs = sum(outcome, na.rm = TRUE),
          pred = sum(pred, na.rm = TRUE),
          .groups = "drop"
        )
      
      temporal_metrics_list[[length(temporal_metrics_list) + 1]] <- tibble(
        endpoint_key = ep_key,
        cv_type = "temporal",
        fold = yy,
        n_test = nrow(test_dat),
        obs_sum = sum(test_dat$outcome, na.rm = TRUE),
        pred_sum = sum(test_dat$pred, na.rm = TRUE),
        mae = mean(abs(test_dat$outcome - test_dat$pred), na.rm = TRUE),
        rmse = sqrt(mean((test_dat$outcome - test_dat$pred)^2, na.rm = TRUE)),
        poisson_deviance = poisson_deviance(test_dat$outcome, test_dat$pred),
        cor_daily_spearman = safe_cor(test_dat$outcome, test_dat$pred, method = "spearman"),
        cor_comm_spearman = safe_cor(by_comm$obs, by_comm$pred, method = "spearman"),
        cor_comm_pearson = safe_cor(by_comm$obs, by_comm$pred, method = "pearson")
      )
      
      temporal_preds_list[[length(temporal_preds_list) + 1]] <- test_dat %>%
        mutate(
          endpoint_key = ep_key,
          cv_type = "temporal",
          split_id = as.character(yy)
        )
    }
  }
  
  # ---- Full fit ----
  dat_ep_full <- dat_ep %>%
    mutate(
      year = factor(as.character(year), levels = levels(dat_ep$year)),
      dow  = factor(as.character(dow), levels = levels(dat_ep$dow))
    ) %>%
    as.data.frame()
  
  fit_full <- tryCatch(
    fit_endpoint_model(
      dat = dat_ep_full,
      form = form_spatial_full,
      family_choice = family_choice
    ),
    error = function(e) NULL
  )
  
  if (is.null(fit_full)) {
    fit_log[[length(fit_log) + 1]] <- tibble(
      endpoint_key = ep_key,
      status = "skipped",
      reason = "full model fit failed",
      n_rows = nrow(dat_ep_full),
      n_events = sum(dat_ep_full$outcome, na.rm = TRUE),
      n_vars = length(z_vars_use),
      vars_used = paste(z_vars_use, collapse = "; ")
    )
    next
  }
  
  model_list[[ep_key]] <- fit_full
  coef_list[[ep_key]] <- extract_coef_table(fit_full, ep_key)
  
  interaction_beta_list[[ep_key]] <- extract_interaction_betas(
    fit_full,
    ep_key,
    heat_var = "heat_dose",
    z_vars = z_vars_use
  )
  
  main_beta_list[[ep_key]] <- extract_main_effect_betas(
    fit_full,
    ep_key,
    z_vars = z_vars_use
  )
  
  fit_log[[length(fit_log) + 1]] <- tibble(
    endpoint_key = ep_key,
    status = "fit",
    reason = NA_character_,
    n_rows = nrow(dat_ep_full),
    n_events = sum(dat_ep_full$outcome, na.rm = TRUE),
    n_vars = length(z_vars_use),
    vars_used = paste(z_vars_use, collapse = "; ")
  )
}

# -----------------------------
# BIND / SUMMARIZE
# -----------------------------
endpoint_var_map <- bind_rows(endpoint_var_map)
fit_log <- bind_rows(fit_log)
spatial_metrics <- bind_rows(spatial_metrics_list)
temporal_metrics <- bind_rows(temporal_metrics_list)
spatial_predictions <- bind_rows(spatial_preds_list)
temporal_predictions <- bind_rows(temporal_preds_list)
endpoint_coefficients <- bind_rows(coef_list)
endpoint_interaction_betas <- bind_rows(interaction_beta_list)
endpoint_main_betas <- bind_rows(main_beta_list)

spatial_summary <- if (nrow(spatial_metrics) > 0) {
  spatial_metrics %>%
    group_by(endpoint_key) %>%
    summarise(
      cv_type = "spatial",
      folds = n(),
      mean_mae = mean(mae, na.rm = TRUE),
      mean_rmse = mean(rmse, na.rm = TRUE),
      mean_poisson_deviance = mean(poisson_deviance, na.rm = TRUE),
      mean_cor_daily_spearman = mean(cor_daily_spearman, na.rm = TRUE),
      mean_cor_comm_spearman = mean(cor_comm_spearman, na.rm = TRUE),
      mean_cor_comm_pearson = mean(cor_comm_pearson, na.rm = TRUE),
      .groups = "drop"
    )
} else {
  tibble()
}

temporal_summary <- if (nrow(temporal_metrics) > 0) {
  temporal_metrics %>%
    group_by(endpoint_key) %>%
    summarise(
      cv_type = "temporal",
      folds = n(),
      mean_mae = mean(mae, na.rm = TRUE),
      mean_rmse = mean(rmse, na.rm = TRUE),
      mean_poisson_deviance = mean(poisson_deviance, na.rm = TRUE),
      mean_cor_daily_spearman = mean(cor_daily_spearman, na.rm = TRUE),
      mean_cor_comm_spearman = mean(cor_comm_spearman, na.rm = TRUE),
      mean_cor_comm_pearson = mean(cor_comm_pearson, na.rm = TRUE),
      .groups = "drop"
    )
} else {
  tibble()
}

endpoint_model_performance <- bind_rows(spatial_summary, temporal_summary) %>%
  left_join(
    hvi_endpoint_metadata %>% select(endpoint_key, outcome_label, source, domain),
    by = "endpoint_key"
  ) %>%
  relocate(outcome_label, source, domain, .after = endpoint_key)

# -----------------------------
# EXPORT
# -----------------------------
write_csv(fit_log, file.path(out_dir, "09a_fit_log.csv"))
write_csv(endpoint_var_map, file.path(out_dir, "09a_endpoint_variable_map.csv"))
write_csv(spatial_metrics, file.path(out_dir, "09a_cv_metrics_spatial_by_split.csv"))
write_csv(temporal_metrics, file.path(out_dir, "09a_cv_metrics_temporal_by_split.csv"))
write_csv(endpoint_model_performance, file.path(out_dir, "09a_endpoint_model_performance.csv"))
write_csv(spatial_predictions, file.path(out_dir, "09a_cv_predictions_spatial.csv"))
write_csv(temporal_predictions, file.path(out_dir, "09a_cv_predictions_temporal.csv"))
write_csv(endpoint_coefficients, file.path(out_dir, "09a_endpoint_coefficients.csv"))
write_csv(endpoint_interaction_betas, file.path(out_dir, "09a_endpoint_interaction_betas.csv"))
write_csv(endpoint_main_betas, file.path(out_dir, "09a_endpoint_main_betas.csv"))

if (save_model_rds) {
  saveRDS(model_list, file.path(out_dir, "09a_endpoint_models.rds"))
}

# -----------------------------
# SAVE TO ENVIRONMENT
# -----------------------------
assign("endpoint_models", model_list, envir = .GlobalEnv)
assign("endpoint_var_map", endpoint_var_map, envir = .GlobalEnv)
assign("endpoint_model_performance", endpoint_model_performance, envir = .GlobalEnv)
assign("endpoint_interaction_betas", endpoint_interaction_betas, envir = .GlobalEnv)
assign("endpoint_main_betas", endpoint_main_betas, envir = .GlobalEnv)
assign("endpoint_fit_log", fit_log, envir = .GlobalEnv)

message("09a complete. Outputs written to: ", out_dir)
message("Endpoints fit successfully: ", sum(fit_log$status == "fit", na.rm = TRUE))
message("Endpoints skipped: ", sum(fit_log$status == "skipped", na.rm = TRUE))