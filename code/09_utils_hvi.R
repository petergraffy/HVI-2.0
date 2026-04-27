# ================================================================================================
# HVI 2.0 | 09_utils_hvi.R
# Shared utilities for endpoint fitting, structural scoring, temperature-grid scoring,
# and operational daily risk scoring.
# ================================================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(mgcv)
  library(janitor)
  library(readr)
  library(stringr)
  library(lubridate)
})

source(file.path(if (dir.exists("code")) "code" else ".", "00_config.R"))

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

poisson_deviance <- function(obs, pred, eps = 1e-10) {
  pred <- pmax(pred, eps)
  obs  <- pmax(obs, 0)
  term <- ifelse(obs == 0, 0, obs * log(obs / pred))
  2 * sum(term - (obs - pred), na.rm = TRUE)
}

make_group_folds <- function(groups, k = 5, seed = 1) {
  set.seed(seed)
  ug <- unique(as.character(groups))
  ug <- sample(ug, length(ug))
  tibble(community = ug, fold = rep(seq_len(k), length.out = length(ug)))
}

infer_selected_var_table <- function(project_dir = NULL, subdir = "variable_selection") {
  if (exists("selected_variables_final", envir = .GlobalEnv)) {
    return(get("selected_variables_final", envir = .GlobalEnv))
  }
  
  if (!is.null(project_dir)) {
    rds_path <- file.path(project_dir, subdir, "selected_variables_final.rds")
    csv_path <- file.path(project_dir, subdir, "selected_variables_final.csv")
    
    if (file.exists(rds_path)) {
      return(readRDS(rds_path))
    }
    
    if (file.exists(csv_path)) {
      return(readr::read_csv(csv_path, show_col_types = FALSE))
    }
  }
  
  NULL
}

infer_global_selected_vars <- function(dat, project_dir = NULL, subdir = "variable_selection") {
  sv <- infer_selected_var_table(project_dir = project_dir, subdir = subdir)
  
  if (!is.null(sv) && all(c("variable", "selected_final") %in% names(sv))) {
    vars <- sv %>%
      filter(selected_final %in% TRUE) %>%
      pull(variable) %>%
      unique()
    vars <- vars[vars %in% names(dat)]
    if (length(vars) > 0) return(vars)
  }
  
  if (exists("final_selected_vulnerability_vars", envir = .GlobalEnv)) {
    vars <- get("final_selected_vulnerability_vars", envir = .GlobalEnv)
    vars <- vars[vars %in% names(dat)]
    if (length(vars) > 0) return(vars)
  }
  
  if (!is.null(project_dir)) {
    rds_path <- file.path(project_dir, subdir, "final_selected_vulnerability_vars.rds")
    csv_path <- file.path(project_dir, subdir, "final_selected_vulnerability_vars.csv")
    
    if (file.exists(rds_path)) {
      vars <- readRDS(rds_path)
      vars <- vars[vars %in% names(dat)]
      if (length(vars) > 0) return(vars)
    }
    
    if (file.exists(csv_path)) {
      vars <- readr::read_csv(csv_path, show_col_types = FALSE) %>%
        pull(variable) %>%
        unique()
      vars <- vars[vars %in% names(dat)]
      if (length(vars) > 0) return(vars)
    }
  }
  
  stop(
    "No selected vulnerability variables found. ",
    "Run 08_variable_selection_hvi.R first or load the saved selected-variable artifacts."
  )
}

infer_endpoint_selected_vars <- function(endpoint_key, dat, fallback_vars = NULL, project_dir = NULL, subdir = "variable_selection") {
  if (is.null(fallback_vars)) {
    fallback_vars <- infer_global_selected_vars(dat, project_dir = project_dir, subdir = subdir)
  }
  
  sv <- infer_selected_var_table(project_dir = project_dir, subdir = subdir)
  if (is.null(sv) || !"variable" %in% names(sv)) return(fallback_vars)
  
  candidate_endpoint_cols <- c("endpoint_key", "outcome", "endpoint")
  endpoint_col <- candidate_endpoint_cols[candidate_endpoint_cols %in% names(sv)][1]
  
  selected_col <- c("selected_final", "selected")
  selected_col <- selected_col[selected_col %in% names(sv)][1]
  
  if (length(endpoint_col) == 1 && !is.na(endpoint_col) &&
      length(selected_col) == 1 && !is.na(selected_col)) {
    vars <- sv %>%
      filter(.data[[endpoint_col]] == endpoint_key, .data[[selected_col]] %in% TRUE) %>%
      pull(variable) %>%
      unique()
    vars <- vars[vars %in% names(dat)]
    if (length(vars) > 0) return(vars)
  }
  
  if (length(selected_col) == 1 && !is.na(selected_col)) {
    vars <- sv %>%
      filter(.data[[selected_col]] %in% TRUE) %>%
      pull(variable) %>%
      unique()
    vars <- vars[vars %in% names(dat)]
    if (length(vars) > 0) return(vars)
  }
  
  fallback_vars
}

get_heat_var <- function(endpoint_key, dat) {
  c(paste0("heat_dose__", endpoint_key), paste0("heat_dose_", endpoint_key)) %>%
    keep(~ .x %in% names(dat)) %>%
    .[1]
}

get_excess_var <- function(endpoint_key, dat) {
  c(paste0("temp_excess__", endpoint_key), paste0("temp_excess_", endpoint_key)) %>%
    keep(~ .x %in% names(dat)) %>%
    .[1]
}

build_formula <- function(outcome_var,
                          heat_var,
                          z_vars,
                          use_offset = FALSE,
                          humidity_var = NULL,
                          include_year = TRUE,
                          include_humidity = FALSE,
                          doy_k = 10) {
  rhs_terms <- c(
    heat_var,
    z_vars,
    paste0(heat_var, ":", z_vars),
    paste0("s(doy, bs = 'cc', k = ", doy_k, ")"),
    "dow"
  )
  
  if (include_year) rhs_terms <- c(rhs_terms, "year")
  if (include_humidity && !is.null(humidity_var)) {
    rhs_terms <- c(rhs_terms, paste0("s(", humidity_var, ", k = 6)"))
  }
  
  rhs_terms <- rhs_terms[!is.na(rhs_terms) & nzchar(rhs_terms)]
  rhs <- paste(rhs_terms, collapse = " + ")
  
  if (use_offset) {
    as.formula(paste0(outcome_var, " ~ ", rhs, " + offset(log(pop_offset))"))
  } else {
    as.formula(paste0(outcome_var, " ~ ", rhs))
  }
}

fit_endpoint_model <- function(dat, form, family_choice = c("nb", "quasipoisson")) {
  family_choice <- match.arg(family_choice)
  dat <- as.data.frame(dat)
  
  if (family_choice == "nb") {
    mgcv::bam(
      formula = form,
      data = dat,
      family = mgcv::nb(),
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

extract_coef_table <- function(fit, endpoint_key) {
  tibble(
    endpoint_key = endpoint_key,
    term = names(coef(fit)),
    estimate = as.numeric(coef(fit))
  )
}

extract_interaction_betas <- function(fit, endpoint_key, heat_var, z_vars) {
  cf <- coef(fit)
  nm <- names(cf)
  tibble(
    endpoint_key = endpoint_key,
    variable = z_vars,
    beta_interaction = map_dbl(z_vars, function(v) {
      cand <- c(paste0(heat_var, ":", v), paste0(v, ":", heat_var))
      hit <- nm[nm %in% cand][1]
      ifelse(is.na(hit), 0, unname(cf[hit]))
    })
  )
}

extract_main_effect_betas <- function(fit, endpoint_key, z_vars) {
  cf <- coef(fit)
  nm <- names(cf)
  tibble(
    endpoint_key = endpoint_key,
    variable = z_vars,
    beta_main = map_dbl(z_vars, function(v) {
      hit <- nm[nm == v][1]
      ifelse(is.na(hit), 0, unname(cf[hit]))
    })
  )
}

score_driver_contributions <- function(base_df, beta_tbl, score_prefix = "interaction") {
  beta_col <- names(beta_tbl)[names(beta_tbl) %in% c("beta_interaction", "beta_main")][1]
  beta_vec <- beta_tbl[[beta_col]]
  names(beta_vec) <- beta_tbl$variable
  score_vars <- intersect(names(beta_vec), names(base_df))

  contrib_long <- base_df %>%
    select(community, year, all_of(score_vars)) %>%
    pivot_longer(cols = all_of(score_vars), names_to = "variable", values_to = "z_value") %>%
    mutate(
      beta = beta_vec[variable],
      contribution = z_value * beta,
      score_component = score_prefix
    )

  contrib_long
}

calc_endpoint_weight_table <- function(perf_tbl, endpoint_meta, default_weight_by_source = c(Mortality = 3, ED = 2, EMS = 2, Other = 1)) {
  perf_base <- perf_tbl %>%
    group_by(endpoint_key) %>%
    summarise(
      spatial_rank = mean(cor_comm_spearman, na.rm = TRUE),
      temporal_rank = mean(cor_comm_spearman, na.rm = TRUE),
      .groups = "drop"
    )

  endpoint_meta %>%
    mutate(source_clean = case_when(
      str_to_upper(source) == "MORTALITY" ~ "Mortality",
      str_to_upper(source) == "ED" ~ "ED",
      str_to_upper(source) == "EMS" ~ "EMS",
      TRUE ~ "Other"
    )) %>%
    left_join(perf_base, by = "endpoint_key") %>%
    mutate(
      source_weight = unname(default_weight_by_source[source_clean]),
      source_weight = ifelse(is.na(source_weight), 1, source_weight),
      performance_weight = pmax(coalesce(spatial_rank, 0), 0),
      endpoint_weight = source_weight * ifelse(performance_weight > 0, performance_weight, 0.25)
    ) %>%
    select(endpoint_key, outcome_label, source, domain, source_clean, endpoint_weight, source_weight, performance_weight)
}

make_top_driver_labels <- function(contrib_long, top_n = 3, positive_only = FALSE) {
  dat <- contrib_long
  if (positive_only) dat <- dat %>% filter(contribution > 0)

  dat %>%
    group_by(community, year, endpoint_key = if ("endpoint_key" %in% names(dat)) endpoint_key else NA_character_) %>%
    arrange(desc(abs(contribution)), .by_group = TRUE) %>%
    slice_head(n = top_n) %>%
    summarise(
      top_drivers = paste(variable, collapse = " | "),
      .groups = "drop"
    )
}

ensure_output_dir <- function(project_dir, subdir = "09_model_outputs") {
  out_dir <- file.path(project_dir, subdir)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  out_dir
}
