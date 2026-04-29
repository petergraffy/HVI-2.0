# ================================================================================================
# HVI 2.0 | 08_variable_selection_hvi.R
# Variable selection pipeline for baseline vulnerability modifiers of heat-related health risk
#
# REQUIRED OBJECTS IN MEMORY:
#   - hvi_model_matrix
#   - hvi_endpoint_metadata
#
# MAIN IDEA:
#   Stage 0: Global vulnerability screening
#       Drop variables with excessive missingness or near-zero variance
#
#   Stage 1: Random forest screening by endpoint
#       outcome ~ heat_dose + vulnerability variables
#
#   Stage 2: Interaction-aware GAM refinement by endpoint
#       outcome ~ heat_dose + selected vulnerability vars + heat_dose:selected vars + seasonal controls
#
#   Stage 3: Combined ranking
#       RF importance + interaction strength + sign consistency + significance
#
#   Stage 4: Correlation pruning
#       Remove redundant variables based on pairwise correlation and rank priority
#
# OUTPUTS:
#   Key tables:
#     - variable_screening_table.csv
#     - variable_drop_log_stage0_screening.csv
#     - rf_importance_by_endpoint.csv
#     - rf_importance_summary.csv
#     - rf_skip_table.csv
#     - rf_candidate_variables.csv
#     - rf_candidate_drop_table.csv
#     - interaction_terms_by_endpoint.csv
#     - interaction_summary.csv
#     - interaction_model_stats.csv
#     - interaction_skip_table.csv
#     - interaction_endpoint_variable_screen.csv
#     - variable_selection_summary.csv
#     - pre_prune_variable_table.csv
#     - pre_prune_drop_table.csv
#     - correlation_pruning_pairs.csv
#     - selected_variables_final.csv
#     - final_selected_vulnerability_vars.csv
#     - final_selected_vulnerability_vars.rds
#     - selected_variables_final.rds
#     - variable_selection_methods_table.csv
#     - selected_variable_correlation_matrix.csv
#     - selected_variable_correlation_long.csv
#
#   Figures:
#     - fig_rf_importance_top15.png
#     - fig_interaction_importance_top15.png
#     - fig_overall_selection_score_lollipop.png
#     - fig_variable_selection_heatmap.png
#     - fig_selected_variable_correlation.png
#     - fig_rf_vs_interaction_scatter.png
#
# NOTES:
#   - Uses z_ variables already present in hvi_model_matrix
#   - Designed for interpretable, publication-ready selection
#   - tmax-based heat_dose variables are assumed already present in hvi_model_matrix
# ================================================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(janitor)
  library(mgcv)
  library(ranger)
  library(glue)
  library(stringr)
  library(readr)
  library(scales)
  library(forcats)
  library(ggrepel)
})

source(file.path(if (dir.exists("code")) "code" else ".", "00_config.R"))

# -----------------------------
# CONFIG
# -----------------------------
project_dir <- HVI_PATHS$private
out_dir <- HVI_PATHS$private_outputs$variable_selection
hvi_dir_create(out_dir)

# objects already in memory
model_matrix_obj  <- "hvi_model_matrix"
endpoint_meta_obj <- "hvi_endpoint_metadata"

# RF settings
rf_num_trees <- 750
rf_min_rows  <- 500
rf_seed      <- 20260401

# GAM settings
gam_k_doy <- 10
gam_family <- quasipoisson(link = "log")
min_total_events <- 50
min_rows_interaction <- 200

# variable screening thresholds
max_missing_prop <- 0.20
min_nonzero_sd   <- 1e-8

# candidate carry-forward
rf_top_n <- 12

# final selection thresholds
min_endpoints_with_rf_signal   <- 3
min_endpoints_with_interaction <- 2
interaction_abs_beta_quantile  <- 0.65
rf_importance_quantile         <- 0.65
max_abs_correlation            <- 0.80

# optional manual drop list after ranking, before final output
manual_drop_vars <- c(
  "z_mean_unemployed"   # prefer SVI unemployment instead
)

# plotting
plot_top_n <- 15
base_size <- 13

# saved artifacts
selected_vars_rds      <- file.path(out_dir, "final_selected_vulnerability_vars.rds")
selected_vars_csv      <- file.path(out_dir, "final_selected_vulnerability_vars.csv")
selection_summary_rds  <- file.path(out_dir, "selected_variables_final.rds")

# -----------------------------
# LABELS / THEMES
# -----------------------------
var_label_dict <- c(
  "z_pop_density_km2" = "Population density (per km²)",
  "z_total_pop" = "Population size",
  "z_mean_age" = "Mean age",
  "z_median_age" = "Median age",
  "z_mean_black" = "Proportion Black",
  "z_mean_hisp" = "Proportion Hispanic",
  "z_mean_white" = "Proportion White",
  "z_mean_asian" = "Proportion Asian",
  "z_mean_income" = "Mean income",
  "z_median_income" = "Median income",
  "z_mean_unemployed" = "Proportion unemployed",
  "z_mean_employed" = "Proportion employed",
  "z_mean_college" = "College educated",
  "z_mean_hs" = "High school educated",
  "z_mean_male" = "Proportion male",
  "z_mean_female" = "Proportion female",
  "z_ndvi" = "NDVI",
  "z_mean_ndvi" = "NDVI",
  "z_ac_prob" = "Air conditioning prevalence",
  "z_ac_cbsa_rank" = "Air conditioning percentile rank",
  "z_no2" = "NO2",
  "z_pm25" = "PM2.5",
  "z_svi_rpl_theme1" = "SVI: Socioeconomic status",
  "z_svi_rpl_theme2" = "SVI: Household composition & disability",
  "z_svi_rpl_theme3" = "SVI: Minority status & language",
  "z_svi_rpl_theme4" = "SVI: Housing type & transportation",
  "z_svi_rpl_themes" = "SVI: Overall vulnerability",
  "z_svi_ep_pov" = "Poverty",
  "z_svi_ep_unemp" = "Unemployment",
  "z_svi_ep_nohsdp" = "No high school diploma",
  "z_svi_ep_age65" = "Age 65+",
  "z_svi_ep_age17" = "Age <18",
  "z_svi_ep_disabl" = "Disability",
  "z_svi_ep_sngpnt" = "Single-parent households",
  "z_svi_ep_limeng" = "Limited English proficiency",
  "z_svi_ep_minrty" = "Minority population",
  "z_svi_ep_munit" = "Multi-unit housing",
  "z_svi_ep_mobile" = "Mobile homes",
  "z_svi_ep_crowd" = "Crowding",
  "z_svi_ep_noveh" = "No vehicle access",
  "z_svi_ep_groupq" = "Group quarters",
  "z_svi_ep_uninsur" = "Uninsured",
  "z_svi_ep_noint" = "No internet access",
  "z_svi_ep_hburd" = "Housing cost burden"
)

pretty_var_label <- function(x) {
  out <- unname(var_label_dict[x])
  fallback <- x %>%
    str_remove("^z_") %>%
    str_replace_all("_", " ") %>%
    str_replace_all("\\bNdvi\\b", "NDVI") %>%
    str_replace_all("\\bNo2\\b", "NO2") %>%
    str_replace_all("\\bPm25\\b", "PM2.5") %>%
    str_replace_all("\\bAc\\b", "AC") %>%
    str_replace_all("\\bSvi\\b", "SVI") %>%
    str_to_title()
  out[is.na(out)] <- fallback[is.na(out)]
  out
}

pretty_endpoint_label <- function(x) {
  case_when(
    x == "deaths" ~ "All-cause mortality",
    x == "death_cvd" ~ "CVD mortality",
    x == "death_respiratory" ~ "Respiratory mortality",
    x == "death_renal" ~ "Renal mortality",
    x == "death_neurologic" ~ "Neurologic mortality",
    x == "death_mental" ~ "Mental health mortality",
    x == "death_gi" ~ "GI mortality",
    x == "death_injury" ~ "Injury mortality",
    x == "death_syncope" ~ "Syncope mortality",
    x == "death_dehydration" ~ "Dehydration mortality",
    x == "death_heat" ~ "Heat-related mortality",
    
    x == "ed_visits" ~ "All-cause ED visits",
    x == "ed_cvd" ~ "CVD ED visits",
    x == "ed_respiratory" ~ "Respiratory ED visits",
    x == "ed_renal" ~ "Renal ED visits",
    x == "ed_neurologic" ~ "Neurologic ED visits",
    x == "ed_mental" ~ "Mental health ED visits",
    x == "ed_gi" ~ "GI ED visits",
    x == "ed_injury" ~ "Injury ED visits",
    x == "ed_syncope" ~ "Syncope ED visits",
    x == "ed_dehydration" ~ "Dehydration ED visits",
    x == "ed_heat" ~ "Heat-related ED visits",
    
    x == "ems_calls" ~ "All-cause EMS calls",
    x == "ems_cvd" ~ "CVD EMS calls",
    x == "ems_respiratory" ~ "Respiratory EMS calls",
    x == "ems_neuro" ~ "Neurologic EMS calls",
    x == "ems_mental" ~ "Mental health EMS calls",
    x == "ems_gi" ~ "GI EMS calls",
    x == "ems_bleeding" ~ "Bleeding EMS calls",
    x == "ems_injury" ~ "Injury EMS calls",
    x == "ems_syncope" ~ "Syncope EMS calls",
    x == "ems_heat" ~ "Heat-related EMS calls",
    TRUE ~ x
  )
}

source_from_endpoint <- function(x) {
  case_when(
    str_starts(x, "death") | x == "deaths" ~ "Mortality",
    str_starts(x, "ed_") | x == "ed_visits" ~ "ED",
    str_starts(x, "ems_") | x == "ems_calls" ~ "EMS",
    TRUE ~ "Other"
  )
}

theme_pub <- function(base_size = 12) {
  theme_minimal(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_line(linewidth = 0.3, color = "grey85"),
      panel.grid.major.x = element_blank(),
      axis.title = element_text(face = "bold"),
      plot.title = element_text(face = "bold"),
      plot.subtitle = element_text(color = "grey30"),
      strip.text = element_text(face = "bold"),
      legend.position = "bottom"
    )
}

# -----------------------------
# HELPER FUNCTIONS
# -----------------------------
rescale_0_100 <- function(x) {
  rng <- range(x, na.rm = TRUE)
  if (anyNA(rng) || diff(rng) == 0) return(rep(50, length(x)))
  100 * (x - rng[1]) / diff(rng)
}

build_rf_formula <- function(outcome_var, heat_var, vuln_vars) {
  as.formula(
    paste0(outcome_var, " ~ ", paste(c(heat_var, vuln_vars), collapse = " + "))
  )
}

build_gam_formula <- function(outcome_var, heat_var, vuln_vars, use_offset = FALSE) {
  rhs_terms <- c(
    heat_var,
    vuln_vars,
    paste0(heat_var, ":", vuln_vars),
    glue("s(doy, bs = 'cc', k = {gam_k_doy})"),
    "factor(dow)",
    "factor(year)"
  )
  rhs <- paste(rhs_terms, collapse = " + ")
  
  if (use_offset) {
    as.formula(paste0(outcome_var, " ~ ", rhs, " + offset(log(pop_offset))"))
  } else {
    as.formula(paste0(outcome_var, " ~ ", rhs))
  }
}

clean_analysis_df <- function(dat, outcome_var, heat_var, vuln_vars, use_offset = TRUE) {
  keep <- c("community", "date", "year", "doy", "dow", "pop_offset", outcome_var, heat_var, vuln_vars)
  keep <- keep[keep %in% names(dat)]
  
  out <- dat %>%
    select(all_of(keep)) %>%
    rename(
      outcome = all_of(outcome_var),
      heat_dose = all_of(heat_var)
    )
  
  if ("pop_offset" %in% names(out) && use_offset) {
    out <- out %>%
      mutate(pop_offset = ifelse(is.na(pop_offset) | pop_offset <= 0, NA_real_, pop_offset))
  }
  
  out
}

median_impute_vars <- function(df, vars) {
  for (v in vars) {
    if (v %in% names(df)) {
      med <- median(df[[v]], na.rm = TRUE)
      if (!is.finite(med)) med <- 0
      df[[v]][is.na(df[[v]])] <- med
    }
  }
  df
}

screen_vulnerability_vars <- function(dat, vars, max_missing_prop = 0.20, min_nonzero_sd = 1e-8) {
  tibble(variable = vars) %>%
    mutate(
      missing_prop = map_dbl(variable, ~ mean(is.na(dat[[.x]]))),
      sd_val       = map_dbl(variable, ~ sd(dat[[.x]], na.rm = TRUE)),
      keep_missing = missing_prop <= max_missing_prop,
      keep_sd      = is.finite(sd_val) & sd_val > min_nonzero_sd,
      keep         = keep_missing & keep_sd,
      drop_reason = case_when(
        !keep_missing & !keep_sd ~ glue("dropped: missing_prop>{max_missing_prop} and near-zero variance"),
        !keep_missing ~ glue("dropped: missing_prop>{max_missing_prop}"),
        !keep_sd ~ glue("dropped: sd<={min_nonzero_sd}"),
        TRUE ~ "kept"
      )
    ) %>%
    arrange(desc(keep), missing_prop, desc(sd_val))
}

run_rf_screen <- function(dat, outcome_var, heat_var, vuln_vars, num_trees = 500, seed = 1) {
  df <- clean_analysis_df(dat, outcome_var, heat_var, vuln_vars, use_offset = FALSE) %>%
    drop_na(outcome, heat_dose)
  
  for (v in vuln_vars) {
    if (v %in% names(df)) {
      med <- median(df[[v]], na.rm = TRUE)
      if (!is.finite(med)) med <- 0
      df[[v]][is.na(df[[v]])] <- med
    }
  }
  
  if (nrow(df) < rf_min_rows) return(NULL)
  
  set.seed(seed)
  
  fit <- ranger(
    formula = build_rf_formula("outcome", "heat_dose", vuln_vars),
    data = df,
    num.trees = num_trees,
    importance = "permutation",
    seed = seed,
    respect.unordered.factors = "order"
  )
  
  tibble(
    variable = names(fit$variable.importance),
    importance = as.numeric(fit$variable.importance)
  ) %>%
    arrange(desc(importance))
}

extract_interaction_terms <- function(fit, heat_var, endpoint_key) {
  s <- summary(fit)
  
  if (is.null(s$p.table)) return(tibble())
  
  pt <- as.data.frame(s$p.table, stringsAsFactors = FALSE)
  pt$term <- rownames(pt)
  rownames(pt) <- NULL
  names(pt) <- make.names(names(pt))
  
  est_col  <- names(pt)[names(pt) %in% c("Estimate")]
  se_col   <- names(pt)[names(pt) %in% c("Std..Error", "Std.Error")]
  stat_col <- names(pt)[names(pt) %in% c("t.value", "z.value")]
  p_col    <- names(pt)[names(pt) %in% c("Pr...t..", "Pr...z..", "p.value")]
  
  if (length(est_col) == 0 || length(se_col) == 0 || length(stat_col) == 0 || length(p_col) == 0) {
    return(tibble())
  }
  
  pt %>%
    transmute(
      term = .data$term,
      estimate = as.numeric(.data[[est_col[1]]]),
      std.error = as.numeric(.data[[se_col[1]]]),
      statistic = as.numeric(.data[[stat_col[1]]]),
      p.value = as.numeric(.data[[p_col[1]]])
    ) %>%
    filter(
      str_detect(term, paste0("^", heat_var, ":")) |
        str_detect(term, paste0(":", heat_var, "$"))
    ) %>%
    mutate(
      variable = case_when(
        str_detect(term, paste0("^", heat_var, ":")) ~ str_remove(term, paste0("^", heat_var, ":")),
        str_detect(term, paste0(":", heat_var, "$")) ~ str_remove(term, paste0(":", heat_var, "$")),
        TRUE ~ term
      ),
      endpoint_key = endpoint_key
    ) %>%
    select(endpoint_key, variable, term, estimate, std.error, statistic, p.value)
}

prune_correlated_vars <- function(dat, vars, ranking_tbl, max_abs_correlation = 0.80) {
  vars <- vars[vars %in% names(dat)]
  if (length(vars) <= 1) {
    return(list(
      keep = vars,
      dropped_pairs = tibble()
    ))
  }
  
  corr_mat <- cor(dat[, vars], use = "pairwise.complete.obs")
  diag(corr_mat) <- 0
  
  ranked <- ranking_tbl %>%
    filter(variable %in% vars) %>%
    arrange(desc(overall_score), desc(mean_abs_beta), desc(mean_rf_importance)) %>%
    pull(variable)
  
  keep <- character(0)
  dropped_pairs <- list()
  
  for (v in ranked) {
    if (length(keep) == 0) {
      keep <- c(keep, v)
    } else {
      cor_with_keep <- abs(corr_mat[v, keep])
      cor_with_keep <- cor_with_keep[is.finite(cor_with_keep)]
      
      if (length(cor_with_keep) == 0 || max(cor_with_keep) < max_abs_correlation) {
        keep <- c(keep, v)
      } else {
        blocking_var <- names(which.max(abs(corr_mat[v, keep])))[1]
        dropped_pairs[[length(dropped_pairs) + 1]] <- tibble(
          dropped_variable = v,
          kept_variable = blocking_var,
          abs_correlation = abs(corr_mat[v, blocking_var]),
          reason = glue("dropped: abs(correlation)>={max_abs_correlation} with higher-ranked variable")
        )
      }
    }
  }
  
  list(
    keep = keep,
    dropped_pairs = bind_rows(dropped_pairs)
  )
}

# -----------------------------
# 1. LOAD OBJECTS
# -----------------------------
if (!exists(model_matrix_obj, envir = .GlobalEnv)) {
  stop("Object not found: ", model_matrix_obj)
}
if (!exists(endpoint_meta_obj, envir = .GlobalEnv)) {
  stop("Object not found: ", endpoint_meta_obj)
}

hvi_model_matrix <- get(model_matrix_obj, envir = .GlobalEnv) %>%
  mutate(
    date = as.Date(date),
    year = as.integer(year),
    doy = as.integer(doy),
    dow = factor(as.character(dow), levels = c("Sun", "Mon", "Tue", "Wed", "Thu", "Fri", "Sat"))
  )

hvi_endpoint_metadata <- get(endpoint_meta_obj, envir = .GlobalEnv)

# -----------------------------
# 2. IDENTIFY AND SCREEN VULNERABILITY VARIABLES
# -----------------------------
z_vuln_vars <- names(hvi_model_matrix)[str_detect(names(hvi_model_matrix), "^z_")]

if (length(z_vuln_vars) == 0) {
  stop("No z_ vulnerability variables found in hvi_model_matrix.")
}

var_screen <- screen_vulnerability_vars(
  dat = hvi_model_matrix,
  vars = z_vuln_vars,
  max_missing_prop = max_missing_prop,
  min_nonzero_sd = min_nonzero_sd
) %>%
  mutate(
    variable_label = pretty_var_label(variable),
    stage = "stage0_global_screening"
  )

vuln_vars_screened <- var_screen %>%
  filter(keep) %>%
  pull(variable)

var_drop_stage0 <- var_screen %>%
  filter(!keep) %>%
  select(stage, variable, variable_label, missing_prop, sd_val, drop_reason)

if (length(vuln_vars_screened) == 0) {
  stop("No vulnerability variables survived missingness / variance screening.")
}

write_csv(var_screen, file.path(out_dir, "variable_screening_table.csv"))
write_csv(var_drop_stage0, file.path(out_dir, "variable_drop_log_stage0_screening.csv"))

# -----------------------------
# 3. RANDOM FOREST SCREENING BY ENDPOINT
# -----------------------------
rf_results <- list()
rf_skip_log <- list()

for (i in seq_len(nrow(hvi_endpoint_metadata))) {
  ep_key      <- hvi_endpoint_metadata$endpoint_key[i]
  outcome_var <- hvi_endpoint_metadata$panel_outcome_col[i]
  
  heat_var_candidates <- c(
    paste0("heat_dose__", ep_key),
    paste0("heat_dose_", ep_key)
  )
  heat_var <- heat_var_candidates[heat_var_candidates %in% names(hvi_model_matrix)][1]
  
  if (!outcome_var %in% names(hvi_model_matrix)) {
    rf_skip_log[[length(rf_skip_log) + 1]] <- tibble(
      endpoint_key = ep_key,
      outcome_var = outcome_var,
      heat_var = NA_character_,
      reason = "outcome_var not found in hvi_model_matrix"
    )
    next
  }
  
  if (length(heat_var) == 0 || is.na(heat_var)) {
    rf_skip_log[[length(rf_skip_log) + 1]] <- tibble(
      endpoint_key = ep_key,
      outcome_var = outcome_var,
      heat_var = paste(heat_var_candidates, collapse = " OR "),
      reason = "heat_var not found in hvi_model_matrix"
    )
    next
  }
  
  df_tmp <- hvi_model_matrix %>%
    select(any_of(c(outcome_var, heat_var, vuln_vars_screened))) %>%
    rename(
      outcome = all_of(outcome_var),
      heat_dose = all_of(heat_var)
    )
  
  n_complete <- df_tmp %>%
    select(outcome, heat_dose) %>%
    drop_na() %>%
    nrow()
  
  total_events <- sum(df_tmp$outcome, na.rm = TRUE)
  outcome_sd   <- sd(df_tmp$outcome, na.rm = TRUE)
  
  if (n_complete < rf_min_rows) {
    rf_skip_log[[length(rf_skip_log) + 1]] <- tibble(
      endpoint_key = ep_key,
      outcome_var = outcome_var,
      heat_var = heat_var,
      reason = paste0("too few complete rows for RF: ", n_complete)
    )
    next
  }
  
  if (!is.finite(total_events) || total_events <= 0) {
    rf_skip_log[[length(rf_skip_log) + 1]] <- tibble(
      endpoint_key = ep_key,
      outcome_var = outcome_var,
      heat_var = heat_var,
      reason = "outcome has no positive events"
    )
    next
  }
  
  if (!is.finite(outcome_sd) || outcome_sd <= 0) {
    rf_skip_log[[length(rf_skip_log) + 1]] <- tibble(
      endpoint_key = ep_key,
      outcome_var = outcome_var,
      heat_var = heat_var,
      reason = "outcome has zero variance"
    )
    next
  }
  
  message("RF screening: ", ep_key)
  
  imp <- tryCatch(
    run_rf_screen(
      dat = hvi_model_matrix,
      outcome_var = outcome_var,
      heat_var = heat_var,
      vuln_vars = vuln_vars_screened,
      num_trees = rf_num_trees,
      seed = rf_seed
    ),
    error = function(e) {
      rf_skip_log[[length(rf_skip_log) + 1]] <<- tibble(
        endpoint_key = ep_key,
        outcome_var = outcome_var,
        heat_var = heat_var,
        reason = paste("run_rf_screen error:", conditionMessage(e))
      )
      NULL
    }
  )
  
  if (is.null(imp) || nrow(imp) == 0) {
    rf_skip_log[[length(rf_skip_log) + 1]] <- tibble(
      endpoint_key = ep_key,
      outcome_var = outcome_var,
      heat_var = heat_var,
      reason = "run_rf_screen returned NULL or empty result"
    )
    next
  }
  
  rf_results[[ep_key]] <- imp %>%
    transmute(
      endpoint_key = ep_key,
      source = source_from_endpoint(ep_key),
      outcome_var = outcome_var,
      heat_var = heat_var,
      variable = as.character(variable),
      importance = as.numeric(importance)
    )
}

rf_skip_table <- bind_rows(rf_skip_log)

if (length(rf_results) == 0) {
  message("\nNo RF results were produced. Skip log:")
  print(rf_skip_table, n = Inf)
  stop("Random forest stage produced no usable endpoint results.")
}

rf_importance_by_endpoint <- bind_rows(rf_results) %>%
  mutate(
    variable_label = pretty_var_label(variable),
    endpoint_label = pretty_endpoint_label(endpoint_key)
  ) %>%
  arrange(endpoint_key, desc(importance))

assign("rf_skip_table", rf_skip_table, envir = .GlobalEnv)
assign("rf_importance_by_endpoint", rf_importance_by_endpoint, envir = .GlobalEnv)

rf_importance_summary <- rf_importance_by_endpoint %>%
  group_by(variable) %>%
  summarise(
    mean_rf_importance = mean(importance, na.rm = TRUE),
    median_rf_importance = median(importance, na.rm = TRUE),
    max_rf_importance = max(importance, na.rm = TRUE),
    sd_rf_importance = sd(importance, na.rm = TRUE),
    n_endpoints_rf = n_distinct(endpoint_key),
    n_positive_importance = sum(importance > 0, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(variable_label = pretty_var_label(variable)) %>%
  arrange(desc(mean_rf_importance))

write_csv(rf_importance_by_endpoint, file.path(out_dir, "rf_importance_by_endpoint.csv"))
write_csv(rf_importance_summary, file.path(out_dir, "rf_importance_summary.csv"))
write_csv(rf_skip_table, file.path(out_dir, "rf_skip_table.csv"))

# RF carry-forward logic
rf_cut_mean <- quantile(rf_importance_summary$mean_rf_importance, rf_importance_quantile, na.rm = TRUE)

rf_candidate_table <- rf_importance_summary %>%
  mutate(
    pass_endpoint_count = n_positive_importance >= min_endpoints_with_rf_signal,
    pass_importance_cut = mean_rf_importance >= rf_cut_mean,
    selected_rf_candidate = pass_endpoint_count & pass_importance_cut
  ) %>%
  arrange(desc(selected_rf_candidate), desc(mean_rf_importance))

rf_candidate_vars <- rf_candidate_table %>%
  filter(selected_rf_candidate) %>%
  slice_head(n = rf_top_n) %>%
  pull(variable)

if (length(rf_candidate_vars) < 3) {
  rf_candidate_vars <- rf_importance_summary %>%
    slice_head(n = min(8, n())) %>%
    pull(variable)
}

rf_candidate_drop_table <- rf_candidate_table %>%
  mutate(
    stage = "stage1_rf_candidate_selection",
    drop_reason = case_when(
      variable %in% rf_candidate_vars ~ "kept for interaction stage",
      !pass_endpoint_count & !pass_importance_cut ~ "dropped: insufficient endpoint recurrence and low RF importance",
      !pass_endpoint_count ~ "dropped: insufficient endpoint recurrence",
      !pass_importance_cut ~ "dropped: below RF importance threshold",
      TRUE ~ "dropped: truncated by top-N carry-forward rule"
    )
  ) %>%
  filter(!(variable %in% rf_candidate_vars)) %>%
  select(stage, variable, variable_label, mean_rf_importance, n_positive_importance, drop_reason)

write_csv(rf_candidate_table, file.path(out_dir, "rf_candidate_variables.csv"))
write_csv(rf_candidate_drop_table, file.path(out_dir, "rf_candidate_drop_table.csv"))

# -----------------------------
# 4. INTERACTION-AWARE GAM REFINEMENT
# -----------------------------
interaction_results <- list()
interaction_fit_stats <- list()
interaction_skip_log <- list()
interaction_var_screen_log <- list()

hvi_model_matrix <- hvi_model_matrix %>%
  mutate(
    date = as.Date(date),
    year = as.integer(year),
    doy = lubridate::yday(date),
    dow = factor(
      lubridate::wday(date, label = TRUE, abbr = TRUE),
      levels = c("Sun", "Mon", "Tue", "Wed", "Thu", "Fri", "Sat")
    )
  )

for (i in seq_len(nrow(hvi_endpoint_metadata))) {
  ep_key      <- hvi_endpoint_metadata$endpoint_key[i]
  outcome_var <- hvi_endpoint_metadata$panel_outcome_col[i]
  
  heat_var_candidates <- c(
    paste0("heat_dose__", ep_key),
    paste0("heat_dose_", ep_key)
  )
  heat_var <- heat_var_candidates[heat_var_candidates %in% names(hvi_model_matrix)][1]
  
  if (!outcome_var %in% names(hvi_model_matrix)) {
    interaction_skip_log[[length(interaction_skip_log) + 1]] <- tibble(
      endpoint_key = ep_key,
      outcome_var = outcome_var,
      heat_var = NA_character_,
      reason = "outcome_var not found in hvi_model_matrix"
    )
    next
  }
  
  if (length(heat_var) == 0 || is.na(heat_var)) {
    interaction_skip_log[[length(interaction_skip_log) + 1]] <- tibble(
      endpoint_key = ep_key,
      outcome_var = outcome_var,
      heat_var = paste(heat_var_candidates, collapse = " OR "),
      reason = "heat_var not found in hvi_model_matrix"
    )
    next
  }
  
  df_tmp <- clean_analysis_df(
    dat = hvi_model_matrix,
    outcome_var = outcome_var,
    heat_var = heat_var,
    vuln_vars = rf_candidate_vars,
    use_offset = TRUE
  ) %>%
    drop_na(outcome, heat_dose, doy, dow, year)
  
  use_offset_tmp <- "pop_offset" %in% names(df_tmp) &&
    sum(is.finite(df_tmp$pop_offset) & df_tmp$pop_offset > 0, na.rm = TRUE) > 0
  
  if (use_offset_tmp) {
    df_tmp <- df_tmp %>%
      mutate(pop_offset = ifelse(is.na(pop_offset) | pop_offset <= 0, NA_real_, pop_offset)) %>%
      drop_na(pop_offset)
  }
  
  vars_here_initial <- rf_candidate_vars[rf_candidate_vars %in% names(df_tmp)]
  
  if (length(vars_here_initial) == 0) {
    interaction_skip_log[[length(interaction_skip_log) + 1]] <- tibble(
      endpoint_key = ep_key,
      outcome_var = outcome_var,
      heat_var = heat_var,
      reason = "no rf_candidate_vars found in analysis data"
    )
    next
  }
  
  for (v in vars_here_initial) {
    df_tmp[[v]] <- suppressWarnings(as.numeric(df_tmp[[v]]))
  }
  
  df_tmp <- median_impute_vars(df_tmp, vars_here_initial)
  
  var_tbl <- tibble(variable = vars_here_initial) %>%
    mutate(
      endpoint_key = ep_key,
      variable_label = pretty_var_label(variable),
      n_unique = map_int(variable, ~ dplyr::n_distinct(df_tmp[[.x]][is.finite(df_tmp[[.x]])])),
      sd_val   = map_dbl(variable, ~ sd(df_tmp[[.x]], na.rm = TRUE)),
      keep     = is.finite(sd_val) & sd_val > 0 & n_unique >= 2,
      drop_reason = case_when(
        keep ~ "kept for endpoint-specific GAM",
        !is.finite(sd_val) ~ "dropped: non-finite SD after imputation",
        sd_val <= 0 ~ "dropped: zero variance after imputation",
        n_unique < 2 ~ "dropped: fewer than 2 unique values after imputation",
        TRUE ~ "dropped: failed endpoint-specific variable screen"
      )
    )
  
  interaction_var_screen_log[[length(interaction_var_screen_log) + 1]] <- var_tbl
  
  vars_here <- var_tbl %>%
    filter(keep) %>%
    pull(variable)
  
  if (length(vars_here) == 0) {
    interaction_skip_log[[length(interaction_skip_log) + 1]] <- tibble(
      endpoint_key = ep_key,
      outcome_var = outcome_var,
      heat_var = heat_var,
      reason = "no candidate variables with non-zero variance"
    )
    next
  }
  
  n_rows <- nrow(df_tmp)
  total_events <- sum(df_tmp$outcome, na.rm = TRUE)
  outcome_sd <- sd(df_tmp$outcome, na.rm = TRUE)
  heat_sd <- sd(df_tmp$heat_dose, na.rm = TRUE)
  
  if (n_rows < min_rows_interaction) {
    interaction_skip_log[[length(interaction_skip_log) + 1]] <- tibble(
      endpoint_key = ep_key,
      outcome_var = outcome_var,
      heat_var = heat_var,
      reason = paste0("too few rows after filtering/imputation: ", n_rows)
    )
    next
  }
  
  if (!is.finite(total_events) || total_events < min_total_events) {
    interaction_skip_log[[length(interaction_skip_log) + 1]] <- tibble(
      endpoint_key = ep_key,
      outcome_var = outcome_var,
      heat_var = heat_var,
      reason = paste0("too few total events: ", total_events)
    )
    next
  }
  
  if (!is.finite(outcome_sd) || outcome_sd <= 0) {
    interaction_skip_log[[length(interaction_skip_log) + 1]] <- tibble(
      endpoint_key = ep_key,
      outcome_var = outcome_var,
      heat_var = heat_var,
      reason = "outcome has zero variance"
    )
    next
  }
  
  if (!is.finite(heat_sd) || heat_sd <= 0) {
    interaction_skip_log[[length(interaction_skip_log) + 1]] <- tibble(
      endpoint_key = ep_key,
      outcome_var = outcome_var,
      heat_var = heat_var,
      reason = "heat_dose has zero variance"
    )
    next
  }
  
  message("Interaction GAM: ", ep_key)
  
  form <- build_gam_formula(
    outcome_var = "outcome",
    heat_var = "heat_dose",
    vuln_vars = vars_here,
    use_offset = use_offset_tmp
  )
  
  fit <- tryCatch(
    mgcv::bam(
      formula = form,
      data = df_tmp,
      family = gam_family,
      method = "fREML",
      discrete = TRUE
    ),
    error = function(e) {
      interaction_skip_log[[length(interaction_skip_log) + 1]] <<- tibble(
        endpoint_key = ep_key,
        outcome_var = outcome_var,
        heat_var = heat_var,
        reason = paste("bam error:", conditionMessage(e))
      )
      NULL
    }
  )
  
  if (is.null(fit)) next
  
  interaction_terms <- tryCatch(
    extract_interaction_terms(
      fit = fit,
      heat_var = "heat_dose",
      endpoint_key = ep_key
    ) %>%
      mutate(
        source = source_from_endpoint(ep_key),
        variable_label = pretty_var_label(variable),
        endpoint_label = pretty_endpoint_label(endpoint_key)
      ),
    error = function(e) {
      interaction_skip_log[[length(interaction_skip_log) + 1]] <<- tibble(
        endpoint_key = ep_key,
        outcome_var = outcome_var,
        heat_var = heat_var,
        reason = paste("extract_interaction_terms error:", conditionMessage(e))
      )
      NULL
    }
  )
  
  if (is.null(interaction_terms) || nrow(interaction_terms) == 0) {
    interaction_skip_log[[length(interaction_skip_log) + 1]] <- tibble(
      endpoint_key = ep_key,
      outcome_var = outcome_var,
      heat_var = heat_var,
      reason = "model fit succeeded but no interaction terms were extracted"
    )
    next
  }
  
  interaction_results[[ep_key]] <- interaction_terms
  
  interaction_fit_stats[[ep_key]] <- tibble(
    endpoint_key = ep_key,
    endpoint_label = pretty_endpoint_label(ep_key),
    source = source_from_endpoint(ep_key),
    n = n_rows,
    total_events = total_events,
    n_vars = length(vars_here),
    vars_used = paste(vars_here, collapse = "; "),
    use_offset = use_offset_tmp,
    deviance_explained = summary(fit)$dev.expl,
    aic = tryCatch(AIC(fit), error = function(e) NA_real_)
  )
}

interaction_terms_by_endpoint <- bind_rows(interaction_results)
interaction_model_stats <- bind_rows(interaction_fit_stats)
interaction_skip_table <- bind_rows(interaction_skip_log)
interaction_endpoint_variable_screen <- bind_rows(interaction_var_screen_log)

assign("interaction_terms_by_endpoint", interaction_terms_by_endpoint, envir = .GlobalEnv)
assign("interaction_model_stats", interaction_model_stats, envir = .GlobalEnv)
assign("interaction_skip_table", interaction_skip_table, envir = .GlobalEnv)

if (nrow(interaction_terms_by_endpoint) == 0) {
  message("\nNo interaction GAM results were produced. Skip log:")
  print(interaction_skip_table, n = Inf)
  stop("Interaction GAM stage produced no results.")
}

interaction_summary <- interaction_terms_by_endpoint %>%
  group_by(variable) %>%
  summarise(
    mean_beta = mean(estimate, na.rm = TRUE),
    mean_abs_beta = mean(abs(estimate), na.rm = TRUE),
    median_abs_beta = median(abs(estimate), na.rm = TRUE),
    max_abs_beta = max(abs(estimate), na.rm = TRUE),
    prop_positive = mean(estimate > 0, na.rm = TRUE),
    sign_consistency = abs(prop_positive - 0.5) * 2,
    n_endpoints_interaction = n_distinct(endpoint_key),
    n_p_lt_0_05 = sum(p.value < 0.05, na.rm = TRUE),
    mean_p = mean(p.value, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(variable_label = pretty_var_label(variable)) %>%
  arrange(desc(mean_abs_beta))

write_csv(interaction_terms_by_endpoint, file.path(out_dir, "interaction_terms_by_endpoint.csv"))
write_csv(interaction_summary, file.path(out_dir, "interaction_summary.csv"))
write_csv(interaction_model_stats, file.path(out_dir, "interaction_model_stats.csv"))
write_csv(interaction_skip_table, file.path(out_dir, "interaction_skip_table.csv"))
write_csv(interaction_endpoint_variable_screen, file.path(out_dir, "interaction_endpoint_variable_screen.csv"))

# -----------------------------
# 5. COMBINE RF + INTERACTION EVIDENCE
# -----------------------------
selection_summary <- rf_importance_summary %>%
  full_join(interaction_summary, by = c("variable", "variable_label")) %>%
  mutate(
    mean_rf_importance = replace_na(mean_rf_importance, 0),
    median_rf_importance = replace_na(median_rf_importance, 0),
    max_rf_importance = replace_na(max_rf_importance, 0),
    sd_rf_importance = replace_na(sd_rf_importance, 0),
    n_endpoints_rf = replace_na(n_endpoints_rf, 0),
    n_positive_importance = replace_na(n_positive_importance, 0),
    
    mean_beta = replace_na(mean_beta, 0),
    mean_abs_beta = replace_na(mean_abs_beta, 0),
    median_abs_beta = replace_na(median_abs_beta, 0),
    max_abs_beta = replace_na(max_abs_beta, 0),
    prop_positive = replace_na(prop_positive, 0.5),
    sign_consistency = replace_na(sign_consistency, 0),
    n_endpoints_interaction = replace_na(n_endpoints_interaction, 0),
    n_p_lt_0_05 = replace_na(n_p_lt_0_05, 0),
    mean_p = replace_na(mean_p, 1),
    
    rf_score = rescale_0_100(mean_rf_importance),
    interaction_score = rescale_0_100(mean_abs_beta),
    consistency_score = 100 * sign_consistency,
    significance_score = rescale_0_100(n_p_lt_0_05),
    
    overall_score = 0.40 * rf_score +
      0.40 * interaction_score +
      0.10 * consistency_score +
      0.10 * significance_score
  ) %>%
  arrange(desc(overall_score))

interaction_beta_cut <- quantile(selection_summary$mean_abs_beta, interaction_abs_beta_quantile, na.rm = TRUE)
rf_importance_cut    <- quantile(selection_summary$mean_rf_importance, rf_importance_quantile, na.rm = TRUE)

pre_prune_table <- selection_summary %>%
  mutate(
    pass_endpoint_recurrence =
      n_endpoints_rf >= min_endpoints_with_rf_signal |
      n_endpoints_interaction >= min_endpoints_with_interaction,
    pass_effect_size =
      mean_abs_beta >= interaction_beta_cut |
      mean_rf_importance >= rf_importance_cut,
    selected_pre_prune = pass_endpoint_recurrence & pass_effect_size,
    stage = "stage3_combined_selection",
    drop_reason = case_when(
      selected_pre_prune ~ "kept for correlation pruning",
      !pass_endpoint_recurrence & !pass_effect_size ~ "dropped: insufficient endpoint recurrence and weak RF/interaction signal",
      !pass_endpoint_recurrence ~ "dropped: insufficient endpoint recurrence",
      !pass_effect_size ~ "dropped: below RF and interaction thresholds",
      TRUE ~ "dropped: did not satisfy combined selection rule"
    )
  ) %>%
  arrange(desc(selected_pre_prune), desc(overall_score))

pre_prune_vars <- pre_prune_table %>%
  filter(selected_pre_prune) %>%
  pull(variable)

if (length(pre_prune_vars) < 5) {
  pre_prune_vars <- selection_summary %>%
    slice_head(n = min(8, n())) %>%
    pull(variable)
  
  pre_prune_table <- pre_prune_table %>%
    mutate(
      selected_pre_prune = variable %in% pre_prune_vars,
      drop_reason = case_when(
        variable %in% pre_prune_vars ~ "kept by fallback top-N rule",
        TRUE ~ drop_reason
      )
    )
}

# continue from:
# if (length(pre_prune_vars) < 5) { ... }

pre_prune_table <- pre_prune_table %>%
  mutate(
    variable_label = pretty_var_label(variable)
  )

pre_prune_drop_table <- pre_prune_table %>%
  filter(!selected_pre_prune) %>%
  select(
    stage, variable, variable_label,
    n_endpoints_rf, n_endpoints_interaction,
    mean_rf_importance, mean_abs_beta, overall_score,
    drop_reason
  )

write_csv(pre_prune_table, file.path(out_dir, "pre_prune_variable_table.csv"))
write_csv(pre_prune_drop_table, file.path(out_dir, "pre_prune_drop_table.csv"))

# -----------------------------
# 6. CORRELATION PRUNING
# -----------------------------
prune_out <- prune_correlated_vars(
  dat = hvi_model_matrix,
  vars = pre_prune_vars,
  ranking_tbl = selection_summary,
  max_abs_correlation = max_abs_correlation
)

vars_after_corr_prune <- prune_out$keep
correlation_pruning_pairs <- prune_out$dropped_pairs %>%
  mutate(
    stage = "stage4_correlation_pruning",
    dropped_variable_label = pretty_var_label(dropped_variable),
    kept_variable_label = pretty_var_label(kept_variable)
  ) %>%
  select(
    stage,
    dropped_variable, dropped_variable_label,
    kept_variable, kept_variable_label,
    abs_correlation, reason
  )

write_csv(
  correlation_pruning_pairs,
  file.path(out_dir, "correlation_pruning_pairs.csv")
)

# -----------------------------
# 7. MANUAL DROP STAGE
# -----------------------------
manual_drop_table <- tibble(
  variable = manual_drop_vars
) %>%
  filter(variable %in% vars_after_corr_prune) %>%
  mutate(
    variable_label = pretty_var_label(variable),
    stage = "stage5_manual_drop",
    drop_reason = "dropped manually after ranking/pruning"
  )

final_selected_vars <- setdiff(vars_after_corr_prune, manual_drop_vars)

# if manual drops remove too many variables, keep correlation-pruned set
if (length(final_selected_vars) == 0) {
  final_selected_vars <- vars_after_corr_prune
  manual_drop_table <- manual_drop_table %>%
    mutate(
      drop_reason = "manual drop requested but ignored because it would remove all final variables"
    )
}

# -----------------------------
# 8. FINAL SELECTION TABLES
# -----------------------------
selected_variables_final <- selection_summary %>%
  mutate(
    variable_label = pretty_var_label(variable),
    selected_pre_prune = variable %in% pre_prune_vars,
    selected_after_corr_prune = variable %in% vars_after_corr_prune,
    selected_final = variable %in% final_selected_vars,
    dropped_stage0_screening = !(variable %in% vuln_vars_screened),
    dropped_stage1_rf = variable %in% rf_candidate_drop_table$variable,
    dropped_stage3_combined = variable %in% pre_prune_drop_table$variable,
    dropped_stage4_correlation = variable %in% correlation_pruning_pairs$dropped_variable,
    dropped_stage5_manual = variable %in% manual_drop_table$variable,
    final_status = case_when(
      selected_final ~ "selected_final",
      dropped_stage5_manual ~ "dropped_stage5_manual",
      dropped_stage4_correlation ~ "dropped_stage4_correlation",
      dropped_stage3_combined ~ "dropped_stage3_combined",
      dropped_stage1_rf ~ "dropped_stage1_rf",
      dropped_stage0_screening ~ "dropped_stage0_screening",
      TRUE ~ "not_selected"
    )
  ) %>%
  arrange(desc(selected_final), desc(overall_score))

write_csv(
  selection_summary,
  file.path(out_dir, "variable_selection_summary.csv")
)

write_csv(
  selected_variables_final,
  file.path(out_dir, "selected_variables_final.csv")
)

readr::write_csv(
  tibble(
    variable = final_selected_vars,
    variable_label = pretty_var_label(final_selected_vars)
  ),
  selected_vars_csv
)

saveRDS(final_selected_vars, selected_vars_rds)
saveRDS(selected_variables_final, selection_summary_rds)

# publication-facing methods table
variable_selection_methods_table <- selected_variables_final %>%
  transmute(
    variable,
    label = variable_label,
    mean_rf_importance = round(mean_rf_importance, 4),
    mean_abs_interaction_beta = round(mean_abs_beta, 4),
    sign_consistency = round(sign_consistency, 3),
    n_endpoints_rf,
    n_endpoints_interaction,
    n_p_lt_0_05,
    overall_score = round(overall_score, 1),
    selected_pre_prune,
    selected_after_corr_prune,
    selected_final,
    final_status
  ) %>%
  arrange(desc(selected_final), desc(overall_score))

write_csv(
  variable_selection_methods_table,
  file.path(out_dir, "variable_selection_methods_table.csv")
)

# -----------------------------
# 9. GLOBAL DROP LOG
# -----------------------------
global_drop_log <- bind_rows(
  var_drop_stage0 %>%
    mutate(
      mean_rf_importance = NA_real_,
      mean_abs_beta = NA_real_,
      overall_score = NA_real_
    ) %>%
    select(stage, variable, variable_label, mean_rf_importance, mean_abs_beta, overall_score, drop_reason),
  
  rf_candidate_drop_table %>%
    mutate(
      mean_abs_beta = NA_real_,
      overall_score = NA_real_
    ) %>%
    select(stage, variable, variable_label, mean_rf_importance, mean_abs_beta, overall_score, drop_reason),
  
  pre_prune_drop_table %>%
    select(stage, variable, variable_label, mean_rf_importance, mean_abs_beta, overall_score, drop_reason),
  
  correlation_pruning_pairs %>%
    transmute(
      stage,
      variable = dropped_variable,
      variable_label = dropped_variable_label,
      mean_rf_importance = NA_real_,
      mean_abs_beta = NA_real_,
      overall_score = NA_real_,
      drop_reason = paste0(
        reason, "; kept variable = ", kept_variable,
        " (", kept_variable_label, "), abs(r) = ",
        round(abs_correlation, 3)
      )
    ),
  
  manual_drop_table %>%
    mutate(
      mean_rf_importance = NA_real_,
      mean_abs_beta = NA_real_,
      overall_score = NA_real_
    ) %>%
    select(stage, variable, variable_label, mean_rf_importance, mean_abs_beta, overall_score, drop_reason)
) %>%
  arrange(stage, variable)

write_csv(
  global_drop_log,
  file.path(out_dir, "variable_drop_log_all_stages.csv")
)

# -----------------------------
# 10. CORRELATION EXPORTS FOR FINAL VARIABLES
# -----------------------------
if (length(final_selected_vars) >= 2) {
  corr_mat <- cor(
    hvi_model_matrix[, final_selected_vars],
    use = "pairwise.complete.obs"
  )
  
  corr_var_order <- selected_variables_final %>%
    filter(variable %in% final_selected_vars) %>%
    arrange(desc(overall_score)) %>%
    pull(variable)
  
  corr_var_order <- corr_var_order[corr_var_order %in% colnames(corr_mat)]
  
  corr_long <- as.data.frame(as.table(corr_mat)) %>%
    as_tibble() %>%
    rename(var1 = Var1, var2 = Var2, correlation = Freq) %>%
    mutate(
      var1 = factor(var1, levels = corr_var_order),
      var2 = factor(var2, levels = rev(corr_var_order)),
      var1_label = pretty_var_label(as.character(var1)),
      var2_label = pretty_var_label(as.character(var2))
    )
  
  write_csv(as.data.frame(corr_mat), file.path(out_dir, "selected_variable_correlation_matrix.csv"))
  write_csv(
    corr_long %>% mutate(var1 = as.character(var1), var2 = as.character(var2)),
    file.path(out_dir, "selected_variable_correlation_long.csv")
  )
}

# -----------------------------
# 11. FIGURES
# -----------------------------
endpoint_focus <- c(
  "deaths", "death_cvd", "death_respiratory",
  "ed_visits", "ed_cvd", "ed_injury", "ed_respiratory",
  "ems_calls", "ems_cvd", "ems_injury", "ems_mental", "ems_syncope"
)

domain_palette <- c(
  "ED" = "#1b9e77",
  "EMS" = "#d95f02",
  "Mortality" = "#7570b3"
)

rf_importance_by_endpoint <- rf_importance_by_endpoint %>%
  mutate(
    domain = factor(source, levels = c("ED", "EMS", "Mortality"))
  )

interaction_terms_by_endpoint <- interaction_terms_by_endpoint %>%
  mutate(
    domain = factor(source, levels = c("ED", "EMS", "Mortality"))
  )

top_rf_vars <- selection_summary %>%
  slice_max(order_by = mean_rf_importance, n = plot_top_n) %>%
  pull(variable)

top_interaction_vars <- selection_summary %>%
  filter(mean_abs_beta > 0) %>%
  slice_max(order_by = mean_abs_beta, n = plot_top_n) %>%
  pull(variable)

top_scatter_vars <- selection_summary %>%
  slice_max(order_by = overall_score, n = min(10, nrow(selection_summary))) %>%
  pull(variable)

# Figure 1: RF importance bar plot
fig_rf <- selection_summary %>%
  slice_max(order_by = mean_rf_importance, n = plot_top_n) %>%
  mutate(variable_label = fct_reorder(pretty_var_label(variable), mean_rf_importance)) %>%
  ggplot(aes(x = mean_rf_importance, y = variable_label)) +
  geom_col(width = 0.75) +
  labs(
    title = "Random forest screening of baseline vulnerability variables",
    subtitle = "Average permutation importance across heat-health endpoints",
    x = "Mean permutation importance",
    y = NULL
  ) +
  theme_pub(base_size = base_size)

ggsave(
  filename = file.path(out_dir, "fig_rf_importance_top15.png"),
  plot = fig_rf,
  width = 10,
  height = 7,
  dpi = 400
)

# Figure 2: Interaction importance bar plot
fig_interaction <- selection_summary %>%
  filter(mean_abs_beta > 0) %>%
  slice_max(order_by = mean_abs_beta, n = plot_top_n) %>%
  mutate(variable_label = fct_reorder(pretty_var_label(variable), mean_abs_beta)) %>%
  ggplot(aes(x = mean_abs_beta, y = variable_label)) +
  geom_col(width = 0.75) +
  labs(
    title = "Interaction-based vulnerability selection",
    subtitle = "Average absolute heat × vulnerability interaction coefficient across endpoints",
    x = "Mean absolute interaction coefficient",
    y = NULL
  ) +
  theme_pub(base_size = base_size)

ggsave(
  filename = file.path(out_dir, "fig_interaction_importance_top15.png"),
  plot = fig_interaction,
  width = 10,
  height = 7,
  dpi = 400
)

# Figure 3: Overall score lollipop
fig_combined_score <- selection_summary %>%
  slice_max(order_by = overall_score, n = plot_top_n) %>%
  mutate(variable_label = fct_reorder(pretty_var_label(variable), overall_score)) %>%
  ggplot(aes(x = overall_score, y = variable_label)) +
  geom_segment(
    aes(x = 0, xend = overall_score, y = variable_label, yend = variable_label),
    linewidth = 1.1,
    color = "grey75"
  ) +
  geom_point(size = 4) +
  labs(
    title = "Overall variable selection ranking",
    subtitle = "Composite score combining RF importance, interaction strength, sign consistency, and significance",
    x = "Overall selection score",
    y = NULL
  ) +
  theme_pub(base_size = base_size)

ggsave(
  filename = file.path(out_dir, "fig_overall_selection_score_lollipop.png"),
  plot = fig_combined_score,
  width = 10,
  height = 7,
  dpi = 400
)

# Figure 4: Heatmap for final selected variables
# Figure 4: Heatmap for final selected variables
if (length(final_selected_vars) >= 1) {
  endpoint_order <- hvi_endpoint_metadata$endpoint_key
  endpoint_order_labels <- pretty_endpoint_label(endpoint_order)
  
  heatmap_dat <- interaction_terms_by_endpoint %>%
    group_by(endpoint_key, endpoint_label, variable_label) %>%
    summarise(
      estimate = mean(estimate, na.rm = TRUE),
      p.value = mean(p.value, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    filter(variable_label %in% pretty_var_label(final_selected_vars)) %>%
    mutate(
      endpoint_label = factor(endpoint_label, levels = endpoint_order_labels),
      variable_label = factor(variable_label, levels = rev(pretty_var_label(final_selected_vars)))
    )
  
  # robust symmetric limits to prevent a few extreme coefficients from dominating
  heat_lim <- quantile(abs(heatmap_dat$estimate), probs = 0.90, na.rm = TRUE)
  if (!is.finite(heat_lim) || heat_lim <= 0) {
    heat_lim <- max(abs(heatmap_dat$estimate), na.rm = TRUE)
  }
  
  fig_heatmap <- ggplot(heatmap_dat, aes(x = endpoint_label, y = variable_label, fill = estimate)) +
    geom_tile(color = "white", linewidth = 0.25) +
    geom_point(
      data = ~ dplyr::filter(.x, p.value < 0.05),
      size = 2.0,
      color = "black"
    ) +
    scale_fill_gradient2(
      low = "#2c7bb6",
      mid = "white",
      high = "#d7191c",
      midpoint = 0,
      limits = c(-heat_lim, heat_lim),
      oob = scales::squish,
      name = "Interaction\nbeta"
    ) +
    guides(
      fill = guide_colorbar(
        barwidth = unit(12, "cm"),
        barheight = unit(0.8, "cm"),
        title.position = "top",
        title.hjust = 0.5
      )
    ) +
    labs(
      title = "Heat-vulnerability interaction patterns across endpoints",
      x = NULL,
      y = NULL
    ) +
    theme_pub(base_size = 11) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "bottom",
      legend.title = element_text(size = 12, face = "bold"),
      legend.text = element_text(size = 11)
    )
  
  ggsave(
    filename = file.path(out_dir, "fig_variable_selection_heatmap.png"),
    plot = fig_heatmap,
    width = 14,
    height = 8,
    dpi = 400
  )
}

# Figure 5: Final selected variable correlation heatmap
if (length(final_selected_vars) >= 2) {
  corr_mat <- cor(
    hvi_model_matrix[, final_selected_vars],
    use = "pairwise.complete.obs"
  )
  
  corr_var_order <- selected_variables_final %>%
    filter(variable %in% final_selected_vars) %>%
    arrange(desc(overall_score)) %>%
    pull(variable)
  
  corr_long <- as.data.frame(as.table(corr_mat)) %>%
    as_tibble() %>%
    rename(var1 = Var1, var2 = Var2, correlation = Freq) %>%
    mutate(
      var1 = factor(var1, levels = corr_var_order),
      var2 = factor(var2, levels = rev(corr_var_order)),
      var1_label = pretty_var_label(as.character(var1)),
      var2_label = pretty_var_label(as.character(var2))
    )
  
  fig_corr <- ggplot(corr_long, aes(x = var1_label, y = var2_label, fill = correlation)) +
    geom_tile(color = "white", linewidth = 0.25) +
    scale_fill_gradient2(
      low = "#2166AC",
      mid = "white",
      high = "#B2182B",
      midpoint = 0,
      limits = c(-1, 1),
      breaks = c(-1, -0.5, 0, 0.5, 1),
      labels = scales::label_number(accuracy = 0.1),
      name = "Correlation"
    ) +
    labs(
      title = "Correlation structure of final selected vulnerability variables",
      subtitle = "Used to prune redundant features before HVI model fitting",
      x = NULL,
      y = NULL
    ) +
    theme_pub(base_size = 11) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "bottom"
    )
  
  ggsave(
    filename = file.path(out_dir, "fig_selected_variable_correlation.png"),
    plot = fig_corr,
    width = 10,
    height = 8,
    dpi = 400
  )
}

# Figure 6: RF vs interaction scatter
fig_rf_vs_interaction <- selection_summary %>%
  mutate(
    variable_label = pretty_var_label(variable),
    label_flag = variable %in% top_scatter_vars
  ) %>%
  ggplot(aes(x = mean_rf_importance, y = mean_abs_beta)) +
  geom_point(aes(size = overall_score, alpha = label_flag)) +
  ggrepel::geom_text_repel(
    data = ~ dplyr::filter(.x, label_flag),
    aes(label = variable_label),
    size = 3.4,
    box.padding = 0.25,
    point.padding = 0.2,
    segment.color = "grey60",
    segment.size = 0.3,
    max.overlaps = Inf,
    min.segment.length = 0
  ) +
  scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = 0.45), guide = "none") +
  labs(
    title = "Concordance between random forest and interaction-based selection",
    subtitle = "Variables high on both axes are strong candidates for final HVI inclusion",
    x = "Mean random forest importance",
    y = "Mean absolute interaction coefficient",
    size = "Overall\nscore"
  ) +
  theme_pub(base_size = base_size)

ggsave(
  filename = file.path(out_dir, "fig_rf_vs_interaction_scatter.png"),
  plot = fig_rf_vs_interaction,
  width = 10,
  height = 8,
  dpi = 400
)

# optional domain-colored endpoint distributions
if (nrow(rf_importance_by_endpoint) > 0 && length(final_selected_vars) > 0) {
  fig_rf_domain <- rf_importance_by_endpoint %>%
    filter(variable %in% final_selected_vars) %>%
    mutate(variable_label = pretty_var_label(variable)) %>%
    ggplot(aes(x = importance, y = fct_reorder(variable_label, importance))) +
    geom_point(
      aes(color = domain),
      position = position_jitter(height = 0.15),
      alpha = 0.8,
      size = 2.5
    ) +
    stat_summary(
      fun = mean,
      geom = "point",
      shape = 18,
      size = 4,
      color = "black"
    ) +
    scale_color_manual(values = domain_palette) +
    labs(
      title = "Random forest importance across endpoints",
      subtitle = "Colored by domain; diamonds represent mean importance",
      x = "Permutation importance",
      y = NULL,
      color = "Endpoint domain"
    ) +
    theme_pub(base_size = 12)
  
  ggsave(
    file.path(out_dir, "fig_rf_importance_domain.png"),
    fig_rf_domain,
    width = 12,
    height = 7,
    dpi = 400
  )
}

if (nrow(interaction_terms_by_endpoint) > 0 && length(final_selected_vars) > 0) {
  fig_violin <- interaction_terms_by_endpoint %>%
    filter(variable %in% final_selected_vars) %>%
    mutate(variable_label = pretty_var_label(variable)) %>%
    ggplot(aes(x = estimate, y = fct_reorder(variable_label, estimate))) +
    geom_violin(
      fill = "grey85",
      color = NA,
      alpha = 0.8
    ) +
    geom_point(
      aes(color = domain, size = -log10(p.value)),
      alpha = 0.9
    ) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
    scale_color_manual(values = domain_palette) +
    labs(
      title = "Heat × vulnerability interaction effects across endpoints",
      subtitle = "Points colored by domain; size reflects statistical significance",
      x = "Interaction coefficient",
      y = NULL,
      color = "Endpoint domain",
      size = expression(-log[10](p))
    ) +
    theme_pub(base_size = 12) +
    facet_grid(~ domain, scales = "free_x", space = "free_x")
  
  ggsave(
    file.path(out_dir, "fig_interaction_violin_domain.png"),
    fig_violin,
    width = 12,
    height = 7,
    dpi = 400
  )
}

# -----------------------------
# 12. SAVE OBJECTS + SUMMARY
# -----------------------------
assign("rf_importance_by_endpoint", rf_importance_by_endpoint, envir = .GlobalEnv)
assign("rf_importance_summary", rf_importance_summary, envir = .GlobalEnv)
assign("interaction_terms_by_endpoint", interaction_terms_by_endpoint, envir = .GlobalEnv)
assign("interaction_summary", interaction_summary, envir = .GlobalEnv)
assign("selected_variables_final", selected_variables_final, envir = .GlobalEnv)
assign("final_selected_vulnerability_vars", final_selected_vars, envir = .GlobalEnv)
assign("global_drop_log", global_drop_log, envir = .GlobalEnv)

message("\nDone.")
message("Objects created in environment:")
message(" - rf_importance_by_endpoint")
message(" - rf_importance_summary")
message(" - interaction_terms_by_endpoint")
message(" - interaction_summary")
message(" - selected_variables_final")
message(" - final_selected_vulnerability_vars")
message(" - global_drop_log")

message("\nSaved variable-selection artifacts:")
message(" - ", selected_vars_csv)
message(" - ", selected_vars_rds)
message(" - ", selection_summary_rds)

message("\nFinal selected variables:")
print(final_selected_vars)

message("\nTop selected variables table:")
print(
  selected_variables_final %>%
    filter(selected_final) %>%
    select(
      variable, variable_label,
      mean_rf_importance, mean_abs_beta,
      n_endpoints_interaction, overall_score
    ) %>%
    arrange(desc(overall_score))
)

message("\nDrop counts by stage:")
print(
  global_drop_log %>%
    count(stage, name = "n_dropped")
)
