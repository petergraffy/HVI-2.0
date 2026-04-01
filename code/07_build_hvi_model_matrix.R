# ================================================================================================
# HVI 2.0 | 06_build_hvi_model_matrix.R
# Build CA x day model matrix for 2019-2022 overlap period
#
# Uses objects already in memory:
#   - panel_overlap (daily CA x day panel with outcomes + climate)
#   - mrt_table_expanded (endpoint MRT summary table)
#
# Uses baseline vulnerability either:
#   - from an object in memory, OR
#   - from disk
#
# Creates:
#   - hvi_endpoint_metadata   (endpoint MRT + lag + panel-column mapping)
#   - hvi_model_matrix        (CA x day model matrix)
#
# Optional CSV exports:
#   - hvi_endpoint_metadata.csv
#   - hvi_model_matrix_2019_2022.csv
#
# Notes:
# - mrt_table_expanded contains MRT but not lag.
# - If a true lag table/object exists in memory, set true_lag_obj below.
# - Otherwise fallback lags are used:
#       no_lag_forced = 0
#       EMS DLNM      = 1
#       Mortality/ED  = 2
# ================================================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(lubridate)
  library(janitor)
  library(slider)
  library(stringr)
  library(readr)
})

# -----------------------------
# CONFIG
# -----------------------------
project_dir <- "C:/Users/Peter Graffy/Box/HVI2.0"

# in-memory objects
daily_panel_obj <- "panel_overlap"
mrt_table_obj   <- "mrt_table_expanded"

# baseline vulnerability: prefer object if it exists, otherwise read from file
baseline_obj  <- "baseline_vulnerability_ca_year"
baseline_path <- file.path(project_dir, "Baseline Burden", "baseline_vulnerability_ca_year.csv")

# optional true lag lookup object in memory
# set to NULL if you do not have one
# expected columns somewhere in it should identify endpoint + max_lag
true_lag_obj <- NULL
# examples later:
# true_lag_obj <- "endpoint_lag_lookup"
# true_lag_obj <- "dlnm_lag_table"

# write outputs to disk?
write_outputs <- TRUE

out_model_matrix  <- file.path(project_dir, "hvi_model_matrix_2019_2022.csv")
out_endpoint_meta <- file.path(project_dir, "hvi_endpoint_metadata.csv")

# temperature / humidity
temp_var_candidates     <- c("tmax", "temp_max", "daily_max_temp")
humidity_var_candidates <- c("humidity", "rh", "relative_humidity")

# column candidates in daily panel
community_col_candidates <- c("community", "community_area", "community_name")
date_col_candidates      <- c("date", "event_date")

# optional population columns
panel_pop_col_candidates    <- c("total_pop", "population", "pop")
baseline_pop_col_candidates <- c("total_pop", "population", "pop")

# recommended endpoint set for first production HVI
# left side = outcome key in mrt_table_expanded$outcome
# panel_candidates = possible matching column names in panel_overlap
endpoint_spec <- tribble(
  ~endpoint_key,      ~panel_candidates,                                 ~include_in_hvi,
  "deaths",           list(c("death_allcause", "deaths")),               TRUE,
  "death_cvd",        list(c("death_cvd")),                              TRUE,
  "death_injury",     list(c("death_injury")),                           TRUE,
  "death_mental",     list(c("death_mental")),                           TRUE,
  "death_renal",      list(c("death_renal")),                            TRUE,
  "death_respiratory",list(c("death_resp", "death_respiratory")),        TRUE,
  
  "ed_visits",        list(c("ed_allcause", "ed_visits")),               TRUE,
  "ed_cvd",           list(c("ed_cvd")),                                 TRUE,
  "ed_injury",        list(c("ed_injury")),                              TRUE,
  "ed_renal",         list(c("ed_renal")),                               TRUE,
  "ed_respiratory",   list(c("ed_resp", "ed_respiratory")),              TRUE,
  "ed_syncope",       list(c("ed_syncope")),                             TRUE,
  "ed_dehydration",   list(c("ed_dehydration")),                         TRUE,
  
  "ems_calls",        list(c("ems_allcause", "ems_calls")),              TRUE,
  "ems_cvd",          list(c("ems_cvd")),                                TRUE,
  "ems_injury",       list(c("ems_injury")),                             TRUE,
  "ems_mental",       list(c("ems_mental")),                             TRUE,
  "ems_neuro",        list(c("ems_neuro", "ems_neurologic")),            TRUE,
  "ems_respiratory",  list(c("ems_resp", "ems_respiratory")),            TRUE,
  "ems_syncope",      list(c("ems_syncope")),                            TRUE,
  "ems_bleeding",     list(c("ems_bleeding")),                           TRUE,
  "ems_gi",           list(c("ems_gi")),                                 TRUE
)

# vulnerability variables to carry into the model matrix
candidate_vuln_vars <- c(
  # ACS / demographics
  "total_pop", "mean_age", "median_age", "mean_black", "mean_hisp", "mean_white",
  "mean_asian", "mean_income", "median_income", "mean_unemployed", "mean_employed",
  "mean_college", "mean_hs", "mean_male", "mean_female",
  
  # NDVI / air pollution / AC
  "ndvi", "mean_ndvi", "ac_prob", "ac_cbsa_rank", "no2", "pm25",
  
  # SVI
  "svi_rpl_theme1", "svi_rpl_theme2", "svi_rpl_theme3", "svi_rpl_theme4", "svi_rpl_themes",
  "svi_ep_pov", "svi_ep_unemp", "svi_ep_nohsdp", "svi_ep_age65", "svi_ep_age17",
  "svi_ep_disabl", "svi_ep_sngpnt", "svi_ep_limeng", "svi_ep_minrty",
  "svi_ep_munit", "svi_ep_mobile", "svi_ep_crowd", "svi_ep_noveh", "svi_ep_groupq",
  "svi_ep_uninsur", "svi_ep_noint", "svi_ep_hburd"
)

center_scale_within_overlap <- TRUE

# -----------------------------
# Helper functions
# -----------------------------
pick_first_existing <- function(df, candidates) {
  hit <- candidates[candidates %in% names(df)][1]
  if (length(hit) == 0 || is.na(hit)) return(NA_character_)
  hit
}

pick_first_existing_value <- function(candidates, available_names) {
  hit <- candidates[candidates %in% available_names][1]
  if (length(hit) == 0 || is.na(hit)) return(NA_character_)
  hit
}

clean_community <- function(x) {
  x %>%
    as.character() %>%
    str_to_upper() %>%
    str_replace_all("&", "AND") %>%
    str_replace_all("[[:punct:]]", " ") %>%
    str_squish()
}

safe_date <- function(x) {
  if (inherits(x, "Date")) return(x)
  as.Date(x)
}

zscore <- function(x) {
  s <- sd(x, na.rm = TRUE)
  m <- mean(x, na.rm = TRUE)
  if (is.na(s) || s == 0) return(rep(0, length(x)))
  as.numeric((x - m) / s)
}

make_lagged_heat_dose <- function(x, max_lag) {
  slider::slide_dbl(
    x,
    .f = ~ sum(.x, na.rm = TRUE),
    .before = max_lag,
    .complete = FALSE
  )
}

resolve_panel_outcome_col <- function(cands, available_names) {
  hit <- cands[cands %in% available_names][1]
  if (length(hit) == 0 || is.na(hit)) return(NA_character_)
  hit
}

# -----------------------------
# 1. Load core objects
# -----------------------------
if (!exists(daily_panel_obj, envir = .GlobalEnv)) {
  stop("Object not found in environment: ", daily_panel_obj)
}
if (!exists(mrt_table_obj, envir = .GlobalEnv)) {
  stop("Object not found in environment: ", mrt_table_obj)
}

daily_panel <- get(daily_panel_obj, envir = .GlobalEnv) %>%
  clean_names()

mrt_table_expanded <- get(mrt_table_obj, envir = .GlobalEnv) %>%
  clean_names()

if (exists(baseline_obj, envir = .GlobalEnv)) {
  baseline <- get(baseline_obj, envir = .GlobalEnv) %>%
    clean_names()
} else {
  baseline <- read_csv(baseline_path, show_col_types = FALSE) %>%
    clean_names()
}

# -----------------------------
# 2. Basic harmonization
# -----------------------------
community_col <- pick_first_existing(daily_panel, community_col_candidates)
date_col      <- pick_first_existing(daily_panel, date_col_candidates)
temp_var      <- pick_first_existing(daily_panel, temp_var_candidates)
humidity_var  <- pick_first_existing(daily_panel, humidity_var_candidates)

if (is.na(community_col)) stop("Could not find community column in daily panel.")
if (is.na(date_col))      stop("Could not find date/event_date column in daily panel.")
if (is.na(temp_var))      stop("Could not find tmax column in daily panel.")

daily_panel <- daily_panel %>%
  mutate(
    community = clean_community(.data[[community_col]]),
    date      = safe_date(.data[[date_col]]),
    year      = year(date),
    month     = month(date),
    doy       = yday(date),
    dow       = wday(date)
  ) %>%
  filter(year >= 2019, year <= 2022)

baseline <- baseline %>%
  mutate(
    community = clean_community(community),
    year = as.integer(year)
  ) %>%
  filter(year >= 2019, year <= 2022)

# -----------------------------
# 3. Vulnerability variables
# -----------------------------
vuln_vars <- candidate_vuln_vars[candidate_vuln_vars %in% names(baseline)]

if (length(vuln_vars) == 0) {
  stop("None of the candidate vulnerability variables were found in baseline.")
}

message("Using vulnerability vars:")
print(vuln_vars)

baseline_pop_col <- pick_first_existing(baseline, baseline_pop_col_candidates)
panel_pop_col    <- pick_first_existing(daily_panel, panel_pop_col_candidates)

baseline_model <- baseline %>%
  select(community, year, any_of(c(vuln_vars, baseline_pop_col))) %>%
  distinct()

if (center_scale_within_overlap) {
  baseline_model <- baseline_model %>%
    mutate(across(all_of(vuln_vars), zscore, .names = "z_{.col}"))
}

z_vuln_vars <- names(baseline_model)[str_detect(names(baseline_model), "^z_")]

baseline_model <- baseline_model %>%
  mutate(
    baseline_vuln_mean_z = rowMeans(select(., all_of(z_vuln_vars)), na.rm = TRUE)
  )

# -----------------------------
# 4. Join baseline to daily panel
# -----------------------------
model_dat <- daily_panel %>%
  left_join(baseline_model, by = c("community", "year"))

if (!is.na(panel_pop_col) && panel_pop_col %in% names(model_dat)) {
  model_dat <- model_dat %>%
    mutate(pop_offset = .data[[panel_pop_col]])
} else if (!is.na(baseline_pop_col) && baseline_pop_col %in% names(model_dat)) {
  model_dat <- model_dat %>%
    mutate(pop_offset = .data[[baseline_pop_col]])
} else {
  model_dat <- model_dat %>%
    mutate(pop_offset = NA_real_)
}

# -----------------------------
# 5. Build endpoint metadata
# -----------------------------
# 5a. Resolve panel outcome columns
endpoint_meta_base <- endpoint_spec %>%
  rowwise() %>%
  mutate(
    panel_outcome_col = resolve_panel_outcome_col(panel_candidates[[1]], names(model_dat))
  ) %>%
  ungroup()

# 5b. Build base metadata from mrt_table_expanded
# mrt_table_expanded fields seen:
# outcome, outcome_label, mrt, n_days, mean_daily_count, model_type, endpoint, domain
endpoint_meta <- mrt_table_expanded %>%
  transmute(
    endpoint_key      = as.character(outcome),
    outcome_label     = outcome_label,
    mrt               = as.numeric(mrt),
    n_days            = as.integer(n_days),
    mean_daily_count  = as.numeric(mean_daily_count),
    model_type        = as.character(model_type),
    source            = as.character(endpoint),
    domain            = as.character(domain)
  ) %>%
  left_join(
    endpoint_meta_base %>%
      select(endpoint_key, panel_outcome_col, include_in_hvi),
    by = "endpoint_key"
  ) %>%
  filter(include_in_hvi %in% TRUE) %>%
  filter(!is.na(panel_outcome_col)) %>%
  distinct(endpoint_key, .keep_all = TRUE)

if (nrow(endpoint_meta) == 0) {
  stop("No endpoints matched between mrt_table_expanded and the daily panel.")
}

# 5c. Attach true lags if available, otherwise use fallback lags
if (!is.null(true_lag_obj) && exists(true_lag_obj, envir = .GlobalEnv)) {
  lag_df <- get(true_lag_obj, envir = .GlobalEnv) %>%
    clean_names()
  
  lag_endpoint_col <- pick_first_existing(lag_df, c("endpoint_key", "outcome", "endpoint"))
  lag_value_col    <- pick_first_existing(lag_df, c("max_lag", "lag", "lag_days", "best_lag"))
  
  if (is.na(lag_endpoint_col) || is.na(lag_value_col)) {
    stop("true_lag_obj exists but no endpoint / lag columns could be identified.")
  }
  
  lag_lookup <- lag_df %>%
    transmute(
      endpoint_key = as.character(.data[[lag_endpoint_col]]),
      max_lag = as.integer(round(as.numeric(.data[[lag_value_col]])))
    ) %>%
    distinct(endpoint_key, .keep_all = TRUE)
  
  endpoint_meta <- endpoint_meta %>%
    left_join(lag_lookup, by = "endpoint_key")
} else {
  endpoint_meta <- endpoint_meta %>%
    mutate(max_lag = NA_integer_)
}

endpoint_meta <- endpoint_meta %>%
  mutate(
    max_lag = case_when(
      !is.na(max_lag)                  ~ max_lag,
      model_type == "no_lag_forced"    ~ 0L,
      str_to_upper(source) == "EMS"    ~ 1L,
      str_to_upper(source) %in% c("ED", "MORTALITY") ~ 2L,
      TRUE                             ~ 0L
    )
  )

message("Endpoint metadata:")
print(endpoint_meta %>% select(endpoint_key, panel_outcome_col, mrt, max_lag, model_type, source, domain))

# -----------------------------
# 6. Build endpoint-specific heat features
# -----------------------------
model_dat <- model_dat %>%
  arrange(community, date)

for (i in seq_len(nrow(endpoint_meta))) {
  ek   <- endpoint_meta$endpoint_key[i]
  mrt  <- endpoint_meta$mrt[i]
  lagk <- endpoint_meta$max_lag[i]
  
  excess_name   <- paste0("temp_excess__", ek)
  heatdose_name <- paste0("heat_dose__", ek)
  mrt_name      <- paste0("mrt__", ek)
  lag_name      <- paste0("max_lag__", ek)
  
  model_dat <- model_dat %>%
    group_by(community) %>%
    arrange(date, .by_group = TRUE) %>%
    mutate(
      !!mrt_name      := mrt,
      !!lag_name      := lagk,
      !!excess_name   := pmax(.data[[temp_var]] - mrt, 0),
      !!heatdose_name := make_lagged_heat_dose(.data[[excess_name]], max_lag = lagk)
    ) %>%
    ungroup()
}

# -----------------------------
# 7. Keep final columns
# -----------------------------
outcome_cols_present <- endpoint_meta$panel_outcome_col[endpoint_meta$panel_outcome_col %in% names(model_dat)]

heat_cols_present <- names(model_dat)[
  str_detect(names(model_dat), "^temp_excess__|^heat_dose__|^mrt__|^max_lag__")
]

keep_cols <- unique(c(
  "community", "date", "year", "month", "doy", "dow",
  temp_var, humidity_var,
  outcome_cols_present,
  "pop_offset",
  vuln_vars,
  z_vuln_vars,
  "baseline_vuln_mean_z",
  heat_cols_present
))

keep_cols <- keep_cols[keep_cols %in% names(model_dat)]

hvi_model_matrix <- model_dat %>%
  select(all_of(keep_cols)) %>%
  arrange(community, date)

hvi_endpoint_metadata <- endpoint_meta %>%
  arrange(source, domain, endpoint_key)

# -----------------------------
# 8. Put objects in environment
# -----------------------------
assign("hvi_model_matrix", hvi_model_matrix, envir = .GlobalEnv)
assign("hvi_endpoint_metadata", hvi_endpoint_metadata, envir = .GlobalEnv)

# -----------------------------
# 9. Optional write to disk
# -----------------------------
if (isTRUE(write_outputs)) {
  write_csv(hvi_model_matrix, out_model_matrix)
  write_csv(hvi_endpoint_metadata, out_endpoint_meta)
}

# -----------------------------
# 10. QA
# -----------------------------
message("\nDone.")
message("Created objects in environment:")
message(" - hvi_model_matrix")
message(" - hvi_endpoint_metadata")

message("\nRows: ", nrow(hvi_model_matrix))
message("Communities: ", n_distinct(hvi_model_matrix$community))
message("Date range: ", min(hvi_model_matrix$date, na.rm = TRUE), " to ", max(hvi_model_matrix$date, na.rm = TRUE))

message("\nEndpoints included:")
print(
  hvi_endpoint_metadata %>%
    select(endpoint_key, panel_outcome_col, outcome_label, mrt, max_lag, source, domain, model_type)
)

message("\nMissingness check for mapped outcome columns:")
print(
  hvi_endpoint_metadata %>%
    mutate(
      n_missing = map_int(panel_outcome_col, ~ sum(is.na(hvi_model_matrix[[.x]]))),
      n_nonzero = map_int(panel_outcome_col, ~ sum(hvi_model_matrix[[.x]] > 0, na.rm = TRUE))
    ) %>%
    select(endpoint_key, panel_outcome_col, n_missing, n_nonzero)
)
