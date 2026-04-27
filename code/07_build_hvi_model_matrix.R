# ================================================================================================
# HVI 2.0 | 07_build_hvi_model_matrix.R
# Build CA x day model matrix for 2019-2022 overlap period
#
# Required objects in memory:
#   - panel_overlap        : daily CA x day panel with outcomes + climate
#   - mrt_table_expanded   : endpoint MRT summary table from 05_dlnm_mrt.R
#
# Baseline vulnerability:
#   - uses object in memory if available, otherwise reads from disk
#
# Creates:
#   - hvi_endpoint_metadata
#   - hvi_model_matrix
#
# Optional CSV exports:
#   - hvi_endpoint_metadata.csv
#   - hvi_model_matrix_2019_2022.csv
#
# Notes:
# - This version uses tmax for all downstream exposure features.
# - This version expects mrt_table_expanded to already contain:
#     outcome, outcome_label, mrt, n_days, mean_daily_count, model_type,
#     endpoint, domain, max_lag, lag_input, study_year_start, study_year_end,
#     study_window_label
# ================================================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(lubridate)
  library(janitor)
  library(slider)
  library(stringr)
  library(readr)
})

source(file.path(if (dir.exists("code")) "code" else ".", "00_config.R"))

# -----------------------------
# CONFIG
# -----------------------------
project_dir <- HVI_PATHS$private

analysis_year_start <- 2019L
analysis_year_end   <- 2022L

daily_panel_obj <- "panel_overlap"
mrt_table_obj   <- "mrt_table_expanded"

baseline_obj  <- "baseline_vulnerability_ca_year"
baseline_path <- file.path(HVI_PATHS$private_outputs$baseline_burden, "baseline_vulnerability_ca_year.csv")

write_outputs <- TRUE
out_model_matrix <- file.path(project_dir, "hvi_model_matrix_2019_2022.csv")
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
  ~endpoint_key,       ~panel_candidates,                                  ~include_in_hvi,
  
  # -----------------------------
  # MORTALITY
  # -----------------------------
  "deaths",            list(c("death_allcause", "deaths")),                TRUE,
  "death_cvd",         list(c("death_cvd")),                               TRUE,
  "death_respiratory", list(c("death_resp", "death_respiratory")),         TRUE,
  "death_renal",       list(c("death_renal")),                             TRUE,
  "death_neurologic",  list(c("death_neuro", "death_neurologic")),         TRUE,
  "death_mental",      list(c("death_mental")),                            TRUE,
  "death_gi",          list(c("death_gi")),                                TRUE,
  "death_injury",      list(c("death_injury")),                            TRUE,
  "death_syncope",     list(c("death_syncope")),                           TRUE,
  "death_dehydration", list(c("death_dehydration")),                       TRUE,
  "death_heat",        list(c("death_heat")),                              TRUE,
  
  # -----------------------------
  # ED
  # -----------------------------
  "ed_visits",         list(c("ed_allcause", "ed_visits")),                TRUE,
  "ed_cvd",            list(c("ed_cvd")),                                  TRUE,
  "ed_respiratory",    list(c("ed_resp", "ed_respiratory")),               TRUE,
  "ed_renal",          list(c("ed_renal")),                                TRUE,
  "ed_neurologic",     list(c("ed_neuro", "ed_neurologic")),               TRUE,
  "ed_mental",         list(c("ed_mental")),                               TRUE,
  "ed_gi",             list(c("ed_gi")),                                   TRUE,
  "ed_injury",         list(c("ed_injury")),                               TRUE,
  "ed_syncope",        list(c("ed_syncope")),                              TRUE,
  "ed_dehydration",    list(c("ed_dehydration")),                          TRUE,
  "ed_heat",           list(c("ed_heat")),                                 TRUE,
  
  # -----------------------------
  # EMS
  # -----------------------------
  "ems_calls",         list(c("ems_allcause", "ems_calls")),               TRUE,
  "ems_cvd",           list(c("ems_cvd")),                                 TRUE,
  "ems_respiratory",   list(c("ems_resp", "ems_respiratory")),             TRUE,
  "ems_neuro",         list(c("ems_neuro", "ems_neurologic")),             TRUE,
  "ems_mental",        list(c("ems_mental")),                              TRUE,
  "ems_gi",            list(c("ems_gi")),                                  TRUE,
  "ems_bleeding",      list(c("ems_bleeding")),                            TRUE,
  "ems_injury",        list(c("ems_injury")),                              TRUE,
  "ems_syncope",       list(c("ems_syncope")),                             TRUE,
  "ems_heat",          list(c("ems_heat")),                                TRUE
)

# vulnerability variables to carry into the model matrix
candidate_vuln_vars <- c(
  "pop_density_km2", "mean_age", "median_age", "mean_black", "mean_hisp", "mean_white",
  "mean_asian", "mean_income", "median_income", "mean_unemployed", "mean_employed",
  "mean_college", "mean_hs", "mean_male", "mean_female",
  "ndvi", "mean_ndvi", "ac_prob", "ac_cbsa_rank", "no2", "pm25",
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
    dow       = factor(wday(date),
                       levels = levels(wday(
                         seq.Date(as.Date("2023-01-01"), by = "day", length.out = 7)
                       )))
  ) %>%
  filter(year >= analysis_year_start, year <= analysis_year_end)

baseline <- baseline %>%
  mutate(
    community = clean_community(community),
    year = as.integer(year)
  ) %>%
  filter(year >= analysis_year_start, year <= analysis_year_end)

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
endpoint_meta_base <- endpoint_spec %>%
  rowwise() %>%
  mutate(
    panel_outcome_col = resolve_panel_outcome_col(panel_candidates[[1]], names(model_dat))
  ) %>%
  ungroup()

required_mrt_cols <- c(
  "outcome", "outcome_label", "mrt", "n_days", "mean_daily_count",
  "model_type", "endpoint", "domain", "max_lag", "lag_input",
  "study_year_start", "study_year_end", "study_window_label"
)

missing_mrt_cols <- setdiff(required_mrt_cols, names(mrt_table_expanded))
if (length(missing_mrt_cols) > 0) {
  stop(
    "mrt_table_expanded is missing required columns: ",
    paste(missing_mrt_cols, collapse = ", "),
    ". Rebuild MRTs from 05_dlnm_mrt.R before running this script."
  )
}

endpoint_meta <- mrt_table_expanded %>%
  transmute(
    endpoint_key       = as.character(outcome),
    outcome_label      = as.character(outcome_label),
    mrt                = as.numeric(mrt),
    n_days             = as.integer(n_days),
    mean_daily_count   = as.numeric(mean_daily_count),
    model_type         = as.character(model_type),
    source             = as.character(endpoint),
    domain             = as.character(domain),
    max_lag            = as.integer(max_lag),
    lag_input          = as.integer(lag_input),
    study_year_start   = as.integer(study_year_start),
    study_year_end     = as.integer(study_year_end),
    study_window_label = as.character(study_window_label)
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

if (any(is.na(endpoint_meta$mrt))) {
  warning("Some endpoints have missing MRT values and will still be carried forward.")
}

if (any(is.na(endpoint_meta$max_lag))) {
  stop("Some endpoints have missing max_lag in mrt_table_expanded.")
}

message("Endpoint metadata:")
print(
  endpoint_meta %>%
    select(endpoint_key, panel_outcome_col, outcome_label, mrt, max_lag, model_type, source, domain)
)

# -----------------------------
# 6. Build endpoint-specific heat features from tmax
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
outcome_cols_present <- endpoint_meta$panel_outcome_col[
  endpoint_meta$panel_outcome_col %in% names(model_dat)
]

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
    select(endpoint_key, panel_outcome_col, outcome_label, mrt, max_lag, source, domain, model_type, study_window_label)
)

message("\nMissingness check for mapped outcome columns:")
print(
  hvi_endpoint_metadata %>%
    mutate(
      n_missing = purrr::map_int(panel_outcome_col, ~ sum(is.na(hvi_model_matrix[[.x]]))),
      n_nonzero = purrr::map_int(panel_outcome_col, ~ sum(hvi_model_matrix[[.x]] > 0, na.rm = TRUE))
    ) %>%
    select(endpoint_key, panel_outcome_col, n_missing, n_nonzero)
)
