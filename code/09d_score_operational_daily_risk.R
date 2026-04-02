# ================================================================================================
# HVI 2.0 | 09d_score_operational_daily_risk.R
# Use observed daily heat-dose values to score endpoint-specific, family-specific,
# and overall operational heat-health risk by community-day.
# ================================================================================================

source("09_utils_hvi.R")

# -----------------------------
# CONFIG
# -----------------------------
project_dir <- getwd()
out_dir <- ensure_output_dir(project_dir, "09_model_outputs")

model_matrix_obj  <- "hvi_model_matrix"
endpoint_meta_obj <- "hvi_endpoint_metadata"
models_obj        <- "endpoint_models"
weights_obj       <- "endpoint_weights"

# Optional alert tiers for the overall 0-100 score
alert_breaks <- c(-Inf, 25, 50, 75, Inf)
alert_labels <- c("Low", "Moderate", "High", "Very High")

# -----------------------------
# LOAD OBJECTS
# -----------------------------
if (!exists(model_matrix_obj, envir = .GlobalEnv)) stop("Object not found: ", model_matrix_obj)
if (!exists(endpoint_meta_obj, envir = .GlobalEnv)) stop("Object not found: ", endpoint_meta_obj)
if (!exists(models_obj, envir = .GlobalEnv)) stop("Run 09a first: missing object ", models_obj)
if (!exists(weights_obj, envir = .GlobalEnv)) stop("Run 09b first: missing object ", weights_obj)

hvi_model_matrix <- get(model_matrix_obj, envir = .GlobalEnv) %>%
  clean_names() %>%
  mutate(
    date = as.Date(date),
    community = as.character(community),
    year = as.integer(year),
    doy = as.integer(doy),
    dow = factor(as.character(dow), levels = c("Sun", "Mon", "Tue", "Wed", "Thu", "Fri", "Sat"))
  )

hvi_endpoint_metadata <- get(endpoint_meta_obj, envir = .GlobalEnv) %>% clean_names()
endpoint_models <- get(models_obj, envir = .GlobalEnv)
endpoint_weights <- get(weights_obj, envir = .GlobalEnv)

# -----------------------------
# SCORE OBSERVED DAILY RISK
# -----------------------------
endpoint_daily_list <- list()

for (ep_key in names(endpoint_models)) {
  fit <- endpoint_models[[ep_key]]
  meta_row <- hvi_endpoint_metadata %>% filter(endpoint_key == ep_key)
  if (nrow(meta_row) == 0) next

  outcome_var <- meta_row$panel_outcome_col[1]
  heat_var <- get_heat_var(ep_key, hvi_model_matrix)
  if (length(heat_var) == 0 || is.na(heat_var)) next

  all_model_vars <- all.vars(formula(fit))
  z_vars_use <- all_model_vars[str_detect(all_model_vars, "^z_")]

  score_df <- hvi_model_matrix %>%
    select(community, date, year, doy, dow, pop_offset, any_of(outcome_var), all_of(heat_var), any_of(z_vars_use)) %>%
    rename(outcome = any_of(outcome_var), heat_dose = all_of(heat_var)) %>%
    drop_na(heat_dose, doy, dow, year)

  if (nrow(score_df) == 0) next

  ref_df <- score_df %>% mutate(heat_dose = 0)
  score_df$predicted_count <- as.numeric(predict(fit, newdata = score_df, type = "response"))
  ref_df$reference_count   <- as.numeric(predict(fit, newdata = ref_df, type = "response"))

  ep_daily <- score_df %>%
    mutate(
      endpoint_key = ep_key,
      outcome_label = meta_row$outcome_label[1],
      source = meta_row$source[1],
      domain = meta_row$domain[1],
      observed_count = outcome,
      reference_count = ref_df$reference_count,
      excess_events = predicted_count - reference_count,
      relative_risk = ifelse(reference_count > 0, predicted_count / reference_count, NA_real_)
    ) %>%
    select(community, date, year, endpoint_key, outcome_label, source, domain,
           observed_count, predicted_count, reference_count, excess_events, relative_risk)

  endpoint_daily_list[[ep_key]] <- ep_daily
}

ca_day_endpoint_risk <- bind_rows(endpoint_daily_list) %>%
  group_by(endpoint_key) %>%
  mutate(endpoint_risk_0_100 = rescale_0_100(excess_events)) %>%
  ungroup() %>%
  left_join(endpoint_weights %>% select(endpoint_key, endpoint_weight, source_weight, performance_weight), by = "endpoint_key")

ca_day_family_risk <- ca_day_endpoint_risk %>%
  group_by(community, date, year, source) %>%
  summarise(
    family_observed_count = sum(observed_count, na.rm = TRUE),
    family_predicted_count = sum(predicted_count, na.rm = TRUE),
    family_reference_count = sum(reference_count, na.rm = TRUE),
    family_excess_events = sum(excess_events, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(source) %>%
  mutate(family_risk_0_100 = rescale_0_100(family_excess_events)) %>%
  ungroup()

# dominant endpoint per day-community
ca_day_dominant_endpoint <- ca_day_endpoint_risk %>%
  group_by(community, date, year) %>%
  slice_max(order_by = excess_events, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  transmute(
    community, date, year,
    dominant_endpoint = endpoint_key,
    dominant_endpoint_label = outcome_label,
    dominant_endpoint_source = source,
    dominant_endpoint_excess = excess_events,
    dominant_endpoint_risk_0_100 = endpoint_risk_0_100
  )

ca_day_overall_operational_hvi <- ca_day_endpoint_risk %>%
  mutate(endpoint_weight = coalesce(endpoint_weight, 1)) %>%
  group_by(community, date, year) %>%
  summarise(
    total_observed_count = sum(observed_count, na.rm = TRUE),
    total_predicted_count = sum(predicted_count, na.rm = TRUE),
    total_reference_count = sum(reference_count, na.rm = TRUE),
    total_excess_events = sum(excess_events, na.rm = TRUE),
    overall_weighted_excess = weighted.mean(excess_events, w = endpoint_weight, na.rm = TRUE),
    n_endpoints = n(),
    .groups = "drop"
  ) %>%
  mutate(
    overall_risk_0_100 = rescale_0_100(overall_weighted_excess),
    alert_tier = cut(overall_risk_0_100, breaks = alert_breaks, labels = alert_labels, include.lowest = TRUE, right = TRUE)
  ) %>%
  left_join(ca_day_dominant_endpoint, by = c("community", "date", "year"))

# -----------------------------
# EXPORT
# -----------------------------
write_csv(ca_day_endpoint_risk, file.path(out_dir, "09d_ca_day_endpoint_risk.csv"))
write_csv(ca_day_family_risk, file.path(out_dir, "09d_ca_day_family_risk.csv"))
write_csv(ca_day_overall_operational_hvi, file.path(out_dir, "09d_ca_day_overall_operational_hvi.csv"))

assign("ca_day_endpoint_risk", ca_day_endpoint_risk, envir = .GlobalEnv)
assign("ca_day_family_risk", ca_day_family_risk, envir = .GlobalEnv)
assign("ca_day_overall_operational_hvi", ca_day_overall_operational_hvi, envir = .GlobalEnv)

message("09d complete. Outputs written to: ", out_dir)
