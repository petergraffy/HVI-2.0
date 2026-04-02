# ================================================================================================
# HVI 2.0 | 09c_score_temperature_grid.R
# Score endpoint-specific and total heat-health risk across a temperature grid.
#
# IMPORTANT:
# The fitted models use lagged heat_dose. This script defaults to a steady-state scenario:
# a temperature value is assumed to persist across the endpoint-specific lag window, so
# heat_dose = pmax(temp_f - MRT, 0) * (max_lag + 1).
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

# temperature grid in Fahrenheit or same native units as tmax in your matrix
# edit if your tmax is Celsius
# temp_grid <- seq(24, 42, by = 1)  # example Celsius grid
if (exists(model_matrix_obj, envir = .GlobalEnv)) {
  tmp <- get(model_matrix_obj, envir = .GlobalEnv) %>% clean_names()
  if ("tmax" %in% names(tmp) && max(tmp$tmax, na.rm = TRUE) > 60) {
    temp_grid <- seq(75, 105, by = 1)
  } else {
    temp_grid <- seq(24, 42, by = 1)
  }
} else {
  temp_grid <- seq(75, 105, by = 1)
}

template_doy <- 200
template_dow <- "Wed"
scenario_label <- "steady_state_lag_window"

# -----------------------------
# LOAD OBJECTS
# -----------------------------
if (!exists(model_matrix_obj, envir = .GlobalEnv)) stop("Object not found: ", model_matrix_obj)
if (!exists(endpoint_meta_obj, envir = .GlobalEnv)) stop("Object not found: ", endpoint_meta_obj)
if (!exists(models_obj, envir = .GlobalEnv)) stop("Run 09a first: missing object ", models_obj)
if (!exists(weights_obj, envir = .GlobalEnv)) stop("Run 09b first: missing object ", weights_obj)

hvi_model_matrix <- get(model_matrix_obj, envir = .GlobalEnv) %>% clean_names()
hvi_endpoint_metadata <- get(endpoint_meta_obj, envir = .GlobalEnv) %>% clean_names()
endpoint_models <- get(models_obj, envir = .GlobalEnv)
endpoint_weights <- get(weights_obj, envir = .GlobalEnv)

all_z_vars <- names(hvi_model_matrix)[str_detect(names(hvi_model_matrix), "^z_")]
community_year_base <- hvi_model_matrix %>%
  select(community, year, any_of(all_z_vars), pop_offset) %>%
  distinct() %>%
  mutate(doy = template_doy, dow = factor(template_dow, levels = c("Sun", "Mon", "Tue", "Wed", "Thu", "Fri", "Sat")))

# -----------------------------
# SCORE GRID
# -----------------------------
endpoint_grid_list <- list()

for (ep_key in names(endpoint_models)) {
  fit <- endpoint_models[[ep_key]]
  meta_row <- hvi_endpoint_metadata %>% filter(endpoint_key == ep_key)
  if (nrow(meta_row) == 0) next

  model_vars <- attr(stats::terms(fit), "term.labels")
  all_model_vars <- all.vars(formula(fit))
  z_vars_use <- all_model_vars[str_detect(all_model_vars, "^z_")]
  z_vars_use <- setdiff(z_vars_use, c("outcome", "heat_dose", "pop_offset", "doy", "year"))

  mrt_val <- meta_row$mrt[1]
  max_lag_val <- meta_row$max_lag[1]
  if (is.na(max_lag_val)) max_lag_val <- 0

  base_ep <- community_year_base %>%
    select(community, year, any_of(z_vars_use), pop_offset, doy, dow)

  grid_ep <- tidyr::expand_grid(
    base_ep,
    temp_value = temp_grid
  ) %>%
    mutate(
      endpoint_key = ep_key,
      outcome_label = meta_row$outcome_label[1],
      source = meta_row$source[1],
      domain = meta_row$domain[1],
      mrt = mrt_val,
      max_lag = max_lag_val,
      excess_above_mrt = pmax(temp_value - mrt_val, 0),
      heat_dose = excess_above_mrt * (max_lag_val + 1),
      scenario = scenario_label
    )

  ref_ep <- grid_ep %>% mutate(heat_dose = 0)

  grid_ep$predicted_count <- as.numeric(predict(fit, newdata = grid_ep, type = "response"))
  ref_ep$reference_count  <- as.numeric(predict(fit, newdata = ref_ep, type = "response"))

  grid_ep <- grid_ep %>%
    mutate(
      reference_count = ref_ep$reference_count,
      excess_events = predicted_count - reference_count,
      relative_risk = ifelse(reference_count > 0, predicted_count / reference_count, NA_real_)
    )

  endpoint_grid_list[[ep_key]] <- grid_ep
}

temp_grid_endpoint_risk <- bind_rows(endpoint_grid_list) %>%
  group_by(endpoint_key) %>%
  mutate(endpoint_risk_0_100 = rescale_0_100(excess_events)) %>%
  ungroup() %>%
  left_join(endpoint_weights %>% select(endpoint_key, endpoint_weight, source_weight, performance_weight), by = "endpoint_key")

# family totals
family_grid <- temp_grid_endpoint_risk %>%
  group_by(community, year, temp_value, source, scenario) %>%
  summarise(
    family_predicted_count = sum(predicted_count, na.rm = TRUE),
    family_reference_count = sum(reference_count, na.rm = TRUE),
    family_excess_events = sum(excess_events, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(source) %>%
  mutate(family_risk_0_100 = rescale_0_100(family_excess_events)) %>%
  ungroup()

# overall operational total
overall_grid <- temp_grid_endpoint_risk %>%
  mutate(endpoint_weight = coalesce(endpoint_weight, 1)) %>%
  group_by(community, year, temp_value, scenario) %>%
  summarise(
    total_predicted_count = sum(predicted_count, na.rm = TRUE),
    total_reference_count = sum(reference_count, na.rm = TRUE),
    total_excess_events = sum(excess_events, na.rm = TRUE),
    overall_weighted_excess = weighted.mean(excess_events, w = endpoint_weight, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(overall_risk_0_100 = rescale_0_100(overall_weighted_excess))

# dominant endpoint at each community-year-temperature
endpoint_dominant <- temp_grid_endpoint_risk %>%
  group_by(community, year, temp_value, scenario) %>%
  slice_max(order_by = excess_events, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  transmute(
    community, year, temp_value, scenario,
    dominant_endpoint = endpoint_key,
    dominant_endpoint_excess = excess_events,
    dominant_endpoint_risk_0_100 = endpoint_risk_0_100
  )

overall_grid <- overall_grid %>% left_join(endpoint_dominant, by = c("community", "year", "temp_value", "scenario"))

# convenience table for the exact 92-degree style query when that value exists
query_temp <- if (any(abs(temp_grid - 92) < 1e-8)) 92 else temp_grid[which.min(abs(temp_grid - 92))]
query_temperature_outputs <- list(
  endpoint = temp_grid_endpoint_risk %>% filter(temp_value == query_temp),
  overall  = overall_grid %>% filter(temp_value == query_temp),
  family   = family_grid %>% filter(temp_value == query_temp)
)

# -----------------------------
# EXPORT
# -----------------------------
write_csv(temp_grid_endpoint_risk, file.path(out_dir, "09c_temperature_grid_endpoint_risk.csv"))
write_csv(family_grid, file.path(out_dir, "09c_temperature_grid_family_risk.csv"))
write_csv(overall_grid, file.path(out_dir, "09c_temperature_grid_overall_risk.csv"))
write_csv(query_temperature_outputs$endpoint, file.path(out_dir, paste0("09c_endpoint_risk_at_", query_temp, ".csv")))
write_csv(query_temperature_outputs$family, file.path(out_dir, paste0("09c_family_risk_at_", query_temp, ".csv")))
write_csv(query_temperature_outputs$overall, file.path(out_dir, paste0("09c_overall_risk_at_", query_temp, ".csv")))

assign("temp_grid_endpoint_risk", temp_grid_endpoint_risk, envir = .GlobalEnv)
assign("temp_grid_family_risk", family_grid, envir = .GlobalEnv)
assign("temp_grid_overall_risk", overall_grid, envir = .GlobalEnv)
assign("query_temperature_outputs", query_temperature_outputs, envir = .GlobalEnv)

message("09c complete. Outputs written to: ", out_dir)
