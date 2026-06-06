# ================================================================================================
# HVI 2.0 | Long-term annual forecast frontend exports
# Builds annual community-endpoint model matrices and 5/10/20-year forecast artifacts.
#
# This is an annual planning layer, not a replacement for the daily operational HVI. It estimates
# how annual endpoint burden changes under a linear annual-temperature trend using aggregated
# community-year temperature and outcome histories.
# ================================================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(janitor)
  library(lubridate)
  library(readr)
  library(stringr)
})

source(file.path(if (dir.exists("code")) "code" else ".", "00_config.R"))
hvi_source("09_utils_hvi.R")

# -----------------------------
# Config
# -----------------------------
project_dir <- HVI_PATHS$private
out_dir <- ensure_output_dir(project_dir, "09_model_outputs")
frontend_dir <- file.path(out_dir, "frontend_exports")
hvi_dir_create(frontend_dir)

horizons_years <- as.integer(strsplit(hvi_env("HVI_FORECAST_HORIZONS_YEARS", "5,10,20"), ",")[[1]])
confidence_level <- suppressWarnings(as.numeric(hvi_env("HVI_FORECAST_CONFIDENCE_LEVEL", "0.95")))
if (!is.finite(confidence_level) || confidence_level <= 0 || confidence_level >= 1) confidence_level <- 0.95

model_matrix_candidates <- c(
  file.path(HVI_PATHS$private, "hvi_model_matrix_2019_2022.csv"),
  file.path(HVI_PATHS$repo, "hvi_model_matrix_2019_2022.csv")
)
endpoint_meta_candidates <- c(
  file.path(HVI_PATHS$private, "hvi_endpoint_metadata.csv"),
  file.path(HVI_PATHS$private_outputs$model_outputs, "hvi_endpoint_metadata.csv")
)
endpoint_weight_candidates <- c(
  file.path(HVI_PATHS$private_outputs$model_outputs, "09b_endpoint_weights.csv"),
  file.path(HVI_PATHS$repo, "code", "09_model_outputs", "09b_endpoint_weights.csv")
)
climate_candidates <- c(
  HVI_PATHS$raw$climate_full,
  HVI_PATHS$raw$climate_hist,
  file.path(HVI_PATHS$private, "Climate", "CA_temps_rh_90-24.csv"),
  file.path(HVI_PATHS$private, "Climate", "CA_temps_rh_90-23.csv")
)

pick_existing <- function(paths, label) {
  hit <- paths[file.exists(paths)][1]
  if (length(hit) == 0 || is.na(hit)) stop("Missing required ", label, ". Checked: ", paste(paths, collapse = "; "))
  hit
}

clean_community <- function(x) {
  x %>%
    as.character() %>%
    str_to_upper() %>%
    str_replace_all("&", "AND") %>%
    str_replace_all("[[:punct:]]", " ") %>%
    str_squish() %>%
    str_replace("^O HARE$", "OHARE")
}

c_to_f <- function(x) (x * 9 / 5) + 32

safe_rate <- function(count, pop) {
  if_else(!is.na(pop) & pop > 0, 100000 * count / pop, NA_real_)
}

fit_endpoint_annual_model <- function(dat) {
  dat <- dat %>%
    filter(
      is.finite(observed_count),
      is.finite(annual_mean_temp_c),
      is.finite(population),
      population > 0
    ) %>%
    mutate(community = factor(community))

  if (nrow(dat) < 20 || n_distinct(dat$community) < 2 || sum(dat$observed_count, na.rm = TRUE) <= 0) {
    return(list(fit = NULL, model_type = "mean_rate_fallback", beta_temp_c = 0, status = "insufficient_annual_history"))
  }

  fit <- tryCatch(
    glm(
      observed_count ~ annual_mean_temp_c + community + offset(log(population)),
      family = quasipoisson(link = "log"),
      data = dat
    ),
    error = function(e) NULL
  )

  if (is.null(fit)) {
    return(list(fit = NULL, model_type = "mean_rate_fallback", beta_temp_c = 0, status = "glm_failed"))
  }

  beta <- unname(coef(fit)[["annual_mean_temp_c"]])
  if (!is.finite(beta)) beta <- 0

  list(fit = fit, model_type = "quasipoisson_community_fe", beta_temp_c = beta, status = "ok")
}

predict_endpoint_annual <- function(model, baseline_dat, forecast_dat) {
  if (!is.null(model$fit)) {
    baseline_pred <- as.numeric(predict(model$fit, newdata = baseline_dat, type = "response"))
    forecast_pred <- as.numeric(predict(model$fit, newdata = forecast_dat, type = "response", se.fit = TRUE)$fit)

    pred_link <- predict(model$fit, newdata = forecast_dat, type = "link", se.fit = TRUE)
    crit <- qnorm(1 - ((1 - confidence_level) / 2))
    forecast_lo <- exp(pred_link$fit - crit * pred_link$se.fit)
    forecast_hi <- exp(pred_link$fit + crit * pred_link$se.fit)

    return(tibble(
      predicted_count_baseline = baseline_pred,
      predicted_count_forecast = forecast_pred,
      predicted_count_forecast_lo = forecast_lo,
      predicted_count_forecast_hi = forecast_hi
    ))
  }

  baseline_rate <- safe_rate(baseline_dat$observed_count, baseline_dat$population) / 100000
  baseline_rate[!is.finite(baseline_rate)] <- mean(baseline_rate, na.rm = TRUE)
  baseline_rate[!is.finite(baseline_rate)] <- 0
  baseline_pred <- pmax(baseline_rate * baseline_dat$population, 0)
  forecast_pred <- baseline_pred * exp(model$beta_temp_c * (forecast_dat$annual_mean_temp_c - baseline_dat$annual_mean_temp_c))

  tibble(
    predicted_count_baseline = baseline_pred,
    predicted_count_forecast = forecast_pred,
    predicted_count_forecast_lo = NA_real_,
    predicted_count_forecast_hi = NA_real_
  )
}

# -----------------------------
# Load inputs
# -----------------------------
model_matrix_path <- pick_existing(model_matrix_candidates, "hvi_model_matrix_2019_2022.csv")
endpoint_meta_path <- pick_existing(endpoint_meta_candidates, "hvi_endpoint_metadata.csv")
endpoint_weight_path <- pick_existing(endpoint_weight_candidates, "09b_endpoint_weights.csv")
climate_path <- pick_existing(climate_candidates, "long climate panel")

hvi_model_matrix <- read_csv(model_matrix_path, show_col_types = FALSE) %>%
  clean_names() %>%
  mutate(
    community = clean_community(community),
    date = as.Date(date),
    year = as.integer(year)
  )

hvi_endpoint_metadata <- read_csv(endpoint_meta_path, show_col_types = FALSE) %>%
  clean_names() %>%
  mutate(endpoint_key = as.character(endpoint_key))

endpoint_weights <- read_csv(endpoint_weight_path, show_col_types = FALSE) %>%
  clean_names() %>%
  mutate(endpoint_key = as.character(endpoint_key))

climate_long <- read_csv(climate_path, show_col_types = FALSE) %>%
  clean_names()

climate_date_col <- c("event_date", "date")[c("event_date", "date") %in% names(climate_long)][1]
climate_temp_col <- c("tmax", "temp_max", "daily_max_temp", "temperature_c", "temperature")[c("tmax", "temp_max", "daily_max_temp", "temperature_c", "temperature") %in% names(climate_long)][1]
if (is.na(climate_date_col) || is.na(climate_temp_col)) {
  stop("Climate panel must include a date column and a daily maximum temperature column.")
}

climate_annual <- climate_long %>%
  transmute(
    community = if ("community" %in% names(climate_long)) clean_community(community) else NA_character_,
    date = as.Date(.data[[climate_date_col]]),
    year = year(date),
    tmax = suppressWarnings(as.numeric(.data[[climate_temp_col]]))
  ) %>%
  filter(!is.na(year), is.finite(tmax)) %>%
  group_by(year) %>%
  summarise(citywide_annual_mean_temp_c = mean(tmax, na.rm = TRUE), .groups = "drop")

trend_override <- suppressWarnings(as.numeric(hvi_env("HVI_FORECAST_TEMP_TREND_C_PER_YEAR", NA_character_)))
if (is.finite(trend_override)) {
  temp_trend_c_per_year <- trend_override
  temp_trend_source <- "HVI_FORECAST_TEMP_TREND_C_PER_YEAR"
} else {
  trend_fit <- lm(citywide_annual_mean_temp_c ~ year, data = climate_annual)
  temp_trend_c_per_year <- unname(coef(trend_fit)[["year"]])
  temp_trend_source <- basename(climate_path)
}

if (!is.finite(temp_trend_c_per_year)) {
  stop("Could not estimate or read a finite annual temperature trend.")
}

# -----------------------------
# Annual historical model matrix
# -----------------------------
temp_col <- c("tmax", "temp_max", "daily_max_temp", "temperature_c", "temperature")
temp_col <- temp_col[temp_col %in% names(hvi_model_matrix)][1]
if (is.na(temp_col)) stop("Could not find a temperature column in hvi_model_matrix.")

endpoint_keys <- hvi_endpoint_metadata$endpoint_key[
  hvi_endpoint_metadata$panel_outcome_col %in% names(hvi_model_matrix)
]

annual_base <- hvi_model_matrix %>%
  group_by(community, year) %>%
  summarise(
    days_observed = n_distinct(date),
    annual_mean_temp_c = mean(.data[[temp_col]], na.rm = TRUE),
    annual_max_temp_c = max(.data[[temp_col]], na.rm = TRUE),
    annual_mean_humidity = if ("humidity" %in% names(hvi_model_matrix)) mean(humidity, na.rm = TRUE) else NA_real_,
    population = if ("pop_offset" %in% names(hvi_model_matrix)) mean(pop_offset, na.rm = TRUE) else NA_real_,
    tree_canopy_pct = if ("tree_canopy_pct" %in% names(hvi_model_matrix)) mean(tree_canopy_pct, na.rm = TRUE) else NA_real_,
    tree_canopy_fraction = if ("tree_canopy_fraction" %in% names(hvi_model_matrix)) mean(tree_canopy_fraction, na.rm = TRUE) else NA_real_,
    ndvi = if ("ndvi" %in% names(hvi_model_matrix)) mean(ndvi, na.rm = TRUE) else NA_real_,
    ac_prob = if ("ac_prob" %in% names(hvi_model_matrix)) mean(ac_prob, na.rm = TRUE) else NA_real_,
    .groups = "drop"
  ) %>%
  mutate(
    annual_mean_temp_f = c_to_f(annual_mean_temp_c),
    annual_max_temp_f = c_to_f(annual_max_temp_c)
  )

outcome_annual <- map_dfr(seq_len(nrow(hvi_endpoint_metadata)), function(i) {
  ep <- hvi_endpoint_metadata$endpoint_key[i]
  outcome_col <- hvi_endpoint_metadata$panel_outcome_col[i]
  if (!outcome_col %in% names(hvi_model_matrix)) return(tibble())

  mrt_val <- suppressWarnings(as.numeric(hvi_endpoint_metadata$mrt[i]))
  heat_dose_col <- paste0("heat_dose_", ep)
  heat_dose_alt <- paste0("heat_dose__", ep)
  heat_dose_col <- c(heat_dose_col, heat_dose_alt)[c(heat_dose_col, heat_dose_alt) %in% names(hvi_model_matrix)][1]

  hvi_model_matrix %>%
    group_by(community, year) %>%
    summarise(
      observed_count = sum(.data[[outcome_col]], na.rm = TRUE),
      annual_heat_days = if (is.finite(mrt_val)) sum(.data[[temp_col]] > mrt_val, na.rm = TRUE) else NA_real_,
      annual_heat_dose = if (!is.na(heat_dose_col)) sum(.data[[heat_dose_col]], na.rm = TRUE) else if (is.finite(mrt_val)) sum(pmax(.data[[temp_col]] - mrt_val, 0), na.rm = TRUE) else NA_real_,
      .groups = "drop"
    ) %>%
    mutate(endpoint_key = ep)
})

frontend_forecast_annual_model_matrix <- outcome_annual %>%
  left_join(annual_base, by = c("community", "year")) %>%
  left_join(
    hvi_endpoint_metadata %>%
      select(endpoint_key, outcome_label, source, domain, mrt, max_lag),
    by = "endpoint_key"
  ) %>%
  mutate(
    event_family = source,
    observed_rate_per_100k = safe_rate(observed_count, population)
  ) %>%
  select(
    community, year, endpoint_key, event_family, outcome_label, domain,
    days_observed, annual_mean_temp_c, annual_mean_temp_f,
    annual_max_temp_c, annual_max_temp_f, annual_mean_humidity,
    annual_heat_days, annual_heat_dose,
    population, observed_count, observed_rate_per_100k,
    tree_canopy_pct, tree_canopy_fraction, ndvi, ac_prob,
    mrt, max_lag
  ) %>%
  arrange(endpoint_key, community, year)

# -----------------------------
# Forecast scenarios
# -----------------------------
baseline_year <- max(frontend_forecast_annual_model_matrix$year, na.rm = TRUE)
baseline_rows <- frontend_forecast_annual_model_matrix %>%
  filter(year == baseline_year)

daily_baseline_temperature <- hvi_model_matrix %>%
  filter(year == baseline_year) %>%
  select(community, date, all_of(temp_col))

forecast_heat_metrics <- map_dfr(seq_len(nrow(hvi_endpoint_metadata)), function(i) {
  ep <- hvi_endpoint_metadata$endpoint_key[i]
  mrt_val <- suppressWarnings(as.numeric(hvi_endpoint_metadata$mrt[i]))
  if (!is.finite(mrt_val) || !ep %in% endpoint_keys) return(tibble())

  map_dfr(horizons_years, function(h) {
    temp_delta <- temp_trend_c_per_year * h
    daily_baseline_temperature %>%
      mutate(forecast_temp_c = .data[[temp_col]] + temp_delta) %>%
      group_by(community) %>%
      summarise(
        endpoint_key = ep,
        horizon_years = h,
        forecast_annual_heat_days = sum(forecast_temp_c > mrt_val, na.rm = TRUE),
        forecast_annual_heat_dose = sum(pmax(forecast_temp_c - mrt_val, 0), na.rm = TRUE),
        .groups = "drop"
      )
  })
})

forecast_endpoint <- map_dfr(endpoint_keys, function(ep) {
  hist_ep <- frontend_forecast_annual_model_matrix %>%
    filter(endpoint_key == ep)
  baseline_ep <- baseline_rows %>%
    filter(endpoint_key == ep)

  if (nrow(hist_ep) == 0 || nrow(baseline_ep) == 0) return(tibble())

  model <- fit_endpoint_annual_model(hist_ep)

  map_dfr(horizons_years, function(h) {
    temp_delta <- temp_trend_c_per_year * h
    forecast_ep <- baseline_ep %>%
      mutate(
        horizon_years = h,
        baseline_year = baseline_year,
        forecast_year = baseline_year + h,
        temp_trend_c_per_year = temp_trend_c_per_year,
        temperature_delta_c = temp_delta,
        temperature_delta_f = temp_delta * 9 / 5,
        annual_mean_temp_c = annual_mean_temp_c + temp_delta,
        annual_mean_temp_f = c_to_f(annual_mean_temp_c),
        annual_max_temp_c = annual_max_temp_c + temp_delta,
        annual_max_temp_f = c_to_f(annual_max_temp_c),
        baseline_annual_heat_days = annual_heat_days,
        baseline_annual_heat_dose = annual_heat_dose,
        model_type = model$model_type,
        model_status = model$status,
        beta_temp_c = model$beta_temp_c
      ) %>%
      left_join(forecast_heat_metrics, by = c("community", "endpoint_key", "horizon_years"))

    pred <- predict_endpoint_annual(model, baseline_ep, forecast_ep)

    bind_cols(
      forecast_ep,
      pred
    ) %>%
      mutate(
        predicted_count_change = predicted_count_forecast - predicted_count_baseline,
        predicted_pct_change = if_else(predicted_count_baseline > 0, 100 * predicted_count_change / predicted_count_baseline, NA_real_),
        predicted_rate_baseline_per_100k = safe_rate(predicted_count_baseline, population),
        predicted_rate_forecast_per_100k = safe_rate(predicted_count_forecast, population),
        predicted_rate_change_per_100k = predicted_rate_forecast_per_100k - predicted_rate_baseline_per_100k
      )
  })
}) %>%
  select(
    community, baseline_year, forecast_year, horizon_years,
    endpoint_key, event_family, outcome_label, domain,
    annual_mean_temp_c, annual_mean_temp_f,
    annual_max_temp_c, annual_max_temp_f,
    temperature_delta_c, temperature_delta_f,
    temp_trend_c_per_year,
    population,
    predicted_count_baseline, predicted_count_forecast,
    predicted_count_forecast_lo, predicted_count_forecast_hi,
    predicted_count_change, predicted_pct_change,
    predicted_rate_baseline_per_100k,
    predicted_rate_forecast_per_100k,
    predicted_rate_change_per_100k,
    baseline_annual_heat_days, baseline_annual_heat_dose,
    forecast_annual_heat_days, forecast_annual_heat_dose,
    model_type, model_status, beta_temp_c
  ) %>%
  arrange(horizon_years, endpoint_key, community)

weight_cols <- endpoint_weights %>%
  select(any_of(c("endpoint_key", "endpoint_weight", "source_weight", "performance_weight"))) %>%
  mutate(endpoint_weight = coalesce(endpoint_weight, 1))

forecast_endpoint_weighted <- forecast_endpoint %>%
  left_join(weight_cols, by = "endpoint_key") %>%
  mutate(endpoint_weight = coalesce(endpoint_weight, 1))

dominant_forecast_endpoint <- forecast_endpoint_weighted %>%
  group_by(community, baseline_year, forecast_year, horizon_years) %>%
  slice_max(order_by = predicted_count_change, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  transmute(
    community, baseline_year, forecast_year, horizon_years,
    dominant_endpoint = endpoint_key,
    dominant_endpoint_label = outcome_label,
    dominant_endpoint_source = event_family,
    dominant_endpoint_count_change = predicted_count_change
  )

frontend_long_term_forecast_overall <- forecast_endpoint_weighted %>%
  group_by(community, baseline_year, forecast_year, horizon_years) %>%
  summarise(
    annual_mean_temp_c = first(annual_mean_temp_c),
    annual_mean_temp_f = first(annual_mean_temp_f),
    temperature_delta_c = first(temperature_delta_c),
    temperature_delta_f = first(temperature_delta_f),
    population = first(population),
    total_predicted_count_baseline = sum(predicted_count_baseline, na.rm = TRUE),
    total_predicted_count_forecast = sum(predicted_count_forecast, na.rm = TRUE),
    total_predicted_count_change = total_predicted_count_forecast - total_predicted_count_baseline,
    overall_weighted_count_change = weighted.mean(predicted_count_change, w = endpoint_weight, na.rm = TRUE),
    n_endpoints = n_distinct(endpoint_key),
    .groups = "drop"
  ) %>%
  group_by(horizon_years) %>%
  mutate(
    forecast_risk_0_100 = rescale_positive_0_100(overall_weighted_count_change),
    forecast_tier = cut(
      forecast_risk_0_100,
      breaks = c(-Inf, 25, 50, 75, Inf),
      labels = c("Low", "Moderate", "High", "Very High"),
      include.lowest = TRUE,
      right = TRUE
    )
  ) %>%
  ungroup() %>%
  left_join(dominant_forecast_endpoint, by = c("community", "baseline_year", "forecast_year", "horizon_years")) %>%
  arrange(horizon_years, desc(forecast_risk_0_100), community)

frontend_long_term_forecast_metadata <- tibble(
  generated_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
  baseline_year = baseline_year,
  horizons_years = paste(horizons_years, collapse = "|"),
  climate_trend_source = temp_trend_source,
  climate_trend_start_year = min(climate_annual$year, na.rm = TRUE),
  climate_trend_end_year = max(climate_annual$year, na.rm = TRUE),
  temp_trend_c_per_year = temp_trend_c_per_year,
  temp_trend_f_per_year = temp_trend_c_per_year * 9 / 5,
  confidence_level = confidence_level,
  model_note = paste(
    "Annual forecasts use quasipoisson endpoint models of annual community counts against annual mean daily maximum temperature,",
    "with community fixed effects and population offsets when estimable. These are planning scenarios, not climate-model projections."
  )
)

# -----------------------------
# Export
# -----------------------------
write_csv(frontend_forecast_annual_model_matrix, file.path(frontend_dir, "frontend_forecast_annual_model_matrix.csv"))
write_csv(forecast_endpoint, file.path(frontend_dir, "frontend_long_term_forecast_endpoint.csv"))
write_csv(frontend_long_term_forecast_overall, file.path(frontend_dir, "frontend_long_term_forecast_overall.csv"))
write_csv(frontend_long_term_forecast_metadata, file.path(frontend_dir, "frontend_long_term_forecast_metadata.csv"))

assign("frontend_forecast_annual_model_matrix", frontend_forecast_annual_model_matrix, envir = .GlobalEnv)
assign("frontend_long_term_forecast_endpoint", forecast_endpoint, envir = .GlobalEnv)
assign("frontend_long_term_forecast_overall", frontend_long_term_forecast_overall, envir = .GlobalEnv)
assign("frontend_long_term_forecast_metadata", frontend_long_term_forecast_metadata, envir = .GlobalEnv)

message("13 complete. Long-term forecast frontend exports written to: ", frontend_dir)
