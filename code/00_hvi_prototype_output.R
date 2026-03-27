# ======================================================================================
# HVI 2.0 | Frontend export bundle
# Assumes required objects already exist in the R environment:
#   panel, panel_causes, panel_overlap, hvi_proto, ca or comm_areas,
#   mrt_table, curve_table
#
# Optional objects used if present:
#   citywide_heat_effects, hvi_ranked, community_month_normals, heat_ref
#
# Writes frontend-ready files to:
#   results/frontend_exports/
# ======================================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(purrr)
  library(lubridate)
  library(data.table)
  library(sf)
  library(jsonlite)
  library(janitor)
})

# --------------------------------------------------------------------------------------
# 0) Setup
# --------------------------------------------------------------------------------------
out_dir <- "results/frontend_exports"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

`%||%` <- function(x, y) if (!is.null(x)) x else y

to_title_safe <- function(x) {
  x %>%
    str_replace_all("_", " ") %>%
    str_squish() %>%
    str_to_title()
}

rescale_01 <- function(x) {
  rng <- range(x, na.rm = TRUE)
  if (!all(is.finite(rng)) || diff(rng) == 0) return(rep(0.5, length(x)))
  (x - rng[1]) / diff(rng)
}

rescale_0100 <- function(x) 100 * rescale_01(x)

safe_ntile_label <- function(x) {
  q <- ntile(x, 4)
  case_when(
    q == 1 ~ "low",
    q == 2 ~ "moderate",
    q == 3 ~ "high",
    q == 4 ~ "severe",
    TRUE ~ NA_character_
  )
}

std_comm <- function(x) {
  x %>%
    as.character() %>%
    str_trim() %>%
    str_to_upper()
}

pick_first_existing <- function(df, choices) {
  nm <- intersect(choices, names(df))
  if (length(nm) == 0) return(NULL)
  nm[[1]]
}

first_nonnull_object <- function(...) {
  xs <- list(...)
  xs[[which(vapply(xs, function(z) !is.null(z), logical(1)))[1]]]
}

# --------------------------------------------------------------------------------------
# 1) Resolve map object
# --------------------------------------------------------------------------------------
map_sf <- NULL

if (exists("ca")) {
  map_sf <- get("ca")
} else if (exists("comm_areas")) {
  map_sf <- get("comm_areas")
}

if (is.null(map_sf)) {
  stop("Need either `ca` or `comm_areas` in the environment.")
}

map_sf <- map_sf %>%
  janitor::clean_names()

comm_col <- pick_first_existing(
  map_sf,
  c("community", "ca_name", "commarea", "community_area", "community_area_name", "name")
)

if (is.null(comm_col)) {
  stop("Could not identify community name field in map object.")
}

map_sf <- map_sf %>%
  mutate(
    community = std_comm(.data[[comm_col]])
  )

# build community_area_id
community_lookup <- tibble(
  community = sort(unique(map_sf$community))
) %>%
  mutate(
    community_area_id = sprintf("CA%02d", row_number()),
    display_name = str_to_title(community)
  )

map_sf <- map_sf %>%
  left_join(community_lookup, by = "community")

# --------------------------------------------------------------------------------------
# 2) Resolve panel objects
# --------------------------------------------------------------------------------------
if (!exists("panel")) stop("`panel` must exist.")
panel_base <- get("panel") %>%
  mutate(community = std_comm(community))

panel_use <- if (exists("panel_causes")) {
  get("panel_causes") %>% mutate(community = std_comm(community))
} else {
  panel_base
}

if (!"event_date" %in% names(panel_use)) {
  stop("Expected `event_date` in panel data.")
}

panel_use <- panel_use %>%
  mutate(event_date = as.Date(event_date)) %>%
  left_join(community_lookup, by = "community")

if (exists("panel_overlap")) {
  panel_overlap_use <- get("panel_overlap") %>%
    mutate(
      community = std_comm(community),
      event_date = as.Date(event_date)
    ) %>%
    left_join(community_lookup, by = "community")
} else {
  panel_overlap_use <- panel_use
}

if (!exists("hvi_proto")) stop("`hvi_proto` must exist.")
hvi_proto_use <- get("hvi_proto") %>%
  mutate(community = std_comm(community)) %>%
  left_join(community_lookup, by = "community")

# --------------------------------------------------------------------------------------
# 3) Community master
# --------------------------------------------------------------------------------------
centroids <- st_centroid(st_geometry(map_sf))
coords <- st_coordinates(centroids)

community_area_master <- map_sf %>%
  st_drop_geometry() %>%
  transmute(
    community_area_id,
    community_area_name = display_name,
    display_name,
    geometry_id = community_area_id
  ) %>%
  distinct() %>%
  mutate(
    population = NA_real_,
    region_group = NA_character_,
    neighborhood_summary = NA_character_
  ) %>%
  bind_cols(
    tibble(
      centroid_lon = coords[, 1],
      centroid_lat = coords[, 2]
    )
  ) %>%
  left_join(
    hvi_proto_use %>%
      transmute(
        community_area_id,
        baseline_vulnerability_score = if ("hvi2_allcause_pct" %in% names(.)) 100 * hvi2_allcause_pct else NA_real_,
        baseline_vulnerability_tier = if ("hvi2_allcause" %in% names(.)) safe_ntile_label(hvi2_allcause) else NA_character_
      ),
    by = "community_area_id"
  ) %>%
  mutate(
    acs_year = 2022L,
    metadata_version = "hvi2_frontend_v1"
  )

# --------------------------------------------------------------------------------------
# 4) Daily risk scores
#    Uses current prototype to create a day-level score:
#    hazard = day-level heat/exposure signal
#    vulnerability = community prototype HVI percentile
#    impact = observed burden that day
# --------------------------------------------------------------------------------------
hazard_components <- intersect(
  c("tmax", "tmean", "humidity", "pm25", "no2", "heat_index"),
  names(panel_use)
)

panel_daily <- panel_use %>%
  mutate(
    heat_index = case_when(
      "heat_index" %in% names(.) ~ heat_index,
      all(c("tmax", "humidity") %in% names(.)) ~ tmax + 0.05 * humidity,
      TRUE ~ NA_real_
    ),
    observed_total_events =
      coalesce(.data[["deaths"]], 0) +
      coalesce(.data[["ed_visits"]], 0) +
      coalesce(.data[["ems_calls"]], 0)
  )

# day-level hazard score
panel_daily <- panel_daily %>%
  group_by(community_area_id) %>%
  mutate(
    z_tmax = if ("tmax" %in% names(.)) as.numeric(scale(tmax)) else NA_real_,
    z_heat_index = if ("heat_index" %in% names(.)) as.numeric(scale(heat_index)) else NA_real_,
    z_pm25 = if ("pm25" %in% names(.)) as.numeric(scale(pm25)) else NA_real_,
    z_no2 = if ("no2" %in% names(.)) as.numeric(scale(no2)) else NA_real_,
    z_obs_events = as.numeric(scale(observed_total_events))
  ) %>%
  ungroup()

panel_daily <- panel_daily %>%
  mutate(
    hazard_raw = rowMeans(
      cbind(z_tmax, z_heat_index, z_pm25, z_no2),
      na.rm = TRUE
    ),
    hazard_raw = ifelse(is.nan(hazard_raw), NA_real_, hazard_raw)
  )

panel_daily <- panel_daily %>%
  left_join(
    hvi_proto_use %>%
      transmute(
        community_area_id,
        vulnerability_raw = if ("hvi2_allcause" %in% names(.)) hvi2_allcause else NA_real_,
        vulnerability_pct = if ("hvi2_allcause_pct" %in% names(.)) 100 * hvi2_allcause_pct else NA_real_,
        mort_excess_per_hot_day = .data[["mort_excess_per_hot_day"]] %||% NA_real_,
        ed_excess_per_hot_day   = .data[["ed_excess_per_hot_day"]] %||% NA_real_,
        ems_excess_per_hot_day  = .data[["ems_excess_per_hot_day"]] %||% NA_real_
      ),
    by = "community_area_id"
  ) %>%
  mutate(
    impact_raw = z_obs_events
  )

panel_daily <- panel_daily %>%
  mutate(
    hazard_score_0_100 = rescale_0100(hazard_raw),
    vulnerability_score_0_100 = ifelse(
      is.na(vulnerability_pct),
      rescale_0100(vulnerability_raw),
      vulnerability_pct
    ),
    impact_score_0_100 = rescale_0100(impact_raw),
    overall_risk_score_0_100 = 0.40 * hazard_score_0_100 +
      0.35 * vulnerability_score_0_100 +
      0.25 * impact_score_0_100,
    percentile_in_city = percent_rank(overall_risk_score_0_100),
    risk_tier = safe_ntile_label(overall_risk_score_0_100),
    confidence_low = pmax(0, overall_risk_score_0_100 - 7.5),
    confidence_high = pmin(100, overall_risk_score_0_100 + 7.5),
    dominant_driver = case_when(
      coalesce(z_tmax, -Inf) >= pmax(coalesce(z_pm25, -Inf), coalesce(z_no2, -Inf), na.rm = TRUE) ~ "heat",
      coalesce(z_pm25, -Inf) >= pmax(coalesce(z_tmax, -Inf), coalesce(z_no2, -Inf), na.rm = TRUE) ~ "pm25",
      coalesce(z_no2, -Inf) >= pmax(coalesce(z_tmax, -Inf), coalesce(z_pm25, -Inf), na.rm = TRUE) ~ "no2",
      TRUE ~ "mixed"
    ),
    dominant_exposure = dominant_driver,
    top_outcome = case_when(
      coalesce(ed_visits, 0) >= pmax(coalesce(deaths, 0), coalesce(ems_calls, 0), na.rm = TRUE) ~ "ed",
      coalesce(ems_calls, 0) >= pmax(coalesce(deaths, 0), coalesce(ed_visits, 0), na.rm = TRUE) ~ "ems",
      TRUE ~ "mortality"
    ),
    alert_flag = risk_tier %in% c("high", "severe"),
    severe_alert_flag = risk_tier %in% c("severe")
  )

daily_risk_scores <- panel_daily %>%
  transmute(
    date = event_date,
    community_area_id,
    overall_risk_score_0_100 = round(overall_risk_score_0_100, 2),
    risk_tier,
    hazard_score_0_100 = round(hazard_score_0_100, 2),
    vulnerability_score_0_100 = round(vulnerability_score_0_100, 2),
    impact_score_0_100 = round(impact_score_0_100, 2),
    dominant_driver,
    dominant_exposure,
    top_outcome,
    percentile_in_city = round(100 * percentile_in_city, 1),
    confidence_low = round(confidence_low, 2),
    confidence_high = round(confidence_high, 2),
    alert_flag,
    severe_alert_flag,
    model_version = "hvi2_prototype_frontend_v1",
    run_timestamp = format(Sys.time(), tz = "UTC", usetz = TRUE)
  )

# --------------------------------------------------------------------------------------
# 5) Daily outcome forecast
#    For observed dates, use actual daily counts and community-level excess/hot-day burden
# --------------------------------------------------------------------------------------
outcomes_long <- panel_daily %>%
  select(
    event_date, community_area_id,
    deaths, ed_visits, ems_calls,
    mort_excess_per_hot_day, ed_excess_per_hot_day, ems_excess_per_hot_day,
    risk_tier
  ) %>%
  pivot_longer(
    cols = c(deaths, ed_visits, ems_calls),
    names_to = "outcome",
    values_to = "expected_events"
  ) %>%
  mutate(
    expected_excess_events = case_when(
      outcome == "deaths" ~ mort_excess_per_hot_day,
      outcome == "ed_visits" ~ ed_excess_per_hot_day,
      outcome == "ems_calls" ~ ems_excess_per_hot_day,
      TRUE ~ NA_real_
    ),
    baseline_events = pmax(expected_events - coalesce(expected_excess_events, 0), 0),
    relative_risk = ifelse(baseline_events > 0, expected_events / baseline_events, NA_real_),
    attributable_fraction = ifelse(!is.na(relative_risk) & relative_risk > 0, (relative_risk - 1) / relative_risk, NA_real_),
    ci_low = pmax(0, coalesce(expected_excess_events, 0) * 0.8),
    ci_high = coalesce(expected_excess_events, 0) * 1.2,
    outcome = recode(
      outcome,
      deaths = "mortality",
      ed_visits = "ed",
      ems_calls = "ems"
    ),
    time_horizon = "observed",
    units = "events/day"
  )

daily_outcome_forecast <- outcomes_long %>%
  transmute(
    date = event_date,
    community_area_id,
    outcome,
    time_horizon,
    expected_events = round(expected_events, 3),
    expected_excess_events = round(expected_excess_events, 3),
    relative_risk = round(relative_risk, 3),
    attributable_fraction = round(attributable_fraction, 3),
    baseline_events = round(baseline_events, 3),
    ci_low = round(ci_low, 3),
    ci_high = round(ci_high, 3),
    units,
    model_version = "hvi2_prototype_frontend_v1"
  )

# --------------------------------------------------------------------------------------
# 6) Thresholds from MRT table
# --------------------------------------------------------------------------------------
if (!exists("mrt_table")) {
  warning("`mrt_table` not found. Building empty thresholds table.")
  thresholds <- tibble(
    community_area_id = character(),
    outcome = character(),
    exposure_metric = character(),
    reference_value = numeric(),
    alert_threshold = numeric(),
    severe_threshold = numeric(),
    lag_structure = character(),
    peak_lag_days = integer(),
    threshold_definition = character(),
    model_version = character()
  )
} else {
  mrt_use <- get("mrt_table")
  
  thresholds <- community_lookup %>%
    crossing(mrt_use) %>%
    transmute(
      community_area_id,
      outcome = outcome,
      exposure_metric = "tmax",
      reference_value = mrt,
      alert_threshold = mrt + 3,
      severe_threshold = mrt + 6,
      lag_structure = paste0("0-", lag, " days"),
      peak_lag_days = lag,
      threshold_definition = "Prototype frontend threshold derived from citywide MRT + fixed offsets",
      model_version = "mrt_warmseason_frontend_v1"
    )
}

# --------------------------------------------------------------------------------------
# 7) Response curves
# --------------------------------------------------------------------------------------
if (!exists("curve_table")) {
  warning("`curve_table` not found. Building empty response_curves table.")
  response_curves <- tibble(
    curve_scope = character(),
    community_area_id = character(),
    outcome = character(),
    exposure_metric = character(),
    exposure_value = numeric(),
    relative_risk = numeric(),
    ci_low = numeric(),
    ci_high = numeric(),
    reference_value = numeric(),
    alert_threshold = numeric(),
    severe_threshold = numeric(),
    model_version = character()
  )
} else {
  curve_use <- get("curve_table")
  
  response_curves <- curve_use %>%
    transmute(
      curve_scope = "citywide",
      community_area_id = NA_character_,
      outcome,
      exposure_metric = "tmax",
      exposure_value = temp,
      relative_risk = rr,
      ci_low = rr_low,
      ci_high = rr_high,
      reference_value = mrt,
      alert_threshold = mrt + 3,
      severe_threshold = mrt + 6,
      model_version = "mrt_warmseason_frontend_v1"
    )
}

# --------------------------------------------------------------------------------------
# 8) Timeseries
# --------------------------------------------------------------------------------------
timeseries <- panel_daily %>%
  left_join(
    daily_risk_scores %>% select(date, community_area_id, overall_risk_score_0_100, risk_tier),
    by = c("event_date" = "date", "community_area_id")
  ) %>%
  transmute(
    date = event_date,
    community_area_id,
    tmax = .data[["tmax"]] %||% NA_real_,
    tmin = .data[["tmin"]] %||% NA_real_,
    tmean = .data[["tmean"]] %||% NA_real_,
    heat_index = .data[["heat_index"]] %||% NA_real_,
    humidity = .data[["humidity"]] %||% NA_real_,
    pm25 = .data[["pm25"]] %||% NA_real_,
    no2 = .data[["no2"]] %||% NA_real_,
    risk_score = round(overall_risk_score_0_100, 2),
    risk_tier,
    expected_excess_events_total = round(
      coalesce(mort_excess_per_hot_day, 0) +
        coalesce(ed_excess_per_hot_day, 0) +
        coalesce(ems_excess_per_hot_day, 0), 3
    ),
    expected_excess_events_ed = round(coalesce(ed_excess_per_hot_day, 0), 3),
    expected_excess_events_ems = round(coalesce(ems_excess_per_hot_day, 0), 3),
    expected_excess_events_mortality = round(coalesce(mort_excess_per_hot_day, 0), 3),
    observed_flag = TRUE,
    model_version = "hvi2_prototype_frontend_v1"
  )

# --------------------------------------------------------------------------------------
# 9) Forecast summary
#    Placeholder v1: repeats latest observed exposures forward 3 days by community.
#    Replace later with a true forecast ingest.
# --------------------------------------------------------------------------------------
latest_day <- max(panel_daily$event_date, na.rm = TRUE)

latest_obs <- panel_daily %>%
  filter(event_date == latest_day) %>%
  select(
    community_area_id, community,
    forecast_tmax = tmax,
    forecast_heat_index = heat_index,
    forecast_pm25 = pm25,
    forecast_no2 = no2,
    vulnerability_score_0_100
  )

forecast_summary <- crossing(
  latest_obs,
  tibble(forecast_day = 1:3)
) %>%
  mutate(
    issue_date = latest_day,
    forecast_date = latest_day + forecast_day,
    # persistence forecast with slight decay
    forecast_tmax = forecast_tmax,
    forecast_heat_index = forecast_heat_index,
    forecast_pm25 = forecast_pm25,
    predicted_risk_score = pmin(
      100,
      round(
        0.55 * rescale_0100(forecast_tmax) +
          0.15 * rescale_0100(forecast_heat_index) +
          0.10 * rescale_0100(forecast_pm25) +
          0.20 * vulnerability_score_0_100,
        2
      )
    ),
    predicted_tier = safe_ntile_label(predicted_risk_score),
    predicted_dominant_driver = case_when(
      coalesce(forecast_pm25, -Inf) > coalesce(forecast_tmax, -Inf) ~ "pm25",
      TRUE ~ "heat"
    ),
    predicted_top_outcome = "mixed",
    delta_vs_today = 0,
    confidence_low = pmax(0, predicted_risk_score - 8),
    confidence_high = pmin(100, predicted_risk_score + 8),
    model_version = "placeholder_persistence_forecast_v1"
  ) %>%
  select(
    issue_date,
    forecast_date,
    forecast_day,
    community_area_id,
    forecast_tmax,
    forecast_heat_index,
    forecast_pm25,
    predicted_risk_score,
    predicted_tier,
    predicted_dominant_driver,
    predicted_top_outcome,
    delta_vs_today,
    confidence_low,
    confidence_high,
    model_version
  )

# --------------------------------------------------------------------------------------
# 10) Map metadata
# --------------------------------------------------------------------------------------
map_metadata <- tibble(
  geography = "Chicago community areas",
  id_field = "community_area_id",
  name_field = "community_area_name",
  geometry_file = "community_areas.geojson",
  default_map_metric = "overall_risk_score_0_100",
  default_join_key = "community_area_id",
  n_areas = n_distinct(community_area_master$community_area_id),
  crs = tryCatch(st_crs(map_sf)$input, error = function(e) NA_character_)
)

# --------------------------------------------------------------------------------------
# 11) Model metadata
# --------------------------------------------------------------------------------------
data_start <- min(panel_use$event_date, na.rm = TRUE)
data_end   <- max(panel_use$event_date, na.rm = TRUE)

model_metadata <- tibble(
  model_version = "hvi2_prototype_frontend_v1",
  run_timestamp = format(Sys.time(), tz = "UTC", usetz = TRUE),
  data_start_date = as.character(data_start),
  data_end_date = as.character(data_end),
  forecast_issue_time = as.character(latest_day),
  weather_source = "Internal pipeline object(s) already loaded in R",
  pollution_source = if (any(c("pm25", "no2") %in% names(panel_use))) "Internal pipeline object(s) already loaded in R" else NA_character_,
  health_data_refresh_date = as.character(data_end),
  notes = "Prototype frontend export built from existing HVI 2.0 environment objects.",
  method_summary = paste(
    "Day-level frontend bundle using observed panel data, community-level prototype HVI,",
    "and MRT response curves. Forecast table is persistence-based placeholder."
  ),
  contact = "Peter Graffy"
)

# --------------------------------------------------------------------------------------
# 12) GeoJSON export
# --------------------------------------------------------------------------------------
map_geo <- map_sf %>%
  select(community_area_id, community_area_name = display_name)

st_write(map_geo, file.path(out_dir, "community_areas.geojson"), delete_dsn = TRUE, quiet = TRUE)

# --------------------------------------------------------------------------------------
# 13) Write CSV/JSON exports
# --------------------------------------------------------------------------------------
fwrite(community_area_master, file.path(out_dir, "community_area_master.csv"))
fwrite(daily_risk_scores, file.path(out_dir, "daily_risk_scores.csv"))
fwrite(daily_outcome_forecast, file.path(out_dir, "daily_outcome_forecast.csv"))
fwrite(thresholds, file.path(out_dir, "thresholds.csv"))
fwrite(response_curves, file.path(out_dir, "response_curves.csv"))
fwrite(timeseries, file.path(out_dir, "timeseries.csv"))
fwrite(forecast_summary, file.path(out_dir, "forecast_summary.csv"))
fwrite(map_metadata, file.path(out_dir, "map_metadata.csv"))
fwrite(model_metadata, file.path(out_dir, "model_metadata.csv"))

write_json(community_area_master, file.path(out_dir, "community_area_master.json"), dataframe = "rows", pretty = TRUE, auto_unbox = TRUE, na = "null")
write_json(daily_risk_scores, file.path(out_dir, "daily_risk_scores.json"), dataframe = "rows", pretty = FALSE, auto_unbox = TRUE, na = "null")
write_json(daily_outcome_forecast, file.path(out_dir, "daily_outcome_forecast.json"), dataframe = "rows", pretty = FALSE, auto_unbox = TRUE, na = "null")
write_json(thresholds, file.path(out_dir, "thresholds.json"), dataframe = "rows", pretty = FALSE, auto_unbox = TRUE, na = "null")
write_json(response_curves, file.path(out_dir, "response_curves.json"), dataframe = "rows", pretty = FALSE, auto_unbox = TRUE, na = "null")
write_json(timeseries, file.path(out_dir, "timeseries.json"), dataframe = "rows", pretty = FALSE, auto_unbox = TRUE, na = "null")
write_json(forecast_summary, file.path(out_dir, "forecast_summary.json"), dataframe = "rows", pretty = FALSE, auto_unbox = TRUE, na = "null")
write_json(map_metadata, file.path(out_dir, "map_metadata.json"), dataframe = "rows", pretty = TRUE, auto_unbox = TRUE, na = "null")
write_json(model_metadata, file.path(out_dir, "model_metadata.json"), dataframe = "rows", pretty = TRUE, auto_unbox = TRUE, na = "null")

# --------------------------------------------------------------------------------------
# 14) Tiny API index manifest for Replit
# --------------------------------------------------------------------------------------
api_manifest <- list(
  generated_at = format(Sys.time(), tz = "UTC", usetz = TRUE),
  files = list(
    community_area_master = "community_area_master.json",
    daily_risk_scores = "daily_risk_scores.json",
    daily_outcome_forecast = "daily_outcome_forecast.json",
    thresholds = "thresholds.json",
    response_curves = "response_curves.json",
    timeseries = "timeseries.json",
    forecast_summary = "forecast_summary.json",
    map_metadata = "map_metadata.json",
    model_metadata = "model_metadata.json",
    community_areas_geojson = "community_areas.geojson"
  )
)

write_json(api_manifest, file.path(out_dir, "manifest.json"), pretty = TRUE, auto_unbox = TRUE)

# --------------------------------------------------------------------------------------
# 15) Console summary
# --------------------------------------------------------------------------------------
cat("\nFrontend export bundle written to:", out_dir, "\n")
cat("Files:\n")
print(list.files(out_dir))
cat("\nRows:\n")
cat("community_area_master :", nrow(community_area_master), "\n")
cat("daily_risk_scores     :", nrow(daily_risk_scores), "\n")
cat("daily_outcome_forecast:", nrow(daily_outcome_forecast), "\n")
cat("thresholds            :", nrow(thresholds), "\n")
cat("response_curves       :", nrow(response_curves), "\n")
cat("timeseries            :", nrow(timeseries), "\n")
cat("forecast_summary      :", nrow(forecast_summary), "\n")