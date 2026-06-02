# ================================================================================================
# HVI 2.0 | 09e_build_frontend_exports.R
# Prepare lightweight long and wide exports for a future interactive map/dashboard.
# ================================================================================================

source(file.path(if (dir.exists("code")) "code" else ".", "00_config.R"))
hvi_source("09_utils_hvi.R")

# -----------------------------
# CONFIG
# -----------------------------
project_dir <- HVI_PATHS$private
out_dir <- ensure_output_dir(project_dir, "09_model_outputs")
frontend_dir <- file.path(out_dir, "frontend_exports")
dir.create(frontend_dir, showWarnings = FALSE, recursive = TRUE)

# -----------------------------
# REQUIRED OBJECTS
# -----------------------------
required_objs <- c(
  "ca_year_overall_structural_hvi",
  "ca_year_endpoint_vulnerability",
  "temp_grid_overall_risk",
  "temp_grid_endpoint_risk",
  "ca_day_overall_operational_hvi",
  "ca_day_endpoint_risk",
  "hvi_model_matrix"
)

missing_objs <- required_objs[!vapply(required_objs, exists, logical(1), envir = .GlobalEnv)]
if (length(missing_objs) > 0) stop("Run 09b-09d first. Missing objects: ", paste(missing_objs, collapse = ", "))

ca_year_overall_structural_hvi <- get("ca_year_overall_structural_hvi", envir = .GlobalEnv) %>%
  mutate(community = as.character(community), year = as.integer(as.character(year)))
ca_year_endpoint_vulnerability <- get("ca_year_endpoint_vulnerability", envir = .GlobalEnv) %>%
  mutate(community = as.character(community), year = as.integer(as.character(year)))
temp_grid_overall_risk <- get("temp_grid_overall_risk", envir = .GlobalEnv)
temp_grid_endpoint_risk <- get("temp_grid_endpoint_risk", envir = .GlobalEnv)
ca_day_overall_operational_hvi <- get("ca_day_overall_operational_hvi", envir = .GlobalEnv) %>%
  mutate(community = as.character(community), date = as.Date(date), year = as.integer(as.character(year)))
ca_day_endpoint_risk <- get("ca_day_endpoint_risk", envir = .GlobalEnv) %>%
  mutate(community = as.character(community), date = as.Date(date), year = as.integer(as.character(year)))
hvi_model_matrix <- get("hvi_model_matrix", envir = .GlobalEnv) %>%
  clean_names() %>%
  mutate(community = as.character(community), date = as.Date(date), year = as.integer(as.character(year)))

# -----------------------------
# PUBLIC-SAFE LABEL HELPERS
# -----------------------------
driver_label_lookup <- c(
  z_ac_prob = "Estimated household air-conditioning access",
  z_ndvi = "Neighborhood greenness",
  z_mean_ndvi = "Neighborhood greenness",
  z_tree_canopy_pct = "Tree canopy coverage",
  z_tree_canopy_fraction = "Tree canopy coverage",
  z_mean_tree_canopy = "Tree canopy coverage",
  z_pm25 = "Fine particulate air pollution",
  z_no2 = "Nitrogen dioxide air pollution",
  z_pop_density_km2 = "Population density",
  z_median_age = "Median age",
  z_mean_female = "Female population share",
  z_mean_employed = "Employment",
  z_mean_hs = "High-school education",
  z_mean_asian = "Asian population share",
  z_svi_ep_age17 = "Children under age 18",
  z_svi_ep_age65 = "Older adults",
  z_svi_ep_crowd = "Crowded housing",
  z_svi_ep_disabl = "Disability",
  z_svi_ep_groupq = "Group quarters residence",
  z_svi_ep_hburd = "Housing cost burden",
  z_svi_ep_limeng = "Limited English proficiency",
  z_svi_ep_minrty = "Race and ethnicity composition",
  z_svi_ep_mobile = "Mobile homes",
  z_svi_ep_munit = "Multi-unit housing",
  z_svi_ep_nohsdp = "No high-school diploma",
  z_svi_ep_noint = "No internet access",
  z_svi_ep_noveh = "No vehicle access",
  z_svi_ep_pov = "Poverty",
  z_svi_ep_sngpnt = "Single-parent households",
  z_svi_ep_unemp = "Unemployment",
  z_svi_ep_uninsur = "Uninsured population",
  z_svi_rpl_theme1 = "Socioeconomic vulnerability",
  z_svi_rpl_theme2 = "Household composition vulnerability",
  z_svi_rpl_theme3 = "Race, ethnicity, and language vulnerability",
  z_svi_rpl_theme4 = "Housing and transportation vulnerability"
)

label_driver <- function(x) {
  x <- as.character(x)
  out <- unname(driver_label_lookup[x])
  missing <- is.na(out) & !is.na(x) & nzchar(x)
  out[missing] <- x[missing] %>%
    str_remove("^z_") %>%
    str_replace_all("^svi_ep_", "") %>%
    str_replace_all("^svi_rpl_", "svi_") %>%
    str_replace_all("_", " ") %>%
    str_to_sentence()
  out
}

impact_tier <- function(x) {
  cut(
    x,
    breaks = c(-Inf, 25, 50, 75, Inf),
    labels = c("Low", "Moderate", "High", "Very High"),
    include.lowest = TRUE,
    right = TRUE
  )
}

c_to_f <- function(x) (x * 9 / 5) + 32

num_env <- function(name, default) {
  val <- suppressWarnings(as.numeric(Sys.getenv(name, unset = NA_character_)))
  ifelse(is.na(val), default, val)
}

cost_assumptions_path <- Sys.getenv("HVI_EXCESS_COST_ASSUMPTIONS", unset = "")
if (nzchar(cost_assumptions_path) && file.exists(cost_assumptions_path)) {
  excess_cost_assumptions <- read_csv(cost_assumptions_path, show_col_types = FALSE) %>%
    clean_names()
} else {
  excess_cost_assumptions <- tibble(
    event_family = c("ED", "EMS", "Mortality"),
    endpoint_key = c("ed_visits", "ems_calls", "deaths"),
    unit_cost_usd = c(
      num_env("HVI_UNIT_COST_ED_VISIT_USD", 2400),
      num_env("HVI_UNIT_COST_EMS_CALL_USD", 3000),
      num_env("HVI_UNIT_COST_DEATH_USD", 11600000)
    ),
    cost_basis = c(
      "Approximate Chicago-area emergency department visit medical cost. Override with HVI_UNIT_COST_ED_VISIT_USD or HVI_EXCESS_COST_ASSUMPTIONS.",
      "Approximate Chicago ambulance/EMS transport medical cost. Override with HVI_UNIT_COST_EMS_CALL_USD or HVI_EXCESS_COST_ASSUMPTIONS.",
      "Mortality risk valuation, not a direct medical bill. Override with HVI_UNIT_COST_DEATH_USD or HVI_EXCESS_COST_ASSUMPTIONS."
    ),
    cost_year = as.integer(Sys.getenv("HVI_UNIT_COST_YEAR", unset = "2026")),
    source_note = c(
      "Default placeholder based on public Chicago-area cost ranges; replace with local payer/hospital average if available.",
      "Default placeholder based on public Chicago ambulance cost ranges; replace with local billing average if available.",
      "Default placeholder aligned with public federal mortality-risk valuation practice; not Chicago-specific medical spending."
    )
  )
}

excess_cost_assumptions <- excess_cost_assumptions %>%
  mutate(
    endpoint_key = as.character(endpoint_key),
    event_family = as.character(event_family),
    unit_cost_usd = suppressWarnings(as.numeric(unit_cost_usd)),
    cost_year = suppressWarnings(as.integer(cost_year))
  )

temp_candidates <- c("tmax", "temp_max", "daily_max_temp", "temperature_c", "temperature")
temp_col <- temp_candidates[temp_candidates %in% names(hvi_model_matrix)][1]
if (is.na(temp_col)) {
  stop("Could not find a daily temperature column in hvi_model_matrix.")
}

historical_temperature <- hvi_model_matrix %>%
  transmute(
    community = as.character(community),
    date = as.Date(date),
    year = as.integer(as.character(year)),
    temperature_c = suppressWarnings(as.numeric(.data[[temp_col]])),
    humidity = if ("humidity" %in% names(hvi_model_matrix)) suppressWarnings(as.numeric(humidity)) else NA_real_,
    population = if ("pop_offset" %in% names(hvi_model_matrix)) suppressWarnings(as.numeric(pop_offset)) else NA_real_
  ) %>%
  distinct(community, date, year, .keep_all = TRUE) %>%
  mutate(temperature_f = c_to_f(temperature_c))

canopy_pct_col <- c("tree_canopy_pct", "mean_tree_canopy", "tree_canopy")
canopy_pct_col <- canopy_pct_col[canopy_pct_col %in% names(hvi_model_matrix)][1]
canopy_fraction_col <- c("tree_canopy_fraction")
canopy_fraction_col <- canopy_fraction_col[canopy_fraction_col %in% names(hvi_model_matrix)][1]

historical_canopy <- hvi_model_matrix %>%
  transmute(
    community = as.character(community),
    year = as.integer(as.character(year)),
    tree_canopy_pct = if (!is.na(canopy_pct_col)) {
      suppressWarnings(as.numeric(.data[[canopy_pct_col]]))
    } else if (!is.na(canopy_fraction_col)) {
      100 * suppressWarnings(as.numeric(.data[[canopy_fraction_col]]))
    } else {
      NA_real_
    },
    tree_canopy_fraction = if (!is.na(canopy_fraction_col)) {
      suppressWarnings(as.numeric(.data[[canopy_fraction_col]]))
    } else {
      tree_canopy_pct / 100
    }
  ) %>%
  distinct(community, year, .keep_all = TRUE)

frontend_tree_canopy_year_ca <- historical_canopy %>%
  arrange(community, year)

# -----------------------------
# FRONTEND-FRIENDLY TABLES
# -----------------------------
frontend_structural_summary <- ca_year_overall_structural_hvi %>%
  select(
    community, year,
    overall_structural_0_100,
    overall_structural_raw,
    dominant_endpoint,
    dominant_endpoint_score,
    n_endpoints
  )

frontend_structural_endpoint_long <- ca_year_endpoint_vulnerability %>%
  mutate(
    top_driver_1_label = label_driver(top_driver_1),
    top_driver_2_label = label_driver(top_driver_2),
    top_driver_3_label = label_driver(top_driver_3)
  ) %>%
  select(
    community, year, endpoint_key, outcome_label, source, domain,
    vulnerability_0_100, vulnerability_raw,
    top_driver_1, top_driver_1_label,
    top_driver_2, top_driver_2_label,
    top_driver_3, top_driver_3_label
  )

frontend_baseline_risk_factors <- frontend_structural_summary %>%
  left_join(
    frontend_structural_endpoint_long,
    by = c("community", "year", "dominant_endpoint" = "endpoint_key")
  ) %>%
  left_join(
    historical_canopy,
    by = c("community", "year")
  ) %>%
  transmute(
    community,
    year,
    tree_canopy_pct,
    tree_canopy_fraction,
    baseline_risk_score_0_100 = overall_structural_0_100,
    baseline_risk_tier = impact_tier(overall_structural_0_100),
    baseline_dominant_endpoint = dominant_endpoint,
    baseline_dominant_endpoint_label = outcome_label,
    baseline_dominant_endpoint_source = source,
    baseline_dominant_endpoint_domain = domain,
    baseline_dominant_endpoint_score_0_100 = dominant_endpoint_score,
    baseline_driver_1 = top_driver_1,
    baseline_driver_1_label = top_driver_1_label,
    baseline_driver_2 = top_driver_2,
    baseline_driver_2_label = top_driver_2_label,
    baseline_driver_3 = top_driver_3,
    baseline_driver_3_label = top_driver_3_label,
    baseline_risk_description = paste0(
      "Baseline heat-health vulnerability reflects structural and environmental factors associated with higher modeled heat-sensitive risk. ",
      "For this community-year, the strongest endpoint signal is ",
      coalesce(outcome_label, dominant_endpoint),
      ". Leading modeled contributors include ",
      coalesce(top_driver_1_label, "available vulnerability covariates"),
      if_else(!is.na(top_driver_2_label), paste0(", ", top_driver_2_label), ""),
      if_else(!is.na(top_driver_3_label), paste0(", and ", top_driver_3_label), ""),
      "."
    )
  )

frontend_temp_query_long <- temp_grid_endpoint_risk %>%
  select(
    community, year, temp_value, endpoint_key, outcome_label, source, domain,
    endpoint_risk_0_100, predicted_count, reference_count, excess_events, relative_risk
  )

frontend_temp_query_overall <- temp_grid_overall_risk %>%
  select(
    community, year, temp_value,
    overall_risk_0_100, total_predicted_count, total_reference_count, total_excess_events,
    dominant_endpoint, dominant_endpoint_risk_0_100
  )

frontend_daily_map <- ca_day_overall_operational_hvi %>%
  left_join(historical_temperature, by = c("community", "date", "year")) %>%
  left_join(historical_canopy, by = c("community", "year")) %>%
  select(
    community, date, year,
    temperature_c, temperature_f, humidity,
    tree_canopy_pct, tree_canopy_fraction,
    overall_risk_0_100, alert_tier,
    total_predicted_count, total_reference_count, total_excess_events,
    dominant_endpoint, dominant_endpoint_label, dominant_endpoint_source,
    dominant_endpoint_risk_0_100
  )

frontend_daily_endpoint_long <- ca_day_endpoint_risk %>%
  select(
    community, date, year, endpoint_key, outcome_label, source, domain,
    endpoint_risk_0_100, observed_count, predicted_count, reference_count, excess_events, relative_risk
  )

all_cause_endpoint_keys <- c("ed_visits", "ems_calls", "deaths")

frontend_excess_endpoint_daily <- ca_day_endpoint_risk %>%
  mutate(excess_events_heat = pmax(excess_events, 0)) %>%
  left_join(
    historical_temperature %>%
      select(community, date, year, temperature_c, temperature_f, population),
    by = c("community", "date", "year")
  ) %>%
  left_join(
    historical_canopy,
    by = c("community", "year")
  ) %>%
  left_join(excess_cost_assumptions, by = "endpoint_key") %>%
  mutate(
    event_family = source,
    excess_rate_per_100k = if_else(
      !is.na(population) & population > 0,
      100000 * excess_events_heat / population,
      NA_real_
    )
  ) %>%
  select(
    community, date, year,
    endpoint_key, event_family, outcome_label, domain,
    temperature_c, temperature_f,
    tree_canopy_pct, tree_canopy_fraction,
    population,
    excess_events_heat, excess_rate_per_100k,
    endpoint_risk_0_100, relative_risk
  )

frontend_excess_allcause_daily <- frontend_excess_endpoint_daily %>%
  filter(endpoint_key %in% all_cause_endpoint_keys) %>%
  left_join(excess_cost_assumptions, by = "endpoint_key") %>%
  mutate(
    event_family = coalesce(event_family.y, event_family.x),
    estimated_cost_usd = excess_events_heat * unit_cost_usd
  ) %>%
  select(
    community, date, year,
    endpoint_key, event_family, outcome_label,
    temperature_c, temperature_f,
    tree_canopy_pct, tree_canopy_fraction,
    population,
    excess_events_heat, excess_rate_per_100k,
    unit_cost_usd, estimated_cost_usd,
    cost_year, cost_basis
  )

frontend_excess_endpoint_annual <- frontend_excess_endpoint_daily %>%
  group_by(community, year, endpoint_key, event_family, outcome_label, domain) %>%
  summarise(
    warm_season_days = n_distinct(date),
    mean_temperature_f = mean(temperature_f, na.rm = TRUE),
    max_temperature_f = max(temperature_f, na.rm = TRUE),
    population = mean(population, na.rm = TRUE),
    tree_canopy_pct = mean(tree_canopy_pct, na.rm = TRUE),
    tree_canopy_fraction = mean(tree_canopy_fraction, na.rm = TRUE),
    excess_events_heat = sum(excess_events_heat, na.rm = TRUE),
    excess_rate_per_100k = if_else(
      !is.na(population) & population > 0,
      100000 * excess_events_heat / population,
      NA_real_
    ),
    .groups = "drop"
  )

frontend_excess_allcause_annual <- frontend_excess_allcause_daily %>%
  group_by(community, year, endpoint_key, event_family, outcome_label) %>%
  summarise(
    warm_season_days = n_distinct(date),
    mean_temperature_f = mean(temperature_f, na.rm = TRUE),
    max_temperature_f = max(temperature_f, na.rm = TRUE),
    population = mean(population, na.rm = TRUE),
    tree_canopy_pct = mean(tree_canopy_pct, na.rm = TRUE),
    tree_canopy_fraction = mean(tree_canopy_fraction, na.rm = TRUE),
    excess_events_heat = sum(excess_events_heat, na.rm = TRUE),
    excess_rate_per_100k = if_else(
      !is.na(population) & population > 0,
      100000 * excess_events_heat / population,
      NA_real_
    ),
    unit_cost_usd = first(unit_cost_usd),
    estimated_cost_usd = sum(estimated_cost_usd, na.rm = TRUE),
    cost_year = first(cost_year),
    cost_basis = first(cost_basis),
    .groups = "drop"
  )

city_population_by_year <- historical_temperature %>%
  distinct(community, year, .keep_all = TRUE) %>%
  group_by(year) %>%
  summarise(population = sum(population, na.rm = TRUE), .groups = "drop")

reference_city_population <- city_population_by_year %>%
  filter(year == max(year, na.rm = TRUE)) %>%
  pull(population)
if (length(reference_city_population) == 0) reference_city_population <- NA_real_

frontend_excess_endpoint_citywide <- frontend_excess_endpoint_daily %>%
  group_by(year, endpoint_key, event_family, outcome_label, domain) %>%
  summarise(
    community_areas = n_distinct(community),
    warm_season_days = n_distinct(date),
    mean_temperature_f = mean(temperature_f, na.rm = TRUE),
    max_temperature_f = max(temperature_f, na.rm = TRUE),
    excess_events_heat = sum(excess_events_heat, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(city_population_by_year, by = "year") %>%
  mutate(
    excess_rate_per_100k = if_else(
      !is.na(population) & population > 0,
      100000 * excess_events_heat / population,
      NA_real_
    )
  ) %>%
  bind_rows(
    frontend_excess_endpoint_daily %>%
      group_by(endpoint_key, event_family, outcome_label, domain) %>%
      summarise(
        year = NA_integer_,
        community_areas = n_distinct(community),
        warm_season_days = n_distinct(date),
        mean_temperature_f = mean(temperature_f, na.rm = TRUE),
        max_temperature_f = max(temperature_f, na.rm = TRUE),
        population = reference_city_population,
        excess_events_heat = sum(excess_events_heat, na.rm = TRUE),
        excess_rate_per_100k = if_else(
          !is.na(population) & population > 0,
          100000 * excess_events_heat / population,
          NA_real_
        ),
        .groups = "drop"
      )
  ) %>%
  arrange(event_family, endpoint_key, year)

frontend_excess_allcause_citywide <- frontend_excess_allcause_daily %>%
  group_by(year, endpoint_key, event_family, outcome_label) %>%
  summarise(
    community_areas = n_distinct(community),
    warm_season_days = n_distinct(date),
    mean_temperature_f = mean(temperature_f, na.rm = TRUE),
    max_temperature_f = max(temperature_f, na.rm = TRUE),
    excess_events_heat = sum(excess_events_heat, na.rm = TRUE),
    unit_cost_usd = first(unit_cost_usd),
    estimated_cost_usd = sum(estimated_cost_usd, na.rm = TRUE),
    cost_year = first(cost_year),
    cost_basis = first(cost_basis),
    .groups = "drop"
  ) %>%
  left_join(city_population_by_year, by = "year") %>%
  mutate(
    excess_rate_per_100k = if_else(
      !is.na(population) & population > 0,
      100000 * excess_events_heat / population,
      NA_real_
    )
  ) %>%
  bind_rows(
    frontend_excess_allcause_daily %>%
      group_by(endpoint_key, event_family, outcome_label) %>%
      summarise(
        year = NA_integer_,
        community_areas = n_distinct(community),
        warm_season_days = n_distinct(date),
        mean_temperature_f = mean(temperature_f, na.rm = TRUE),
        max_temperature_f = max(temperature_f, na.rm = TRUE),
        population = reference_city_population,
        excess_events_heat = sum(excess_events_heat, na.rm = TRUE),
        excess_rate_per_100k = if_else(
          !is.na(population) & population > 0,
          100000 * excess_events_heat / population,
          NA_real_
        ),
        unit_cost_usd = first(unit_cost_usd),
        estimated_cost_usd = sum(estimated_cost_usd, na.rm = TRUE),
        cost_year = first(cost_year),
        cost_basis = first(cost_basis),
        .groups = "drop"
      )
  ) %>%
  arrange(endpoint_key, year)

frontend_historical_daily_summary <- frontend_daily_map %>%
  transmute(
    community, date, year,
    temperature_c, temperature_f, humidity,
    tree_canopy_pct, tree_canopy_fraction,
    daily_risk_score_0_100 = overall_risk_0_100,
    daily_risk_tier = alert_tier,
    daily_dominant_endpoint = dominant_endpoint,
    daily_dominant_endpoint_label = dominant_endpoint_label,
    daily_dominant_endpoint_source = dominant_endpoint_source,
    daily_dominant_endpoint_risk_0_100 = dominant_endpoint_risk_0_100
  )

frontend_historical_daily_endpoint_ranked <- ca_day_endpoint_risk %>%
  group_by(community, date, year) %>%
  mutate(
    positive_excess = pmax(excess_events, 0),
    positive_excess_total = sum(positive_excess, na.rm = TRUE),
    impact_share_pct = if_else(
      positive_excess_total > 0,
      100 * positive_excess / positive_excess_total,
      NA_real_
    ),
    impact_rank = min_rank(desc(excess_events)),
    impact_tier = impact_tier(endpoint_risk_0_100)
  ) %>%
  ungroup() %>%
  arrange(community, date, impact_rank, endpoint_key) %>%
  select(
    community, date, year,
    endpoint_key, outcome_label, source, domain,
    impact_rank, endpoint_risk_0_100, impact_tier,
    relative_risk, impact_share_pct
  )

# optional wide versions for easier JS consumption
frontend_temp_query_wide <- frontend_temp_query_long %>%
  select(community, year, temp_value, endpoint_key, endpoint_risk_0_100) %>%
  mutate(metric_name = paste0("risk_", endpoint_key)) %>%
  select(-endpoint_key) %>%
  pivot_wider(names_from = metric_name, values_from = endpoint_risk_0_100)

frontend_daily_wide <- frontend_daily_endpoint_long %>%
  select(community, date, endpoint_key, endpoint_risk_0_100) %>%
  mutate(metric_name = paste0("risk_", endpoint_key)) %>%
  select(-endpoint_key) %>%
  pivot_wider(names_from = metric_name, values_from = endpoint_risk_0_100)

# -----------------------------
# EXPORT
# -----------------------------
write_csv(frontend_structural_summary, file.path(frontend_dir, "frontend_structural_summary.csv"))
write_csv(frontend_structural_endpoint_long, file.path(frontend_dir, "frontend_structural_endpoint_long.csv"))
write_csv(frontend_baseline_risk_factors, file.path(frontend_dir, "frontend_baseline_risk_factors.csv"))
write_csv(frontend_tree_canopy_year_ca, file.path(frontend_dir, "frontend_tree_canopy_year_ca.csv"))
write_csv(frontend_temp_query_long, file.path(frontend_dir, "frontend_temperature_query_endpoint_long.csv"))
write_csv(frontend_temp_query_overall, file.path(frontend_dir, "frontend_temperature_query_overall.csv"))
write_csv(frontend_temp_query_wide, file.path(frontend_dir, "frontend_temperature_query_wide.csv"))
write_csv(frontend_daily_map, file.path(frontend_dir, "frontend_daily_map.csv"))
write_csv(frontend_daily_endpoint_long, file.path(frontend_dir, "frontend_daily_endpoint_long.csv"))
write_csv(frontend_daily_wide, file.path(frontend_dir, "frontend_daily_wide.csv"))
write_csv(frontend_historical_daily_summary, file.path(frontend_dir, "frontend_historical_daily_summary.csv"))
write_csv(frontend_historical_daily_endpoint_ranked, file.path(frontend_dir, "frontend_historical_daily_endpoint_ranked.csv"))
write_csv(frontend_excess_endpoint_daily, file.path(frontend_dir, "frontend_excess_endpoint_daily.csv"))
write_csv(frontend_excess_endpoint_annual, file.path(frontend_dir, "frontend_excess_endpoint_annual.csv"))
write_csv(frontend_excess_endpoint_citywide, file.path(frontend_dir, "frontend_excess_endpoint_citywide.csv"))
write_csv(frontend_excess_allcause_daily, file.path(frontend_dir, "frontend_excess_allcause_daily.csv"))
write_csv(frontend_excess_allcause_annual, file.path(frontend_dir, "frontend_excess_allcause_annual.csv"))
write_csv(frontend_excess_allcause_citywide, file.path(frontend_dir, "frontend_excess_allcause_citywide.csv"))
write_csv(excess_cost_assumptions, file.path(frontend_dir, "frontend_excess_cost_assumptions.csv"))

assign("frontend_structural_summary", frontend_structural_summary, envir = .GlobalEnv)
assign("frontend_baseline_risk_factors", frontend_baseline_risk_factors, envir = .GlobalEnv)
assign("frontend_tree_canopy_year_ca", frontend_tree_canopy_year_ca, envir = .GlobalEnv)
assign("frontend_temp_query_long", frontend_temp_query_long, envir = .GlobalEnv)
assign("frontend_temp_query_overall", frontend_temp_query_overall, envir = .GlobalEnv)
assign("frontend_daily_map", frontend_daily_map, envir = .GlobalEnv)
assign("frontend_daily_endpoint_long", frontend_daily_endpoint_long, envir = .GlobalEnv)
assign("frontend_historical_daily_summary", frontend_historical_daily_summary, envir = .GlobalEnv)
assign("frontend_historical_daily_endpoint_ranked", frontend_historical_daily_endpoint_ranked, envir = .GlobalEnv)
assign("frontend_excess_endpoint_daily", frontend_excess_endpoint_daily, envir = .GlobalEnv)
assign("frontend_excess_endpoint_annual", frontend_excess_endpoint_annual, envir = .GlobalEnv)
assign("frontend_excess_endpoint_citywide", frontend_excess_endpoint_citywide, envir = .GlobalEnv)
assign("frontend_excess_allcause_daily", frontend_excess_allcause_daily, envir = .GlobalEnv)
assign("frontend_excess_allcause_annual", frontend_excess_allcause_annual, envir = .GlobalEnv)
assign("frontend_excess_allcause_citywide", frontend_excess_allcause_citywide, envir = .GlobalEnv)
assign("excess_cost_assumptions", excess_cost_assumptions, envir = .GlobalEnv)

message("09e complete. Frontend exports written to: ", frontend_dir)
