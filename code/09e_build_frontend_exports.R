# ================================================================================================
# HVI 2.0 | 09e_build_frontend_exports.R
# Prepare lightweight long and wide exports for a future interactive map/dashboard.
# ================================================================================================

source("09_utils_hvi.R")

# -----------------------------
# CONFIG
# -----------------------------
project_dir <- getwd()
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
  "ca_day_endpoint_risk"
)

missing_objs <- required_objs[!vapply(required_objs, exists, logical(1), envir = .GlobalEnv)]
if (length(missing_objs) > 0) stop("Run 09b-09d first. Missing objects: ", paste(missing_objs, collapse = ", "))

ca_year_overall_structural_hvi <- get("ca_year_overall_structural_hvi", envir = .GlobalEnv)
ca_year_endpoint_vulnerability <- get("ca_year_endpoint_vulnerability", envir = .GlobalEnv)
temp_grid_overall_risk <- get("temp_grid_overall_risk", envir = .GlobalEnv)
temp_grid_endpoint_risk <- get("temp_grid_endpoint_risk", envir = .GlobalEnv)
ca_day_overall_operational_hvi <- get("ca_day_overall_operational_hvi", envir = .GlobalEnv)
ca_day_endpoint_risk <- get("ca_day_endpoint_risk", envir = .GlobalEnv)

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
  select(
    community, year, endpoint_key, outcome_label, source, domain,
    vulnerability_0_100, vulnerability_raw,
    top_driver_1, top_driver_2, top_driver_3
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
  select(
    community, date, year,
    overall_risk_0_100, alert_tier,
    total_predicted_count, total_reference_count, total_excess_events,
    dominant_endpoint, dominant_endpoint_label, dominant_endpoint_source
  )

frontend_daily_endpoint_long <- ca_day_endpoint_risk %>%
  select(
    community, date, year, endpoint_key, outcome_label, source, domain,
    endpoint_risk_0_100, observed_count, predicted_count, reference_count, excess_events, relative_risk
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
write_csv(frontend_temp_query_long, file.path(frontend_dir, "frontend_temperature_query_endpoint_long.csv"))
write_csv(frontend_temp_query_overall, file.path(frontend_dir, "frontend_temperature_query_overall.csv"))
write_csv(frontend_temp_query_wide, file.path(frontend_dir, "frontend_temperature_query_wide.csv"))
write_csv(frontend_daily_map, file.path(frontend_dir, "frontend_daily_map.csv"))
write_csv(frontend_daily_endpoint_long, file.path(frontend_dir, "frontend_daily_endpoint_long.csv"))
write_csv(frontend_daily_wide, file.path(frontend_dir, "frontend_daily_wide.csv"))

assign("frontend_structural_summary", frontend_structural_summary, envir = .GlobalEnv)
assign("frontend_temp_query_long", frontend_temp_query_long, envir = .GlobalEnv)
assign("frontend_temp_query_overall", frontend_temp_query_overall, envir = .GlobalEnv)
assign("frontend_daily_map", frontend_daily_map, envir = .GlobalEnv)
assign("frontend_daily_endpoint_long", frontend_daily_endpoint_long, envir = .GlobalEnv)

message("09e complete. Frontend exports written to: ", frontend_dir)
