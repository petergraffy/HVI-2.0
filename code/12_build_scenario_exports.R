# ================================================================================================
# HVI 2.0 | Scenario slider exports
# Backend artifact for public dashboard controls such as temperature and NDVI sliders.
#
# This script precomputes a public-safe scenario grid. It intentionally uses base R plus mgcv
# so the export can run in lightweight automation environments.
#
# Temperature changes are scored through endpoint-specific MRT/lag heat-dose terms.
# NDVI changes are association-based perturbations of z_ndvi and should not be interpreted
# as a causal intervention without a separate causal design.
# ================================================================================================

source(file.path(if (dir.exists("code")) "code" else ".", "00_config.R"))

rescale_0_100_local <- function(x) {
  rng <- range(x, na.rm = TRUE)
  if (anyNA(rng) || diff(rng) == 0) return(rep(50, length(x)))
  100 * (x - rng[1]) / diff(rng)
}

read_csv_base <- function(path) {
  utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
}

clean_names_base <- function(x) {
  x <- tolower(x)
  x <- gsub("[^a-z0-9]+", "_", x)
  x <- gsub("^_|_$", "", x)
  x
}

coalesce_num <- function(x, fallback) {
  x[is.na(x)] <- fallback
  x
}

# ------------------------------------------------------------------------------------------------
# Config
# ------------------------------------------------------------------------------------------------
model_output_dirs <- c(
  HVI_PATHS$private_outputs$model_outputs,
  file.path(HVI_PATHS$repo, "code", "09_model_outputs")
)
model_output_dirs <- model_output_dirs[dir.exists(model_output_dirs)]
if (length(model_output_dirs) == 0) stop("No model output directory found.")
required_model_files <- c("09a_endpoint_models.rds", "09b_endpoint_weights.csv")
complete_dir <- vapply(
  model_output_dirs,
  function(path) all(file.exists(file.path(path, required_model_files))),
  logical(1)
)
if (!any(complete_dir)) {
  stop("No model output directory contains required files: ", paste(required_model_files, collapse = ", "))
}
model_output_dir <- model_output_dirs[which(complete_dir)[1]]

scenario_out_dir <- file.path(model_output_dir, "scenario_exports")
hvi_dir_create(scenario_out_dir)

scenario_year <- as.integer(Sys.getenv("HVI_SCENARIO_YEAR", unset = "2022"))
template_doy <- as.integer(Sys.getenv("HVI_SCENARIO_DOY", unset = "200"))
template_dow <- Sys.getenv("HVI_SCENARIO_DOW", unset = "Wed")

temperature_grid_f <- seq(75, 105, by = 1)
ndvi_delta_grid <- c(-0.10, -0.05, 0, 0.05, 0.10, 0.15, 0.20)
humidity_grid <- as.numeric(strsplit(Sys.getenv("HVI_SCENARIO_HUMIDITY_GRID", unset = "50"), ",")[[1]])

scenario_label <- "temperature_ndvi_slider_grid"

# ------------------------------------------------------------------------------------------------
# Load artifacts
# ------------------------------------------------------------------------------------------------
model_rds <- file.path(model_output_dir, "09a_endpoint_models.rds")
if (!file.exists(model_rds)) stop("Missing fitted model artifact: ", model_rds)
endpoint_models <- readRDS(model_rds)

matrix_candidates <- c(
  file.path(HVI_PATHS$private, "hvi_model_matrix_2019_2022.csv"),
  file.path(HVI_PATHS$repo, "hvi_model_matrix_2019_2022.csv")
)
matrix_path <- matrix_candidates[file.exists(matrix_candidates)][1]
if (is.na(matrix_path)) stop("Missing hvi_model_matrix_2019_2022.csv.")
hvi_model_matrix <- read_csv_base(matrix_path)
names(hvi_model_matrix) <- clean_names_base(names(hvi_model_matrix))

meta_candidates <- c(
  file.path(HVI_PATHS$private, "hvi_endpoint_metadata.csv"),
  file.path(model_output_dir, "hvi_endpoint_metadata.csv")
)
meta_path <- meta_candidates[file.exists(meta_candidates)][1]
if (is.na(meta_path)) stop("Missing hvi_endpoint_metadata.csv.")
hvi_endpoint_metadata <- read_csv_base(meta_path)
names(hvi_endpoint_metadata) <- clean_names_base(names(hvi_endpoint_metadata))

weights_path <- file.path(model_output_dir, "09b_endpoint_weights.csv")
if (!file.exists(weights_path)) stop("Missing endpoint weights artifact: ", weights_path)
endpoint_weights <- read_csv_base(weights_path)
names(endpoint_weights) <- clean_names_base(names(endpoint_weights))

# ------------------------------------------------------------------------------------------------
# Model support checks and base table
# ------------------------------------------------------------------------------------------------
model_terms <- unique(unlist(lapply(endpoint_models, function(fit) all.vars(stats::formula(fit)))))
humidity_active <- "humidity" %in% model_terms
ndvi_active <- any(vapply(endpoint_models, function(fit) "z_ndvi" %in% all.vars(stats::formula(fit)), logical(1)))

if (!"ndvi" %in% names(hvi_model_matrix) || !"z_ndvi" %in% names(hvi_model_matrix)) {
  stop("Scenario builder requires ndvi and z_ndvi in hvi_model_matrix.")
}

raw_ndvi_sd <- stats::sd(hvi_model_matrix$ndvi, na.rm = TRUE)
if (!is.finite(raw_ndvi_sd) || raw_ndvi_sd <= 0) {
  stop("Cannot build NDVI slider grid because raw ndvi has zero or missing variance.")
}

hvi_model_matrix$year_int <- suppressWarnings(as.integer(as.character(hvi_model_matrix$year)))
base_year <- hvi_model_matrix[hvi_model_matrix$year_int == scenario_year, , drop = FALSE]
if (nrow(base_year) == 0) {
  scenario_year <- max(hvi_model_matrix$year_int, na.rm = TRUE)
  base_year <- hvi_model_matrix[hvi_model_matrix$year_int == scenario_year, , drop = FALSE]
}

z_vars <- grep("^z_", names(base_year), value = TRUE)
summary_vars <- unique(c("pop_offset", "ndvi", "z_ndvi", "humidity", z_vars))
summary_vars <- summary_vars[summary_vars %in% names(base_year)]

community_keys <- unique(base_year[, c("community", "year_int"), drop = FALSE])
community_year_base <- community_keys
for (v in summary_vars) {
  community_year_base[[v]] <- ave(base_year[[v]], base_year$community, base_year$year_int, FUN = function(x) mean(x, na.rm = TRUE))[match(
    paste(community_year_base$community, community_year_base$year_int),
    paste(base_year$community, base_year$year_int)
  )]
}
community_year_base <- unique(community_year_base)
names(community_year_base)[names(community_year_base) == "year_int"] <- "year"
community_year_base$doy <- template_doy
community_year_base$dow <- factor(template_dow, levels = c("Sun", "Mon", "Tue", "Wed", "Thu", "Fri", "Sat"))

# ------------------------------------------------------------------------------------------------
# Score scenario grid
# ------------------------------------------------------------------------------------------------
endpoint_scenario_list <- list()

for (ep_key in names(endpoint_models)) {
  fit <- endpoint_models[[ep_key]]
  meta_row <- hvi_endpoint_metadata[hvi_endpoint_metadata$endpoint_key == ep_key, , drop = FALSE]
  if (nrow(meta_row) == 0) next

  all_model_vars <- all.vars(stats::formula(fit))
  z_vars_use <- grep("^z_", all_model_vars, value = TRUE)
  z_vars_use <- setdiff(z_vars_use, c("outcome", "heat_dose", "pop_offset", "doy", "year"))

  mrt_val <- suppressWarnings(as.numeric(meta_row$mrt[1]))
  max_lag_val <- suppressWarnings(as.numeric(meta_row$max_lag[1]))
  if (is.na(max_lag_val)) max_lag_val <- 0
  outcome_col <- meta_row$panel_outcome_col[1]
  observed_max_count <- if (outcome_col %in% names(hvi_model_matrix)) {
    max(hvi_model_matrix[[outcome_col]], na.rm = TRUE)
  } else {
    NA_real_
  }
  display_count_cap <- if (is.finite(observed_max_count)) max(observed_max_count * 2, 1) else Inf

  base_cols <- unique(c("community", "year", z_vars_use, "pop_offset", "ndvi", "z_ndvi", "humidity", "doy", "dow"))
  base_cols <- base_cols[base_cols %in% names(community_year_base)]
  base_ep <- community_year_base[, base_cols, drop = FALSE]

  scenario_index <- expand.grid(
    row_id = seq_len(nrow(base_ep)),
    temperature_f = temperature_grid_f,
    ndvi_delta = ndvi_delta_grid,
    humidity_scenario = humidity_grid,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  grid_ep <- cbind(base_ep[scenario_index$row_id, , drop = FALSE], scenario_index[, -1, drop = FALSE])

  grid_ep$endpoint_key <- ep_key
  grid_ep$outcome_label <- meta_row$outcome_label[1]
  grid_ep$source <- meta_row$source[1]
  grid_ep$domain <- meta_row$domain[1]
  grid_ep$scenario <- scenario_label
  grid_ep$mrt <- mrt_val
  grid_ep$max_lag <- max_lag_val
  grid_ep$temperature_c <- (grid_ep$temperature_f - 32) * 5 / 9
  grid_ep$temp_value <- grid_ep$temperature_c
  grid_ep$temp_value_model_units <- grid_ep$temperature_c
  grid_ep$excess_above_mrt <- pmax(grid_ep$temp_value_model_units - mrt_val, 0)
  grid_ep$heat_dose_uncapped <- grid_ep$excess_above_mrt * (max_lag_val + 1)
  heat_dose_col <- paste0("heat_dose__", ep_key)
  max_supported_heat_dose <- if (heat_dose_col %in% names(hvi_model_matrix)) {
    max(hvi_model_matrix[[heat_dose_col]], na.rm = TRUE)
  } else {
    max(grid_ep$heat_dose_uncapped, na.rm = TRUE)
  }
  grid_ep$heat_dose_capped <- is.finite(max_supported_heat_dose) & grid_ep$heat_dose_uncapped > max_supported_heat_dose
  grid_ep$max_supported_heat_dose <- max_supported_heat_dose
  grid_ep$heat_dose <- pmin(grid_ep$heat_dose_uncapped, max_supported_heat_dose)
  grid_ep$ndvi_baseline <- grid_ep$ndvi
  grid_ep$ndvi_scenario <- pmin(pmax(grid_ep$ndvi + grid_ep$ndvi_delta, -1), 1)
  grid_ep$z_ndvi <- grid_ep$z_ndvi + (grid_ep$ndvi_delta / raw_ndvi_sd)
  if (humidity_active) grid_ep$humidity <- grid_ep$humidity_scenario

  if (!is.null(fit$xlevels) && "dow" %in% names(fit$xlevels)) {
    grid_ep$dow <- factor(as.character(grid_ep$dow), levels = fit$xlevels[["dow"]])
  }
  if (!is.null(fit$xlevels) && "year" %in% names(fit$xlevels)) {
    grid_ep$year <- factor(as.character(grid_ep$year), levels = fit$xlevels[["year"]])
  }

  ref_ep <- grid_ep
  ref_ep$heat_dose <- 0

  vars_needed <- intersect(all_model_vars, names(grid_ep))
  keep <- stats::complete.cases(grid_ep[, vars_needed, drop = FALSE])
  grid_ep <- grid_ep[keep, , drop = FALSE]
  ref_ep <- ref_ep[keep, , drop = FALSE]
  if (nrow(grid_ep) == 0 || nrow(ref_ep) == 0) next

  raw_predicted_count <- as.numeric(predict(fit, newdata = grid_ep, type = "response"))
  ref_ep$reference_count <- as.numeric(predict(fit, newdata = ref_ep, type = "response"))
  grid_ep$prediction_capped_for_display <- is.finite(display_count_cap) & raw_predicted_count > display_count_cap
  grid_ep$predicted_count <- pmin(raw_predicted_count, display_count_cap)
  grid_ep$reference_count <- pmin(ref_ep$reference_count, display_count_cap)
  grid_ep$excess_events <- grid_ep$predicted_count - grid_ep$reference_count
  grid_ep$relative_risk <- ifelse(grid_ep$reference_count > 0, grid_ep$predicted_count / grid_ep$reference_count, NA_real_)
  grid_ep$observed_max_count <- observed_max_count
  grid_ep$display_count_cap <- display_count_cap

  keep_cols <- c(
    "community", "year", "scenario", "temperature_f", "temperature_c", "humidity_scenario", "ndvi_delta",
    "ndvi_baseline", "ndvi_scenario", "endpoint_key", "outcome_label", "source", "domain",
    "heat_dose_uncapped", "heat_dose", "heat_dose_capped", "max_supported_heat_dose",
    "prediction_capped_for_display", "observed_max_count", "display_count_cap",
    "predicted_count", "reference_count", "excess_events", "relative_risk"
  )
  endpoint_scenario_list[[ep_key]] <- grid_ep[, keep_cols, drop = FALSE]
}

scenario_endpoint <- do.call(rbind, endpoint_scenario_list)
scenario_endpoint$endpoint_risk_0_100 <- unlist(tapply(
  scenario_endpoint$excess_events,
  scenario_endpoint$endpoint_key,
  rescale_0_100_local
), use.names = FALSE)

weight_cols <- c("endpoint_key", "endpoint_weight", "source_weight", "performance_weight")
weight_cols <- weight_cols[weight_cols %in% names(endpoint_weights)]
scenario_endpoint <- merge(scenario_endpoint, endpoint_weights[, weight_cols, drop = FALSE], by = "endpoint_key", all.x = TRUE)
scenario_endpoint$endpoint_weight <- coalesce_num(scenario_endpoint$endpoint_weight, 1)

group_cols <- c("community", "year", "scenario", "temperature_f", "humidity_scenario", "ndvi_delta")
group_key <- do.call(paste, c(scenario_endpoint[group_cols], sep = "\r"))
overall_split <- split(seq_len(nrow(scenario_endpoint)), group_key)

scenario_overall <- do.call(rbind, lapply(overall_split, function(idx) {
  dat <- scenario_endpoint[idx, , drop = FALSE]
  w <- dat$endpoint_weight
  data.frame(
    community = dat$community[1],
    year = dat$year[1],
    scenario = dat$scenario[1],
    temperature_f = dat$temperature_f[1],
    humidity_scenario = dat$humidity_scenario[1],
    ndvi_delta = dat$ndvi_delta[1],
    total_predicted_count = sum(dat$predicted_count, na.rm = TRUE),
    total_reference_count = sum(dat$reference_count, na.rm = TRUE),
    total_excess_events = sum(dat$excess_events, na.rm = TRUE),
    overall_weighted_excess = stats::weighted.mean(dat$excess_events, w = w, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}))
scenario_overall$overall_risk_0_100 <- rescale_0_100_local(scenario_overall$overall_weighted_excess)

dominant <- do.call(rbind, lapply(overall_split, function(idx) {
  dat <- scenario_endpoint[idx, , drop = FALSE]
  hit <- which.max(dat$excess_events)
  data.frame(
    community = dat$community[hit],
    year = dat$year[hit],
    scenario = dat$scenario[hit],
    temperature_f = dat$temperature_f[hit],
    humidity_scenario = dat$humidity_scenario[hit],
    ndvi_delta = dat$ndvi_delta[hit],
    dominant_endpoint = dat$endpoint_key[hit],
    dominant_endpoint_label = dat$outcome_label[hit],
    dominant_endpoint_source = dat$source[hit],
    dominant_endpoint_excess = dat$excess_events[hit],
    dominant_endpoint_risk_0_100 = dat$endpoint_risk_0_100[hit],
    stringsAsFactors = FALSE
  )
}))
scenario_overall <- merge(scenario_overall, dominant, by = group_cols, all.x = TRUE)

scenario_variable_metadata <- data.frame(
  variable = c("temperature_f", "humidity_scenario", "ndvi_delta"),
  label = c("Temperature", "Relative humidity", "NDVI change"),
  min = c(min(temperature_grid_f), min(humidity_grid), min(ndvi_delta_grid)),
  max = c(max(temperature_grid_f), max(humidity_grid), max(ndvi_delta_grid)),
  step = c(1, NA, 0.05),
  units = c("deg F", "percent", "index delta"),
  mode = c("absolute", "absolute", "delta"),
  active_in_current_model = c(TRUE, humidity_active, ndvi_active),
  interpretation = c(
    "Heat exposure scenario displayed in Fahrenheit and converted to Celsius for MRT-based model scoring.",
    "Only changes predictions when endpoint models include humidity terms.",
    "Association-based perturbation of neighborhood greenness; not a causal intervention estimate."
  ),
  stringsAsFactors = FALSE
)

scenario_baseline_values <- community_year_base[, c("community", "year", "ndvi", "humidity", "pop_offset"), drop = FALSE]
names(scenario_baseline_values) <- c("community", "year", "ndvi_baseline", "humidity_baseline", "pop_offset")

# ------------------------------------------------------------------------------------------------
# Export private model outputs and public dashboard copies
# ------------------------------------------------------------------------------------------------
utils::write.csv(scenario_endpoint, file.path(scenario_out_dir, "scenario_grid_endpoint.csv"), row.names = FALSE, na = "")
utils::write.csv(scenario_overall, file.path(scenario_out_dir, "scenario_grid_overall.csv"), row.names = FALSE, na = "")
utils::write.csv(scenario_variable_metadata, file.path(scenario_out_dir, "scenario_variable_metadata.csv"), row.names = FALSE, na = "")
utils::write.csv(scenario_baseline_values, file.path(scenario_out_dir, "scenario_baseline_values.csv"), row.names = FALSE, na = "")

public_scenario_dir <- file.path(HVI_PATHS$public$dashboard, "scenarios")
hvi_dir_create(public_scenario_dir)
hvi_write_public_csv(scenario_endpoint, file.path(public_scenario_dir, "scenario_grid_endpoint.csv"))
hvi_write_public_csv(scenario_overall, file.path(public_scenario_dir, "scenario_grid_overall.csv"))
hvi_write_public_csv(scenario_variable_metadata, file.path(public_scenario_dir, "scenario_variable_metadata.csv"))
hvi_write_public_csv(scenario_baseline_values, file.path(public_scenario_dir, "scenario_baseline_values.csv"))

assign("scenario_endpoint", scenario_endpoint, envir = .GlobalEnv)
assign("scenario_overall", scenario_overall, envir = .GlobalEnv)
assign("scenario_variable_metadata", scenario_variable_metadata, envir = .GlobalEnv)
assign("scenario_baseline_values", scenario_baseline_values, envir = .GlobalEnv)

message("Scenario exports written to: ", scenario_out_dir)
message("Public scenario exports written to: ", public_scenario_dir)
