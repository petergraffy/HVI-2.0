# ================================================================================================
# HVI 2.0 | 10_publication_output.R
# Publication-ready figures, maps, and tables for HVI 2.0
#
# Expected objects in memory (as available):
#   - hvi_model_matrix
#   - hvi_endpoint_metadata
#   - endpoint_model_performance
#   - endpoint_weights
#   - endpoint_interaction_betas
#   - ca_year_endpoint_vulnerability
#   - ca_year_overall_structural_hvi
#   - ca_day_endpoint_risk
#   - ca_day_overall_operational_hvi
#   - temp_grid_endpoint_risk   OR temperature_grid_endpoint_risk
#   - temp_grid_overall_risk    OR temperature_grid_overall_risk
#   - temp_grid_family_risk     OR temperature_grid_family_risk
#   - cas OR commareas OR other sf object with Chicago community-area geometry
#
# Outputs:
#   {project_dir}/publication_outputs/
#       figures/
#       tables/
#       supplements/
# ================================================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(sf)
  library(scales)
  library(glue)
  library(janitor)
  library(patchwork)
  library(gt)
  library(viridis)
  library(lubridate)
  library(forcats)
  library(stringr)
  library(ggrepel)
})

source(file.path(if (dir.exists("code")) "code" else ".", "00_config.R"))

# ------------------------------------------------------------------------------------------------
# CONFIG
# ------------------------------------------------------------------------------------------------
project_dir <- HVI_PATHS$private

pub_dir <- HVI_PATHS$private_outputs$publication_outputs
fig_dir <- file.path(pub_dir, "figures")
tab_dir <- file.path(pub_dir, "tables")
sup_dir <- file.path(pub_dir, "supplements")

dir.create(pub_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(sup_dir, recursive = TRUE, showWarnings = FALSE)

target_year <- 2022

# publication temperatures in Fahrenheit; will be converted if grid is Celsius
target_temps_f <- c(85, 91, 97)

# choose a small endpoint panel for faceted maps if needed
selected_endpoint_maps <- c(
  "deaths", "death_cvd", "death_renal", "death_injury",
  "ed_visits", "ed_cvd", "ed_respiratory", "ed_injury",
  "ems_calls", "ems_cvd", "ems_respiratory", "ems_injury"
)

base_size <- 11

# ------------------------------------------------------------------------------------------------
# THEME
# ------------------------------------------------------------------------------------------------
theme_pub <- function(base_size = 11) {
  theme_minimal(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      axis.title = element_text(face = "bold"),
      plot.title = element_text(face = "bold", size = base_size + 2),
      plot.subtitle = element_text(size = base_size, color = "grey25"),
      plot.caption = element_text(size = base_size - 1, color = "grey35"),
      strip.text = element_text(face = "bold"),
      legend.title = element_text(face = "bold"),
      legend.position = "right"
    )
}

# ------------------------------------------------------------------------------------------------
# HELPERS
# ------------------------------------------------------------------------------------------------
save_pub_plot <- function(plot_obj, filename, width = 9, height = 7, dpi = 400, bg = "white") {
  ggsave(
    filename = file.path(fig_dir, filename),
    plot = plot_obj,
    width = width,
    height = height,
    dpi = dpi,
    bg = bg
  )
}

safe_write_csv <- function(df, path) {
  readr::write_csv(df, path, na = "")
}

rescale_0_100 <- function(x) {
  if (all(is.na(x))) return(rep(NA_real_, length(x)))
  rng <- range(x, na.rm = TRUE)
  if (!all(is.finite(rng)) || diff(rng) == 0) return(rep(50, length(x)))
  scales::rescale(x, to = c(0, 100), from = rng)
}

f_to_c <- function(x) (x - 32) * 5 / 9
c_to_f <- function(x) x * 9 / 5 + 32

detect_temperature_units <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return("unknown")
  if (max(x, na.rm = TRUE) > 60) "F" else "C"
}

closest_available_temperature <- function(x, target) {
  x[which.min(abs(x - target))]
}

harmonize_year_values <- function(x) {
  x_chr <- as.character(x)
  if (all(na.omit(unique(x_chr)) %in% c("1", "2", "3", "4"))) {
    dplyr::recode(
      x_chr,
      "1" = "2019",
      "2" = "2020",
      "3" = "2021",
      "4" = "2022",
      .default = x_chr
    )
  } else {
    x_chr
  }
}

normalize_community_key <- function(x) {
  out <- x %>%
    as.character() %>%
    stringr::str_to_upper() %>%
    stringr::str_replace_all("&", "AND") %>%
    stringr::str_replace_all("[[:punct:]]", " ") %>%
    stringr::str_squish()
  
  out <- dplyr::case_when(
    out %in% c("O HARE", "OHARE") ~ "OHARE",
    out %in% c("MCKINLEY PARK") ~ "MCKINLEY PARK",
    out %in% c("NEAR WEST SIDE") ~ "NEAR WEST SIDE",
    out %in% c("NEAR NORTH SIDE") ~ "NEAR NORTH SIDE",
    out %in% c("NEAR SOUTH SIDE") ~ "NEAR SOUTH SIDE",
    TRUE ~ out
  )
  
  out
}

label_community_pretty <- function(x) {
  x %>%
    as.character() %>%
    str_replace_all("_", " ") %>%
    str_squish() %>%
    str_to_title() %>%
    str_replace("^Ohare$", "O'Hare")
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

get_first_match <- function(df, patterns) {
  hits <- names(df)[str_detect(names(df), patterns)]
  if (length(hits) == 0) return(NA_character_)
  hits[1]
}

get_temp_col <- function(df) {
  if ("temperature" %in% names(df)) return("temperature")
  get_first_match(df, "temp_value|temp_f|temperature")
}

get_endpoint_risk_col <- function(df) {
  if ("endpoint_risk_0_100" %in% names(df)) return("endpoint_risk_0_100")
  stop("Expected endpoint_risk_0_100 in dataframe, but it was not found.")
}

get_overall_risk_col <- function(df) {
  if ("overall_risk_0_100" %in% names(df)) return("overall_risk_0_100")
  stop("Expected overall_risk_0_100 in dataframe, but it was not found.")
}

get_structural_col <- function(df) {
  if ("overall_structural_0_100" %in% names(df)) return("overall_structural_0_100")
  stop("Expected overall_structural_0_100 in dataframe, but it was not found.")
}

get_excess_col <- function(df) {
  get_first_match(df, "excess_events|total_excess_events|overall_weighted_excess")
}

safe_entropy <- function(x) {
  x <- x[is.finite(x) & x > 0]
  if (length(x) == 0) return(NA_real_)
  p <- x / sum(x)
  -sum(p * log(p))
}

safe_weighted_mean <- function(x, w) {
  ok <- is.finite(x) & is.finite(w)
  if (sum(ok) == 0) return(NA_real_)
  if (sum(w[ok], na.rm = TRUE) <= 0) return(mean(x[ok], na.rm = TRUE))
  weighted.mean(x[ok], w[ok], na.rm = TRUE)
}

infer_map_sf <- function() {
  if (exists("ca", inherits = TRUE)) {
    obj <- get("ca", inherits = TRUE)
    if (inherits(obj, "sf")) return(obj)
  }
  if (exists("commareas", inherits = TRUE)) {
    obj <- get("commareas", inherits = TRUE)
    if (inherits(obj, "sf")) return(obj)
  }
  stop("No sf map object found. Expected something like `cas` or `commareas` in the environment.")
}

standardize_map_key <- function(sf_obj) {
  nm <- names(sf_obj)
  candidate_cols <- c("community", "ca_name", "commarea", "community_area", "area_name", "community_name")
  candidate_cols <- candidate_cols[candidate_cols %in% nm]
  
  if (length(candidate_cols) == 0) {
    stop("Map object does not have a recognizable community identifier column.")
  }
  
  score_col <- purrr::map_dbl(candidate_cols, function(col) {
    x <- sf_obj[[col]]
    if (is.factor(x)) x <- as.character(x)
    if (!is.character(x)) return(-Inf)
    x2 <- str_to_upper(str_squish(x))
    mean(str_detect(x2, "[A-Z]"), na.rm = TRUE)
  })
  
  best_col <- candidate_cols[which.max(score_col)]
  
  sf_obj %>%
    mutate(
      community = normalize_community_key(.data[[best_col]])
    )
}

make_choropleth <- function(sf_df, value_col, title, subtitle = NULL, fill_label = NULL,
                            na_fill = "grey90") {
  ggplot(sf_df) +
    geom_sf(aes(fill = .data[[value_col]]), color = "white", linewidth = 0.15) +
    coord_sf(expand = FALSE) +
    scale_fill_hvi_0_100(fill_label %||% value_col) +
    wide_colorbar +
    labs(
      title = title,
      subtitle = subtitle,
      caption = "Chicago community areas"
    ) +
    theme_pub_map(base_size = base_size)
}

make_categorical_map <- function(sf_df, value_col, title, subtitle = NULL) {
  ggplot(sf_df) +
    geom_sf(aes(fill = .data[[value_col]]), color = "white", linewidth = 0.15) +
    coord_sf(expand = FALSE) +
    scale_fill_viridis_d(na.value = "grey90", begin = 0.1, end = 0.9, name = NULL) +
    labs(
      title = title,
      subtitle = subtitle,
      caption = "Chicago community areas"
    ) +
    theme_pub(base_size = base_size) +
    theme(
      axis.text = element_blank(),
      axis.title = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank(),
      legend.position = "bottom"
    )
}

scale_fill_hvi_0_100 <- function(name) {
  scale_fill_viridis_c(
    option = "C",
    limits = c(0, 100),
    breaks = c(0, 25, 50, 75, 100),
    labels = scales::label_number(accuracy = 1),
    oob = scales::squish,
    name = name
  )
}

theme_pub_map_shared_legend <- function(base_size = 11) {
  theme_pub(base_size = base_size) +
    theme(
      axis.text = element_blank(),
      axis.title = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank(),
      legend.position = "bottom",
      legend.title = element_text(face = "bold", size = base_size + 1),
      legend.text = element_text(size = base_size),
      legend.key.width = unit(1.8, "cm"),
      legend.key.height = unit(0.5, "cm")
    )
}

wide_colorbar <- guides(
  fill = guide_colorbar(
    barwidth = unit(12, "cm"),
    barheight = unit(0.7, "cm"),
    title.position = "top",
    title.hjust = 0.5
  )
)

theme_pub_map <- function(base_size = 11) {
  theme_pub(base_size = base_size) +
    theme(
      axis.text = element_blank(),
      axis.title = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank(),
      legend.position = "bottom",
      legend.title = element_text(face = "bold", size = base_size + 1),
      legend.text = element_text(size = base_size),
      legend.key.width = unit(1.5, "cm"),
      legend.key.height = unit(0.6, "cm")
    )
}

wide_colorbar <- guides(
  fill = guide_colorbar(
    barwidth = unit(10, "cm"),
    barheight = unit(0.7, "cm"),
    title.position = "top",
    title.hjust = 0.5
  )
)

# ------------------------------------------------------------------------------------------------
# REQUIRED OBJECTS
# ------------------------------------------------------------------------------------------------
required_any <- c(
  "ca_year_overall_structural_hvi",
  "ca_day_overall_operational_hvi",
  "temp_grid_overall_risk",
  "temperature_grid_overall_risk"
)

if (!any(vapply(required_any, exists, logical(1), inherits = TRUE))) {
  stop("No core HVI output objects found in the environment.")
}

have_model_matrix  <- exists("hvi_model_matrix", inherits = TRUE)
have_meta          <- exists("hvi_endpoint_metadata", inherits = TRUE)
have_perf          <- exists("endpoint_model_performance", inherits = TRUE)
have_weights       <- exists("endpoint_weights", inherits = TRUE)
have_interactions  <- exists("endpoint_interaction_betas", inherits = TRUE)
have_structural_ep <- exists("ca_year_endpoint_vulnerability", inherits = TRUE)
have_structural_ov <- exists("ca_year_overall_structural_hvi", inherits = TRUE)
have_daily_ep      <- exists("ca_day_endpoint_risk", inherits = TRUE)
have_daily_ov      <- exists("ca_day_overall_operational_hvi", inherits = TRUE)

temp_endpoint_obj <- case_when(
  exists("temperature_grid_endpoint_risk", inherits = TRUE) ~ "temperature_grid_endpoint_risk",
  exists("temp_grid_endpoint_risk", inherits = TRUE) ~ "temp_grid_endpoint_risk",
  TRUE ~ NA_character_
)

temp_overall_obj <- case_when(
  exists("temperature_grid_overall_risk", inherits = TRUE) ~ "temperature_grid_overall_risk",
  exists("temp_grid_overall_risk", inherits = TRUE) ~ "temp_grid_overall_risk",
  TRUE ~ NA_character_
)

temp_family_obj <- case_when(
  exists("temperature_grid_family_risk", inherits = TRUE) ~ "temperature_grid_family_risk",
  exists("temp_grid_family_risk", inherits = TRUE) ~ "temp_grid_family_risk",
  TRUE ~ NA_character_
)

have_temp_endpoint <- !is.na(temp_endpoint_obj)
have_temp_overall  <- !is.na(temp_overall_obj)
have_temp_family   <- !is.na(temp_family_obj)

# ------------------------------------------------------------------------------------------------
# LOAD OBJECTS
# ------------------------------------------------------------------------------------------------
if (have_model_matrix)  hvi_model_matrix <- get("hvi_model_matrix", inherits = TRUE) %>% clean_names()
if (have_meta)          ep_meta <- get("hvi_endpoint_metadata", inherits = TRUE) %>% clean_names()
if (have_perf)          perf <- get("endpoint_model_performance", inherits = TRUE) %>% clean_names()
if (have_weights)       ep_weights <- get("endpoint_weights", inherits = TRUE) %>% clean_names()
if (have_interactions)  ep_betas <- get("endpoint_interaction_betas", inherits = TRUE) %>% clean_names()
if (have_structural_ep) structural_ep <- get("ca_year_endpoint_vulnerability", inherits = TRUE) %>% clean_names()
if (have_structural_ov) structural_overall <- get("ca_year_overall_structural_hvi", inherits = TRUE) %>% clean_names()
if (have_daily_ep)      daily_ep <- get("ca_day_endpoint_risk", inherits = TRUE) %>% clean_names()
if (have_daily_ov)      daily_overall <- get("ca_day_overall_operational_hvi", inherits = TRUE) %>% clean_names()
if (have_temp_endpoint) temp_ep <- get(temp_endpoint_obj, inherits = TRUE) %>% clean_names()
if (have_temp_overall)  temp_overall <- get(temp_overall_obj, inherits = TRUE) %>% clean_names()
if (have_temp_family)   temp_family <- get(temp_family_obj, inherits = TRUE) %>% clean_names()

chi_sf <- infer_map_sf() %>%
  standardize_map_key() %>%
  st_as_sf() %>%
  st_make_valid()

message("O'Hare in map object: ", "OHARE" %in% chi_sf$community)
if (have_temp_overall) message("O'Hare in temp_overall: ", "OHARE" %in% temp_overall$community)
if (have_temp_endpoint) message("O'Hare in temp_ep: ", "OHARE" %in% temp_ep$community)

# ------------------------------------------------------------------------------------------------
# STANDARDIZE KEYS AND YEARS
# ------------------------------------------------------------------------------------------------
if (have_model_matrix) {
  hvi_model_matrix <- hvi_model_matrix %>%
    mutate(
      community = normalize_community_key(community),
      year = harmonize_year_values(year)
    )
}

if (have_structural_ep) {
  structural_ep <- structural_ep %>%
    mutate(
      community = normalize_community_key(community),
      year = harmonize_year_values(year)
    )
}

if (have_structural_ov) {
  structural_overall <- structural_overall %>%
    mutate(
      community = normalize_community_key(community),
      year = harmonize_year_values(year)
    )
}

if (have_daily_ep) {
  daily_ep <- daily_ep %>%
    mutate(
      community = normalize_community_key(community),
      date = as.Date(date),
      year = harmonize_year_values(year)
    )
}

if (have_daily_ov) {
  daily_overall <- daily_overall %>%
    mutate(
      community = normalize_community_key(community),
      date = as.Date(date),
      year = harmonize_year_values(year)
    )
}

if (have_temp_endpoint) {
  temp_ep <- temp_ep %>%
    mutate(
      community = normalize_community_key(community),
      year = harmonize_year_values(year)
    )
  temp_col_tmp <- get_temp_col(temp_ep)
  temp_ep <- temp_ep %>% mutate(temperature = as.numeric(.data[[temp_col_tmp]]))
}

if (have_temp_overall) {
  temp_overall <- temp_overall %>%
    mutate(
      community = normalize_community_key(community),
      year = harmonize_year_values(year)
    )
  temp_col_tmp <- get_temp_col(temp_overall)
  temp_overall <- temp_overall %>% mutate(temperature = as.numeric(.data[[temp_col_tmp]]))
}

if (have_temp_family) {
  temp_family <- temp_family %>%
    mutate(
      community = normalize_community_key(community),
      year = harmonize_year_values(year)
    )
  temp_col_tmp <- get_temp_col(temp_family)
  temp_family <- temp_family %>% mutate(temperature = as.numeric(.data[[temp_col_tmp]]))
}

if (have_temp_endpoint && have_meta && !"outcome_label" %in% names(temp_ep)) {
  temp_ep <- temp_ep %>%
    left_join(
      ep_meta %>% select(any_of(c("endpoint_key", "outcome_label", "source", "domain"))),
      by = "endpoint_key"
    )
}

if (have_daily_ep && have_meta && !"outcome_label" %in% names(daily_ep)) {
  daily_ep <- daily_ep %>%
    left_join(
      ep_meta %>% select(any_of(c("endpoint_key", "outcome_label", "source", "domain"))),
      by = "endpoint_key"
    )
}

# ------------------------------------------------------------------------------------------------
# DIAGNOSTICS: O'Hare and temperature coverage
# ------------------------------------------------------------------------------------------------
diag_ohare <- tibble(
  dataset = c(
    "chi_sf",
    if (exists("hvi_model_matrix")) "hvi_model_matrix" else NULL,
    if (exists("temp_overall")) "temp_overall" else NULL,
    if (exists("temp_ep")) "temp_ep" else NULL
  ),
  has_ohare = c(
    "OHARE" %in% chi_sf$community,
    if (exists("hvi_model_matrix")) "OHARE" %in% hvi_model_matrix$community else NA,
    if (exists("temp_overall")) "OHARE" %in% temp_overall$community else NA,
    if (exists("temp_ep")) "OHARE" %in% temp_ep$community else NA
  )
)

safe_write_csv(diag_ohare, file.path(tab_dir, "diagnostic_ohare_presence.csv"))

if (have_temp_overall) {
  diag_temp_overall <- temp_overall %>%
    count(temperature, name = "n_rows") %>%
    arrange(temperature)
  safe_write_csv(diag_temp_overall, file.path(tab_dir, "diagnostic_temp_overall_temperature_counts.csv"))
}

if (have_temp_endpoint) {
  diag_temp_ep <- temp_ep %>%
    count(temperature, name = "n_rows") %>%
    arrange(temperature)
  safe_write_csv(diag_temp_ep, file.path(tab_dir, "diagnostic_temp_endpoint_temperature_counts.csv"))
}

if (have_temp_overall) {
  diag_ohare_temp_overall <- temp_overall %>%
    filter(community == "OHARE") %>%
    count(temperature, name = "n_rows") %>%
    arrange(temperature)
  safe_write_csv(diag_ohare_temp_overall, file.path(tab_dir, "diagnostic_ohare_temp_overall.csv"))
}

if (have_temp_endpoint) {
  diag_ohare_temp_ep <- temp_ep %>%
    filter(community == "OHARE") %>%
    count(temperature, endpoint_key, name = "n_rows") %>%
    arrange(temperature, endpoint_key)
  safe_write_csv(diag_ohare_temp_ep, file.path(tab_dir, "diagnostic_ohare_temp_endpoint.csv"))
}

# ------------------------------------------------------------------------------------------------
# TEMPERATURE SETUP
# ------------------------------------------------------------------------------------------------
temp_units <- case_when(
  have_temp_overall ~ detect_temperature_units(temp_overall$temperature),
  have_temp_endpoint ~ detect_temperature_units(temp_ep$temperature),
  TRUE ~ "unknown"
)

target_temps_native <- if (temp_units == "C") f_to_c(target_temps_f) else target_temps_f

target_temp_lookup <- tibble(
  target_f = target_temps_f,
  target_native = target_temps_native
)

if (have_temp_overall) {
  available_temps <- sort(unique(temp_overall$temperature))
} else if (have_temp_endpoint) {
  available_temps <- sort(unique(temp_ep$temperature))
} else {
  available_temps <- numeric(0)
}

target_temp_lookup <- target_temp_lookup %>%
  mutate(
    matched_native = purrr::map_dbl(target_native, ~ closest_available_temperature(available_temps, .x)),
    matched_f = if (temp_units == "C") c_to_f(matched_native) else matched_native
  )

safe_write_csv(target_temp_lookup, file.path(tab_dir, "temperature_targets_lookup.csv"))

# ------------------------------------------------------------------------------------------------
# COMMUNITY COVERAGE QA
# ------------------------------------------------------------------------------------------------
coverage_tbl <- tibble(
  source = c(
    "map_sf",
    if (have_structural_ov) "structural_overall" else NULL,
    if (have_temp_overall) "temp_overall" else NULL,
    if (have_temp_endpoint) "temp_endpoint" else NULL,
    if (have_daily_ov) "daily_overall" else NULL
  ),
  n_communities = c(
    n_distinct(chi_sf$community),
    if (have_structural_ov) n_distinct(structural_overall$community) else NULL,
    if (have_temp_overall) n_distinct(temp_overall$community) else NULL,
    if (have_temp_endpoint) n_distinct(temp_ep$community) else NULL,
    if (have_daily_ov) n_distinct(daily_overall$community) else NULL
  )
)

safe_write_csv(coverage_tbl, file.path(tab_dir, "community_coverage_check.csv"))

o_hare_check <- tibble(
  source = c(
    "map_sf",
    if (have_structural_ov) "structural_overall" else NULL,
    if (have_temp_overall) "temp_overall" else NULL,
    if (have_temp_endpoint) "temp_endpoint" else NULL,
    if (have_daily_ov) "daily_overall" else NULL
  ),
  has_ohare = c(
    "OHARE" %in% chi_sf$community,
    if (have_structural_ov) "OHARE" %in% structural_overall$community else NULL,
    if (have_temp_overall) "OHARE" %in% temp_overall$community else NULL,
    if (have_temp_endpoint) "OHARE" %in% temp_ep$community else NULL,
    if (have_daily_ov) "OHARE" %in% daily_overall$community else NULL
  )
)

safe_write_csv(o_hare_check, file.path(tab_dir, "ohare_presence_check.csv"))

# ------------------------------------------------------------------------------------------------
# TABLE 1: Endpoint model performance
# ------------------------------------------------------------------------------------------------
if (have_perf) {
  perf_summary <- perf
  
  if (have_meta) {
    perf_summary <- perf_summary %>%
      left_join(
        ep_meta %>% select(any_of(c("endpoint_key", "outcome_label", "source", "domain"))),
        by = "endpoint_key",
        suffix = c("", "_meta")
      ) %>%
      mutate(
        outcome_label = coalesce(outcome_label, outcome_label_meta),
        source = coalesce(source, source_meta),
        domain = coalesce(domain, domain_meta)
      ) %>%
      select(-any_of(c("outcome_label_meta", "source_meta", "domain_meta")))
  }
  
  safe_write_csv(perf_summary, file.path(tab_dir, "table_model_performance_summary.csv"))
  
  perf_gt <- perf_summary %>%
    mutate(across(where(is.numeric), ~ round(.x, 3))) %>%
    gt(groupname_col = "cv_type") %>%
    tab_header(
      title = md("**Table. Cross-validated model performance by endpoint**"),
      subtitle = "Spatial and temporal validation summaries"
    ) %>%
    fmt_number(columns = where(is.numeric), decimals = 3)
  
  gtsave(perf_gt, file.path(tab_dir, "table_model_performance_summary.html"))
}

# ------------------------------------------------------------------------------------------------
# TABLE 2: Endpoint weights
# ------------------------------------------------------------------------------------------------
if (have_weights) {
  ep_weights_tbl <- ep_weights %>%
    mutate(across(where(is.numeric), ~ round(.x, 3)))
  
  safe_write_csv(ep_weights_tbl, file.path(tab_dir, "table_endpoint_weights.csv"))
  
  ep_weights_gt <- ep_weights_tbl %>%
    arrange(desc(endpoint_weight)) %>%
    gt() %>%
    tab_header(
      title = md("**Table. Endpoint weights**"),
      subtitle = "Performance-based endpoint contributions to composite risk"
    ) %>%
    fmt_number(columns = where(is.numeric), decimals = 3)
  
  gtsave(ep_weights_gt, file.path(tab_dir, "table_endpoint_weights.html"))
}

# ------------------------------------------------------------------------------------------------
# TABLE 3: Community rankings at multiple temperatures
# ------------------------------------------------------------------------------------------------
if (have_temp_overall) {
  overall_risk_col <- get_overall_risk_col(temp_overall)
  
  top_rankings_multi <- purrr::map_dfr(seq_len(nrow(target_temp_lookup)), function(i) {
    t_match <- target_temp_lookup$matched_native[i]
    t_f <- target_temp_lookup$target_f[i]
    
    temp_overall %>%
      filter(abs(temperature - t_match) < 1e-8) %>%
      arrange(desc(.data[[overall_risk_col]])) %>%
      mutate(
        rank = row_number(),
        target_temp_f = t_f,
        target_temp_native = t_match
      ) %>%
      select(
        target_temp_f, target_temp_native, rank, community,
        all_of(overall_risk_col),
        any_of(c("dominant_endpoint", "dominant_endpoint_risk_0_100"))
      )
  }) %>%
    rename(overall_risk_0_100 = all_of(overall_risk_col)) %>%
    mutate(community_pretty = label_community_pretty(community))
  
  safe_write_csv(top_rankings_multi, file.path(tab_dir, "table_top_communities_multiple_temperatures.csv"))
}

# ------------------------------------------------------------------------------------------------
# DRIVER CONTRIBUTIONS BY TEMPERATURE
# Community-temp-variable weighted structural contribution
# ------------------------------------------------------------------------------------------------
have_driver_temp <- have_model_matrix && have_interactions && have_temp_endpoint && have_weights

if (have_driver_temp) {
  all_z_vars <- names(hvi_model_matrix)[str_detect(names(hvi_model_matrix), "^z_")]
  
  community_year_z <- hvi_model_matrix %>%
    select(community, year, any_of(all_z_vars)) %>%
    distinct()
  
  community_year_z_long <- community_year_z %>%
    pivot_longer(
      cols = all_of(all_z_vars),
      names_to = "variable",
      values_to = "z_value"
    ) %>%
    distinct(community, year, variable, .keep_all = TRUE)
  
  endpoint_driver_base <- ep_betas %>%
    filter(variable %in% all_z_vars) %>%
    distinct(endpoint_key, variable, beta_interaction, .keep_all = TRUE) %>%
    left_join(
      community_year_z_long,
      by = "variable",
      relationship = "many-to-many"
    ) %>%
    mutate(
      base_driver_contribution = z_value * beta_interaction
    )
  
  endpoint_excess_temp <- temp_ep %>%
    select(
      community, year, temperature, endpoint_key,
      excess_events, endpoint_risk_0_100,
      outcome_label, source, domain
    ) %>%
    distinct(community, year, temperature, endpoint_key, .keep_all = TRUE) %>%
    left_join(
      ep_weights %>% select(endpoint_key, endpoint_weight) %>% distinct(),
      by = "endpoint_key"
    ) %>%
    mutate(
      endpoint_weight = coalesce(endpoint_weight, 1)
    )
  
  driver_temp_long <- endpoint_driver_base %>%
    left_join(
      endpoint_excess_temp,
      by = c("community", "year", "endpoint_key"),
      relationship = "many-to-many"
    ) %>%
    mutate(
      weighted_driver_contribution = base_driver_contribution * pmax(excess_events, 0) * endpoint_weight
    )
  
  city_driver_by_temp <- driver_temp_long %>%
    group_by(temperature, variable) %>%
    summarise(
      driver_importance_raw = sum(weighted_driver_contribution, na.rm = TRUE),
      abs_driver_importance = sum(abs(weighted_driver_contribution), na.rm = TRUE),
      .groups = "drop"
    ) %>%
    group_by(temperature) %>%
    mutate(
      driver_importance_0_100 = rescale_0_100(abs_driver_importance)
    ) %>%
    ungroup() %>%
    mutate(
      variable_label = case_when(
        str_detect(variable, "^z_") ~ str_remove(variable, "^z_") %>%
          str_replace_all("_", " ") %>%
          str_to_title(),
        TRUE ~ variable
      )
    )
  
  safe_write_csv(
    driver_temp_long,
    file.path(sup_dir, "supplement_driver_contributions_by_temperature_long.csv")
  )
  
  safe_write_csv(
    city_driver_by_temp,
    file.path(tab_dir, "table_city_driver_importance_by_temperature.csv")
  )
}

# ------------------------------------------------------------------------------------------------
# FIGURE 1: Overall structural vulnerability map
# ------------------------------------------------------------------------------------------------
if (have_structural_ov) {
  structural_col <- get_structural_col(structural_overall)
  
  structural_map_df <- chi_sf %>%
    left_join(
      structural_overall %>% filter(year == as.character(target_year)),
      by = "community"
    ) %>%
    st_as_sf() %>%
    st_make_valid()
  
  p_structural <- make_choropleth(
    structural_map_df,
    value_col = structural_col,
    title = "Structural heat vulnerability across Chicago",
    subtitle = glue("Community-level structural HVI, {target_year}"),
    fill_label = "Structural HVI\n(0-100)"
  )
  
  save_pub_plot(p_structural, "fig1_structural_overall_hvi_map.png", width = 7.6, height = 7.2)
}

# ------------------------------------------------------------------------------------------------
# FIGURE 2: Overall dynamic risk maps at multiple temperatures
# Uses focal-temperature raw burden, rescaled within each temperature slice to 0-100
# ------------------------------------------------------------------------------------------------
if (have_temp_overall) {
  dynamic_source_col <- case_when(
    "overall_weighted_excess" %in% names(temp_overall) ~ "overall_weighted_excess",
    "total_excess_events" %in% names(temp_overall) ~ "total_excess_events",
    "overall_risk_0_100" %in% names(temp_overall) ~ "overall_risk_0_100",
    TRUE ~ NA_character_
  )
  
  if (is.na(dynamic_source_col)) {
    stop("No usable dynamic source column found in temp_overall.")
  }
  
  map_list <- purrr::map(seq_len(nrow(target_temp_lookup)), function(i) {
    t_match <- target_temp_lookup$matched_native[i]
    t_f <- target_temp_lookup$target_f[i]
    
    dat_i <- temp_overall %>%
      filter(abs(temperature - t_match) < 1e-8) %>%
      group_by(community) %>%
      summarise(
        dynamic_raw = mean(.data[[dynamic_source_col]], na.rm = TRUE),
        .groups = "drop"
      ) %>%
      mutate(
        dynamic_plot_0_100 = rescale_0_100(dynamic_raw)
      )
    
    map_i <- tibble(community = unique(chi_sf$community)) %>%
      left_join(dat_i, by = "community") %>%
      mutate(dynamic_plot_0_100 = replace_na(dynamic_plot_0_100, 0)) %>%
      left_join(
        chi_sf %>% select(community, geometry),
        by = "community"
      ) %>%
      st_as_sf() %>%
      st_make_valid()
    
    ggplot(map_i) +
      geom_sf(aes(fill = dynamic_plot_0_100), color = "white", linewidth = 0.15) +
      coord_sf(expand = FALSE) +
      scale_fill_hvi_0_100("Overall risk\n(0-100)") +
      wide_colorbar +
      labs(
        title = glue("{t_f}°F"),
        subtitle = NULL,
      ) +
      theme_pub_map_shared_legend(base_size = base_size) +
      theme(
        plot.title = element_text(hjust = 0.5),
        legend.position = "bottom"
      )
  })
  
  fig2 <- wrap_plots(map_list, ncol = length(map_list), guides = "collect") &
    theme(legend.position = "bottom")
  
  fig2 <- fig2 +
    plot_annotation(
      title = "Dynamic heat-health risk across multiple temperatures",
      theme = theme(
        plot.title = element_text(face = "bold", size = 14),
        plot.subtitle = element_text(size = 11)
      )
    )
  
  save_pub_plot(fig2, "fig2_overall_dynamic_risk_multi_temp.png", width = 15, height = 6.5)
}

# ------------------------------------------------------------------------------------------------
# FIGURE 3: Endpoint composition by temperature
# ------------------------------------------------------------------------------------------------
if (have_temp_endpoint) {
  exclude_endpoints <- c("deaths", "ed_visits", "ems_calls")
  
  endpoint_comp <- purrr::map_dfr(seq_len(nrow(target_temp_lookup)), function(i) {
    t_match <- target_temp_lookup$matched_native[i]
    t_f <- target_temp_lookup$target_f[i]
    
    temp_ep %>%
      filter(
        abs(temperature - t_match) < 1e-8,
        !endpoint_key %in% exclude_endpoints
      ) %>%
      group_by(endpoint_key, outcome_label, source, domain) %>%
      summarise(
        total_excess = sum(pmax(excess_events, 0), na.rm = TRUE),
        .groups = "drop"
      ) %>%
      mutate(
        target_temp_f = t_f,
        endpoint_label = if_else(!is.na(outcome_label), outcome_label, pretty_endpoint_label(endpoint_key))
      )
  })
  
  if (nrow(endpoint_comp) > 0) {
    all_endpoint_labels <- endpoint_comp %>%
      distinct(endpoint_label) %>%
      arrange(endpoint_label) %>%
      pull(endpoint_label)
    
    endpoint_comp <- endpoint_comp %>%
      select(target_temp_f, endpoint_label, total_excess) %>%
      tidyr::complete(
        target_temp_f = target_temps_f,
        endpoint_label = all_endpoint_labels,
        fill = list(total_excess = 0)
      ) %>%
      group_by(target_temp_f) %>%
      mutate(
        total_temp_excess = sum(total_excess, na.rm = TRUE),
        share = if_else(total_temp_excess > 0, total_excess / total_temp_excess, 0)
      ) %>%
      ungroup()
    
    p_endpoint_comp <- endpoint_comp %>%
      mutate(
        endpoint_label = fct_reorder(endpoint_label, share, .fun = max, .desc = TRUE),
        target_temp_f = factor(target_temp_f, levels = target_temps_f)
      ) %>%
      ggplot(aes(x = target_temp_f, y = share, fill = endpoint_label)) +
      geom_col(position = "fill", width = 0.75) +
      scale_y_continuous(labels = percent_format()) +
      labs(
        title = "Cause-specific endpoint composition changes with temperature",
        subtitle = "Stacked proportions of total modeled excess events across endpoints",
        x = "Temperature (°F)",
        y = "Share of modeled excess events",
        fill = "Endpoint"
      ) +
      theme_pub(base_size = base_size) +
      theme(legend.position = "right")
    
    save_pub_plot(p_endpoint_comp, "fig3_endpoint_composition_by_temperature.png", width = 14, height = 9)
    
    safe_write_csv(endpoint_comp, file.path(tab_dir, "table_endpoint_composition_by_temperature.csv"))
  }
}

# ------------------------------------------------------------------------------------------------
# FIGURE 4: Outcome entropy maps at multiple temperatures
# ------------------------------------------------------------------------------------------------
if (have_temp_endpoint) {
  entropy_all <- purrr::map_dfr(seq_len(nrow(target_temp_lookup)), function(i) {
    t_match <- target_temp_lookup$matched_native[i]
    t_f <- target_temp_lookup$target_f[i]
    
    temp_ep %>%
      filter(abs(temperature - t_match) < 1e-8) %>%
      group_by(community, year) %>%
      summarise(
        endpoint_entropy_raw = safe_entropy(pmax(excess_events, 0)),
        .groups = "drop"
      ) %>%
      mutate(target_temp_f = t_f)
  }) %>%
    group_by(target_temp_f) %>%
    mutate(endpoint_entropy_0_100 = rescale_0_100(endpoint_entropy_raw)) %>%
    ungroup()
  
  entropy_maps <- purrr::map(seq_len(nrow(target_temp_lookup)), function(i) {
    t_f <- target_temp_lookup$target_f[i]
    
    entropy_dat <- entropy_all %>%
      filter(target_temp_f == t_f) %>%
      group_by(community) %>%
      summarise(
        endpoint_entropy_plot = mean(endpoint_entropy_0_100, na.rm = TRUE),
        .groups = "drop"
      )
    
    entropy_sf <- tibble(community = unique(chi_sf$community)) %>%
      left_join(entropy_dat, by = "community") %>%
      mutate(endpoint_entropy_plot = replace_na(endpoint_entropy_plot, 0)) %>%
      left_join(
        chi_sf %>% select(community, geometry),
        by = "community"
      ) %>%
      st_as_sf() %>%
      st_make_valid()
    
    ggplot(entropy_sf) +
      geom_sf(aes(fill = endpoint_entropy_plot), color = "white", linewidth = 0.15) +
      coord_sf(expand = FALSE) +
      scale_fill_hvi_0_100("Endpoint\nentropy\n(0-100)") +
      wide_colorbar +
      labs(
        title = glue("{t_f}°F"),
        subtitle = NULL,
        caption = "Chicago community areas"
      ) +
      theme_pub_map_shared_legend(base_size = base_size) +
      theme(
        plot.title = element_text(hjust = 0.5),
        legend.position = "bottom"
      )
  })
  
  fig4 <- wrap_plots(entropy_maps, ncol = length(entropy_maps), guides = "collect") &
    theme(legend.position = "bottom")
  
  fig4 <- fig4 +
    plot_annotation(
      title = "Outcome entropy across Chicago community areas",
      subtitle = "Higher entropy indicates a more diffuse endpoint burden profile; lower entropy indicates domination by fewer endpoints",
      theme = theme(
        plot.title = element_text(face = "bold", size = 14),
        plot.subtitle = element_text(size = 11)
      )
    )
  
  save_pub_plot(fig4, "fig4_entropy_maps_multi_temp.png", width = 15, height = 6.5)
}

# ------------------------------------------------------------------------------------------------
# FIGURE 5: Structural driver importance by temperature
# ------------------------------------------------------------------------------------------------
if (have_driver_temp) {
  fig5_dat <- city_driver_by_temp %>%
    filter(temperature %in% target_temp_lookup$matched_native) %>%
    group_by(temperature) %>%
    slice_max(order_by = abs_driver_importance, n = 12, with_ties = FALSE) %>%
    ungroup() %>%
    left_join(
      target_temp_lookup %>% select(target_f, matched_native),
      by = c("temperature" = "matched_native")
    ) %>%
    mutate(
      target_f = factor(target_f, levels = target_temps_f),
      variable_label = fct_reorder(variable_label, abs_driver_importance, .desc = TRUE)
    )
  
  p_driver_temp <- ggplot(fig5_dat, aes(x = abs_driver_importance, y = variable_label)) +
    geom_col() +
    facet_wrap(~ target_f, scales = "free_x") +
    labs(
      title = "Structural drivers of modeled heat-health risk across temperatures",
      subtitle = "Citywide weighted driver influence across endpoints",
      x = "Absolute weighted driver contribution",
      y = NULL
    ) +
    theme_pub(base_size = base_size)
  
  save_pub_plot(p_driver_temp, "fig5_structural_driver_importance_by_temperature.png", width = 13, height = 9)
}

# ------------------------------------------------------------------------------------------------
# FIGURE 6: All-community endpoint heatmap at focal temperature (cause-specific only)
# ------------------------------------------------------------------------------------------------
if (have_temp_endpoint && have_temp_overall) {
  focal_f <- 91
  focal_match <- target_temp_lookup %>%
    filter(target_f == focal_f) %>%
    pull(matched_native)
  
  endpoint_risk_col <- get_endpoint_risk_col(temp_ep)
  overall_risk_col  <- get_overall_risk_col(temp_overall)
  
  exclude_endpoints <- c("deaths", "ed_visits", "ems_calls")
  
  comm_order <- temp_overall %>%
    filter(abs(temperature - focal_match) < 1e-8) %>%
    group_by(community) %>%
    summarise(
      overall_risk = mean(.data[[overall_risk_col]], na.rm = TRUE),
      .groups = "drop"
    ) %>%
    right_join(
      tibble(community = unique(chi_sf$community)),
      by = "community"
    ) %>%
    mutate(overall_risk = replace_na(overall_risk, 0)) %>%
    arrange(desc(overall_risk)) %>%
    pull(community) %>%
    unique()
  
  endpoint_levels <- temp_ep %>%
    filter(
      abs(temperature - focal_match) < 1e-8,
      !endpoint_key %in% exclude_endpoints
    ) %>%
    group_by(endpoint_key) %>%
    summarise(
      mean_risk = mean(.data[[endpoint_risk_col]], na.rm = TRUE),
      outcome_label = dplyr::first(outcome_label),
      .groups = "drop"
    ) %>%
    mutate(
      endpoint_label = if_else(!is.na(outcome_label), outcome_label, pretty_endpoint_label(endpoint_key))
    ) %>%
    arrange(desc(mean_risk)) %>%
    pull(endpoint_label)
  
  heatmap_full <- temp_ep %>%
    filter(
      abs(temperature - focal_match) < 1e-8,
      !endpoint_key %in% exclude_endpoints
    ) %>%
    group_by(community, endpoint_key) %>%
    summarise(
      endpoint_risk_plot = mean(.data[[endpoint_risk_col]], na.rm = TRUE),
      outcome_label = dplyr::first(outcome_label),
      .groups = "drop"
    ) %>%
    mutate(
      endpoint_label = if_else(!is.na(outcome_label), outcome_label, pretty_endpoint_label(endpoint_key))
    ) %>%
    select(community, endpoint_label, endpoint_risk_plot) %>%
    tidyr::complete(
      community = unique(chi_sf$community),
      endpoint_label = endpoint_levels,
      fill = list(endpoint_risk_plot = 0)
    ) %>%
    mutate(
      community = factor(community, levels = comm_order),
      endpoint_label = factor(endpoint_label, levels = endpoint_levels)
    )
  
  p_heatmap_full <- ggplot(
    heatmap_full,
    aes(x = endpoint_label, y = community, fill = endpoint_risk_plot)
  ) +
    geom_tile(color = "white", linewidth = 0.12) +
    scale_fill_hvi_0_100("Endpoint risk\n(0-100)") +
    wide_colorbar +
    labs(
      title = "Cause-specific endpoint risk architecture across all Chicago community areas",
      subtitle = glue("All mapped community areas at ~{focal_f}°F; all-cause endpoints excluded"),
      x = NULL,
      y = NULL
    ) +
    theme_pub(base_size = 13) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
      axis.text.y = element_text(size = 10),
      panel.grid = element_blank(),
      legend.position = "bottom",
      legend.title = element_text(size = 12, face = "bold"),
      legend.text = element_text(size = 11)
    )
  
  save_pub_plot(
    p_heatmap_full,
    "fig6_all_community_endpoint_heatmap_92f.png",
    width = 13,
    height = 18
  )
}

# ------------------------------------------------------------------------------------------------
# FIGURE 7: Structural vs dynamic risk scatter
# ------------------------------------------------------------------------------------------------
if (have_structural_ov && have_temp_overall) {
  structural_col <- get_structural_col(structural_overall)
  
  focal_f <- 91
  focal_match <- target_temp_lookup %>%
    filter(target_f == focal_f) %>%
    pull(matched_native)
  
  # use raw focal-temperature burden, then rescale within the focal slice
  dynamic_source_col <- if ("overall_weighted_excess" %in% names(temp_overall)) {
    "overall_weighted_excess"
  } else if ("total_excess_events" %in% names(temp_overall)) {
    "total_excess_events"
  } else if ("overall_risk_0_100" %in% names(temp_overall)) {
    "overall_risk_0_100"
  } else {
    stop("No usable dynamic-risk column found in temp_overall.")
  }
  
  dynamic_dat <- temp_overall %>%
    filter(abs(temperature - focal_match) < 1e-8) %>%
    group_by(community) %>%
    summarise(
      dynamic_raw = mean(.data[[dynamic_source_col]], na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      dynamic_score = rescale_0_100(dynamic_raw)
    )
  
  scatter_dat <- structural_overall %>%
    filter(year == as.character(target_year)) %>%
    select(community, structural_score = all_of(structural_col)) %>%
    left_join(dynamic_dat, by = "community") %>%
    drop_na() %>%
    mutate(
      community_label = label_community_pretty(community)
    )
  
  top_labs <- scatter_dat %>%
    mutate(label_score = structural_score + dynamic_score) %>%
    slice_max(label_score, n = 6)
  
  p_scatter <- ggplot(scatter_dat, aes(x = structural_score, y = dynamic_score)) +
    geom_point(size = 2.7, alpha = 0.85) +
    geom_smooth(method = "lm", se = FALSE, linewidth = 0.8, linetype = 2, color = "#3B6CF6") +
    ggrepel::geom_text_repel(
      data = top_labs,
      aes(label = community_label),
      size = 3.2,
      max.overlaps = Inf,
      box.padding = 0.2,
      point.padding = 0.2
    ) +
    scale_x_continuous(limits = c(0, 100), breaks = c(0, 25, 50, 75, 100)) +
    scale_y_continuous(limits = c(0, 100), breaks = c(0, 25, 50, 75, 100)) +
    labs(
      title = "Structural vulnerability versus dynamic risk",
      subtitle = glue("{target_year} structural HVI versus focal-temperature dynamic risk at ~{focal_f}°F"),
      x = "Structural HVI (0-100)",
      y = "Dynamic heat-health risk (0-100)"
    ) +
    theme_pub(base_size = base_size)
  
  save_pub_plot(p_scatter, "fig7_structural_vs_dynamic_scatter.png", width = 8.5, height = 6.5)
}

# ------------------------------------------------------------------------------------------------
# FIGURE 8: Model performance
# ------------------------------------------------------------------------------------------------
if (have_perf) {
  perf_col <- case_when(
    "mean_cor_comm_spearman" %in% names(perf) ~ "mean_cor_comm_spearman",
    "mean_cor_comm_pearson" %in% names(perf) ~ "mean_cor_comm_pearson",
    TRUE ~ NA_character_
  )
  
  if (!is.na(perf_col)) {
    perf_fig_dat <- perf
    
    if (have_meta) {
      perf_fig_dat <- perf_fig_dat %>%
        left_join(
          ep_meta %>% select(any_of(c("endpoint_key", "outcome_label", "source"))),
          by = "endpoint_key",
          suffix = c("", "_meta")
        ) %>%
        mutate(
          outcome_label = coalesce(outcome_label, outcome_label_meta),
          source = coalesce(source, source_meta)
        ) %>%
        select(-any_of(c("outcome_label_meta", "source_meta")))
    }
    
    perf_fig_dat <- perf_fig_dat %>%
      mutate(
        endpoint_label = if_else(!is.na(outcome_label), outcome_label, endpoint_key),
        endpoint_label = fct_reorder(endpoint_label, .data[[perf_col]])
      )
    
    p_perf <- ggplot(perf_fig_dat, aes(x = .data[[perf_col]], y = endpoint_label)) +
      geom_col() +
      facet_wrap(~ cv_type, ncol = 2, scales = "free_y") +
      labs(
        title = "Cross-validated model performance across endpoints",
        subtitle = "Community-level correlation between observed and predicted burden",
        x = ifelse(
          perf_col == "mean_cor_comm_spearman",
          "Mean community-level Spearman correlation",
          "Mean community-level Pearson correlation"
        ),
        y = NULL
      ) +
      theme_pub(base_size = base_size)
    
    save_pub_plot(p_perf, "fig8_model_performance.png", width = 12, height = 9)
  }
}

# ------------------------------------------------------------------------------------------------
# SUPPLEMENTAL EXPORTS
# ------------------------------------------------------------------------------------------------
if (have_daily_ep) safe_write_csv(daily_ep, file.path(sup_dir, "supplement_daily_endpoint_risk.csv"))
if (have_daily_ov) safe_write_csv(daily_overall, file.path(sup_dir, "supplement_daily_overall_operational_hvi.csv"))
if (have_temp_endpoint) safe_write_csv(temp_ep, file.path(sup_dir, "supplement_temperature_grid_endpoint_risk.csv"))
if (have_temp_overall) safe_write_csv(temp_overall, file.path(sup_dir, "supplement_temperature_grid_overall_risk.csv"))
if (have_temp_family) safe_write_csv(temp_family, file.path(sup_dir, "supplement_temperature_grid_family_risk.csv"))
if (have_structural_ep) safe_write_csv(structural_ep, file.path(sup_dir, "supplement_structural_endpoint_vulnerability.csv"))
if (have_structural_ov) safe_write_csv(structural_overall, file.path(sup_dir, "supplement_structural_overall_hvi.csv"))
if (have_interactions) safe_write_csv(ep_betas, file.path(sup_dir, "supplement_endpoint_interaction_betas.csv"))

# ------------------------------------------------------------------------------------------------
# README
# ------------------------------------------------------------------------------------------------
summary_lines <- c(
  glue("Publication outputs generated on {Sys.Date()}."),
  glue("Structural map year: {target_year}."),
  glue("Temperatures requested (F): {paste(target_temps_f, collapse = ', ')}."),
  glue("Detected temperature units in scoring grid: {temp_units}."),
  glue("Matched native temperatures: {paste(round(target_temp_lookup$matched_native, 2), collapse = ', ')}."),
  glue("Matched Fahrenheit temperatures: {paste(round(target_temp_lookup$matched_f, 1), collapse = ', ')}."),
  "All-community endpoint heatmap includes every community area found after key harmonization.",
  "Community keys were normalized to preserve O'Hare as OHARE during map joins."
)

writeLines(summary_lines, con = file.path(pub_dir, "README_publication_outputs.txt"))

message("Done. Publication outputs written to: ", normalizePath(pub_dir))
message("Communities in map sf: ", n_distinct(chi_sf$community))
if (have_structural_ov) message("Communities in structural overall: ", n_distinct(structural_overall$community))
if (have_temp_overall) message("Communities in temp overall: ", n_distinct(temp_overall$community))
if (have_temp_endpoint) message("Communities in temp endpoint: ", n_distinct(temp_ep$community))
message("O'Hare present in map sf: ", "OHARE" %in% chi_sf$community)
if (have_structural_ov) message("O'Hare present in structural overall: ", "OHARE" %in% structural_overall$community)
if (have_temp_overall) message("O'Hare present in temp overall: ", "OHARE" %in% temp_overall$community)
if (have_temp_endpoint) message("O'Hare present in temp endpoint: ", "OHARE" %in% temp_ep$community)

heatmap_ohare_check <- if (exists("heatmap_full")) {
  tibble(
    figure = "fig6_heatmap",
    has_ohare = "OHARE" %in% unique(as.character(heatmap_full$community))
  )
} else {
  tibble(figure = "fig6_heatmap", has_ohare = NA)
}

map_ohare_check <- tibble(
  figure = c("fig2_dynamic_maps", "fig4_entropy_maps"),
  has_ohare = c(
    "OHARE" %in% chi_sf$community,
    "OHARE" %in% chi_sf$community
  )
)

safe_write_csv(
  bind_rows(heatmap_ohare_check, map_ohare_check),
  file.path(tab_dir, "publication_ohare_check.csv")
)


# ------------------------------------------------------------------------------------------------
# FIGURE 9: Community-specific temperature response curves
# ------------------------------------------------------------------------------------------------
if (have_temp_endpoint) {
  selected_communities_for_curves <- c(
    "AUSTIN",
    "NEAR NORTH SIDE",
    "ENGLEWOOD",
    "SOUTH SHORE",
    "NEAR WEST SIDE",
    "HYDE PARK",
    "AVONDALE",
    "GARFIELD RIDGE"
    
  )
  
  selected_endpoints_for_curves <- c(
    "deaths",
    "death_cvd",
    "ed_visits",
    "ed_cvd",
    "ems_calls",
    "ems_cvd"
  )
  
  endpoint_risk_col <- get_endpoint_risk_col(temp_ep)
  
  curve_dat <- temp_ep %>%
    filter(
      community %in% selected_communities_for_curves,
      endpoint_key %in% selected_endpoints_for_curves
    ) %>%
    group_by(community, endpoint_key, temperature) %>%
    summarise(
      endpoint_risk_plot = mean(.data[[endpoint_risk_col]], na.rm = TRUE),
      outcome_label = dplyr::first(outcome_label),
      .groups = "drop"
    ) %>%
    mutate(
      endpoint_label = if_else(!is.na(outcome_label), outcome_label, pretty_endpoint_label(endpoint_key)),
      community_label = label_community_pretty(community),
      temperature_f_plot = if (temp_units == "C") c_to_f(temperature) else temperature
    )
  
  if (nrow(curve_dat) > 0) {
    p_curves <- ggplot(
      curve_dat,
      aes(x = temperature_f_plot, y = endpoint_risk_plot, color = community_label)
    ) +
      geom_line(linewidth = 1) +
      facet_wrap(~ endpoint_label, ncol = 2, scales = "fixed") +
      scale_y_continuous(limits = c(0, 100), breaks = c(0, 25, 50, 75, 100)) +
      labs(
        title = "Community-specific heat-effect modification curves",
        subtitle = "Selected communities across all-cause and cardiovascular endpoints",
        x = "Temperature (°F)",
        y = "Endpoint risk (0-100)",
        color = "Community"
      ) +
      theme_pub(base_size = base_size) +
      theme(legend.position = "bottom")
    
    save_pub_plot(p_curves, "fig9_community_effect_modification_curves.png", width = 12, height = 9)
  }
}


# ------------------------------------------------------------------------------------------------
# FIGURE 10: Community-specific curves by source family
# ------------------------------------------------------------------------------------------------
if (have_temp_endpoint) {
  selected_communities_for_curves <- c(
    "AUSTIN",
    "NEAR NORTH SIDE",
    "ENGLEWOOD",
    "SOUTH SHORE",
    "NEAR WEST SIDE",
    "HYDE PARK",
    "AVONDALE",
    "GARFIELD RIDGE"
  )
  
  endpoint_risk_col <- get_endpoint_risk_col(temp_ep)
  
  curve_dat_all <- temp_ep %>%
    filter(
      community %in% selected_communities_for_curves
    ) %>%
    group_by(community, endpoint_key, temperature, source) %>%
    summarise(
      endpoint_risk_plot = mean(.data[[endpoint_risk_col]], na.rm = TRUE),
      outcome_label = dplyr::first(outcome_label),
      .groups = "drop"
    ) %>%
    mutate(
      endpoint_label = if_else(!is.na(outcome_label), outcome_label, pretty_endpoint_label(endpoint_key)),
      community_label = label_community_pretty(community),
      temperature_f_plot = if (temp_units == "C") c_to_f(temperature) else temperature
    )
  
  for (src in c("Mortality", "ED", "EMS")) {
    curve_dat_src <- curve_dat_all %>%
      filter(source == src)
    
    if (nrow(curve_dat_src) == 0) next
    
    top_endpoints_src <- curve_dat_src %>%
      group_by(endpoint_key, endpoint_label) %>%
      summarise(mean_risk = mean(endpoint_risk_plot, na.rm = TRUE), .groups = "drop") %>%
      slice_max(mean_risk, n = 6, with_ties = FALSE) %>%
      pull(endpoint_key)
    
    curve_dat_src <- curve_dat_src %>%
      filter(endpoint_key %in% top_endpoints_src)
    
    p_src <- ggplot(
      curve_dat_src,
      aes(x = temperature_f_plot, y = endpoint_risk_plot, color = community_label)
    ) +
      geom_line(linewidth = 1) +
      facet_wrap(~ endpoint_label, ncol = 2, scales = "fixed") +
      scale_y_continuous(limits = c(0, 100), breaks = c(0, 25, 50, 75, 100)) +
      labs(
        title = glue("{src} community-specific temperature response curves"),
        subtitle = "Selected communities and top-risk endpoints",
        x = "Temperature (°F)",
        y = "Endpoint risk (0-100)",
        color = "Community"
      ) +
      theme_pub(base_size = base_size) +
      theme(legend.position = "bottom")
    
    save_pub_plot(
      p_src,
      glue("fig10_{tolower(src)}_community_curves.png"),
      width = 12,
      height = 9
    )
  }
}

# ------------------------------------------------------------------------------------------------
# FIGURE X: Faceted endpoint maps at ~91°F
# Cause-specific endpoints only; all mapped communities retained
# ------------------------------------------------------------------------------------------------

# --- Build endpoint label lookup FIRST ---
endpoint_label_lookup <- temp_ep %>%
  filter(endpoint_key %in% endpoint_panel_keys) %>%
  group_by(endpoint_key) %>%
  summarise(
    outcome_label = dplyr::first(outcome_label),
    .groups = "drop"
  ) %>%
  mutate(
    endpoint_label = if_else(!is.na(outcome_label),
                             outcome_label,
                             pretty_endpoint_label(endpoint_key))
  ) %>%
  select(endpoint_key, endpoint_label)

# --- Build main dataset ---
endpoint_panel_dat <- temp_ep %>%
  filter(
    abs(temperature - focal_match_map) < 1e-8,
    endpoint_key %in% endpoint_panel_keys
  ) %>%
  group_by(community, endpoint_key) %>%
  summarise(
    endpoint_risk_plot = mean(.data[[endpoint_risk_col]], na.rm = TRUE),
    .groups = "drop"
  ) %>%
  tidyr::complete(
    community = unique(chi_sf$community),
    endpoint_key = endpoint_panel_keys,
    fill = list(endpoint_risk_plot = 0)
  ) %>%
  left_join(endpoint_label_lookup, by = "endpoint_key") %>%  # <-- critical fix
  mutate(
    endpoint_label = factor(endpoint_label, levels = endpoint_panel_levels)
  ) %>%
  left_join(
    chi_sf %>% select(community, geometry),
    by = "community"
  ) %>%
  st_as_sf() %>%
  st_make_valid()

# ------------------------------------------------------------------------------------------------
# FIGURE X: Faceted endpoint maps at ~91°F
# Cause-specific endpoints only; all mapped communities retained
# ------------------------------------------------------------------------------------------------
if (have_temp_endpoint) {
  focal_f_map <- 91
  
  focal_match_map <- target_temp_lookup %>%
    mutate(abs_diff = abs(target_f - focal_f_map)) %>%
    arrange(abs_diff) %>%
    slice(1) %>%
    pull(matched_native)
  
  endpoint_risk_col <- get_endpoint_risk_col(temp_ep)
  
  exclude_endpoints <- c("deaths", "ed_visits", "ems_calls")
  
  manual_endpoint_panel <- c(
    "death_gi", "death_renal", "death_injury",
    "ed_cvd", "ed_respiratory", "ed_renal",
    "ems_cvd", "ems_respiratory", "ems_injury"
  )
  
  endpoint_panel_keys <- manual_endpoint_panel[manual_endpoint_panel %in% unique(temp_ep$endpoint_key)]
  
  # build a named lookup vector instead of relying on join-created columns
  endpoint_label_lookup <- temp_ep %>%
    filter(endpoint_key %in% endpoint_panel_keys) %>%
    group_by(endpoint_key) %>%
    summarise(
      outcome_label = dplyr::first(outcome_label),
      .groups = "drop"
    ) %>%
    mutate(
      endpoint_label_clean = dplyr::if_else(
        !is.na(outcome_label),
        outcome_label,
        pretty_endpoint_label(endpoint_key)
      )
    )
  
  endpoint_label_vec <- endpoint_label_lookup$endpoint_label_clean
  names(endpoint_label_vec) <- endpoint_label_lookup$endpoint_key
  
  endpoint_panel_levels <- endpoint_label_lookup %>%
    filter(endpoint_key %in% endpoint_panel_keys) %>%
    pull(endpoint_label_clean) %>%
    unique()
  
  endpoint_panel_dat <- temp_ep %>%
    filter(
      abs(temperature - focal_match_map) < 1e-8,
      endpoint_key %in% endpoint_panel_keys
    ) %>%
    group_by(community, endpoint_key) %>%
    summarise(
      endpoint_risk_plot = mean(.data[[endpoint_risk_col]], na.rm = TRUE),
      .groups = "drop"
    ) %>%
    tidyr::complete(
      community = unique(chi_sf$community),
      endpoint_key = endpoint_panel_keys,
      fill = list(endpoint_risk_plot = 0)
    ) %>%
    mutate(
      endpoint_label = unname(endpoint_label_vec[endpoint_key]),
      endpoint_label = factor(endpoint_label, levels = endpoint_panel_levels)
    ) %>%
    left_join(
      chi_sf %>% st_drop_geometry() %>% select(community),
      by = "community"
    ) %>%
    left_join(
      chi_sf %>% select(community, geometry),
      by = "community"
    ) %>%
    st_as_sf() %>%
    st_make_valid()
  
  p_endpoint_map_panel <- ggplot(endpoint_panel_dat) +
    geom_sf(aes(fill = endpoint_risk_plot), color = "white", linewidth = 0.12) +
    coord_sf(expand = FALSE) +
    facet_wrap(~ endpoint_label, ncol = 3) +
    scale_fill_hvi_0_100("Endpoint risk\n(0-100)") +
    wide_colorbar +
    labs(
      title = "Cause-specific endpoint risk surfaces at ~91°F"
    ) +
    theme_pub_map_shared_legend(base_size = 10) +
    theme(
      strip.text = element_text(face = "bold", size = 10),
      legend.position = "bottom"
    )
  
  save_pub_plot(
    p_endpoint_map_panel,
    "fig_endpoint_map_panel_91f.png",
    width = 12.5,
    height = 10
  )
}
