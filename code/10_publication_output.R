# ================================================================================================
# HVI 2.0 | Publication Outputs
# PI: Peter Graffy
#
# Purpose:
#   Build publication-ready maps, figures, and tables from the HVI 2.0 pipeline.
#
# Expected objects already in environment (as available):
#   ca_year_endpoint_structural_hvi
#   ca_year_overall_structural_hvi
#   ca_day_endpoint_risk
#   ca_day_overall_operational_hvi
#   temperature_grid_endpoint_risk
#   temperature_grid_overall_risk
#   endpoint_model_performance
#   endpoint_weights
#   hvi_endpoint_metadata
#   cas   OR   commareas   OR other sf object with community geometry
#
# Outputs:
#   outputs/publication/
#     figures/
#     tables/
#     supplements/
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
})

source(file.path(if (dir.exists("code")) "code" else ".", "00_config.R"))

# ------------------------------------------------------------------------------------------------
# Directories
# ------------------------------------------------------------------------------------------------
pub_dir <- HVI_PATHS$private_outputs$publication_outputs
fig_dir <- file.path(pub_dir, "figures")
tab_dir <- file.path(pub_dir, "tables")
sup_dir <- file.path(pub_dir, "supplements")

dir.create(pub_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(sup_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------------------------------------
# User-settable options
# ------------------------------------------------------------------------------------------------
target_temp_f <- 92
target_year   <- 2022
top_n_labels  <- 8

selected_endpoint_maps <- c(
  "deaths", "death_cvd", "death_injury", "ed_renal", "ems_injury", "ems_calls", "ems_cvd"
)

selected_endpoint_curves <- c(
  "deaths", "death_cvd", "death_injury", "ed_renal", "ems_injury", "ems_calls", "ems_cvd"
)

selected_communities_for_curves <- c(
  "AUSTIN", "WEST TOWN", "NEAR NORTH SIDE", "ENGLEWOOD", "SOUTH SHORE", "NORTH LAWNDALE"
)

# detect temperature units from temp grid once available later
target_temp_native <- NULL

# ------------------------------------------------------------------------------------------------
# Theme
# ------------------------------------------------------------------------------------------------
theme_pub <- function(base_size = 11) {
  theme_minimal(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      axis.title = element_text(face = "bold"),
      plot.title = element_text(face = "bold", size = base_size + 2),
      plot.subtitle = element_text(size = base_size),
      plot.caption = element_text(size = base_size - 1, color = "grey35"),
      strip.text = element_text(face = "bold"),
      legend.title = element_text(face = "bold"),
      legend.position = "right"
    )
}

# ------------------------------------------------------------------------------------------------
# Helpers
# ------------------------------------------------------------------------------------------------
normalize_community_key <- function(x) {
  x %>%
    as.character() %>%
    stringr::str_squish() %>%
    stringr::str_to_upper()
}

label_community_pretty <- function(x) {
  x %>%
    as.character() %>%
    stringr::str_replace_all("_", " ") %>%
    stringr::str_squish() %>%
    stringr::str_to_title()
}

rescale_0_100 <- function(x) {
  if (all(is.na(x))) return(rep(NA_real_, length(x)))
  rng <- range(x, na.rm = TRUE)
  if (diff(rng) == 0) return(rep(50, length(x)))
  scales::rescale(x, to = c(0, 100), from = rng)
}

safe_write_csv <- function(df, path) {
  readr::write_csv(df, path, na = "")
}

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

fmt_num <- function(x, digits = 2) formatC(x, digits = digits, format = "f", big.mark = ",")

infer_map_sf <- function() {
  if (exists("cas", inherits = TRUE)) {
    obj <- get("cas", inherits = TRUE)
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
  
  candidate_cols <- c("community", "CA_NAME", "ca_name", "commarea", "community_area", "area_name")
  candidate_cols <- candidate_cols[candidate_cols %in% nm]
  
  if (length(candidate_cols) == 0) {
    stop("Map object does not have a recognizable community identifier column.")
  }
  
  # Prefer columns that look like character community names rather than numeric IDs
  score_col <- purrr::map_dbl(candidate_cols, function(col) {
    x <- sf_obj[[col]]
    if (is.factor(x)) x <- as.character(x)
    if (!is.character(x)) return(-Inf)
    
    x2 <- stringr::str_to_upper(stringr::str_squish(x))
    mean(stringr::str_detect(x2, "[A-Z]"), na.rm = TRUE)
  })
  
  best_col <- candidate_cols[which.max(score_col)]
  
  sf_obj %>%
    mutate(community = as.character(.data[[best_col]]))
}

get_source_clean <- function(x) {
  case_when(
    str_to_lower(x) %in% c("mortality", "death", "deaths") ~ "Mortality",
    str_to_lower(x) == "ed" ~ "ED",
    str_to_lower(x) == "ems" ~ "EMS",
    TRUE ~ as.character(x)
  )
}

make_choropleth <- function(sf_df, value_col, title, subtitle = NULL, fill_label = NULL,
                            palette_option = "C", na_fill = "grey90") {
  sf_df <- sf_df %>%
    sf::st_as_sf() %>%
    sf::st_make_valid()
  
  ggplot(sf_df) +
    geom_sf(aes(fill = .data[[value_col]]), color = "white", linewidth = 0.15) +
    coord_sf(default_crs = NULL, expand = FALSE) +
    scale_fill_viridis_c(
      option = palette_option,
      na.value = na_fill,
      labels = scales::label_number(accuracy = 1),
      name = fill_label %||% value_col
    ) +
    labs(
      title = title,
      subtitle = subtitle,
      caption = "Chicago community areas"
    ) +
    theme_pub() +
    theme(
      axis.text = element_blank(),
      axis.title = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank()
    )
}

make_categorical_map <- function(sf_df, value_col, title, subtitle = NULL) {
  sf_df <- sf_df %>%
    sf::st_as_sf() %>%
    sf::st_make_valid()
  
  ggplot(sf_df) +
    geom_sf(aes(fill = .data[[value_col]]), color = "white", linewidth = 0.15) +
    coord_sf(default_crs = NULL, expand = FALSE) +
    scale_fill_viridis_d(na.value = "grey90", begin = 0.1, end = 0.9, name = NULL) +
    labs(
      title = title,
      subtitle = subtitle,
      caption = "Chicago community areas"
    ) +
    theme_pub() +
    theme(
      axis.text = element_blank(),
      axis.title = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank(),
      legend.position = "bottom"
    )
}

get_first_match <- function(x, patterns) {
  hits <- names(x)[stringr::str_detect(names(x), patterns)]
  if (length(hits) == 0) return(NA_character_)
  hits[1]
}

get_temp_col <- function(df) {
  if ("temperature" %in% names(df)) return("temperature")
  get_first_match(df, "temp_value|temp_f|temperature")
}

get_endpoint_risk_col <- function(df) {
  get_first_match(df, "^endpoint_risk_0_100$|endpoint.*risk.*0_100|risk_0_100")
}

get_overall_risk_col <- function(df) {
  get_first_match(df, "^overall_risk_0_100$|overall.*risk.*0_100|overall_risk|risk_0_100")
}

get_structural_col <- function(df) {
  get_first_match(df, "overall.*0_100|hvi_0_100|vulnerability_0_100|structural.*0_100")
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

# ------------------------------------------------------------------------------------------------
# Validate required objects
# ------------------------------------------------------------------------------------------------
have_structural_endpoint <- exists("ca_year_endpoint_structural_hvi", inherits = TRUE)
have_structural_overall  <- exists("ca_year_overall_structural_hvi", inherits = TRUE)
have_daily_endpoint      <- exists("ca_day_endpoint_risk", inherits = TRUE)
have_daily_overall       <- exists("ca_day_overall_operational_hvi", inherits = TRUE)

# support both naming conventions from 09c
temp_endpoint_obj <- dplyr::case_when(
  exists("temperature_grid_endpoint_risk", inherits = TRUE) ~ "temperature_grid_endpoint_risk",
  exists("temp_grid_endpoint_risk", inherits = TRUE) ~ "temp_grid_endpoint_risk",
  TRUE ~ NA_character_
)

temp_overall_obj <- dplyr::case_when(
  exists("temperature_grid_overall_risk", inherits = TRUE) ~ "temperature_grid_overall_risk",
  exists("temp_grid_overall_risk", inherits = TRUE) ~ "temp_grid_overall_risk",
  TRUE ~ NA_character_
)

temp_family_obj <- dplyr::case_when(
  exists("temperature_grid_family_risk", inherits = TRUE) ~ "temperature_grid_family_risk",
  exists("temp_grid_family_risk", inherits = TRUE) ~ "temp_grid_family_risk",
  TRUE ~ NA_character_
)

have_temp_endpoint <- !is.na(temp_endpoint_obj)
have_temp_overall  <- !is.na(temp_overall_obj)
have_temp_family   <- !is.na(temp_family_obj)

have_perf    <- exists("endpoint_model_performance", inherits = TRUE)
have_weights <- exists("endpoint_weights", inherits = TRUE)
have_meta    <- exists("hvi_endpoint_metadata", inherits = TRUE)

if (!have_structural_endpoint && !have_daily_endpoint && !have_temp_endpoint) {
  stop("No publication-ready model outputs found in the environment.")
}

chi_sf <- infer_map_sf() %>%
  standardize_map_key() %>%
  mutate(
    community = normalize_community_key(community)
  ) %>%
  sf::st_as_sf() %>%
  sf::st_make_valid()

# ------------------------------------------------------------------------------------------------
# Pull objects
# ------------------------------------------------------------------------------------------------
if (have_structural_endpoint) structural_ep <- get("ca_year_endpoint_structural_hvi", inherits = TRUE)
if (have_structural_overall)  structural_overall <- get("ca_year_overall_structural_hvi", inherits = TRUE)
if (have_daily_endpoint)      daily_ep <- get("ca_day_endpoint_risk", inherits = TRUE)
if (have_daily_overall)       daily_overall <- get("ca_day_overall_operational_hvi", inherits = TRUE)

if (have_temp_endpoint) temp_ep <- get(temp_endpoint_obj, inherits = TRUE)
if (have_temp_overall)  temp_overall <- get(temp_overall_obj, inherits = TRUE)
if (have_temp_family)   temp_family <- get(temp_family_obj, inherits = TRUE)

if (have_perf)    perf <- get("endpoint_model_performance", inherits = TRUE)
if (have_weights) ep_weights <- get("endpoint_weights", inherits = TRUE)
if (have_meta)    ep_meta <- get("hvi_endpoint_metadata", inherits = TRUE)

# ------------------------------------------------------------------------------------------------
# Standardize common fields
# ------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------
# Standardize common fields
# ------------------------------------------------------------------------------------------------
if (have_structural_endpoint) {
  structural_ep <- structural_ep %>%
    mutate(
      community = normalize_community_key(community),
      year = harmonize_year_values(year)
    )
}

if (have_structural_overall) {
  structural_overall <- structural_overall %>%
    mutate(
      community = normalize_community_key(community),
      year = harmonize_year_values(year)
    )
}

if (have_daily_endpoint) {
  daily_ep <- daily_ep %>%
    mutate(
      community = normalize_community_key(community),
      date = as.Date(date)
    )
}

if (have_daily_overall) {
  daily_overall <- daily_overall %>%
    mutate(
      community = normalize_community_key(community),
      date = as.Date(date)
    )
}

if (have_temp_endpoint) {
  temp_ep <- temp_ep %>%
    mutate(
      community = normalize_community_key(community),
      year = harmonize_year_values(year)
    )
  
  if ("temperature" %in% names(temp_ep)) {
    temp_ep <- temp_ep %>% mutate(temperature = as.numeric(temperature))
  } else if ("temp_value" %in% names(temp_ep)) {
    temp_ep <- temp_ep %>% mutate(temperature = as.numeric(temp_value))
  } else if ("temp_f" %in% names(temp_ep)) {
    temp_ep <- temp_ep %>% mutate(temperature = as.numeric(temp_f))
  } else {
    stop("temp_ep is missing a temperature column. Expected one of: temperature, temp_value, temp_f")
  }
}

if (have_temp_overall) {
  temp_overall <- temp_overall %>%
    mutate(
      community = normalize_community_key(community),
      year = harmonize_year_values(year)
    )
  
  if ("temperature" %in% names(temp_overall)) {
    temp_overall <- temp_overall %>% mutate(temperature = as.numeric(temperature))
  } else if ("temp_value" %in% names(temp_overall)) {
    temp_overall <- temp_overall %>% mutate(temperature = as.numeric(temp_value))
  } else if ("temp_f" %in% names(temp_overall)) {
    temp_overall <- temp_overall %>% mutate(temperature = as.numeric(temp_f))
  } else {
    stop("temp_overall is missing a temperature column. Expected one of: temperature, temp_value, temp_f")
  }
}

if (have_temp_family) {
  temp_family <- temp_family %>%
    mutate(
      community = normalize_community_key(community),
      year = harmonize_year_values(year)
    )
  
  if ("temperature" %in% names(temp_family)) {
    temp_family <- temp_family %>% mutate(temperature = as.numeric(temperature))
  } else if ("temp_value" %in% names(temp_family)) {
    temp_family <- temp_family %>% mutate(temperature = as.numeric(temp_value))
  } else if ("temp_f" %in% names(temp_family)) {
    temp_family <- temp_family %>% mutate(temperature = as.numeric(temp_f))
  } else {
    stop("temp_family is missing a temperature column. Expected one of: temperature, temp_value, temp_f")
  }
}


# ------------------------------------------------------------------------------------------------
# Temperature target in native units
# ------------------------------------------------------------------------------------------------
if (have_temp_overall) {
  temp_units <- detect_temperature_units(temp_overall$temperature)
} else if (have_temp_endpoint) {
  temp_units <- detect_temperature_units(temp_ep$temperature)
} else {
  temp_units <- "unknown"
}

if (temp_units == "C") {
  target_temp_native <- f_to_c(target_temp_f)
} else {
  target_temp_native <- target_temp_f
}

message("Detected temperature units in grid: ", temp_units)
message("Requested focal temperature: ", target_temp_f, " F")
message("Native focal temperature used for filtering: ", round(target_temp_native, 2))
# ------------------------------------------------------------------------------------------------
# Join endpoint metadata if useful
# ------------------------------------------------------------------------------------------------
if (have_meta && have_daily_endpoint && !"outcome_label" %in% names(daily_ep)) {
  daily_ep <- daily_ep %>%
    left_join(
      ep_meta %>%
        select(any_of(c("endpoint_key", "outcome_label", "source", "domain"))),
      by = "endpoint_key"
    )
}

if (have_meta && have_temp_endpoint && !"outcome_label" %in% names(temp_ep)) {
  temp_ep <- temp_ep %>%
    left_join(
      ep_meta %>%
        select(any_of(c("endpoint_key", "outcome_label", "source", "domain"))),
      by = "endpoint_key"
    )
}

# ------------------------------------------------------------------------------------------------
# TABLE 1: Endpoint model performance
# ------------------------------------------------------------------------------------------------
if (have_perf) {
  perf_summary <- perf %>%
    clean_names() %>%
    mutate(
      endpoint_key = as.character(endpoint_key),
      cv_type = as.character(cv_type)
    )
  
  if (have_meta) {
    perf_summary <- perf_summary %>%
      left_join(
        ep_meta %>%
          select(any_of(c("endpoint_key", "outcome_label", "source", "domain"))),
        by = "endpoint_key",
        suffix = c("", "_meta")
      ) %>%
      mutate(
        outcome_label = coalesce(outcome_label, outcome_label_meta),
        source = coalesce(source, source_meta),
        domain = coalesce(domain, domain_meta)
      ) %>%
      select(-any_of(c("outcome_label_meta", "source_meta", "domain_meta"))) %>%
      mutate(
        source = if ("source" %in% names(.)) get_source_clean(source) else NA_character_
      )
  }
  
  safe_write_csv(
    perf_summary,
    file.path(tab_dir, "table_model_performance_summary.csv")
  )
  
  sort_col <- dplyr::case_when(
    "mean_cor_comm_spearman" %in% names(perf_summary) ~ "mean_cor_comm_spearman",
    "mean_cor_comm_pearson" %in% names(perf_summary) ~ "mean_cor_comm_pearson",
    TRUE ~ NA_character_
  )
  
  perf_gt_dat <- perf_summary
  if (!is.na(sort_col)) {
    perf_gt_dat <- perf_gt_dat %>%
      arrange(cv_type, desc(.data[[sort_col]]))
  } else {
    perf_gt_dat <- perf_gt_dat %>%
      arrange(cv_type, endpoint_key)
  }
  
  perf_gt <- perf_gt_dat %>%
    mutate(across(where(is.numeric), ~ round(.x, 3))) %>%
    gt(groupname_col = "cv_type") %>%
    tab_header(
      title = md("**Table. Cross-validated model performance by endpoint**"),
      subtitle = "Spatial and temporal validation summaries"
    )
  
  perf_label_map <- c(
    endpoint_key = "Endpoint key",
    outcome_label = "Outcome",
    source = "Source",
    domain = "Domain",
    folds = "Folds",
    mean_rmse = "RMSE",
    mean_mae = "MAE",
    mean_poisson_deviance = "Poisson deviance",
    mean_cor_daily_spearman = "Daily Spearman r",
    mean_cor_comm_spearman = "Community Spearman r",
    mean_cor_comm_pearson = "Community Pearson r"
  )
  
  perf_gt <- perf_gt %>%
    cols_label(.list = perf_label_map[names(perf_label_map) %in% names(perf_summary)]) %>%
    fmt_number(columns = where(is.numeric), decimals = 3)
  
  gtsave(
    perf_gt,
    file.path(tab_dir, "table_model_performance_summary.html")
  )
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
      title = md("**Table. Operational endpoint weights**"),
      subtitle = "Source-priority and performance-weighted contributions"
    ) %>%
    fmt_number(columns = where(is.numeric), decimals = 3)
  
  gtsave(ep_weights_gt, file.path(tab_dir, "table_endpoint_weights.html"))
}

# ------------------------------------------------------------------------------------------------
# TABLE 3: Community rankings at target temperature
# ------------------------------------------------------------------------------------------------
if (have_temp_overall) {
  temp_col <- get_temp_col(temp_overall)
  risk_col <- get_overall_risk_col(temp_overall)
  
  if (!is.na(temp_col) && !is.na(risk_col)) {
    target_temp_match <- closest_available_temperature(temp_overall[[temp_col]], target_temp_native)
    
    top_communities_92 <- temp_overall %>%
      filter(abs(.data[[temp_col]] - target_temp_match) < 1e-8) %>%
      arrange(desc(.data[[risk_col]])) %>%
      mutate(rank = row_number()) %>%
      select(rank, community, all_of(temp_col), all_of(risk_col), any_of(c("dominant_endpoint", "dominant_driver"))) %>%
      rename(temperature_native = all_of(temp_col), overall_risk_0_100 = all_of(risk_col)) %>%
      mutate(
        temperature_f = if (temp_units == "C") c_to_f(temperature_native) else temperature_native
      )
    
    safe_write_csv(top_communities_92, file.path(tab_dir, glue("table_top_communities_{target_temp_f}f.csv")))
    
    top_communities_gt <- top_communities_92 %>%
      slice_head(n = 20) %>%
      gt() %>%
      tab_header(
        title = md(glue("**Table. Community rankings at ~{target_temp_f}°F**")),
        subtitle = glue("Overall heat-health risk by community area (native grid value = {round(target_temp_match, 1)}°{temp_units})")
      ) %>%
      fmt_number(columns = c(temperature_native, temperature_f, overall_risk_0_100), decimals = 1)
    
    gtsave(top_communities_gt, file.path(tab_dir, glue("table_top_communities_{target_temp_f}f.html")))
  }
}

# ------------------------------------------------------------------------------------------------
# FIGURE 1: Overall structural vulnerability map
# ------------------------------------------------------------------------------------------------
if (have_structural_overall) {
  structural_overall_map_df <- chi_sf %>%
    left_join(
      structural_overall %>%
        filter(year == as.character(target_year)),
      by = "community"
    ) %>%
    sf::st_as_sf() %>%
    sf::st_make_valid()
  
  score_col <- get_structural_col(structural_overall_map_df)
  
  if (!is.na(score_col)) {
    message("Figure 1 matched non-missing rows: ", sum(!is.na(structural_overall_map_df[[score_col]])))
    
    p_structural <- make_choropleth(
      structural_overall_map_df,
      value_col = score_col,
      title = "Structural heat vulnerability across Chicago",
      subtitle = glue("Community-level overall structural HVI, {target_year}"),
      fill_label = "Structural HVI\n(0-100)"
    )
    
    save_pub_plot(p_structural, "fig1_structural_overall_hvi_map.png", width = 7.5, height = 7.2)
  } else {
    message("No structural overall 0-100 score column detected for map figure.")
  }
}

# ------------------------------------------------------------------------------------------------
# FIGURE 2: Overall dynamic risk map at target temperature
# ------------------------------------------------------------------------------------------------
if (have_temp_overall) {
  temp_col <- get_temp_col(temp_overall)
  risk_col <- get_overall_risk_col(temp_overall)
  dominant_col <- get_first_match(temp_overall, "dominant_endpoint")
  
  if (!is.na(temp_col) && !is.na(risk_col)) {
    target_temp_match <- closest_available_temperature(temp_overall[[temp_col]], target_temp_native)
    
    overall_92_dat <- temp_overall %>%
      filter(abs(.data[[temp_col]] - target_temp_match) < 1e-8) %>%
      group_by(community) %>%
      summarise(
        across(any_of(c(risk_col, "dominant_endpoint", "dominant_endpoint_excess", "dominant_endpoint_risk_0_100")), ~ dplyr::first(.x)),
        .groups = "drop"
      )
    
    overall_92_map_df <- chi_sf %>%
      left_join(overall_92_dat, by = "community") %>%
      sf::st_as_sf() %>%
      sf::st_make_valid()
    
    message("Figure 2 matched non-missing rows: ", sum(!is.na(overall_92_map_df[[risk_col]])))
    
    p_dynamic <- make_choropleth(
      overall_92_map_df,
      value_col = risk_col,
      title = glue("Overall heat-health risk at ~{target_temp_f}°F"),
      subtitle = glue("Temperature-conditional dynamic HVI across community areas (grid value = {round(target_temp_match, 1)}°{temp_units})"),
      fill_label = "Overall risk\n(0-100)"
    )
    
    if (!is.na(dominant_col)) {
      p_dominant <- make_categorical_map(
        overall_92_map_df,
        value_col = dominant_col,
        title = "Dominant endpoint at high heat",
        subtitle = glue("Leading contributor to overall risk at ~{target_temp_f}°F")
      )
      
      fig2 <- p_dynamic + p_dominant + patchwork::plot_layout(ncol = 2)
      save_pub_plot(fig2, glue("fig2_overall_dynamic_risk_{target_temp_f}f.png"), width = 13, height = 6.8)
    } else {
      save_pub_plot(p_dynamic, glue("fig2_overall_dynamic_risk_{target_temp_f}f.png"), width = 7.5, height = 7.2)
    }
  } else {
    message("Figure 2 skipped: missing temperature or overall risk column in temp_overall.")
  }
}

# ------------------------------------------------------------------------------------------------
# FIGURE 3: Endpoint-specific map panel at target temperature
# ------------------------------------------------------------------------------------------------
if (have_temp_endpoint) {
  temp_col <- get_temp_col(temp_ep)
  risk_col <- get_endpoint_risk_col(temp_ep)
  
  if (!is.na(temp_col) && !is.na(risk_col)) {
    target_temp_match <- closest_available_temperature(temp_ep[[temp_col]], target_temp_native)
    
    endpoint_map_dat <- chi_sf %>%
      left_join(
        temp_ep %>%
          filter(abs(.data[[temp_col]] - target_temp_match) < 1e-8, endpoint_key %in% selected_endpoint_maps) %>%
          mutate(community = as.character(community)),
        by = "community"
      ) %>%
      sf::st_as_sf() %>%
      sf::st_make_valid() %>%
      mutate(
        facet_label = if ("outcome_label" %in% names(.)) outcome_label else endpoint_key,
        facet_label = forcats::fct_inorder(facet_label)
      ) %>%
      filter(!is.na(endpoint_key))
    
    if (nrow(endpoint_map_dat) > 0) {
      p_endpoint_panel <- ggplot(endpoint_map_dat) +
        geom_sf(aes(fill = .data[[risk_col]]), color = "white", linewidth = 0.12) +
        coord_sf(default_crs = NULL, expand = FALSE) +
        scale_fill_viridis_c(option = "C", name = "Risk\n(0-100)", labels = scales::label_number(accuracy = 1)) +
        facet_wrap(~ facet_label, ncol = 3) +
        labs(
          title = glue("Endpoint-specific heat-health risk at ~{target_temp_f}°F"),
          subtitle = glue("Illustrative endpoint panel across Chicago community areas (grid value = {round(target_temp_match, 1)}°{temp_units})"),
          caption = "Each panel shows a separate modeled endpoint"
        ) +
        theme_pub() +
        theme(
          axis.text = element_blank(),
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          panel.grid = element_blank()
        )
      
      save_pub_plot(p_endpoint_panel, glue("fig3_endpoint_map_panel_{target_temp_f}f.png"), width = 12.5, height = 10)
    } else {
      message("Figure 3 skipped: no endpoint rows remained after filtering.")
    }
  }
}

# ------------------------------------------------------------------------------------------------
# FIGURE 4: Temperature-response curves for selected communities and endpoints
# ------------------------------------------------------------------------------------------------
if (have_temp_endpoint) {
  temp_col <- get_temp_col(temp_ep)
  risk_col <- get_endpoint_risk_col(temp_ep)
  
  if (!is.na(temp_col) && !is.na(risk_col)) {
    curve_dat <- temp_ep %>%
      filter(
        endpoint_key %in% selected_endpoint_curves,
        community %in% selected_communities_for_curves
      ) %>%
      group_by(community, endpoint_key, outcome_label, .data[[temp_col]]) %>%
      summarise(
        endpoint_risk_plot = mean(.data[[risk_col]], na.rm = TRUE),
        .groups = "drop"
      ) %>%
      mutate(
        endpoint_label = if ("outcome_label" %in% names(.)) outcome_label else endpoint_key,
        endpoint_label = forcats::fct_inorder(endpoint_label),
        community_label = label_community_pretty(community),
        temperature_f_plot = if (temp_units == "C") c_to_f(.data[[temp_col]]) else .data[[temp_col]]
      )
    
    if (nrow(curve_dat) > 0) {
      p_curves <- ggplot(
        curve_dat,
        aes(x = temperature_f_plot, y = endpoint_risk_plot, group = community_label)
      ) +
        geom_line(aes(color = community_label), linewidth = 0.9, alpha = 0.95) +
        facet_wrap(~ endpoint_label, ncol = 2, scales = "fixed") +
        labs(
          title = "Temperature-response risk curves by community",
          subtitle = "Selected communities and selected endpoints",
          x = "Temperature (°F)",
          y = "Endpoint risk (0-100)",
          color = "Community"
        ) +
        theme_pub() +
        theme(legend.position = "bottom")
      
      save_pub_plot(p_curves, "fig4_temperature_response_curves.png", width = 12, height = 9)
    }
  }
}

# ------------------------------------------------------------------------------------------------
# FIGURE 5: Heatmap of endpoint risk at target temperature
# ------------------------------------------------------------------------------------------------
if (have_temp_endpoint) {
  temp_col <- get_temp_col(temp_ep)
  risk_col <- get_endpoint_risk_col(temp_ep)
  
  if (!is.na(temp_col) && !is.na(risk_col)) {
    target_temp_match <- closest_available_temperature(temp_ep[[temp_col]], target_temp_native)
    
    heatmap_dat <- temp_ep %>%
      filter(abs(.data[[temp_col]] - target_temp_match) < 1e-8) %>%
      group_by(community) %>%
      summarise(total_risk = sum(.data[[risk_col]], na.rm = TRUE), .groups = "drop") %>%
      arrange(desc(total_risk)) %>%
      slice_head(n = 77) %>%
      select(community)
    
    heatmap_plot_dat <- temp_ep %>%
      filter(abs(.data[[temp_col]] - target_temp_match) < 1e-8, community %in% heatmap_dat$community) %>%
      mutate(
        endpoint_label = if ("outcome_label" %in% names(.)) outcome_label else endpoint_key
      ) %>%
      group_by(community) %>%
      mutate(comm_total = sum(.data[[risk_col]], na.rm = TRUE)) %>%
      ungroup() %>%
      mutate(
        community = forcats::fct_reorder(community, comm_total, .desc = TRUE),
        endpoint_label = forcats::fct_reorder(endpoint_label, .data[[risk_col]], .fun = mean, .desc = TRUE)
      )
    
    if (nrow(heatmap_plot_dat) > 0) {
      p_heatmap <- ggplot(heatmap_plot_dat, aes(x = endpoint_label, y = community, fill = .data[[risk_col]])) +
        geom_tile(color = "white", linewidth = 0.15) +
        scale_fill_viridis_c(option = "C", name = "Risk\n(0-100)") +
        labs(
          title = glue("Endpoint risk architecture at ~{target_temp_f}°F"),
          subtitle = glue("Top-risk communities across all modeled endpoints (grid value = {round(target_temp_match, 1)}°{temp_units})"),
          x = NULL,
          y = NULL
        ) +
        theme_pub() +
        theme(
          axis.text.x = element_text(angle = 45, hjust = 1),
          panel.grid = element_blank()
        )
      
      save_pub_plot(p_heatmap, glue("fig5_endpoint_heatmap_{target_temp_f}f.png"), width = 12, height = 8)
    }
  }
}
# ------------------------------------------------------------------------------------------------
# FIGURE 6: Structural vs dynamic risk scatter
# ------------------------------------------------------------------------------------------------
if (have_structural_overall && have_temp_overall) {
  structural_col <- get_structural_col(structural_overall)
  temp_col <- get_temp_col(temp_overall)
  risk_col <- get_overall_risk_col(temp_overall)
  
  if (!is.na(structural_col) && !is.na(temp_col) && !is.na(risk_col)) {
    target_temp_match <- closest_available_temperature(temp_overall[[temp_col]], target_temp_native)
    
    dynamic_dat <- temp_overall %>%
      filter(abs(.data[[temp_col]] - target_temp_match) < 1e-8) %>%
      group_by(community) %>%
      summarise(
        dynamic_score = mean(.data[[risk_col]], na.rm = TRUE),
        .groups = "drop"
      )
    
    scatter_dat <- structural_overall %>%
      filter(year == as.character(target_year)) %>%
      select(community, structural_score = all_of(structural_col)) %>%
      left_join(dynamic_dat, by = "community") %>%
      tidyr::drop_na() %>%
      mutate(
        community_label = label_community_pretty(community)
      )
    
    if (nrow(scatter_dat) > 0) {
      top_labs <- scatter_dat %>%
        mutate(label_score = structural_score + dynamic_score) %>%
        slice_max(label_score, n = 5)
      
      p_scatter <- ggplot(scatter_dat, aes(x = structural_score, y = dynamic_score)) +
        geom_point(size = 2.8, alpha = 0.85) +
        geom_smooth(method = "lm", se = FALSE, linewidth = 0.8, linetype = 2, color = "#3B6CF6") +
        geom_text(
          data = top_labs,
          aes(label = community_label),
          size = 3.2,
          nudge_y = 0.35,
          check_overlap = TRUE
        ) +
        labs(
          title = "Structural vulnerability versus dynamic risk",
          subtitle = glue("{target_year} structural HVI vs overall heat-health risk at ~{target_temp_f}°F"),
          x = "Structural HVI (0-100)",
          y = "Dynamic heat-health risk (0-100)"
        ) +
        theme_pub()
      
      save_pub_plot(p_scatter, "fig6_structural_vs_dynamic_scatter.png", width = 8.5, height = 6.5)
    } else {
      message("Figure 6 skipped: no rows after joining structural and dynamic data.")
    }
  }
}

# ------------------------------------------------------------------------------------------------
# FIGURE 7: Model performance figure
# ------------------------------------------------------------------------------------------------
if (have_perf) {
  perf_fig_dat <- perf %>%
    clean_names() %>%
    mutate(
      endpoint_key = as.character(endpoint_key),
      cv_type = as.character(cv_type)
    )
  
  # If metadata is available, merge carefully without duplicating existing columns
  if (have_meta) {
    perf_fig_dat <- perf_fig_dat %>%
      left_join(
        ep_meta %>%
          select(any_of(c("endpoint_key", "outcome_label", "source"))),
        by = "endpoint_key",
        suffix = c("", "_meta")
      ) %>%
      mutate(
        outcome_label = coalesce(outcome_label, outcome_label_meta),
        source = coalesce(source, source_meta)
      ) %>%
      select(-any_of(c("outcome_label_meta", "source_meta"))) %>%
      mutate(
        endpoint_label = if_else(!is.na(outcome_label), outcome_label, endpoint_key),
        source = if ("source" %in% names(.)) get_source_clean(source) else "Other"
      )
  } else {
    perf_fig_dat <- perf_fig_dat %>%
      mutate(
        endpoint_label = endpoint_key,
        source = "Other"
      )
  }
  
  # Choose the best available performance column
  perf_col <- case_when(
    "mean_cor_comm_spearman" %in% names(perf_fig_dat) ~ "mean_cor_comm_spearman",
    "mean_cor_comm_pearson" %in% names(perf_fig_dat) ~ "mean_cor_comm_pearson",
    TRUE ~ NA_character_
  )
  
  if (is.na(perf_col)) {
    stop(
      "No community-level performance column found in `perf`. Available columns are: ",
      paste(names(perf_fig_dat), collapse = ", ")
    )
  }
  
  perf_fig_dat <- perf_fig_dat %>%
    mutate(
      endpoint_label = forcats::fct_reorder(endpoint_label, .data[[perf_col]])
    )
  
  p_perf <- ggplot(perf_fig_dat, aes(x = .data[[perf_col]], y = endpoint_label)) +
    geom_col() +
    facet_wrap(~ cv_type, ncol = 2, scales = "free_y") +
    labs(
      title = "Cross-validated performance across endpoints",
      subtitle = "Community-level correlation between observed and predicted burden",
      x = ifelse(
        perf_col == "mean_cor_comm_spearman",
        "Mean community-level Spearman correlation",
        "Mean community-level Pearson correlation"
      ),
      y = NULL
    ) +
    theme_pub()
  
  save_pub_plot(
    p_perf,
    "fig7_model_performance.png",
    width = 12,
    height = 9
  )
}

# ------------------------------------------------------------------------------------------------
# Supplemental tables
# ------------------------------------------------------------------------------------------------
if (have_daily_endpoint) {
  daily_summary <- daily_ep %>%
    group_by(endpoint_key, community) %>%
    summarise(
      mean_predicted = mean(predicted_count, na.rm = TRUE),
      mean_excess = mean(excess_events, na.rm = TRUE),
      mean_rr = mean(relative_risk, na.rm = TRUE),
      .groups = "drop"
    )
  
  safe_write_csv(daily_summary, file.path(sup_dir, "supplement_daily_endpoint_summary.csv"))
}

if (have_temp_endpoint) {
  safe_write_csv(temp_ep, file.path(sup_dir, "supplement_temperature_grid_endpoint_risk.csv"))
}

if (have_temp_overall) {
  safe_write_csv(temp_overall, file.path(sup_dir, "supplement_temperature_grid_overall_risk.csv"))
}

if (have_temp_family) {
  safe_write_csv(temp_family, file.path(sup_dir, "supplement_temperature_grid_family_risk.csv"))
}

if (have_structural_endpoint) {
  safe_write_csv(structural_ep, file.path(sup_dir, "supplement_structural_endpoint_hvi.csv"))
}

if (have_structural_overall) {
  safe_write_csv(structural_overall, file.path(sup_dir, "supplement_structural_overall_hvi.csv"))
}

# ------------------------------------------------------------------------------------------------
# Manuscript-ready summary text snippets
# ------------------------------------------------------------------------------------------------
summary_lines <- c(
  glue("Publication outputs generated on {Sys.Date()}."),
  glue("Target temperature for focal figures: {target_temp_f}°F."),
  glue("Target year for structural maps: {target_year}."),
  if (have_perf) glue("Model performance summary table written to: {file.path(tab_dir, 'table_model_performance_summary.csv')}") else NULL,
  if (have_temp_overall) glue("Overall temperature-grid risk table written to: {file.path(sup_dir, 'supplement_temperature_grid_overall_risk.csv')}") else NULL,
  if (have_temp_endpoint) glue("Endpoint temperature-grid risk table written to: {file.path(sup_dir, 'supplement_temperature_grid_endpoint_risk.csv')}") else NULL,
  if (have_temp_family) glue("Family temperature-grid risk table written to: {file.path(sup_dir, 'supplement_temperature_grid_family_risk.csv')}") else NULL
)

writeLines(summary_lines, con = file.path(pub_dir, "README_publication_outputs.txt"))

message("Done. Publication outputs written to: ", normalizePath(pub_dir))


message("Communities in chi_sf: ", dplyr::n_distinct(chi_sf$community))
if (have_structural_overall) message("Communities in structural_overall: ", dplyr::n_distinct(structural_overall$community))
if (have_temp_overall) message("Communities in temp_overall: ", dplyr::n_distinct(temp_overall$community))
if (have_temp_endpoint) message("Communities in temp_ep: ", dplyr::n_distinct(temp_ep$community))

if (have_structural_overall) {
  message("Structural join overlap: ",
          sum(unique(chi_sf$community) %in% unique(structural_overall$community)))
}

if (have_temp_overall) {
  message("Temp overall join overlap: ",
          sum(unique(chi_sf$community) %in% unique(temp_overall$community)))
}



# ======================================================================================
# TABLE 1: Community characteristics and outcome burden
# ======================================================================================

library(dplyr)
library(gt)

# ---------------------------
# 1. Structural variables
# ---------------------------

# structural variables only
struct_vars <- names(hvi_model_matrix) %>%
  stringr::str_subset("^z_")

struct_summary <- hvi_model_matrix %>%
  select(community, all_of(struct_vars)) %>%
  distinct() %>%
  summarise(
    across(
      all_of(struct_vars),
      list(
        mean = ~mean(.x, na.rm = TRUE),
        sd   = ~sd(.x, na.rm = TRUE),
        min  = ~min(.x, na.rm = TRUE),
        max  = ~max(.x, na.rm = TRUE)
      ),
      .names = "{.col}_{.fn}"
    )
  ) %>%
  tidyr::pivot_longer(
    everything(),
    names_to = c("variable", "stat"),
    names_sep = "_(?=[^_]+$)"
  ) %>%
  tidyr::pivot_wider(names_from = stat, values_from = value) %>%
  mutate(
    variable_label = unname(var_label_dict[variable])
  ) %>%
  select(variable, variable_label, mean, sd, min, max)

# Add labels
struct_summary <- struct_summary %>%
  mutate(
    variable_label = var_label_dict[variable] %>% as.character()
  ) %>%
  select(variable, variable_label, mean, sd, min, max)

# ---------------------------
# 2. Outcome summaries
# ---------------------------

outcome_vars <- names(hvi_model_matrix) %>%
  stringr::str_subset("^(death|ed_|ems_)")

outcome_summary <- hvi_model_matrix %>%
  select(community, date, all_of(outcome_vars)) %>%
  summarise(across(
    all_of(outcome_vars),
    list(
      total = ~sum(.x, na.rm = TRUE),
      mean_daily = ~mean(.x, na.rm = TRUE)
    ),
    .names = "{.col}_{.fn}"
  )) %>%
  tidyr::pivot_longer(
    everything(),
    names_to = c("endpoint", "stat"),
    names_sep = "_(?=[^_]+$)"
  ) %>%
  tidyr::pivot_wider(names_from = stat, values_from = value)

# ---------------------------
# 3. Population / time frame
# ---------------------------

n_communities <- dplyr::n_distinct(hvi_model_matrix$community)
n_days        <- dplyr::n_distinct(hvi_model_matrix$date)
date_range    <- range(hvi_model_matrix$date, na.rm = TRUE)

# ---------------------------
# 4. Save outputs
# ---------------------------

write.csv(struct_summary, "table1_structural_characteristics.csv", row.names = FALSE)
write.csv(outcome_summary, "table1_outcomes_summary.csv", row.names = FALSE)

# ---------------------------
# 5. Publication-style table (structural)
# ---------------------------

table1_gt <- struct_summary %>%
  mutate(across(c(mean, sd, min, max), ~round(.x, 2))) %>%
  gt() %>%
  cols_label(
    variable_label = "Variable",
    mean = "Mean",
    sd = "SD",
    min = "Min",
    max = "Max"
  ) %>%
  tab_header(
    title = md("**Table 1. Community-level structural characteristics**"),
    subtitle = glue::glue("{n_communities} community areas; study period {date_range[1]} to {date_range[2]}")
  )

gtsave(table1_gt, "table1_structural.html")




# ======================================================================================
# RESULTS PARAGRAPH METRICS
# ======================================================================================

total_population <- hvi_model_matrix %>%
  group_by(community) %>%
  summarise(pop = mean(pop_offset, na.rm = TRUE)) %>%
  summarise(total = sum(pop, na.rm = TRUE)) %>%
  pull(total)

total_events <- hvi_model_matrix %>%
  summarise(across(starts_with(c("death", "ed_", "ems_")), ~sum(.x, na.rm = TRUE))) %>%
  pivot_longer(everything()) %>%
  summarise(total = sum(value)) %>%
  pull(total)

total_deaths <- hvi_model_matrix %>%
  summarise(across(starts_with("death"), ~sum(.x, na.rm = TRUE))) %>%
  rowSums()

total_ed <- hvi_model_matrix %>%
  summarise(across(starts_with("ed_"), ~sum(.x, na.rm = TRUE))) %>%
  rowSums()

total_ems <- hvi_model_matrix %>%
  summarise(across(starts_with("ems_"), ~sum(.x, na.rm = TRUE))) %>%
  rowSums()



library(dplyr)
library(tidyr)
library(stringr)
library(gt)
library(scales)

# ------------------------------------------------------------------------------
# TABLE 1: Source-level panel characteristics
# ------------------------------------------------------------------------------

# Primary all-cause source variables
source_defs <- tibble::tribble(
  ~source,      ~outcome_var,
  "ED",         "ed_visits",
  "EMS",        "ems_calls",
  "Mortality",  "deaths"
)

# Community-level characteristics to show
table1_vars <- c(
  "total_pop",
  "median_age",
  "mean_black",
  "mean_hisp",
  "median_income",
  "ac_prob",
  "ndvi",
  "pm25",
  "no2",
  "svi_rpl_themes"
)

# Match to actual column names in hvi_model_matrix
find_col <- function(x) {
  candidates <- c(x, paste0("z_", x))
  hit <- candidates[candidates %in% names(hvi_model_matrix)][1]
  if (is.na(hit)) NA_character_ else hit
}

table1_var_map <- tibble(
  display = c(
    "Median population",
    "Median age",
    "% Black",
    "% Hispanic",
    "Median income",
    "Air conditioning prevalence",
    "NDVI",
    "PM2.5",
    "NO2",
    "SVI overall"
  ),
  raw = table1_vars,
  col = vapply(table1_vars, find_col, character(1))
) %>%
  filter(!is.na(col))

# helper to back-transform z-scores if both raw and z_ do not exist
# if only z_ exists, we will still summarize it, but label clearly
pretty_value <- function(x, digits = 2) {
  if (all(is.na(x))) return(NA_character_)
  sprintf(paste0("%.", digits, "f"), median(x, na.rm = TRUE))
}

summarize_source_panel <- function(source_name, outcome_var, dat) {
  if (!outcome_var %in% names(dat)) {
    return(NULL)
  }
  
  df <- dat %>%
    filter(!is.na(.data[[outcome_var]]))
  
  if (nrow(df) == 0) return(NULL)
  
  # basic panel stats
  out <- tibble(
    characteristic = c(
      "Community areas, n",
      "Study days, n",
      "Community-days, n",
      "Total events, n",
      "Mean daily events per community-day",
      "Median daily events per community-day"
    ),
    value = c(
      n_distinct(df$community),
      n_distinct(df$date),
      nrow(df),
      sum(df[[outcome_var]], na.rm = TRUE),
      mean(df[[outcome_var]], na.rm = TRUE),
      median(df[[outcome_var]], na.rm = TRUE)
    )
  )
  
  # community characteristic rows
  char_rows <- purrr::map_dfr(seq_len(nrow(table1_var_map)), function(i) {
    col_i <- table1_var_map$col[i]
    disp_i <- table1_var_map$display[i]
    
    tibble(
      characteristic = disp_i,
      value = median(df[[col_i]], na.rm = TRUE)
    )
  })
  
  bind_rows(out, char_rows) %>%
    mutate(source = source_name)
}

table1_long <- purrr::map2_dfr(
  source_defs$source,
  source_defs$outcome_var,
  ~ summarize_source_panel(.x, .y, hvi_model_matrix)
)

table1_wide <- table1_long %>%
  select(source, characteristic, value) %>%
  pivot_wider(names_from = source, values_from = value)

# Format rows differently
count_rows <- c(
  "Community areas, n",
  "Study days, n",
  "Community-days, n",
  "Total events, n"
)

table1_display <- table1_wide %>%
  mutate(
    ED = ifelse(characteristic %in% count_rows, comma(round(ED)), sprintf("%.2f", ED)),
    EMS = ifelse(characteristic %in% count_rows, comma(round(EMS)), sprintf("%.2f", EMS)),
    Mortality = ifelse(characteristic %in% count_rows, comma(round(Mortality)), sprintf("%.2f", Mortality))
  )

write.csv(table1_display, "table1_source_panels.csv", row.names = FALSE)

table1_gt <- table1_display %>%
  gt(rowname_col = "characteristic") %>%
  tab_header(
    title = md("**Table 1. Source-specific analytic panel characteristics**"),
    subtitle = "Community-day panels for emergency department visits, EMS calls, and mortality"
  ) %>%
  cols_label(
    ED = "ED",
    EMS = "EMS",
    Mortality = "Mortality"
  )

gtsave(table1_gt, "table1_source_panels.html")



clean_label <- function(x) {
  
  x <- as.character(x)
  
  case_when(
    
    # -----------------------------
    # SVI THEMES
    # -----------------------------
    x == "svi_rpl_theme1" ~ "Socioeconomic vulnerability (SVI)",
    x == "svi_rpl_theme2" ~ "Household composition & disability (SVI)",
    x == "svi_rpl_theme3" ~ "Minority status & language (SVI)",
    x == "svi_rpl_theme4" ~ "Housing & transportation (SVI)",
    x == "svi_rpl_themes" ~ "Overall social vulnerability index",
    
    # -----------------------------
    # SVI COMPONENTS
    # -----------------------------
    x == "svi_ep_pov" ~ "Population below poverty level",
    x == "svi_ep_unemp" ~ "Unemployment rate",
    x == "svi_ep_nohsdp" ~ "No high school diploma",
    x == "svi_ep_age65" ~ "Population aged ≥65 years",
    x == "svi_ep_age17" ~ "Population aged <18 years",
    x == "svi_ep_disabl" ~ "Population with disability",
    x == "svi_ep_sngpnt" ~ "Single-parent households",
    x == "svi_ep_limeng" ~ "Limited English proficiency",
    x == "svi_ep_minrty" ~ "Minority population",
    x == "svi_ep_munit" ~ "Multi-unit housing",
    x == "svi_ep_mobile" ~ "Mobile homes",
    x == "svi_ep_crowd" ~ "Crowded housing",
    x == "svi_ep_noveh" ~ "No vehicle access",
    x == "svi_ep_groupq" ~ "Group quarters population",
    x == "svi_ep_uninsur" ~ "Uninsured population",
    x == "svi_ep_noint" ~ "No internet access",
    x == "svi_ep_hburd" ~ "Housing cost burden",
    
    # -----------------------------
    # OTHER VARIABLES
    # -----------------------------
    x == "ac_prob" ~ "Air conditioning prevalence",
    x == "ndvi" ~ "Greenness (NDVI)",
    x == "mean_ndvi" ~ "Greenness (NDVI)",
    x == "pm25" ~ "PM2.5",
    x == "no2" ~ "NO2",
    x == "pop_density_km2" ~ "Population density (per km²)",
    x == "median_age" ~ "Median age",
    x == "mean_age" ~ "Mean age",
    x == "mean_black" ~ "Black population (%)",
    x == "mean_hisp" ~ "Hispanic population (%)",
    x == "mean_white" ~ "White population (%)",
    x == "mean_asian" ~ "Asian population (%)",
    x == "median_income" ~ "Median household income",
    x == "mean_income" ~ "Mean income",
    x == "mean_unemployed" ~ "Unemployment rate",
    x == "mean_employed" ~ "Employment rate",
    x == "mean_college" ~ "College education (%)",
    x == "mean_hs" ~ "High school education (%)",
    x == "mean_male" ~ "Male population (%)",
    x == "mean_female" ~ "Female population (%)",
    
    TRUE ~ x
  )
}


# final selected variables
final_vars <- selected_variables_final$variable

# keep only structural (z_) variables
struct_vars_z <- final_vars[str_detect(final_vars, "^z_")]

# map to unscaled names
struct_vars_raw <- str_remove(struct_vars_z, "^z_")

# find actual columns available (prefer raw, fallback to z_)
get_best_var <- function(v) {
  if (v %in% names(hvi_model_matrix)) return(v)
  z_v <- paste0("z_", v)
  if (z_v %in% names(hvi_model_matrix)) return(z_v)
  return(NA_character_)
}

struct_vars_use <- unique(na.omit(vapply(struct_vars_raw, get_best_var, character(1))))

total_counts <- tibble(
  section = "Total counts",
  row = c("Community areas", "Study days", "Community-days"),
  value = c(
    n_distinct(hvi_model_matrix$community),
    n_distinct(hvi_model_matrix$date),
    nrow(hvi_model_matrix)
  )
)

source_counts <- tibble(
  section = "Total events by source",
  row = c("ED visits", "EMS calls", "Deaths"),
  value = c(
    sum(hvi_model_matrix$ed_visits, na.rm = TRUE),
    sum(hvi_model_matrix$ems_calls, na.rm = TRUE),
    sum(hvi_model_matrix$deaths, na.rm = TRUE)
  )
)


endpoint_counts <- ca_day_endpoint_risk %>%
  group_by(endpoint_key) %>%
  summarise(
    total_events = sum(observed_count, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(
    ep_meta %>% select(endpoint_key, outcome_label, source),
    by = "endpoint_key"
  ) %>%
  mutate(
    row = paste0(outcome_label, " (", source, ")"),
    section = "Events by endpoint",
    value = total_events
  ) %>%
  select(section, row, value)

struct_summary <- hvi_model_matrix %>%
  select(community, all_of(struct_vars_use)) %>%
  distinct() %>%
  summarise(
    across(
      -community,
      list(
        mean = ~mean(.x, na.rm = TRUE),
        sd   = ~sd(.x, na.rm = TRUE),
        p25  = ~quantile(.x, 0.25, na.rm = TRUE),
        median = ~median(.x, na.rm = TRUE),
        p75  = ~quantile(.x, 0.75, na.rm = TRUE)
      ),
      .names = "{.col}_{.fn}"
    )
  ) %>%
  pivot_longer(
    everything(),
    names_to = c("variable", "stat"),
    names_sep = "_(?=[^_]+$)"
  ) %>%
  pivot_wider(names_from = stat, values_from = value) %>%
  mutate(
    section = "Structural covariates",
    row = variable,
    value = sprintf(
      "%.2f (%.2f–%.2f)",
      median, p25, p75
    )
  ) %>%
  select(section, row, value)


total_counts <- total_counts %>%
  mutate(value = as.character(value))

source_counts <- source_counts %>%
  mutate(value = as.character(value))

endpoint_counts <- endpoint_counts %>%
  mutate(value = as.character(value))

struct_summary <- struct_summary %>%
  mutate(value = as.character(value))

table1_final <- bind_rows(
  total_counts,
  source_counts,
  endpoint_counts,
  struct_summary
)

table1_final <- table1_final %>%
  mutate(
    row = ifelse(row %in% names(var_label_dict),
                 var_label_dict[row],
                 row)
  )

table1_final <- table1_final %>%
  mutate(
    row = clean_label(row)
  )

write.csv(table1_final, "table1_full.csv", row.names = FALSE)

library(gt)

table1_gt <- table1_final %>%
  gt(groupname_col = "section") %>%
  tab_header(
    title = md("**Table 1. Study population, endpoint burden, and structural characteristics**")
  ) %>%
  cols_label(
    row = "Characteristic",
    value = "Value"
  )

gtsave(table1_gt, "table1_full.html")

library(dplyr)
library(lubridate)
library(tidyr)

# -------------------------------------------------------
# Helper: filter to analytic period
# -------------------------------------------------------
filter_warm_season <- function(df, date_col = "event_date") {
  
  df %>%
    mutate(date_tmp = as.Date(.data[[date_col]])) %>%
    filter(
      year(date_tmp) %in% 2019:2022,
      month(date_tmp) %in% 5:9
    ) %>%
    select(-date_tmp)
}

# Apply to DLNM input datasets
city_ed_filt     <- filter_warm_season(city_ed_full)
city_ems_filt    <- filter_warm_season(city_ems_full)
city_deaths_filt <- filter_warm_season(city_deaths_full)

# -------------------------------------------------------
# Endpoint lists (same as before)
# -------------------------------------------------------
ed_endpoints <- c(
  "ed_visits", "ed_cvd", "ed_dehydration", "ed_injury",
  "ed_renal", "ed_respiratory", "ed_syncope",
  "ed_mental", "ed_neuro", "ed_gi"
)

ems_endpoints <- c(
  "ems_calls", "ems_bleeding", "ems_cvd", "ems_gi",
  "ems_injury", "ems_mental", "ems_neuro",
  "ems_respiratory", "ems_syncope"
)

death_endpoints <- c(
  "deaths", "death_cvd", "death_injury", "death_mental",
  "death_renal", "death_respiratory",
  "death_neuro", "death_gi"
)

# -------------------------------------------------------
# Count helper
# -------------------------------------------------------
count_endpoints <- function(df, endpoints, source_name) {
  
  endpoints <- endpoints[endpoints %in% names(df)]
  
  df %>%
    summarise(across(
      all_of(endpoints),
      ~sum(.x, na.rm = TRUE)
    )) %>%
    pivot_longer(
      everything(),
      names_to = "endpoint_key",
      values_to = "total_events"
    ) %>%
    mutate(source = source_name)
}

# -------------------------------------------------------
# Compute counts (FILTERED DATA)
# -------------------------------------------------------
ed_counts    <- count_endpoints(city_ed_filt, ed_endpoints, "ED")
ems_counts   <- count_endpoints(city_ems_filt, ems_endpoints, "EMS")
death_counts <- count_endpoints(city_deaths_filt, death_endpoints, "Mortality")

endpoint_counts_raw <- bind_rows(ed_counts, ems_counts, death_counts)

endpoint_label_map <- c(
  
  # ED
  ed_visits       = "All-cause ED visits",
  ed_cvd          = "Cardiovascular ED visits",
  ed_dehydration  = "Dehydration ED visits",
  ed_injury       = "Injury-related ED visits",
  ed_renal        = "Renal ED visits",
  ed_respiratory  = "Respiratory ED visits",
  ed_syncope      = "Syncope ED visits",
  ed_mental       = "Mental health ED visits",
  ed_neuro        = "Neurologic ED visits",
  ed_gi           = "Gastrointestinal ED visits",
  
  # EMS
  ems_calls       = "All-cause EMS encounters",
  ems_bleeding    = "Bleeding EMS encounters",
  ems_cvd         = "Cardiovascular EMS encounters",
  ems_gi          = "Gastrointestinal EMS encounters",
  ems_injury      = "Injury-related EMS encounters",
  ems_mental      = "Mental health EMS encounters",
  ems_neuro       = "Neurologic EMS encounters",
  ems_respiratory = "Respiratory EMS encounters",
  ems_syncope     = "Syncope EMS encounters",
  
  # Mortality
  deaths             = "All-cause mortality",
  death_cvd          = "Cardiovascular mortality",
  death_injury       = "Injury-related mortality",
  death_mental       = "Mental health–related mortality",
  death_renal        = "Renal mortality",
  death_respiratory  = "Respiratory mortality",
  death_neuro        = "Neurologic mortality",
  death_gi           = "Gastrointestinal mortality"
)

endpoint_counts_clean <- endpoint_counts_raw %>%
  mutate(
    row = endpoint_label_map[endpoint_key],
    section = case_when(
      source == "ED" ~ "Health events by endpoint — ED",
      source == "EMS" ~ "Health events by endpoint — EMS",
      source == "Mortality" ~ "Health events by endpoint — Mortality"
    ),
    value = scales::comma(total_events)
  ) %>%
  arrange(section, desc(total_events)) %>%
  select(section, row, value)


sum(city_ed_filt$ed_neurologic)
