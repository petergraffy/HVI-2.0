# =====================================================================================
# HVI 2.0 | Heat-response analysis + prototype HVI 2.0
# Warm season only: April-October
# Shared overlap window for 3-endpoint HVI: 2019-2022
# =====================================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(lubridate)
  library(stringr)
  library(purrr)
  library(broom)
  library(ggplot2)
  library(sf)
  library(scales)
  library(splines)
  library(data.table)
  library(patchwork)
})

# -------------------------------------------------------------------------------------
# Output folders
# -------------------------------------------------------------------------------------

dir.create("results", showWarnings = FALSE, recursive = TRUE)
dir.create("figures", showWarnings = FALSE, recursive = TRUE)

# -------------------------------------------------------------------------------------
# Helpers
# -------------------------------------------------------------------------------------

# -----------------------------
# SHARED STUDY WINDOW CONFIG
# -----------------------------
analysis_year_start <- 2019L
analysis_year_end   <- 2022L
analysis_month_start <- 5L
analysis_month_end   <- 9L

in_analysis_window <- function(date) {
  yr <- lubridate::year(date)
  mo <- lubridate::month(date)
  yr >= analysis_year_start &
    yr <= analysis_year_end &
    mo >= analysis_month_start &
    mo <= analysis_month_end
}

standardize_community <- function(x) {
  x <- as.character(x)
  x <- str_squish(x)
  x <- str_to_title(x)
  
  x <- case_when(
    x %in% c("Lakeview") ~ "Lake View",
    x %in% c("Ohare", "O'hare", "O Hare") ~ "O'Hare",
    TRUE ~ x
  )
  
  x
}

safe_scale <- function(x) {
  if (all(is.na(x)) || sd(x, na.rm = TRUE) == 0) return(rep(NA_real_, length(x)))
  as.numeric(scale(x))
}

save_plot <- function(plot_obj, filename, width = 9, height = 6, dpi = 300) {
  ggsave(
    filename = file.path("figures", filename),
    plot = plot_obj,
    width = width,
    height = height,
    dpi = dpi
  )
}

# -------------------------------------------------------------------------------------
# 0) Start from panel_causes
# Assumes panel_causes and comm_areas already exist in memory
# -------------------------------------------------------------------------------------

panel_heat <- panel_causes %>%
  mutate(
    community = standardize_community(community),
    event_date = as.Date(event_date),
    year = year(event_date),
    month = month(event_date),
    dow = wday(event_date),
    warm_season = month %in% 5:9,
    overlap_3src = year >= 2019 & year <= 2022
  ) %>%
  filter(warm_season)

# Community areas for mapping
ca <- comm_areas %>%
  janitor::clean_names()

if (!"community" %in% names(ca)) {
  nm <- intersect(c("community", "ca_name", "commarea", "community_area"), names(ca))[1]
  if (is.na(nm)) stop("Could not find community field in comm_areas.")
  ca <- ca %>% mutate(community = .data[[nm]])
}

ca <- ca %>%
  mutate(community = standardize_community(community))

# -------------------------------------------------------------------------------------
# 1) Heat metrics
# Build these on the shared overlap window for 3-endpoint comparisons
# -------------------------------------------------------------------------------------

heat_ref <- panel_heat %>%
  filter(overlap_3src)

# citywide thresholds
city_tmax_p90 <- quantile(heat_ref$tmax, 0.90, na.rm = TRUE)
city_tmax_p95 <- quantile(heat_ref$tmax, 0.95, na.rm = TRUE)
city_tmean_p90 <- quantile(heat_ref$tmean, 0.90, na.rm = TRUE)
city_humidity_p90 <- quantile(heat_ref$humidity, 0.90, na.rm = TRUE)

# community-month climatology for anomaly-based metrics
community_month_normals <- panel_heat %>%
  group_by(community, month) %>%
  summarise(
    tmax_month_norm = mean(tmax, na.rm = TRUE),
    tmean_month_norm = mean(tmean, na.rm = TRUE),
    humidity_month_norm = mean(humidity, na.rm = TRUE),
    .groups = "drop"
  )

heat_ref <- heat_ref %>%
  left_join(community_month_normals, by = c("community", "month")) %>%
  arrange(community, event_date) %>%
  group_by(community) %>%
  mutate(
    tmax_anom = tmax - tmax_month_norm,
    tmean_anom = tmean - tmean_month_norm,
    humidity_anom = humidity - humidity_month_norm,
    
    hot_day90 = as.integer(tmax >= city_tmax_p90),
    extreme_heat95 = as.integer(tmax >= city_tmax_p95),
    high_tmean90 = as.integer(tmean >= city_tmean_p90),
    humid_day90 = as.integer(humidity >= city_humidity_p90),
    compound_heat = as.integer(tmax >= city_tmax_p90 & humidity >= city_humidity_p90),
    anomaly_hot = as.integer(tmax_anom >= quantile(tmax_anom, 0.90, na.rm = TRUE)),
    
    early_season = month %in% 4:6,
    early_heat = as.integer(early_season & extreme_heat95),
    
    lag_hot1 = lag(hot_day90, 1, default = 0),
    lag_hot2 = lag(hot_day90, 2, default = 0),
    heatwave_2d = as.integer(hot_day90 == 1 & lag_hot1 == 1),
    heatwave_3d = as.integer(hot_day90 == 1 & lag_hot1 == 1 & lag_hot2 == 1)
  ) %>%
  ungroup()

# -------------------------------------------------------------------------------------
# 2) Quick descriptive "what happens on hot days?"
# Use shared overlap window 2019-2022 for all three endpoints
# -------------------------------------------------------------------------------------

heat_compare <- heat_ref %>%
  mutate(
    heat_group = case_when(
      extreme_heat95 == 1 ~ "Extreme heat (>=95th pct)",
      TRUE ~ "Non-extreme"
    )
  ) %>%
  group_by(heat_group) %>%
  summarise(
    days = n(),
    mean_deaths = mean(deaths, na.rm = TRUE),
    mean_ed = mean(ed_visits, na.rm = TRUE),
    mean_ems = mean(ems_calls, na.rm = TRUE),
    
    mean_death_heat = mean(death_heat, na.rm = TRUE),
    mean_ed_heat = mean(ed_heat, na.rm = TRUE),
    mean_ems_heat = mean(ems_heat, na.rm = TRUE),
    
    mean_death_cvd = mean(death_cvd, na.rm = TRUE),
    mean_ed_cvd = mean(ed_cvd, na.rm = TRUE),
    mean_ems_cvd = mean(ems_cvd, na.rm = TRUE),
    
    mean_death_respiratory = mean(death_respiratory, na.rm = TRUE),
    mean_ed_respiratory = mean(ed_respiratory, na.rm = TRUE),
    mean_ems_respiratory = mean(ems_respiratory, na.rm = TRUE),
    .groups = "drop"
  )

fwrite(heat_compare, "results/heat_day_vs_nonheat_summary_2019_2022.csv")

# Plot: hot vs non-hot mean all-cause burden
plot_heat_contrast_dat <- heat_ref %>%
  mutate(heat_group = case_when(
      extreme_heat95 == 1 ~ "Extreme heat (>=95th pct)",
      TRUE ~ "Non-extreme"
    )) %>% 
  filter(overlap_3src) %>%
  mutate(
    heat_group = ifelse(extreme_heat95 == 1, "Extreme heat", "Non-extreme")
  ) %>%
  group_by(heat_group) %>%
  summarise(
    Mortality = mean(deaths, na.rm = TRUE),
    `ED Visits` = mean(ed_visits, na.rm = TRUE),
    `EMS Calls` = mean(ems_calls, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_longer(-heat_group, names_to = "outcome", values_to = "mean_count")

plot_heat_contrast <- ggplot(plot_heat_contrast_dat, aes(x = heat_group, y = mean_count, fill = outcome)) +
  geom_col(position = "dodge") +
  labs(
    title = "Average Warm-Season Daily Counts on Extreme vs Non-Extreme Heat Days",
    subtitle = "May-September, 2019-2022",
    x = NULL,
    y = "Mean daily count",
    fill = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"))

save_plot(plot_heat_contrast, "heat_vs_nonheat_allcause_2019_2022.png", width = 10, height = 6)
save_plot(plot_heat_contrast, "heat_vs_nonheat_allcause_2019_2022.pdf", width = 10, height = 6)

# -------------------------------------------------------------------------------------
# 3) Community-specific heat effect models
# Prototype approach:
# community-stratified quasi-Poisson models with:
# outcome ~ extreme_heat95 + humidity + factor(year) + factor(month) + factor(dow)
# -------------------------------------------------------------------------------------

fit_heat_model_by_community <- function(df, outcome_var, heat_var = "extreme_heat95") {
  
  dat <- df %>%
    filter(!is.na(.data[[outcome_var]]),
           !is.na(.data[[heat_var]]),
           !is.na(humidity)) %>%
    arrange(event_date)
  
  # need both exposed and unexposed days
  if (nrow(dat) < 50) {
    return(tibble(
      rr = NA_real_,
      rr_low = NA_real_,
      rr_high = NA_real_,
      beta = NA_real_,
      se = NA_real_,
      excess_total = NA_real_,
      excess_per_hot_day = NA_real_,
      n_days = nrow(dat),
      n_hot_days = sum(dat[[heat_var]] == 1, na.rm = TRUE)
    ))
  }
  
  if (length(unique(dat[[heat_var]])) < 2) {
    return(tibble(
      rr = NA_real_,
      rr_low = NA_real_,
      rr_high = NA_real_,
      beta = NA_real_,
      se = NA_real_,
      excess_total = NA_real_,
      excess_per_hot_day = NA_real_,
      n_days = nrow(dat),
      n_hot_days = sum(dat[[heat_var]] == 1, na.rm = TRUE)
    ))
  }
  
  form <- as.formula(
    paste0(
      outcome_var,
      " ~ ", heat_var,
      " + scale(humidity) + factor(year) + factor(month) + factor(dow)"
    )
  )
  
  mod <- tryCatch(
    glm(form, family = quasipoisson(link = "log"), data = dat),
    error = function(e) NULL
  )
  
  if (is.null(mod)) {
    return(tibble(
      rr = NA_real_,
      rr_low = NA_real_,
      rr_high = NA_real_,
      beta = NA_real_,
      se = NA_real_,
      excess_total = NA_real_,
      excess_per_hot_day = NA_real_,
      n_days = nrow(dat),
      n_hot_days = sum(dat[[heat_var]] == 1, na.rm = TRUE)
    ))
  }
  
  coef_tab <- broom::tidy(mod)
  hit <- coef_tab %>% filter(term == heat_var)
  
  if (nrow(hit) == 0) {
    return(tibble(
      rr = NA_real_,
      rr_low = NA_real_,
      rr_high = NA_real_,
      beta = NA_real_,
      se = NA_real_,
      excess_total = NA_real_,
      excess_per_hot_day = NA_real_,
      n_days = nrow(dat),
      n_hot_days = sum(dat[[heat_var]] == 1, na.rm = TRUE)
    ))
  }
  
  beta <- hit$estimate[1]
  se <- hit$std.error[1]
  
  rr <- exp(beta)
  rr_low <- exp(beta - 1.96 * se)
  rr_high <- exp(beta + 1.96 * se)
  
  hot_dat <- dat %>% filter(.data[[heat_var]] == 1)
  
  if (nrow(hot_dat) == 0) {
    excess_total <- NA_real_
    excess_per_hot_day <- NA_real_
  } else {
    pred_obs <- predict(mod, newdata = hot_dat, type = "response")
    
    cf_dat <- hot_dat
    cf_dat[[heat_var]] <- 0
    pred_cf <- predict(mod, newdata = cf_dat, type = "response")
    
    excess_total <- sum(pred_obs - pred_cf, na.rm = TRUE)
    excess_per_hot_day <- excess_total / nrow(hot_dat)
  }
  
  tibble(
    rr = rr,
    rr_low = rr_low,
    rr_high = rr_high,
    beta = beta,
    se = se,
    excess_total = excess_total,
    excess_per_hot_day = excess_per_hot_day,
    n_days = nrow(dat),
    n_hot_days = sum(dat[[heat_var]] == 1, na.rm = TRUE)
  )
}

# Shared overlap window for 3-endpoint HVI
panel_overlap <- heat_ref %>%
  filter(overlap_3src)

# All-cause endpoints
hvi_mortality <- panel_overlap %>%
  group_by(community) %>%
  group_modify(~ fit_heat_model_by_community(.x, outcome_var = "deaths")) %>%
  ungroup() %>%
  rename_with(~ paste0("mort_", .x), -community)

hvi_ed <- panel_overlap %>%
  group_by(community) %>%
  group_modify(~ fit_heat_model_by_community(.x, outcome_var = "ed_visits")) %>%
  ungroup() %>%
  rename_with(~ paste0("ed_", .x), -community)

hvi_ems <- panel_overlap %>%
  group_by(community) %>%
  group_modify(~ fit_heat_model_by_community(.x, outcome_var = "ems_calls")) %>%
  ungroup() %>%
  rename_with(~ paste0("ems_", .x), -community)

# Optional: heat-specific endpoint models
heat_specific_mort <- panel_overlap %>%
  group_by(community) %>%
  group_modify(~ fit_heat_model_by_community(.x, outcome_var = "death_cvd")) %>%
  ungroup() %>%
  rename_with(~ paste0("mort_heat_", .x), -community)

heat_specific_ed <- panel_overlap %>%
  group_by(community) %>%
  group_modify(~ fit_heat_model_by_community(.x, outcome_var = "ed_cvd")) %>%
  ungroup() %>%
  rename_with(~ paste0("ed_heat_", .x), -community)

heat_specific_ems <- panel_overlap %>%
  group_by(community) %>%
  group_modify(~ fit_heat_model_by_community(.x, outcome_var = "ems_cvd")) %>%
  ungroup() %>%
  rename_with(~ paste0("ems_heat_", .x), -community)

# -------------------------------------------------------------------------------------
# 4) Build prototype HVI 2.0
# Main version: all-cause excess burden per hot day, 2019-2022 warm season
# -------------------------------------------------------------------------------------

hvi_proto <- hvi_mortality %>%
  left_join(hvi_ed, by = "community") %>%
  left_join(hvi_ems, by = "community") %>%
  left_join(heat_specific_mort, by = "community") %>%
  left_join(heat_specific_ed, by = "community") %>%
  left_join(heat_specific_ems, by = "community") %>%
  mutate(
    # Primary HVI components: all-cause excess burden per hot day
    z_mort_excess = safe_scale(mort_excess_per_hot_day),
    z_ed_excess   = safe_scale(ed_excess_per_hot_day),
    z_ems_excess  = safe_scale(ems_excess_per_hot_day),
    
    # Optional supportive metrics: heat-specific observed burden
    z_mort_heat_excess = safe_scale(mort_heat_excess_per_hot_day),
    z_ed_heat_excess   = safe_scale(ed_heat_excess_per_hot_day),
    z_ems_heat_excess  = safe_scale(ems_heat_excess_per_hot_day),
    
    # Transparent primary composite: equal weights
    hvi2_allcause = rowMeans(
      cbind(z_mort_excess, z_ed_excess, z_ems_excess),
      na.rm = TRUE
    ),
    
    # Broader heat-sensitive composite if desired
    hvi2_heat_sensitive = rowMeans(
      cbind(z_mort_excess, z_ed_excess, z_ems_excess,
            z_mort_heat_excess, z_ed_heat_excess, z_ems_heat_excess),
      na.rm = TRUE
    ),
    
    # Percentile ranks
    hvi2_allcause_pct = percent_rank(hvi2_allcause),
    hvi2_heat_sensitive_pct = percent_rank(hvi2_heat_sensitive)
  )

fwrite(hvi_proto, "results/hvi2_prototype_community_metrics_2019_2022.csv")

# -------------------------------------------------------------------------------------
# 5) Maps
# -------------------------------------------------------------------------------------

map_hvi <- ca %>%
  left_join(hvi_proto, by = "community")

map_mort_excess <- ggplot(map_hvi) +
  geom_sf(aes(fill = mort_excess_per_hot_day), color = NA) +
  scale_fill_viridis_c(option = "plasma", labels = label_number(accuracy = 0.01)) +
  theme_void() +
  labs(
    title = "Excess Mortality per Extreme Heat Day",
    subtitle = "Warm season, 2019-2022",
    fill = "Excess/day"
  )

map_ed_excess <- ggplot(map_hvi) +
  geom_sf(aes(fill = ed_excess_per_hot_day), color = NA) +
  scale_fill_viridis_c(option = "plasma", labels = label_number(accuracy = 0.01)) +
  theme_void() +
  labs(
    title = "Excess ED Visits per Extreme Heat Day",
    subtitle = "Warm season, 2019-2022",
    fill = "Excess/day"
  )

map_ems_excess <- ggplot(map_hvi) +
  geom_sf(aes(fill = ems_excess_per_hot_day), color = NA) +
  scale_fill_viridis_c(option = "plasma", labels = label_number(accuracy = 0.01)) +
  theme_void() +
  labs(
    title = "Excess EMS Calls per Extreme Heat Day",
    subtitle = "Warm season, 2019-2022",
    fill = "Excess/day"
  )

map_hvi_allcause <- ggplot(map_hvi) +
  geom_sf(aes(fill = hvi2_allcause), color = NA) +
  scale_fill_viridis_c(option = "magma", labels = label_number(accuracy = 0.1)) +
  theme_void() +
  labs(
    title = "HVI 2.0 Prototype",
    subtitle = "Composite all-cause excess burden on extreme heat days, 2019-2022",
    fill = "HVI 2.0"
  )

map_hvi_heat_sensitive <- ggplot(map_hvi) +
  geom_sf(aes(fill = hvi2_heat_sensitive), color = NA) +
  scale_fill_viridis_c(option = "magma", labels = label_number(accuracy = 0.1)) +
  theme_void() +
  labs(
    title = "HVI 2.0 Heat-Sensitive Prototype",
    subtitle = "Composite all-cause + heat-specific excess burden, 2019-2022",
    fill = "HVI 2.0"
  )

map_excess_combo <- map_mort_excess | map_ed_excess | map_ems_excess

save_plot(map_mort_excess, "map_excess_mortality_per_hot_day_2019_2022.png", width = 7, height = 7)
save_plot(map_ed_excess, "map_excess_ed_per_hot_day_2019_2022.png", width = 7, height = 7)
save_plot(map_ems_excess, "map_excess_ems_per_hot_day_2019_2022.png", width = 7, height = 7)
save_plot(map_excess_combo, "map_excess_three_panel_2019_2022.png", width = 16, height = 6)
save_plot(map_excess_combo, "map_excess_three_panel_2019_2022.pdf", width = 16, height = 6)

save_plot(map_hvi_allcause, "map_hvi2_allcause_prototype_2019_2022.png", width = 8, height = 7)
save_plot(map_hvi_allcause, "map_hvi2_allcause_prototype_2019_2022.pdf", width = 8, height = 7)
save_plot(map_hvi_heat_sensitive, "map_hvi2_heat_sensitive_prototype_2019_2022.png", width = 8, height = 7)
save_plot(map_hvi_heat_sensitive, "map_hvi2_heat_sensitive_prototype_2019_2022.pdf", width = 8, height = 7)

# -------------------------------------------------------------------------------------
# 6) Community ranking table
# -------------------------------------------------------------------------------------

hvi_ranked <- hvi_proto %>%
  arrange(desc(hvi2_allcause)) %>%
  select(
    community,
    mort_excess_per_hot_day,
    ed_excess_per_hot_day,
    ems_excess_per_hot_day,
    hvi2_allcause,
    hvi2_allcause_pct,
    hvi2_heat_sensitive,
    hvi2_heat_sensitive_pct
  )

fwrite(hvi_ranked, "results/hvi2_prototype_ranked_communities_2019_2022.csv")

# -------------------------------------------------------------------------------------
# 7) Optional: citywide model for quick overall signal
# -------------------------------------------------------------------------------------

city_mod_deaths <- glm(
  deaths ~ extreme_heat95 + scale(humidity) + factor(community) + factor(year) + factor(month) + factor(dow),
  family = quasipoisson(link = "log"),
  data = panel_overlap
)

city_mod_ed <- glm(
  ed_visits ~ extreme_heat95 + scale(humidity) + factor(community) + factor(year) + factor(month) + factor(dow),
  family = quasipoisson(link = "log"),
  data = panel_overlap
)

city_mod_ems <- glm(
  ems_calls ~ extreme_heat95 + scale(humidity) + factor(community) + factor(year) + factor(month) + factor(dow),
  family = quasipoisson(link = "log"),
  data = panel_overlap
)

citywide_heat_effects <- bind_rows(
  tidy(city_mod_deaths) %>% filter(term == "extreme_heat95") %>% mutate(outcome = "Mortality"),
  tidy(city_mod_ed)     %>% filter(term == "extreme_heat95") %>% mutate(outcome = "ED Visits"),
  tidy(city_mod_ems)    %>% filter(term == "extreme_heat95") %>% mutate(outcome = "EMS Calls")
) %>%
  mutate(
    rr = exp(estimate),
    rr_low = exp(estimate - 1.96 * std.error),
    rr_high = exp(estimate + 1.96 * std.error)
  )

fwrite(citywide_heat_effects, "results/citywide_extreme_heat_effects_2019_2022.csv")

# -------------------------------------------------------------------------------------
# 8) Sanity checks
# -------------------------------------------------------------------------------------

cat("\nCitywide heat thresholds:\n")
print(list(
  city_tmax_p90 = city_tmax_p90,
  city_tmax_p95 = city_tmax_p95,
  city_tmean_p90 = city_tmean_p90,
  city_humidity_p90 = city_humidity_p90
))

cat("\nPrototype HVI summary:\n")
print(summary(hvi_proto))

cat("\nCitywide extreme heat effects:\n")
print(citywide_heat_effects)


# -------------------------------------------------------------------------------------
# 7) Optional: citywide model for quick overall signal
# -------------------------------------------------------------------------------------

city_mod_deaths90 <- glm(
  deaths ~ hot_day90 + scale(humidity) + factor(community) + factor(year) + factor(month) + factor(dow),
  family = quasipoisson(link = "log"),
  data = panel_overlap
)

city_mod_ed90 <- glm(
  ed_visits ~ hot_day90 + scale(humidity) + factor(community) + factor(year) + factor(month) + factor(dow),
  family = quasipoisson(link = "log"),
  data = panel_overlap
)

city_mod_ems90 <- glm(
  ems_calls ~ hot_day90 + scale(humidity) + factor(community) + factor(year) + factor(month) + factor(dow),
  family = quasipoisson(link = "log"),
  data = panel_overlap
)

citywide_heat_effects <- bind_rows(
  tidy(city_mod_deaths90) %>% filter(term == "hot_day90") %>% mutate(outcome = "Mortality"),
  tidy(city_mod_ed90)     %>% filter(term == "hot_day90") %>% mutate(outcome = "ED Visits"),
  tidy(city_mod_ems90)    %>% filter(term == "hot_day90") %>% mutate(outcome = "EMS Calls")
) %>%
  mutate(
    rr = exp(estimate),
    rr_low = exp(estimate - 1.96 * std.error),
    rr_high = exp(estimate + 1.96 * std.error)
  )

fwrite(citywide_heat_effects, "results/citywide_extreme_heat_effects90_2019_2022.csv")

cat("\nCitywide extreme heat effects:\n")
print(citywide_heat_effects)


