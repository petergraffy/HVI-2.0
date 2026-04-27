# =====================================================================================
# HVI 2.0 | Cause-specific aggregation for Mortality, ED, and EMS
# Uses pre-filtered analytic inputs: deaths_use, ed_use, ems_use
# =====================================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(stringr)
  library(lubridate)
  library(janitor)
  library(tidyr)
})

source(file.path(if (dir.exists("code")) "code" else ".", "00_config.R"))

# -------------------------------------------------------------------------------------
# File paths
# -------------------------------------------------------------------------------------

path_panel  <- file.path(HVI_PATHS$private_outputs$derived, "community_day_panel_with_climate.csv")
out_panel   <- file.path(HVI_PATHS$private_outputs$derived, "community_day_panel_with_climate_causes.csv")

# -------------------------------------------------------------------------------------
# Helpers
# -------------------------------------------------------------------------------------

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

clean_code <- function(x) {
  x <- as.character(x)
  x <- str_squish(x)
  x[x %in% c("", "NA", "N/A", "NULL", "null")] <- NA_character_
  x <- str_to_upper(x)
  x
}

normalize_icd <- function(x) {
  x <- clean_code(x)
  x <- str_replace_all(x, "[^A-Z0-9]", "")
  x
}

icd3 <- function(x) {
  x <- normalize_icd(x)
  ifelse(is.na(x), NA_character_, str_sub(x, 1, 3))
}

is_icd10 <- function(x) {
  x <- normalize_icd(x)
  !is.na(x) & str_detect(x, "^[A-TV-Z]")
}

is_icd9_numeric <- function(x) {
  x <- normalize_icd(x)
  !is.na(x) & str_detect(x, "^[0-9]")
}

icd9_num <- function(x) {
  x <- normalize_icd(x)
  suppressWarnings(as.numeric(str_extract(x, "^[0-9]{1,3}")))
}

classify_icd_broad <- function(code) {
  code_norm <- normalize_icd(code)
  n9 <- icd9_num(code_norm)
  
  tibble(
    cvd = case_when(
      is_icd10(code_norm) & str_detect(code_norm, "^I") ~ 1L,
      is_icd9_numeric(code_norm) & !is.na(n9) & n9 >= 390 & n9 <= 459 ~ 1L,
      TRUE ~ 0L
    ),
    
    respiratory = case_when(
      is_icd10(code_norm) & str_detect(code_norm, "^J") ~ 1L,
      is_icd9_numeric(code_norm) & !is.na(n9) & n9 >= 460 & n9 <= 519 ~ 1L,
      TRUE ~ 0L
    ),
    
    renal = case_when(
      is_icd10(code_norm) & str_detect(code_norm, "^N1[7-9]") ~ 1L,
      is_icd10(code_norm) & str_detect(code_norm, "^N2[0-9]") ~ 1L,
      is_icd9_numeric(code_norm) & !is.na(n9) & n9 >= 580 & n9 <= 589 ~ 1L,
      TRUE ~ 0L
    ),
    
    heat = case_when(
      is_icd10(code_norm) & str_detect(code_norm, "^T67") ~ 1L,
      is_icd10(code_norm) & str_detect(code_norm, "^X30") ~ 1L,
      is_icd9_numeric(code_norm) & !is.na(n9) & n9 == 992 ~ 1L,
      code_norm %in% c("E9000") ~ 1L,
      TRUE ~ 0L
    ),
    
    dehydration = case_when(
      is_icd10(code_norm) & str_detect(code_norm, "^E86") ~ 1L,
      is_icd9_numeric(code_norm) & !is.na(n9) & n9 == 276 ~ 1L,
      TRUE ~ 0L
    ),
    
    syncope = case_when(
      is_icd10(code_norm) & str_detect(code_norm, "^R55") ~ 1L,
      is_icd9_numeric(code_norm) & !is.na(n9) & n9 == 780 ~ 1L,
      TRUE ~ 0L
    ),
    
    mental = case_when(
      is_icd10(code_norm) & str_detect(code_norm, "^F") ~ 1L,
      is_icd9_numeric(code_norm) & !is.na(n9) & n9 >= 290 & n9 <= 319 ~ 1L,
      TRUE ~ 0L
    ),
    
    neurologic = case_when(
      is_icd10(code_norm) & str_detect(code_norm, "^G") ~ 1L,
      is_icd10(code_norm) & str_detect(code_norm, "^R41") ~ 1L,
      is_icd9_numeric(code_norm) & !is.na(n9) & n9 >= 320 & n9 <= 389 ~ 1L,
      TRUE ~ 0L
    ),
    
    injury = case_when(
      is_icd10(code_norm) & str_detect(code_norm, "^(S|T|V|W|X|Y)") ~ 1L,
      is_icd9_numeric(code_norm) & !is.na(n9) & n9 >= 800 & n9 <= 999 ~ 1L,
      TRUE ~ 0L
    ),
    
    gi = case_when(
      is_icd10(code_norm) & str_detect(code_norm, "^R11") ~ 1L,
      is_icd10(code_norm) & str_detect(code_norm, "^K") ~ 1L,
      is_icd9_numeric(code_norm) & !is.na(n9) & n9 >= 520 & n9 <= 579 ~ 1L,
      TRUE ~ 0L
    )
  )
}

classify_ems_symptom <- function(x) {
  y <- clean_code(x)
  
  tibble(
    ems_respiratory = ifelse(!is.na(y) & str_detect(y, "DYSPNEA|SHORTNESS OF BREATH|RESPIR"), 1L, 0L),
    ems_syncope     = ifelse(!is.na(y) & str_detect(y, "SYNCOPE|NEAR SYNCOPE|FAINT"), 1L, 0L),
    ems_neuro       = ifelse(!is.na(y) & str_detect(y, "ALTERED MENTAL STATUS|DIZZINESS|VERTIGO|DYSPHASIA|SEIZURE|HEADACHE"), 1L, 0L),
    ems_mental      = ifelse(!is.na(y) & str_detect(y, "BEHAVIORAL|PSYCH|ANXIETY|WORRIES|CRISIS"), 1L, 0L),
    ems_injury      = ifelse(!is.na(y) & str_detect(y, "INJURY|TRAUMA|HEAD INJURY|LEG INJURY|PENETRATING"), 1L, 0L),
    ems_gi          = ifelse(!is.na(y) & str_detect(y, "NAUSEA|VOMITING|ABDOM"), 1L, 0L),
    ems_bleeding    = ifelse(!is.na(y) & str_detect(y, "BLEEDING|EPISTAXIS|NOSEBLEED"), 1L, 0L),
    ems_cvd         = ifelse(!is.na(y) & str_detect(y, "HYPERTENSION|CHEST PAIN|CARDIAC|PALPIT"), 1L, 0L),
    ems_heat        = ifelse(!is.na(y) & str_detect(y, "HEAT|HYPERTHERMIA|DEHYDRAT"), 1L, 0L)
  )
}

aggregate_daily_counts <- function(df, prefix) {
  indicator_cols <- names(df)[str_detect(names(df), paste0("^", prefix))]
  df %>%
    group_by(community, event_date) %>%
    summarise(across(all_of(indicator_cols), ~ sum(.x, na.rm = TRUE)), .groups = "drop")
}

# -------------------------------------------------------------------------------------
# Inputs: use already-filtered analytic datasets
# Assumes deaths_use, ed_use, ems_use already exist in memory
# -------------------------------------------------------------------------------------

if (!exists("deaths_use")) stop("deaths_use not found in environment.")
if (!exists("ed_use"))     stop("ed_use not found in environment.")
if (!exists("ems_use"))    stop("ems_use not found in environment.")

panel <- fread(path_panel) %>% clean_names()

deaths_use <- deaths_use %>%
  mutate(
    community = standardize_community(community),
    event_date = as.Date(event_date)
  )

ed_use <- ed_use %>%
  mutate(
    community = standardize_community(community),
    event_date = as.Date(event_date)
  )

ems_use <- ems_use %>%
  mutate(
    community = standardize_community(community),
    event_date = as.Date(event_date)
  )

panel <- panel %>%
  mutate(
    community = standardize_community(community),
    event_date = as.Date(event_date)
  )

# -------------------------------------------------------------------------------------
# Deaths
# -------------------------------------------------------------------------------------

if (!"undcause" %in% names(deaths_use)) {
  stop("deaths_use does not contain 'undcause'.")
}

deaths_class <- bind_cols(
  deaths_use,
  classify_icd_broad(deaths_use$undcause) %>%
    transmute(
      death_cvd         = cvd,
      death_respiratory = respiratory,
      death_renal       = renal,
      death_heat        = heat,
      death_dehydration = dehydration,
      death_syncope     = syncope,
      death_mental      = mental,
      death_neurologic  = neurologic,
      death_injury      = injury,
      death_gi          = gi
    )
)

deaths_daily_cause <- aggregate_daily_counts(deaths_class, "death_")

# -------------------------------------------------------------------------------------
# ED
# -------------------------------------------------------------------------------------

if (!"diagnosis_diagnosis_code_1" %in% names(ed_use)) {
  stop("ed_use does not contain 'diagnosis_diagnosis_code_1'.")
}

ed_class <- bind_cols(
  ed_use,
  classify_icd_broad(ed_use$diagnosis_diagnosis_code_1) %>%
    transmute(
      ed_cvd         = cvd,
      ed_respiratory = respiratory,
      ed_renal       = renal,
      ed_heat        = heat,
      ed_dehydration = dehydration,
      ed_syncope     = syncope,
      ed_mental      = mental,
      ed_neurologic  = neurologic,
      ed_injury      = injury,
      ed_gi          = gi
    )
)

ed_daily_cause <- aggregate_daily_counts(ed_class, "ed_")

# -------------------------------------------------------------------------------------
# EMS
# -------------------------------------------------------------------------------------

aggregate_daily_counts <- function(df, prefix) {
  indicator_cols <- names(df)[str_detect(names(df), paste0("^", prefix))]
  indicator_cols <- indicator_cols[sapply(df[indicator_cols], is.numeric)]
  
  df %>%
    group_by(community, event_date) %>%
    summarise(across(all_of(indicator_cols), ~ sum(.x, na.rm = TRUE)), .groups = "drop")
}

if (!"situation_primary_symptom_esituation_09" %in% names(ems_use)) {
  stop("ems_use does not contain 'situation_primary_symptom_esituation_09'.")
}

ems_class <- bind_cols(
  ems_use,
  classify_ems_symptom(ems_use$situation_primary_symptom_esituation_09)
)

ems_daily_cause <- aggregate_daily_counts(ems_class, "ems_")

# -------------------------------------------------------------------------------------
# Join to panel
# -------------------------------------------------------------------------------------

panel_causes <- panel %>%
  left_join(deaths_daily_cause, by = c("community", "event_date")) %>%
  left_join(ed_daily_cause,     by = c("community", "event_date")) %>%
  left_join(ems_daily_cause,    by = c("community", "event_date")) %>%
  mutate(
    across(
      c(starts_with("death_"), starts_with("ed_"), starts_with("ems_")),
      ~ coalesce(.x, 0L)
    )
  )

# -------------------------------------------------------------------------------------
# QA
# -------------------------------------------------------------------------------------

cat("\nDeaths cause totals:\n")
print(colSums(select(deaths_daily_cause, starts_with("death_")), na.rm = TRUE))

cat("\nED cause totals:\n")
print(colSums(select(ed_daily_cause, starts_with("ed_")), na.rm = TRUE))

cat("\nEMS cause totals:\n")
print(colSums(select(ems_daily_cause, starts_with("ems_")), na.rm = TRUE))

cat("\nPanel with causes dimensions:\n")
print(dim(panel_causes))

# -------------------------------------------------------------------------------------
# Save
# -------------------------------------------------------------------------------------

fwrite(panel_causes, out_panel)
cat("\nSaved to:\n", out_panel, "\n")

# =====================================================================================
# HVI 2.0 | Exploratory plots and maps for all-cause + cause-specific outcomes
# =====================================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(lubridate)
  library(stringr)
  library(sf)
  library(janitor)
  library(patchwork)
  library(scales)
})

# -------------------------------------------------------------------------------------
# Output directory
# -------------------------------------------------------------------------------------

fig_dir <- file.path(HVI_PATHS$private_outputs$publication_outputs, "exploratory_figures")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# -------------------------------------------------------------------------------------
# Helpers
# -------------------------------------------------------------------------------------

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

save_plot <- function(plot_obj, filename, width = 9, height = 6, dpi = 300) {
  ggsave(
    filename = file.path(fig_dir, filename),
    plot = plot_obj,
    width = width,
    height = height,
    dpi = dpi
  )
}

base_theme <- theme_minimal(base_size = 12) +
  theme(
    legend.position = "right",
    plot.title = element_text(face = "bold"),
    strip.text = element_text(face = "bold")
  )

theme_set(base_theme)

# -------------------------------------------------------------------------------------
# Standardize panel and define endpoint availability windows
# -------------------------------------------------------------------------------------

panel_causes <- panel_causes %>%
  mutate(
    community = standardize_community(community),
    event_date = as.Date(event_date),
    year = year(event_date),
    has_deaths = year >= 1993 & year < 2023 & month %in% 5:9,
    has_ed     = year >= 2011 & year < 2024 & month %in% 5:9,
    has_ems    = year >= 2019 & month %in% 5:9
  )

# Standardize community areas object
ca <- comm_areas %>%
  clean_names()

if (!"community" %in% names(ca)) {
  nm <- intersect(c("community", "ca_name", "commarea", "community_area"), names(ca))[1]
  if (is.na(nm)) stop("Could not find community field in comm_areas.")
  ca <- ca %>% mutate(community = .data[[nm]])
}

ca <- ca %>%
  mutate(community = standardize_community(community))

# -------------------------------------------------------------------------------------
# 1) Distribution of all-cause daily counts
# Use endpoint-specific observation windows to avoid structural zeros
# -------------------------------------------------------------------------------------

plot_dist_all_dat <- bind_rows(
  panel_causes %>%
    filter(has_deaths) %>%
    transmute(outcome = "Mortality", count = deaths),
  
  panel_causes %>%
    filter(has_ed) %>%
    transmute(outcome = "ED Visits", count = ed_visits),
  
  panel_causes %>%
    filter(has_ems) %>%
    transmute(outcome = "EMS Calls", count = ems_calls)
)

plot_dist_all <- ggplot(plot_dist_all_dat, aes(x = count)) +
  geom_histogram(bins = 50, fill = "steelblue", alpha = 0.8) +
  facet_wrap(~ outcome, scales = "free") +
  labs(
    title = "Distribution of Daily Counts by Outcome",
    x = "Daily Count",
    y = "Frequency"
  )

save_plot(plot_dist_all, "distribution_all_cause_daily_counts.png", width = 10, height = 6)
save_plot(plot_dist_all, "distribution_all_cause_daily_counts.pdf", width = 10, height = 6)

# -------------------------------------------------------------------------------------
# 2) Distribution of selected cause-specific daily counts
# Again use endpoint-specific windows
# -------------------------------------------------------------------------------------

cause_dist_dat <- bind_rows(
  panel_causes %>%
    filter(has_deaths) %>%
    transmute(
      death_cvd,
      death_respiratory
    ) %>%
    pivot_longer(everything(), names_to = "cause", values_to = "count"),
  
  panel_causes %>%
    filter(has_ed) %>%
    transmute(
      ed_cvd,
      ed_respiratory
    ) %>%
    pivot_longer(everything(), names_to = "cause", values_to = "count"),
  
  panel_causes %>%
    filter(has_ems) %>%
    transmute(
      ems_respiratory,
      ems_syncope,
      ems_neuro
    ) %>%
    pivot_longer(everything(), names_to = "cause", values_to = "count")
)

plot_dist_cause <- ggplot(cause_dist_dat, aes(x = count)) +
  geom_histogram(bins = 40, fill = "darkred", alpha = 0.8) +
  facet_wrap(~ cause, scales = "free") +
  labs(
    title = "Distribution of Cause-Specific Daily Counts",
    x = "Daily Count",
    y = "Frequency"
  )

save_plot(plot_dist_cause, "distribution_cause_specific_daily_counts.png", width = 12, height = 8)
save_plot(plot_dist_cause, "distribution_cause_specific_daily_counts.pdf", width = 12, height = 8)

# -------------------------------------------------------------------------------------
# 3) Time series by outcome
# Use daily totals across community areas and restrict each outcome to its true window
# -------------------------------------------------------------------------------------

ts_deaths <- panel_causes %>%
  filter(has_deaths) %>%
  group_by(event_date) %>%
  summarise(count = sum(deaths, na.rm = TRUE), outcome = "Mortality", .groups = "drop")

ts_ed <- panel_causes %>%
  filter(has_ed) %>%
  group_by(event_date) %>%
  summarise(count = sum(ed_visits, na.rm = TRUE), outcome = "ED Visits", .groups = "drop")

ts_ems <- panel_causes %>%
  filter(has_ems) %>%
  group_by(event_date) %>%
  summarise(count = sum(ems_calls, na.rm = TRUE), outcome = "EMS Calls", .groups = "drop")

plot_ts_dat <- bind_rows(ts_deaths, ts_ed, ts_ems) %>%
  mutate(plot_year = year(event_date))

plot_time_series <- ggplot(plot_ts_dat, aes(x = event_date, y = count, color = outcome)) +
  geom_line(alpha = 0.7, linewidth = 0.35) +
  facet_wrap(~ plot_year, scales = "free_x") +
  labs(
    title = "Daily Time Series by Outcome",
    x = "Date",
    y = "Total Daily Count",
    color = NULL
  ) +
  theme(
    axis.text.x = element_text(size = 7)
  )

save_plot(plot_time_series, "daily_time_series_by_outcome.png", width = 16, height = 10)
save_plot(plot_time_series, "daily_time_series_by_outcome.pdf", width = 16, height = 10)

# -------------------------------------------------------------------------------------
# 4) All-cause maps for all three endpoints
# For direct comparison, use the common overlap window: 2019+
# -------------------------------------------------------------------------------------

map_dat_allcause <- panel_causes %>%
  filter(year >= 2019 & year < 2023 & month %in% 5:9) %>%
  group_by(community) %>%
  summarise(
    deaths = mean(deaths, na.rm = TRUE),
    ed = mean(ed_visits, na.rm = TRUE),
    ems = mean(ems_calls, na.rm = TRUE),
    .groups = "drop"
  )

ca_map_allcause <- ca %>%
  left_join(map_dat_allcause, by = "community")

map_deaths <- ggplot(ca_map_allcause) +
  geom_sf(aes(fill = deaths), color = NA) +
  scale_fill_viridis_c(option = "plasma", labels = label_number()) +
  theme_void() +
  labs(
    title = "Average Daily Mortality",
    subtitle = "2019–2022 overlap window",
    fill = "Count"
  )

map_ed <- ggplot(ca_map_allcause) +
  geom_sf(aes(fill = ed), color = NA) +
  scale_fill_viridis_c(option = "plasma", labels = label_number()) +
  theme_void() +
  labs(
    title = "Average Daily ED Visits",
    subtitle = "2019–2022 overlap window",
    fill = "Count"
  )

map_ems <- ggplot(ca_map_allcause) +
  geom_sf(aes(fill = ems), color = NA) +
  scale_fill_viridis_c(option = "plasma", labels = label_number()) +
  theme_void() +
  labs(
    title = "Average Daily EMS Calls",
    subtitle = "2019–2022 overlap window",
    fill = "Count"
  )

map_allcause_combo <- map_deaths | map_ed | map_ems

save_plot(map_deaths, "map_allcause_mortality.png", width = 7, height = 7)
save_plot(map_ed, "map_allcause_ed.png", width = 7, height = 7)
save_plot(map_ems, "map_allcause_ems.png", width = 7, height = 7)
save_plot(map_allcause_combo, "map_allcause_three_panel.png", width = 16, height = 6)
save_plot(map_allcause_combo, "map_allcause_three_panel.pdf", width = 16, height = 6)

# -------------------------------------------------------------------------------------
# 5) Cause-specific map example: CVD
# Use common overlap window for cross-endpoint comparability
# -------------------------------------------------------------------------------------

map_dat_cvd <- panel_causes %>%
  filter(year >= 2019 & year < 2023 & month %in% 5:9) %>%
  group_by(community) %>%
  summarise(
    death_cvd = mean(death_cvd, na.rm = TRUE),
    ed_cvd    = mean(ed_cvd, na.rm = TRUE),
    ems_cvd   = mean(ems_cvd, na.rm = TRUE),
    .groups = "drop"
  )

ca_map_cvd <- ca %>%
  left_join(map_dat_cvd, by = "community")

map_death_cvd <- ggplot(ca_map_cvd) +
  geom_sf(aes(fill = death_cvd), color = NA) +
  scale_fill_viridis_c(labels = label_number()) +
  theme_void() +
  labs(
    title = "CVD Mortality",
    subtitle = "Mean daily count, 2019–2022",
    fill = "Count"
  )

map_ed_cvd <- ggplot(ca_map_cvd) +
  geom_sf(aes(fill = ed_cvd), color = NA) +
  scale_fill_viridis_c(labels = label_number()) +
  theme_void() +
  labs(
    title = "CVD ED Visits",
    subtitle = "Mean daily count, 2019–2022",
    fill = "Count"
  )

map_ems_cvd <- ggplot(ca_map_cvd) +
  geom_sf(aes(fill = ems_cvd), color = NA) +
  scale_fill_viridis_c(labels = label_number()) +
  theme_void() +
  labs(
    title = "CVD EMS Calls",
    subtitle = "Mean daily count, 2019–2022",
    fill = "Count"
  )

map_cvd_combo <- map_death_cvd | map_ed_cvd | map_ems_cvd

save_plot(map_death_cvd, "map_cvd_mortality.png", width = 7, height = 7)
save_plot(map_ed_cvd, "map_cvd_ed.png", width = 7, height = 7)
save_plot(map_ems_cvd, "map_cvd_ems.png", width = 7, height = 7)
save_plot(map_cvd_combo, "map_cvd_three_panel.png", width = 16, height = 6)
save_plot(map_cvd_combo, "map_cvd_three_panel.pdf", width = 16, height = 6)

# -------------------------------------------------------------------------------------
# 6) Composite health burden map
# Must use common overlap window
# -------------------------------------------------------------------------------------

map_dat_composite <- panel_causes %>%
  filter(year >= 2019 & year < 2023 & month %in% 5:9) %>%
  group_by(community) %>%
  summarise(
    deaths = mean(deaths, na.rm = TRUE),
    ed = mean(ed_visits, na.rm = TRUE),
    ems = mean(ems_calls, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    z_deaths = as.numeric(scale(deaths)),
    z_ed     = as.numeric(scale(ed)),
    z_ems    = as.numeric(scale(ems)),
    composite = z_deaths + z_ed + z_ems
  )

ca_map_composite <- ca %>%
  left_join(map_dat_composite, by = "community")

map_composite <- ggplot(ca_map_composite) +
  geom_sf(aes(fill = composite), color = NA) +
  scale_fill_viridis_c(labels = label_number(accuracy = 0.1)) +
  theme_void() +
  labs(
    title = "Composite Health Burden (Prototype HVI)",
    subtitle = "Sum of z-scored daily mean mortality, ED, and EMS counts (2019–2022)",
    fill = "Composite"
  )

save_plot(map_composite, "map_composite_health_burden.png", width = 8, height = 7)
save_plot(map_composite, "map_composite_health_burden.pdf", width = 8, height = 7)

# -------------------------------------------------------------------------------------
# Optional QA prints
# -------------------------------------------------------------------------------------

cat("\nSaved figures to:", fig_dir, "\n")

cat("\nComposite map summary:\n")
print(summary(map_dat_composite))

cat("\nAll-cause map summary:\n")
print(summary(map_dat_allcause))

# =====================================================================================
# Warm-season citywide time series, shared overlap window only, with secondary y-axis
# =====================================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(lubridate)
  library(scales)
})

dir.create("figures", showWarnings = FALSE, recursive = TRUE)

# -----------------------------------------------------------------------------
# Restrict to true common analytic window
# -----------------------------------------------------------------------------

ts_dat <- panel_causes %>%
  mutate(
    event_date = as.Date(event_date),
    year = year(event_date),
    month = month(event_date),
    warm_season = month %in% 5:9,
    common_overlap = year >= 2019 & year <= 2022
  ) %>%
  filter(common_overlap, warm_season) %>%
  group_by(event_date) %>%
  summarise(
    deaths = sum(deaths, na.rm = TRUE),
    ed_visits = sum(ed_visits, na.rm = TRUE),
    ems_calls = sum(ems_calls, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(plot_year = year(event_date))

# -----------------------------------------------------------------------------
# Choose a scaling factor for mortality so it is visible on the same plot
# Secondary axis will map back to original mortality counts
# -----------------------------------------------------------------------------

scale_factor <- max(c(ts_dat$ed_visits, ts_dat$ems_calls), na.rm = TRUE) /
  max(ts_dat$deaths, na.rm = TRUE)

# -----------------------------------------------------------------------------
# Plot
# -----------------------------------------------------------------------------

plot_time_series_dual <- ggplot(ts_dat, aes(x = event_date)) +
  geom_line(aes(y = ed_visits, color = "ED Visits"), alpha = 0.8, linewidth = 0.4) +
  geom_line(aes(y = ems_calls, color = "EMS Calls"), alpha = 0.8, linewidth = 0.4) +
  geom_line(aes(y = deaths * scale_factor, color = "Mortality"), alpha = 0.9, linewidth = 0.5) +
  facet_wrap(~ plot_year, scales = "free_x", ncol = 2) +
  scale_color_manual(
    values = c(
      "ED Visits" = "#F8766D",
      "EMS Calls" = "#00BA38",
      "Mortality" = "#619CFF"
    )
  ) +
  scale_y_continuous(
    name = "Total Daily Count (ED / EMS)",
    labels = label_number(),
    sec.axis = sec_axis(
      ~ . / scale_factor,
      name = "Total Daily Count (Mortality)",
      labels = label_number()
    )
  ) +
  scale_x_date(
    date_breaks = "2 months",
    date_labels = "%b %Y"
  ) +
  labs(
    title = "Warm-Season Daily Time Series by Outcome",
    subtitle = "April–October only; shared overlap window 2019–2022",
    x = "Date",
    color = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "right",
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(size = 7, angle = 45, hjust = 1)
  )

plot_time_series_dual

ggsave(
  filename = "figures/daily_time_series_by_outcome_warmseason_2019_2022_dual_axis.png",
  plot = plot_time_series_dual,
  width = 14,
  height = 10,
  dpi = 300
)

ggsave(
  filename = "figures/daily_time_series_by_outcome_warmseason_2019_2022_dual_axis.pdf",
  plot = plot_time_series_dual,
  width = 14,
  height = 10
)

suppressPackageStartupMessages({
  library(zoo)
})

ts_dat_smooth <- ts_dat %>%
  arrange(event_date) %>%
  mutate(
    deaths_ma7 = zoo::rollmean(deaths, 7, fill = NA, align = "center"),
    ed_visits_ma7 = zoo::rollmean(ed_visits, 7, fill = NA, align = "center"),
    ems_calls_ma7 = zoo::rollmean(ems_calls, 7, fill = NA, align = "center")
  )

scale_factor_ma <- max(c(ts_dat_smooth$ed_visits_ma7, ts_dat_smooth$ems_calls_ma7), na.rm = TRUE) /
  max(ts_dat_smooth$deaths_ma7, na.rm = TRUE)

plot_time_series_dual_ma <- ggplot(ts_dat_smooth, aes(x = event_date)) +
  geom_line(aes(y = ed_visits_ma7, color = "ED Visits"), alpha = 0.9, linewidth = 0.5) +
  geom_line(aes(y = ems_calls_ma7, color = "EMS Calls"), alpha = 0.9, linewidth = 0.5) +
  geom_line(aes(y = deaths_ma7 * scale_factor_ma, color = "Mortality"), alpha = 0.95, linewidth = 0.6) +
  facet_wrap(~ plot_year, scales = "free_x", ncol = 2) +
  scale_color_manual(
    values = c(
      "ED Visits" = "#F8766D",
      "EMS Calls" = "#00BA38",
      "Mortality" = "#619CFF"
    )
  ) +
  scale_y_continuous(
    name = "7-day Mean Daily Count (ED / EMS)",
    sec.axis = sec_axis(
      ~ . / scale_factor_ma,
      name = "7-day Mean Daily Count (Mortality)"
    )
  ) +
  scale_x_date(
    date_breaks = "2 months",
    date_labels = "%b %Y"
  ) +
  labs(
    title = "Warm-Season Daily Time Series by Outcome",
    subtitle = "7-day moving average; April–October only; shared overlap window 2019–2022",
    x = "Date",
    color = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "right",
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(size = 7, angle = 45, hjust = 1)
  )

plot_time_series_dual_ma

ggsave(
  filename = "figures/daily_time_series_by_outcome_warmseason_2019_2022_dual_axis_ma7.png",
  plot = plot_time_series_dual_ma,
  width = 14,
  height = 10,
  dpi = 300
)





