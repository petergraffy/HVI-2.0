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

# -------------------------------------------------------------------------------------
# File paths
# -------------------------------------------------------------------------------------

path_panel  <- "data/derived/community_day_panel_with_climate.csv"
out_panel   <- "data/derived/community_day_panel_with_climate_causes.csv"

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


library(ggplot2)
library(dplyr)

plot_dat <- panel_causes %>%
  filter(year(event_date) >= 2010)

plot_dist_all <- plot_dat %>%
  select(deaths, ed_visits, ems_calls) %>%
  pivot_longer(everything(), names_to = "outcome", values_to = "count")

ggplot(plot_dist_all, aes(x = count)) +
  geom_histogram(bins = 50, fill = "steelblue", alpha = 0.7) +
  facet_wrap(~ outcome, scales = "free") +
  labs(
    title = "Distribution of Daily Counts by Outcome",
    x = "Daily Count",
    y = "Frequency"
  ) +
  theme_minimal()

cause_vars <- c(
  "death_cvd", "death_respiratory",
  "ed_cvd", "ed_respiratory",
  "ems_respiratory", "ems_syncope", "ems_neuro"
)

plot_dat <- panel_causes %>%
  select(all_of(cause_vars)) %>%
  pivot_longer(everything(), names_to = "cause", values_to = "count")

ggplot(plot_dat, aes(x = count)) +
  geom_histogram(bins = 40, fill = "darkred", alpha = 0.7) +
  facet_wrap(~ cause, scales = "free") +
  theme_minimal() +
  labs(title = "Distribution of Cause-Specific Daily Counts")

ggplot(panel_causes, aes(x = event_date)) +
  geom_line(aes(y = deaths, color = "Mortality"), alpha = 0.7) +
  geom_line(aes(y = ed_visits, color = "ED"), alpha = 0.7) +
  geom_line(aes(y = ems_calls, color = "EMS"), alpha = 0.7) +
  facet_wrap(~ year(event_date), scales = "free_x") +
  theme_minimal() +
  labs(title = "Daily Time Series by Outcome")

map_dat <- panel_causes %>%
  filter(year(event_date) >= 2015) %>%
  group_by(community) %>%
  summarise(
    deaths = mean(deaths, na.rm = TRUE),
    ed     = mean(ed_visits, na.rm = TRUE),
    ems    = mean(ems_calls, na.rm = TRUE),
    .groups = "drop"
  )

ca <- comm_areas %>% select(community)

ca_map <- ca %>%
  mutate(community = str_to_title(community)) %>%
  left_join(map_dat, by = "community")

library(ggplot2)

map_deaths <- ggplot(ca_map) +
  geom_sf(aes(fill = deaths), color = NA) +
  scale_fill_viridis_c(option = "plasma") +
  theme_void() +
  labs(title = "Average Daily Mortality")

map_ed <- ggplot(ca_map) +
  geom_sf(aes(fill = ed), color = NA) +
  scale_fill_viridis_c(option = "plasma") +
  theme_void() +
  labs(title = "Average Daily ED Visits")

map_ems <- ggplot(ca_map) +
  geom_sf(aes(fill = ems), color = NA) +
  scale_fill_viridis_c(option = "plasma") +
  theme_void() +
  labs(title = "Average Daily EMS Calls")

library(patchwork)

map_deaths | map_ed | map_ems

map_dat <- panel_causes %>%
  filter(year(event_date) >= 2015) %>%
  group_by(community) %>%
  summarise(
    death_cvd = mean(death_cvd, na.rm = TRUE),
    ed_cvd    = mean(ed_cvd, na.rm = TRUE),
    ems_cvd   = mean(ems_cvd, na.rm = TRUE),
    .groups = "drop"
  )

ca_map <- ca %>%
  left_join(map_dat, by = "community")

ggplot(ca_map) +
  geom_sf(aes(fill = death_cvd), color = NA) +
  scale_fill_viridis_c() +
  theme_void() +
  labs(title = "CVD Mortality (Mean Daily)")

map_dat <- panel_causes %>%
  filter(year(event_date) >= 2015) %>%
  group_by(community) %>%
  summarise(
    z_deaths = scale(mean(deaths)),
    z_ed     = scale(mean(ed_visits)),
    z_ems    = scale(mean(ems_calls)),
    .groups = "drop"
  ) %>%
  mutate(
    composite = z_deaths + z_ed + z_ems
  )

ca_map <- ca %>%
  left_join(map_dat, by = "community")

ggplot(ca_map) +
  geom_sf(aes(fill = composite), color = NA) +
  scale_fill_viridis_c() +
  theme_void() +
  labs(title = "Composite Health Burden (Prototype HVI)")

























