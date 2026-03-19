# =====================================================================================
# HVI 2.0 | Standardize time and space across mortality, ED, and EMS
# =====================================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(lubridate)
  library(sf)
  library(janitor)
  library(stringr)
})

# -------------------------------------------------------------------------------------
# File paths
# -------------------------------------------------------------------------------------

path_deaths <- "C:/Users/Peter Graffy/Box/HVI2.0/Mortality/all_deaths.csv"
path_ed     <- "C:/Users/Peter Graffy/Box/HVI2.0/ED/ed_outcomes_complete.csv"
path_ems    <- "C:/Users/Peter Graffy/Box/HVI2.0/EMS/new_files_deduplicated_clean.csv"
path_ca     <- "C:/Users/Peter Graffy/Box/HVI2.0/comm_areas.geojson"

# -------------------------------------------------------------------------------------
# Load data
# -------------------------------------------------------------------------------------

deaths <- fread(path_deaths) %>% clean_names()
ed     <- fread(path_ed) %>% clean_names()
ems    <- fread(path_ems) %>% clean_names()

comm_areas <- st_read(path_ca, quiet = TRUE) %>% clean_names()

# -------------------------------------------------------------------------------------
# Inspect community area columns
# -------------------------------------------------------------------------------------

cat("\nCommunity area shapefile columns:\n")
print(names(comm_areas))

cat("\nCRS of community areas:\n")
print(st_crs(comm_areas))

# -------------------------------------------------------------------------------------
# 1. Standardize mortality geography
# -------------------------------------------------------------------------------------
# Assumes mortality already has ca_name
# Create a harmonized community variable

deaths <- deaths %>%
  mutate(
    community = ca_name,
    community = str_squish(as.character(community)),
    source = "mortality"
  )

# -------------------------------------------------------------------------------------
# 2. Standardize ED geography
# -------------------------------------------------------------------------------------
# Assumes ED already has community

ed <- ed %>%
  mutate(
    community = str_squish(as.character(community)),
    source = "ed"
  )

# -------------------------------------------------------------------------------------
# 3. EMS spatial join to community areas
# -------------------------------------------------------------------------------------
# EMS has:
#   scene_gps_latitude_escene_11
#   scene_gps_longitude_escene_11

# Convert to numeric safely
ems <- ems %>%
  mutate(
    scene_gps_latitude_escene_11  = suppressWarnings(as.numeric(scene_gps_latitude_escene_11)),
    scene_gps_longitude_escene_11 = suppressWarnings(as.numeric(scene_gps_longitude_escene_11))
  )

# Check missing / impossible coordinates
cat("\nEMS coordinate summary:\n")
print(summary(ems$scene_gps_latitude_escene_11))
print(summary(ems$scene_gps_longitude_escene_11))

ems_valid <- ems %>%
  filter(
    !is.na(scene_gps_latitude_escene_11),
    !is.na(scene_gps_longitude_escene_11),
    between(scene_gps_latitude_escene_11,  41, 43),
    between(scene_gps_longitude_escene_11, -89, -87)
  )

cat("\nEMS rows before coordinate filtering:", nrow(ems), "\n")
cat("EMS rows after coordinate filtering:", nrow(ems_valid), "\n")

# Convert EMS to sf points using WGS84
ems_sf <- st_as_sf(
  ems_valid,
  coords = c("scene_gps_longitude_escene_11", "scene_gps_latitude_escene_11"),
  crs = 4326,
  remove = FALSE
)

# Transform EMS CRS to match community areas if needed
if (st_crs(ems_sf) != st_crs(comm_areas)) {
  ems_sf <- st_transform(ems_sf, st_crs(comm_areas))
}

# Spatial join
ems_joined <- st_join(ems_sf, comm_areas, join = st_within, left = TRUE)

cat("\nEMS rows with matched community area:\n")
print(sum(!is.na(ems_joined$community)))
cat("EMS rows without matched community area:\n")
print(sum(is.na(ems_joined$community)))

# Optional fallback: nearest feature for unmatched points
unmatched_idx <- which(is.na(ems_joined$community))

if (length(unmatched_idx) > 0) {
  nearest_idx <- st_nearest_feature(ems_joined[unmatched_idx, ], comm_areas)
  ems_joined$community[unmatched_idx] <- comm_areas$community[nearest_idx]
}

cat("\nEMS rows still unmatched after nearest-feature fallback:\n")
print(sum(is.na(ems_joined$community)))

ems_joined <- ems_joined %>%
  mutate(
    community = str_squish(as.character(community)),
    source = "ems"
  )

# Drop geometry for downstream tabular work
ems_final <- ems_joined %>% st_drop_geometry()

# ------------------------------------------------------------------------------
# Source-specific parsing helpers
# ------------------------------------------------------------------------------

clean_chr_datetime <- function(x) {
  x <- as.character(x)
  x <- str_squish(x)
  x[x %in% c("", "NA", "N/A", "NULL", "null")] <- NA_character_
  x
}

parse_date_ymd_safe <- function(x) {
  x <- clean_chr_datetime(x)
  ymd(x, quiet = TRUE)
}

parse_datetime_mdy_hms_ampm_safe <- function(x, tz = "America/Chicago") {
  x <- clean_chr_datetime(x)
  
  # Primary parse for strings like "4/28/2025 8:24:12 PM"
  out <- mdy_hms(x, tz = tz, quiet = TRUE)
  
  # Fallback for rows missing seconds, e.g. "4/28/2025 8:24 PM"
  need_fallback <- is.na(out) & !is.na(x)
  if (any(need_fallback)) {
    out[need_fallback] <- mdy_hm(x[need_fallback], tz = tz, quiet = TRUE)
  }
  
  out
}

clean_chr_datetime <- function(x) {
  x <- as.character(x)
  x <- str_squish(x)
  x[x %in% c("", "NA", "N/A", "NULL", "null")] <- NA_character_
  x
}

parse_ems_datetime <- function(x, tz = "America/Chicago") {
  x <- clean_chr_datetime(x)
  
  out <- rep(as.POSIXct(NA, tz = tz), length(x))
  
  # 1) m/d/Y H:M:S AM/PM
  idx <- !is.na(x)
  out[idx] <- mdy_hms(x[idx], tz = tz, quiet = TRUE)
  
  # 2) m/d/Y H:M AM/PM
  idx <- is.na(out) & !is.na(x)
  if (any(idx)) {
    out[idx] <- mdy_hm(x[idx], tz = tz, quiet = TRUE)
  }
  
  # 3) Y-m-d H:M:S
  idx <- is.na(out) & !is.na(x)
  if (any(idx)) {
    out[idx] <- ymd_hms(x[idx], tz = tz, quiet = TRUE)
  }
  
  # 4) Y-m-d H:M
  idx <- is.na(out) & !is.na(x)
  if (any(idx)) {
    out[idx] <- ymd_hm(x[idx], tz = tz, quiet = TRUE)
  }
  
  out
}

extract_ems_date_fallback <- function(x) {
  x <- clean_chr_datetime(x)
  
  out <- rep(as.Date(NA), length(x))
  
  # format like 2019-04-20 22:32:00 -> first 10 chars is date
  idx_ymd <- !is.na(x) & str_detect(x, "^\\d{4}-\\d{2}-\\d{2}")
  if (any(idx_ymd)) {
    out[idx_ymd] <- ymd(str_sub(x[idx_ymd], 1, 10), quiet = TRUE)
  }
  
  # format like 4/28/2025 8:24:12 PM -> extract m/d/Y before space
  idx_mdy <- !is.na(x) & str_detect(x, "^\\d{1,2}/\\d{1,2}/\\d{4}")
  if (any(idx_mdy)) {
    first_token <- str_extract(x[idx_mdy], "^\\d{1,2}/\\d{1,2}/\\d{4}")
    out[idx_mdy] <- mdy(first_token, quiet = TRUE)
  }
  
  out
}

# -------------------------------------------------------------------------------------
# 5. Standardize date fields
# -------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# Mortality
# dod example: 1993-12-16
# ------------------------------------------------------------------------------

deaths <- deaths %>%
  mutate(
    event_date = parse_date_ymd_safe(dod),
    event_datetime = as.POSIXct(event_date, tz = "America/Chicago")
  )

# ------------------------------------------------------------------------------
# ED
# enc_admit_date example: 2014-02-22
# ------------------------------------------------------------------------------

ed <- ed %>%
  mutate(
    event_date = parse_date_ymd_safe(enc_admit_date),
    event_datetime = as.POSIXct(event_date, tz = "America/Chicago")
  )

# ------------------------------------------------------------------------------
# EMS
# patient_assessment_date_time_exam_03 example: 4/28/2025 8:24:12 PM
# ------------------------------------------------------------------------------

# EMS
ems_final <- ems_final %>%
  mutate(
    ems_raw_datetime = clean_chr_datetime(patient_assessment_date_time_eexam_03),
    event_datetime = parse_ems_datetime(ems_raw_datetime),
    event_date = as.Date(event_datetime),
    event_date = coalesce(event_date, extract_ems_date_fallback(ems_raw_datetime))
  )

# -------------------------------------------------------------------------------------
# 6. Quick QA
# -------------------------------------------------------------------------------------

qa_summary <- function(df, name) {
  cat("\n=================================================\n")
  cat("QA:", name, "\n")
  cat("=================================================\n")
  cat("Rows:", nrow(df), "\n")
  
  if ("community" %in% names(df)) {
    cat("Missing community:", sum(is.na(df$community)), "\n")
    cat("Unique communities:", length(unique(df$community)), "\n")
  }
  
  if ("event_date" %in% names(df)) {
    cat("Missing event_date:", sum(is.na(df$event_date)), "\n")
    print(summary(df$event_date))
  }
  
  if ("event_datetime" %in% names(df)) {
    cat("Missing event_datetime:", sum(is.na(df$event_datetime)), "\n")
    print(head(df$event_datetime, 5))
  }
}

qa_summary(deaths, "Mortality")
qa_summary(ed, "ED")
qa_summary(ems_final, "EMS")

cat("Missing EMS datetime:", sum(is.na(ems_final$event_datetime)), "\n")
cat("Missing EMS date:", sum(is.na(ems_final$event_date)), "\n")

summary(ems_final$event_date)

ems_final %>%
  mutate(
    raw_type = case_when(
      is.na(ems_raw_datetime) ~ "missing_raw",
      str_detect(ems_raw_datetime, "^\\d{1,2}/\\d{1,2}/\\d{4}") ~ "mdy_format",
      str_detect(ems_raw_datetime, "^\\d{4}-\\d{2}-\\d{2}") ~ "ymd_format",
      TRUE ~ "other"
    )
  ) %>%
  count(raw_type)

raw_col <- "patient_assessment_date_time_eexam_03"

# Inspect raw values exactly as read
ems %>%
  mutate(raw = as.character(.data[[raw_col]])) %>%
  summarise(
    n_total = n(),
    n_na = sum(is.na(raw)),
    n_blank = sum(raw == "", na.rm = TRUE),
    n_whitespace = sum(trimws(raw) == "", na.rm = TRUE),
    n_literal_na = sum(raw %in% c("NA", "N/A", "NULL", "null"), na.rm = TRUE)
  )

# -------------------------------------------------------------------------------------
# 7. Save standardized files
# -------------------------------------------------------------------------------------

dir.create("data/derived", recursive = TRUE, showWarnings = FALSE)

fwrite(deaths,   "data/derived/mortality_standardized.csv")
fwrite(ed,       "data/derived/ed_standardized.csv")
fwrite(ems_final,"data/derived/ems_standardized.csv")

cat("\nSaved standardized files to data/derived/\n")