# =====================================================================================
# HVI 2.0 | Climate QA + 2024 Daymet extraction by Chicago community area
# =====================================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(lubridate)
  library(stringr)
  library(janitor)
  library(sf)
})

# -------------------------------------------------------------------------------------
# File paths
# -------------------------------------------------------------------------------------

path_climate_hist <- "C:/Users/Peter Graffy/Box/HVI2.0/Climate/CA_temps_rh_90-23.csv"
path_ca           <- "C:/Users/Peter Graffy/Box/HVI2.0/comm_areas.geojson"

# output folders
dir.create("data/derived/climate", recursive = TRUE, showWarnings = FALSE)
dir.create("data/raw/appeears", recursive = TRUE, showWarnings = FALSE)

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

safe_date_parse <- function(x) {
  x <- as.character(x)
  x <- str_squish(x)
  x[x %in% c("", "NA", "N/A", "NULL", "null")] <- NA_character_
  
  # try a few common date forms
  out <- suppressWarnings(ymd(x, quiet = TRUE))
  idx <- is.na(out) & !is.na(x)
  if (any(idx)) out[idx] <- suppressWarnings(mdy(x[idx], quiet = TRUE))
  idx <- is.na(out) & !is.na(x)
  if (any(idx)) out[idx] <- suppressWarnings(dmy(x[idx], quiet = TRUE))
  out
}

qa_missing_table <- function(df) {
  tibble(
    variable = names(df),
    n_missing = sapply(df, function(x) sum(is.na(x))),
    pct_missing = sapply(df, function(x) mean(is.na(x))) * 100
  ) %>%
    arrange(desc(pct_missing))
}

# -------------------------------------------------------------------------------------
# 1) Read historical community-area climate file
# -------------------------------------------------------------------------------------

clim_hist <- fread(path_climate_hist) %>%
  clean_names()

cat("\nHistorical climate file dimensions:\n")
print(dim(clim_hist))

cat("\nColumn names:\n")
print(names(clim_hist))

cat("\nColumn classes:\n")
print(sapply(clim_hist, class))

cat("\nMissingness summary:\n")
print(qa_missing_table(clim_hist))

cat("\nPreview:\n")
print(head(clim_hist, 5))

# -------------------------------------------------------------------------------------
# 2) Detect likely key columns
# -------------------------------------------------------------------------------------

cat("\nPotential date columns:\n")
print(names(clim_hist)[grepl("date|day|time", names(clim_hist), ignore.case = TRUE)])

cat("\nPotential community columns:\n")
print(names(clim_hist)[grepl("community|ca|comm", names(clim_hist), ignore.case = TRUE)])

cat("\nPotential temperature / humidity columns:\n")
print(names(clim_hist)[grepl("tmax|tmin|tmean|temp|rh|humid", names(clim_hist), ignore.case = TRUE)])

# -------------------------------------------------------------------------------------
# 3) Standardize obvious core fields
#    EDIT THESE IF YOUR COLUMN NAMES DIFFER
# -------------------------------------------------------------------------------------
# Assumed likely columns:
#   - community
#   - date
#   - tmax
#   - tmin
#   - tmean
#   - rh
# Adjust as needed after inspecting names(clim_hist).

# detect a likely date column
date_candidates <- c("date", "day", "observation_date")
date_col <- date_candidates[date_candidates %in% names(clim_hist)][1]

if (is.na(date_col)) {
  stop("Could not auto-detect a date column in CA_temps_rh_90-23.csv. Inspect names(clim_hist) and set date_col manually.")
}

# detect a likely community column
community_candidates <- c("community", "ca_name", "ca", "commarea", "community_area")
community_col <- community_candidates[community_candidates %in% names(clim_hist)][1]

if (is.na(community_col)) {
  stop("Could not auto-detect a community column in CA_temps_rh_90-23.csv. Inspect names(clim_hist) and set community_col manually.")
}

clim_hist <- clim_hist %>%
  mutate(
    community = standardize_community(.data[[community_col]]),
    event_date = safe_date_parse(.data[[date_col]])
  )

# -------------------------------------------------------------------------------------
# 4) Historical climate QA
# -------------------------------------------------------------------------------------

cat("\nHistorical climate QA:\n")
cat("Missing event_date:", sum(is.na(clim_hist$event_date)), "\n")
cat("Missing community:", sum(is.na(clim_hist$community)), "\n")
cat("Unique communities:", length(unique(clim_hist$community)), "\n")
print(summary(clim_hist$event_date))

# numeric coercion for common climate fields if present
for (v in c("tmax", "tmin", "tmean", "humidity")) {
  if (v %in% names(clim_hist)) {
    clim_hist[[v]] <- suppressWarnings(as.numeric(clim_hist[[v]]))
  }
}

# # quick range checks
# range_check <- tibble(
#   variable = intersect(c("tmax", "tmin", "tmean", "humidity"), names(clim_hist)),
#   min_value = sapply(clim_hist[intersect(c("tmax", "tmin", "tmean", "humidity"), names(clim_hist))], min, na.rm = TRUE),
#   p01 = sapply(clim_hist[intersect(c("tmax", "tmin", "tmean", "humidity"), names(clim_hist))], quantile, probs = 0.01, na.rm = TRUE),
#   median = sapply(clim_hist[intersect(c("tmax", "tmin", "tmean", "humidity"), names(clim_hist))], median, na.rm = TRUE),
#   p99 = sapply(clim_hist[intersect(c("tmax", "tmin", "tmean", "humidity"), names(clim_hist))], quantile, probs = 0.99, na.rm = TRUE),
#   max_value = sapply(clim_hist[intersect(c("tmax", "tmin", "tmean", "humidity"), names(clim_hist))], max, na.rm = TRUE)
# )
# 
# cat("\nRange checks:\n")
# print(range_check)

# duplicate check at community-day
dup_hist <- clim_hist %>%
  count(community, event_date, name = "n") %>%
  filter(n > 1)

cat("\nHistorical duplicate community-day rows:", nrow(dup_hist), "\n")

# annual completeness check
hist_completeness <- clim_hist %>%
  mutate(year = year(event_date)) %>%
  count(year, community, name = "n_days") %>%
  group_by(year) %>%
  summarise(
    communities = n_distinct(community),
    min_days = min(n_days, na.rm = TRUE),
    median_days = median(n_days, na.rm = TRUE),
    max_days = max(n_days, na.rm = TRUE),
    .groups = "drop"
  )

cat("\nHistorical annual completeness summary:\n")
print(hist_completeness)

fwrite(clim_hist, "data/derived/climate/ca_temps_rh_1990_2023_standardized.csv")
fwrite(hist_completeness, "data/derived/climate/ca_temps_rh_1990_2023_completeness.csv")

# -------------------------------------------------------------------------------------
# 5) Prepare community polygons for AppEEARS 2024 request
# -------------------------------------------------------------------------------------

comm_areas <- st_read(path_ca, quiet = TRUE) %>%
  clean_names()

cat("\nCommunity area shapefile columns:\n")
print(names(comm_areas))

# standardize a community field
if (!"community" %in% names(comm_areas)) {
  # try to infer
  comm_candidates <- c("community", "ca_name", "commarea", "community_area")
  shp_comm_col <- comm_candidates[comm_candidates %in% names(comm_areas)][1]
  
  if (is.na(shp_comm_col)) {
    stop("Could not find a community-name field in the community-area GeoJSON.")
  }
  
  comm_areas <- comm_areas %>%
    mutate(community = .data[[shp_comm_col]])
}

comm_areas <- comm_areas %>%
  mutate(community = standardize_community(community)) %>%
  st_make_valid()

# AppEEARS area requests accept polygon/vector inputs; keeping a clean GeoJSON is useful.
# Save a request-ready copy in EPSG:4326
comm_areas_wgs84 <- st_transform(comm_areas, 4326)

st_write(
  comm_areas_wgs84,
  "data/raw/appeears/chicago_community_areas_wgs84.geojson",
  delete_dsn = TRUE,
  quiet = TRUE
)

# -------------------------------------------------------------------------------------
# 6) AppEEARS scaffolding for 2024 Daymet area request
# -------------------------------------------------------------------------------------
# Notes:
# - AppEEARS supports area requests by polygon and time range.
# - The CRAN {appeears} package provides helpers such as rs_products(),
#   rs_layers(), and rs_request().
# - You will need valid Earthdata credentials available in your session/environment.
# -------------------------------------------------------------------------------------

# # install.packages("appeears")  # if needed
#  library(appeears)
# 
# cat("\n--- AppEEARS scaffold below ---\n")
# cat("After installing/loading {appeears}, first query products and layers.\n")
# 
# options(keyring_backend = "file")
# 
# rs_set_key(
#   user = "grafpe01",
#   password = "A%98,Syes-,3hcb"
# )
# 
# login <- rs_login(user = "grafpe01")
# rs_products()
# 
# prods <- rs_products()
# prods %>% filter(grepl("DAYMET", ProductAndVersion))
# 
# layers <- rs_layers(product = "DAYMET.004")
# print(layers$Layer)
# 
# # --- 1) clean polygons (as you had) ---
# comm_areas <- st_read("C:/Users/Peter Graffy/Box/HVI2.0/comm_areas.geojson", quiet = TRUE)
# 
# if (!"community" %in% names(comm_areas)) {
#   nm <- intersect(c("community","ca_name","commarea","community_area"), names(comm_areas))[1]
#   comm_areas <- comm_areas %>% mutate(community = .data[[nm]])
# }
# 
# comm_areas_clean <- comm_areas %>%
#   st_make_valid() %>%
#   st_transform(4326) %>%
#   st_collection_extract("POLYGON", warn = FALSE) %>%
#   st_cast("MULTIPOLYGON", warn = FALSE) %>%
#   filter(!st_is_empty(geometry)) %>%
#   filter(st_is_valid(geometry)) %>%
#   mutate(id = community)
# 
# # sanity checks
# stopifnot(st_crs(comm_areas_clean)$epsg == 4326)
# stopifnot(all(st_is_valid(comm_areas_clean)))
# stopifnot(!any(st_is_empty(comm_areas_clean)))
# 
# # --- 2) convert sf -> GeoJSON list for AppEEARS ---
# geojson_txt <- geojsonsf::sf_geojson(comm_areas_clean)
# geo_obj <- jsonlite::fromJSON(geojson_txt, simplifyVector = FALSE)
# 
# # --- 3) build task ---
# task <- list(
#   task_type = "area",
#   task_name = "hvi2_daymet_2024_chicago_community_areas",
#   params = list(
#     dates = list(list(startDate = "01-01-2024", endDate = "12-31-2024")),
#     layers = list(
#       list(product = "DAYMET.004", layer = "tmax"),
#       list(product = "DAYMET.004", layer = "tmin"),
#       list(product = "DAYMET.004", layer = "vp")
#     ),
#     geo = geo_obj,
#     output = list(
#       format = list(type = "netcdf4"),
#       projection = "geographic"
#     )
#   )
# )
# 
# # --- 4) submit ---
# req <- rs_request(
#   request = task,
#   user = "grafpe01",
#   transfer = TRUE,
#   path = "data/raw/appeears"
# )
# 
# req

# =====================================================================================
# HVI 2.0 | Aggregate 2024 Daymet NetCDF by community area using terra
# =====================================================================================

suppressPackageStartupMessages({
  library(terra)
  library(sf)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(janitor)
  library(lubridate)
  library(data.table)
})

# -------------------------------------------------------------------------------------
# File paths
# -------------------------------------------------------------------------------------

nc_path  <- "C:/Users/Peter Graffy/Box/HVI2.0/Climate/DAYMET.004_1km_aid0001.nc"
ca_path  <- "C:/Users/Peter Graffy/Box/HVI2.0/comm_areas.geojson"
out_csv  <- "C:/Users/Peter Graffy/Box/HVI2.0/Climate/ca_daymet_2024_aggregated.csv"

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

to_long_climate <- function(extract_df, value_name, dates) {
  extract_df %>%
    rename(ID = 1) %>%
    pivot_longer(-c(ID, community), names_to = "layer", values_to = value_name) %>%
    group_by(ID, community) %>%
    mutate(event_date = dates[row_number()]) %>%
    ungroup() %>%
    select(community, event_date, all_of(value_name))
}

# -------------------------------------------------------------------------------------
# Read polygons
# -------------------------------------------------------------------------------------

comm_areas <- st_read(ca_path, quiet = TRUE) %>%
  clean_names()

if (!"community" %in% names(comm_areas)) {
  nm <- intersect(c("community", "ca_name", "commarea", "community_area"), names(comm_areas))[1]
  if (is.na(nm)) stop("Could not find a community field in the shapefile.")
  comm_areas <- comm_areas %>% mutate(community = .data[[nm]])
}

comm_areas <- comm_areas %>%
  mutate(community = standardize_community(community)) %>%
  st_make_valid() %>%
  st_transform(4326)

# convert to terra vector
comm_vect <- vect(comm_areas)

# -------------------------------------------------------------------------------------
# Read NetCDF variables
# -------------------------------------------------------------------------------------

r_tmax <- rast(nc_path, subds = "tmax")
r_tmin <- rast(nc_path, subds = "tmin")
r_vp   <- rast(nc_path, subds = "vp")

# IMPORTANT:
# AppEEARS exported this with geographic-looking extent but malformed CRS metadata.
# Force CRS to EPSG:4326 so polygons and raster align.
crs(r_tmax) <- "EPSG:4326"
crs(r_tmin) <- "EPSG:4326"
crs(r_vp)   <- "EPSG:4326"

# attach dates
tmax_dates <- as.Date(time(r_tmax))
tmin_dates <- as.Date(time(r_tmin))
vp_dates   <- as.Date(time(r_vp))

# QA
cat("\nRaster QA:\n")
print(r_tmax)
print(global(r_tmax[[1]], fun = c("min", "max", "mean"), na.rm = TRUE))
print(ext(r_tmax))
print(crs(r_tmax))

cat("\nPolygon QA:\n")
print(ext(comm_vect))
print(crs(comm_vect))

# -------------------------------------------------------------------------------------
# Extract area-weighted means by polygon
# -------------------------------------------------------------------------------------

cat("\nExtracting tmax...\n")
tmax_ext <- terra::extract(
  r_tmax,
  comm_vect,
  fun = mean,
  na.rm = TRUE,
  exact = TRUE
)

cat("\nExtracting tmin...\n")
tmin_ext <- terra::extract(
  r_tmin,
  comm_vect,
  fun = mean,
  na.rm = TRUE,
  exact = TRUE
)

cat("\nExtracting vp...\n")
vp_ext <- terra::extract(
  r_vp,
  comm_vect,
  fun = mean,
  na.rm = TRUE,
  exact = TRUE
)

# add community names
tmax_ext$community <- comm_areas$community
tmin_ext$community <- comm_areas$community
vp_ext$community   <- comm_areas$community

# -------------------------------------------------------------------------------------
# Reshape to long
# -------------------------------------------------------------------------------------

tmax_long <- to_long_climate(tmax_ext, "tmax", tmax_dates)
tmin_long <- to_long_climate(tmin_ext, "tmin", tmin_dates)
vp_long   <- to_long_climate(vp_ext,   "vp",   vp_dates)

# -------------------------------------------------------------------------------------
# Join and derive project-specific climate variables
# -------------------------------------------------------------------------------------

b1 <- 610.78
b2 <- 17.269
b3 <- 237.3

climate_2024 <- tmax_long %>%
  left_join(tmin_long, by = c("community", "event_date")) %>%
  left_join(vp_long,   by = c("community", "event_date")) %>%
  mutate(
    tmean = 0.606 * tmax + 0.394 * tmin,
    humidity = 100 * (vp / (b1 * exp((b2 * tmean) / (b3 + tmean))))
  ) %>%
  arrange(community, event_date)

# -------------------------------------------------------------------------------------
# QA
# -------------------------------------------------------------------------------------

cat("\nAggregated climate QA:\n")
cat("Rows:", nrow(climate_2024), "\n")
cat("Unique communities:", n_distinct(climate_2024$community), "\n")
print(summary(climate_2024$event_date))

print(
  climate_2024 %>%
    summarise(
      missing_tmax = sum(is.na(tmax)),
      missing_tmin = sum(is.na(tmin)),
      missing_vp = sum(is.na(vp)),
      missing_tmean = sum(is.na(tmean)),
      missing_humidity = sum(is.na(humidity))
    )
)

print(
  climate_2024 %>%
    summarise(
      tmax_min = min(tmax, na.rm = TRUE),
      tmax_max = max(tmax, na.rm = TRUE),
      tmin_min = min(tmin, na.rm = TRUE),
      tmin_max = max(tmin, na.rm = TRUE),
      vp_min = min(vp, na.rm = TRUE),
      vp_max = max(vp, na.rm = TRUE),
      humidity_p01 = quantile(humidity, 0.01, na.rm = TRUE),
      humidity_p50 = quantile(humidity, 0.50, na.rm = TRUE),
      humidity_p99 = quantile(humidity, 0.99, na.rm = TRUE)
    )
)

dup_check <- climate_2024 %>%
  count(community, event_date) %>%
  filter(n > 1)

cat("Duplicate community-day rows:", nrow(dup_check), "\n")

# -------------------------------------------------------------------------------------
# Save
# -------------------------------------------------------------------------------------

fwrite(climate_2024, out_csv)
cat("\nSaved to:\n", out_csv, "\n")


# =====================================================================================
# HVI 2.0 | Combine historical + 2024 climate and merge into community-day panel
# =====================================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(stringr)
  library(lubridate)
  library(janitor)
})

# -------------------------------------------------------------------------------------
# File paths
# -------------------------------------------------------------------------------------

path_clim_hist  <- "C:/Users/Peter Graffy/Box/HVI2.0/Climate/CA_temps_rh_90-23.csv"
path_clim_2024  <- "C:/Users/Peter Graffy/Box/HVI2.0/Climate/ca_daymet_2024_aggregated.csv"

# update this if your panel is saved somewhere else
path_panel      <- "data/derived/community_day_panel.csv"

out_clim_full   <- "C:/Users/Peter Graffy/Box/HVI2.0/Climate/CA_temps_rh_90-24.csv"
out_panel_clim  <- "data/derived/community_day_panel_with_climate.csv"

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

clean_missing_chr <- function(x) {
  x <- as.character(x)
  x <- str_squish(x)
  x[x %in% c("", "NA", "N/A", "NULL", "null")] <- NA_character_
  x
}

parse_date_safe <- function(x) {
  x <- clean_missing_chr(x)
  
  out <- suppressWarnings(ymd(x, quiet = TRUE))
  idx <- is.na(out) & !is.na(x)
  if (any(idx)) out[idx] <- suppressWarnings(mdy(x[idx], quiet = TRUE))
  idx <- is.na(out) & !is.na(x)
  if (any(idx)) out[idx] <- suppressWarnings(dmy(x[idx], quiet = TRUE))
  
  as.Date(out)
}

# -------------------------------------------------------------------------------------
# Read files
# -------------------------------------------------------------------------------------

clim_hist <- fread(path_clim_hist) %>% clean_names()
clim_2024 <- fread(path_clim_2024) %>% clean_names()
panel     <- fread(path_panel) %>% clean_names()

cat("\nHistorical climate columns:\n")
print(names(clim_hist))

cat("\n2024 climate columns:\n")
print(names(clim_2024))

cat("\nPanel columns:\n")
print(names(panel))

# -------------------------------------------------------------------------------------
# Standardize historical climate file
# -------------------------------------------------------------------------------------
# This assumes the historical file already contains something like:
# community, date, tmax, tmin, tmean, rh or humidity
# Adjust detection below if needed.

# detect community column
hist_comm_col <- c("community", "ca_name", "ca", "commarea", "community_area")
hist_comm_col <- hist_comm_col[hist_comm_col %in% names(clim_hist)][1]

if (is.na(hist_comm_col)) {
  stop("Could not detect community column in historical climate file.")
}

# detect date column
hist_date_col <- c("event_date", "date", "day")
hist_date_col <- hist_date_col[hist_date_col %in% names(clim_hist)][1]

if (is.na(hist_date_col)) {
  stop("Could not detect date column in historical climate file.")
}

clim_hist <- clim_hist %>%
  mutate(
    community = standardize_community(.data[[hist_comm_col]]),
    event_date = parse_date_safe(.data[[hist_date_col]])
  )

# standardize climate variable names
if ("rh" %in% names(clim_hist) && !("humidity" %in% names(clim_hist))) {
  clim_hist <- clim_hist %>% rename(humidity = rh)
}
if ("mean_tmax" %in% names(clim_hist) && !("tmax" %in% names(clim_hist))) {
  clim_hist <- clim_hist %>% rename(tmax = mean_tmax)
}
if ("mean_tmin" %in% names(clim_hist) && !("tmin" %in% names(clim_hist))) {
  clim_hist <- clim_hist %>% rename(tmin = mean_tmin)
}
if ("mean_tmean" %in% names(clim_hist) && !("tmean" %in% names(clim_hist))) {
  clim_hist <- clim_hist %>% rename(tmean = mean_tmean)
}

# coerce numeric
for (v in c("tmax", "tmin", "tmean", "humidity", "rh", "vp")) {
  if (v %in% names(clim_hist)) {
    clim_hist[[v]] <- suppressWarnings(as.numeric(clim_hist[[v]]))
  }
}

# keep only core columns if present
hist_keep <- intersect(c("community", "event_date", "tmax", "tmin", "tmean", "humidity"), names(clim_hist))
clim_hist_std <- clim_hist %>%
  select(all_of(hist_keep))

# -------------------------------------------------------------------------------------
# Standardize 2024 climate file
# -------------------------------------------------------------------------------------

clim_2024 <- clim_2024 %>%
  mutate(
    community = standardize_community(community),
    event_date = parse_date_safe(event_date)
  )

for (v in c("tmax", "tmin", "tmean", "humidity", "vp")) {
  if (v %in% names(clim_2024)) {
    clim_2024[[v]] <- suppressWarnings(as.numeric(clim_2024[[v]]))
  }
}

clim_2024_std <- clim_2024 %>%
  select(any_of(c("community", "event_date", "tmax", "tmin", "tmean", "humidity")))

# -------------------------------------------------------------------------------------
# QA before append
# -------------------------------------------------------------------------------------

cat("\nHistorical climate QA:\n")
cat("Rows:", nrow(clim_hist_std), "\n")
cat("Communities:", n_distinct(clim_hist_std$community), "\n")
print(summary(clim_hist_std$event_date))

cat("\n2024 climate QA:\n")
cat("Rows:", nrow(clim_2024_std), "\n")
cat("Communities:", n_distinct(clim_2024_std$community), "\n")
print(summary(clim_2024_std$event_date))

# duplicates
hist_dups <- clim_hist_std %>%
  count(community, event_date) %>%
  filter(n > 1)

new_dups <- clim_2024_std %>%
  count(community, event_date) %>%
  filter(n > 1)

cat("\nHistorical duplicate community-day rows:", nrow(hist_dups), "\n")
cat("2024 duplicate community-day rows:", nrow(new_dups), "\n")

# -------------------------------------------------------------------------------------
# Remove any accidental overlap before bind
# -------------------------------------------------------------------------------------

clim_hist_std <- clim_hist_std %>%
  filter(event_date < as.Date("2024-01-01"))

clim_2024_std <- clim_2024_std %>%
  filter(event_date >= as.Date("2024-01-01"),
         event_date <= as.Date("2024-12-31"))

# -------------------------------------------------------------------------------------
# Combine full climate series
# -------------------------------------------------------------------------------------

clim_full <- bind_rows(clim_hist_std, clim_2024_std) %>%
  arrange(community, event_date)

# final duplicate check
dup_full <- clim_full %>%
  count(community, event_date) %>%
  filter(n > 1)

cat("\nCombined duplicate community-day rows:", nrow(dup_full), "\n")

# -------------------------------------------------------------------------------------
# Save combined climate file
# -------------------------------------------------------------------------------------

fwrite(clim_full, out_clim_full)

cat("\nSaved combined climate file to:\n", out_clim_full, "\n")

# -------------------------------------------------------------------------------------
# Standardize panel before join
# -------------------------------------------------------------------------------------

panel <- panel %>%
  mutate(
    community = standardize_community(community),
    event_date = parse_date_safe(event_date)
  )

# -------------------------------------------------------------------------------------
# Join climate into panel
# -------------------------------------------------------------------------------------

panel_clim <- panel %>%
  left_join(clim_full, by = c("community", "event_date")) %>%
  arrange(community, event_date)

# -------------------------------------------------------------------------------------
# QA merged panel
# -------------------------------------------------------------------------------------

cat("\nMerged panel QA:\n")
cat("Rows:", nrow(panel_clim), "\n")
cat("Communities:", n_distinct(panel_clim$community), "\n")
print(summary(panel_clim$event_date))

print(
  panel_clim %>%
    summarise(
      missing_tmax = sum(is.na(tmax)),
      missing_tmin = sum(is.na(tmin)),
      missing_tmean = sum(is.na(tmean)),
      missing_humidity = sum(is.na(humidity))
    )
)

# Missing climate by year
print(
  panel_clim %>%
    mutate(year = year(event_date)) %>%
    group_by(year) %>%
    summarise(
      n = n(),
      pct_missing_tmax = mean(is.na(tmax)) * 100,
      pct_missing_tmin = mean(is.na(tmin)) * 100,
      pct_missing_tmean = mean(is.na(tmean)) * 100,
      pct_missing_humidity = mean(is.na(humidity)) * 100,
      .groups = "drop"
    )
)

# Spot-check date coverage
print(
  panel_clim %>%
    group_by(community) %>%
    summarise(
      min_date = min(event_date, na.rm = TRUE),
      max_date = max(event_date, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    head()
)

library(dplyr)
library(zoo)

panel_clim <- panel_clim %>%
  arrange(community, event_date) %>%
  group_by(community) %>%
  mutate(
    tmax = na.approx(tmax, x = event_date, na.rm = FALSE, maxgap = 3),
    tmin = na.approx(tmin, x = event_date, na.rm = FALSE, maxgap = 3),
    tmean = na.approx(tmean, x = event_date, na.rm = FALSE, maxgap = 3),
    humidity = na.approx(humidity, x = event_date, na.rm = FALSE, maxgap = 3)
  ) %>%
  ungroup()


# -------------------------------------------------------------------------------------
# Save merged panel
# -------------------------------------------------------------------------------------

fwrite(panel_clim, out_panel_clim)

cat("\nSaved climate-merged panel to:\n", out_panel_clim, "\n")




