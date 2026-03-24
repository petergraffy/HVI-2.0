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

# install.packages("appeears")  # if needed
 library(appeears)

cat("\n--- AppEEARS scaffold below ---\n")
cat("After installing/loading {appeears}, first query products and layers.\n")

options(keyring_backend = "file")

rs_set_key(
  user = "grafpe01",
  password = "A%98,Syes-,3hcb"
)

login <- rs_login(user = "grafpe01")
rs_products()

prods <- rs_products()
prods %>% filter(grepl("DAYMET", ProductAndVersion))

layers <- rs_layers(product = "DAYMET.004")
print(layers$Layer)

# --- 1) clean polygons (as you had) ---
comm_areas <- st_read("C:/Users/Peter Graffy/Box/HVI2.0/comm_areas.geojson", quiet = TRUE)

if (!"community" %in% names(comm_areas)) {
  nm <- intersect(c("community","ca_name","commarea","community_area"), names(comm_areas))[1]
  comm_areas <- comm_areas %>% mutate(community = .data[[nm]])
}

comm_areas_clean <- comm_areas %>%
  st_make_valid() %>%
  st_transform(4326) %>%
  st_collection_extract("POLYGON", warn = FALSE) %>%
  st_cast("MULTIPOLYGON", warn = FALSE) %>%
  filter(!st_is_empty(geometry)) %>%
  filter(st_is_valid(geometry)) %>%
  mutate(id = community)

# sanity checks
stopifnot(st_crs(comm_areas_clean)$epsg == 4326)
stopifnot(all(st_is_valid(comm_areas_clean)))
stopifnot(!any(st_is_empty(comm_areas_clean)))

# --- 2) convert sf -> GeoJSON list for AppEEARS ---
geojson_txt <- geojsonsf::sf_geojson(comm_areas_clean)
geo_obj <- jsonlite::fromJSON(geojson_txt, simplifyVector = FALSE)

# --- 3) build task ---
task <- list(
  task_type = "area",
  task_name = "hvi2_daymet_2024_chicago_community_areas",
  params = list(
    dates = list(list(startDate = "01-01-2024", endDate = "12-31-2024")),
    layers = list(
      list(product = "DAYMET.004", layer = "tmax"),
      list(product = "DAYMET.004", layer = "tmin"),
      list(product = "DAYMET.004", layer = "vp")
    ),
    geo = geo_obj,
    output = list(
      format = list(type = "netcdf4"),
      projection = "geographic"
    )
  )
)

# --- 4) submit ---
req <- rs_request(
  request = task,
  user = "grafpe01",
  transfer = TRUE,
  path = "data/raw/appeears"
)

req


















