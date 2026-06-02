# ================================================================================================
# HVI 2.0 | Download and aggregate annual tree canopy cover by community area
#
# Primary source:
#   USDA Forest Service / MRLC NLCD Tree Canopy Cover (TCC), CONUS, 30 m annual rasters.
#   https://data.fs.usda.gov/geodata/rastergateway/treecanopycover/
#
# Outputs:
#   - tree_canopy_year_ca.csv
#   - tree_canopy_year_ca_qa.csv
#   - tree_canopy_download_manifest.csv
#
# Runtime options:
#   HVI_TCC_YEARS          comma/range years; default "2019:2022"
#   HVI_TCC_LOCAL_DIR      optional directory of already-downloaded TCC rasters/zips
#   HVI_TCC_URL_TEMPLATE   optional URL template with {year}, used if scraping cannot resolve links
#   HVI_TCC_PAGE_URL       optional source page override
#   HVI_TCC_OVERWRITE      true/false; default false
# ================================================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(janitor)
  library(readr)
  library(sf)
  library(stringr)
  library(tibble)
  library(tidyr)
})

source(file.path(if (dir.exists("code")) "code" else ".", "00_config.R"))

if (!requireNamespace("terra", quietly = TRUE)) {
  stop("The terra package is required to aggregate tree canopy rasters.")
}

options(timeout = max(getOption("timeout"), 7200))

base_dir <- HVI_PATHS$private_outputs$baseline_burden
ca_shp_path <- HVI_PATHS$raw$community_areas
tcc_dir <- file.path(base_dir, "tree_canopy")
raw_dir <- file.path(tcc_dir, "raw")
hvi_dir_create(base_dir)
hvi_dir_create(raw_dir)

tcc_page_url <- hvi_env(
  "HVI_TCC_PAGE_URL",
  "https://data.fs.usda.gov/geodata/rastergateway/treecanopycover/"
)
tcc_source_preference <- tolower(hvi_env("HVI_TCC_SOURCE", "nlcd"))
overwrite_downloads <- tolower(hvi_env("HVI_TCC_OVERWRITE", "false")) %in% c("true", "1", "yes")

parse_years <- function(x) {
  x <- str_squish(x)
  if (str_detect(x, "^[0-9]{4}\\s*:\\s*[0-9]{4}$")) {
    bounds <- as.integer(str_split(x, "\\s*:\\s*", simplify = TRUE))
    return(seq(bounds[1], bounds[2]))
  }
  as.integer(str_split(x, "\\s*,\\s*", simplify = FALSE)[[1]])
}

study_years <- parse_years(hvi_env("HVI_TCC_YEARS", "2019:2022"))
study_years <- study_years[!is.na(study_years)]
if (length(study_years) == 0) stop("No valid HVI_TCC_YEARS were supplied.")

clean_community <- function(x) {
  x %>%
    as.character() %>%
    str_to_upper() %>%
    str_replace_all("&", "AND") %>%
    str_replace_all("[[:punct:]]", " ") %>%
    str_squish() %>%
    str_replace("^O HARE$", "OHARE")
}

pick_first_existing <- function(x) {
  x <- x[file.exists(x)]
  if (length(x) == 0) return(NA_character_)
  x[1]
}

guess_year <- function(path) {
  hit <- str_extract(basename(path), "(?<![0-9])(19|20)[0-9]{2}(?![0-9])")
  suppressWarnings(as.integer(hit))
}

download_url <- function(url, dest) {
  if (file.exists(dest) && !overwrite_downloads && file.info(dest)$size > 0) {
    return("skipped_exists")
  }
  status <- tryCatch(
    {
      utils::download.file(url, dest, mode = "wb", quiet = FALSE, method = "libcurl")
      "downloaded"
    },
    error = function(e) {
      warning("Download failed for ", url, ": ", conditionMessage(e))
      "failed"
    }
  )
  if (!file.exists(dest) || file.info(dest)$size == 0) status <- "failed"
  status
}

is_valid_zip <- function(path) {
  if (!file.exists(path) || tolower(tools::file_ext(path)) != "zip") return(TRUE)
  ok <- tryCatch({
    utils::unzip(path, list = TRUE)
    TRUE
  }, error = function(e) FALSE)
  ok
}

scrape_tcc_links <- function(page_url) {
  page <- tryCatch(readLines(page_url, warn = FALSE), error = function(e) character())
  if (length(page) == 0) return(tibble(year = integer(), url = character(), label = character()))
  html <- paste(page, collapse = "\n")

  option_matches <- str_match_all(
    html,
    "<option[^>]+value=[\"']([^\"']+)[\"'][^>]*>([^<]+)</option>"
  )[[1]]
  if (nrow(option_matches) == 0) {
    link_matches <- str_match_all(
      html,
      "<a[^>]+href=[\"']([^\"']+)[\"'][^>]*>([^<]+)</a>"
    )[[1]]
    option_matches <- link_matches
  }
  if (nrow(option_matches) == 0) return(tibble(year = integer(), url = character(), label = character()))

  root <- str_match(page_url, "^(https?://[^/]+)")[, 2]

  links <- tibble(
    url = option_matches[, 2],
    label = str_squish(option_matches[, 3])
  ) %>%
    mutate(
      url = case_when(
        str_detect(url, "^https?://") ~ url,
        str_detect(url, "^/") & !is.na(root) ~ paste0(root, url),
        TRUE ~ paste0(dirname(page_url), "/", url)
      ),
      year = suppressWarnings(as.integer(str_extract(label, "(?<![0-9])(19|20)[0-9]{2}(?![0-9])"))),
      source_rank = case_when(
        tcc_source_preference == "science" & !str_detect(str_to_lower(label), "\\bnlcd\\b") ~ 1L,
        tcc_source_preference == "science" & str_detect(str_to_lower(label), "\\bnlcd\\b") ~ 2L,
        tcc_source_preference != "science" & str_detect(str_to_lower(label), "\\bnlcd\\b") ~ 1L,
        str_detect(str_to_lower(label), "\\btcc\\b") ~ 2L,
        TRUE ~ 9L
      )
    ) %>%
    filter(!is.na(year), str_detect(str_to_lower(label), "tcc|tree canopy")) %>%
    arrange(source_rank, label) %>%
    distinct(year, .keep_all = TRUE) %>%
    select(year, url, label)

  links
}

local_tcc_files <- function(local_dir, years) {
  if (is.null(local_dir) || !nzchar(local_dir) || !dir.exists(local_dir)) {
    return(tibble(year = integer(), path = character(), status = character(), url = character()))
  }
  files <- list.files(
    local_dir,
    pattern = "\\.(tif|tiff|img|vrt|zip)$",
    full.names = TRUE,
    recursive = TRUE,
    ignore.case = TRUE
  )
  tibble(path = files) %>%
    mutate(year = vapply(path, guess_year, integer(1))) %>%
    filter(year %in% years) %>%
    arrange(year, path) %>%
    distinct(year, .keep_all = TRUE) %>%
    mutate(status = "local_file", url = NA_character_)
}

resolve_tcc_files <- function(years) {
  local_dir <- hvi_env("HVI_TCC_LOCAL_DIR", "")
  local_files <- local_tcc_files(local_dir, years)
  missing_years <- setdiff(years, local_files$year)

  manifest <- local_files
  if (length(missing_years) == 0) return(manifest)

  url_template <- hvi_env("HVI_TCC_URL_TEMPLATE", "")
  if (nzchar(url_template)) {
    url_tbl <- tibble(
      year = missing_years,
      url = str_replace_all(url_template, fixed("{year}"), as.character(missing_years)),
      label = paste("Template TCC", missing_years)
    )
  } else {
    scraped <- scrape_tcc_links(tcc_page_url)
    url_tbl <- scraped %>%
      filter(year %in% missing_years)
  }

  if (nrow(url_tbl) == 0) {
    stop(
      "Could not resolve tree canopy download URLs for years: ",
      paste(missing_years, collapse = ", "),
      ". Set HVI_TCC_LOCAL_DIR to a directory of NLCD TCC rasters/zips, or set ",
      "HVI_TCC_URL_TEMPLATE to a direct URL template containing {year}."
    )
  }

  downloaded <- lapply(seq_len(nrow(url_tbl)), function(i) {
    yr <- url_tbl$year[i]
    url <- url_tbl$url[i]
    ext <- tools::file_ext(strsplit(url, "\\?")[[1]][1])
    if (!nzchar(ext)) ext <- "zip"
    dest <- file.path(raw_dir, paste0("tree_canopy_", yr, ".", ext))
    status <- download_url(url, dest)
    tibble(year = yr, path = dest, status = status, url = url)
  }) %>%
    bind_rows()

  bind_rows(manifest, downloaded) %>%
    arrange(year)
}

unpack_tcc_path <- function(path, year) {
  ext <- tolower(tools::file_ext(path))
  if (ext != "zip") return(path)
  if (!is_valid_zip(path)) {
    stop("Invalid or incomplete ZIP file: ", path)
  }

  exdir <- file.path(raw_dir, paste0("tree_canopy_", year))
  hvi_dir_create(exdir)
  utils::unzip(path, exdir = exdir)
  candidates <- list.files(
    exdir,
    pattern = "\\.(tif|tiff|img|vrt)$",
    full.names = TRUE,
    recursive = TRUE,
    ignore.case = TRUE
  )
  if (length(candidates) == 0) stop("No raster file found inside ", path)
  candidates <- candidates[
    str_detect(str_to_lower(basename(candidates)), "tcc|tree|canopy|nlcd")
  ]
  if (length(candidates) == 0) {
    candidates <- list.files(
      exdir,
      pattern = "\\.(tif|tiff|img|vrt)$",
      full.names = TRUE,
      recursive = TRUE,
      ignore.case = TRUE
    )
  }
  candidates[1]
}

read_community_areas <- function(path) {
  ca_sf <- st_read(path, quiet = TRUE) %>%
    clean_names()
  name_candidates <- c("community", "community_area", "community_name", "pri_neigh", "commarea", "area_numbe")
  ca_name_col <- intersect(name_candidates, names(ca_sf))[1]
  if (is.na(ca_name_col)) {
    stop("Could not identify the community area name column in the CA boundary file.")
  }
  ca_sf %>%
    transmute(
      community = clean_community(.data[[ca_name_col]]),
      geometry = geometry
    ) %>%
    st_make_valid() %>%
    group_by(community) %>%
    summarise(geometry = st_union(geometry), .groups = "drop")
}

aggregate_tcc <- function(raster_path, year, ca_sf) {
  r <- terra::rast(raster_path)
  if (terra::nlyr(r) > 1) {
    r <- r[[1]]
  }

  ca_proj <- st_transform(ca_sf, terra::crs(r))
  ca_vect <- terra::vect(ca_proj)
  r_chicago <- terra::crop(r, terra::ext(ca_vect), snap = "out")
  r_chicago <- terra::ifel(r_chicago < 0 | r_chicago > 100, NA, r_chicago)
  extracted <- terra::extract(r_chicago, ca_vect, fun = mean, na.rm = TRUE)

  value_col <- setdiff(names(extracted), "ID")[1]
  tibble(
    community = ca_proj$community,
    year = as.integer(year),
    tree_canopy_pct = suppressWarnings(as.numeric(extracted[[value_col]])),
    tree_canopy_fraction = tree_canopy_pct / 100,
    tree_canopy_source = if_else(tcc_source_preference == "science", "USFS Science TCC", "USFS/MRLC NLCD TCC"),
    tree_canopy_units = "percent"
  )
}

tcc_manifest <- resolve_tcc_files(study_years) %>%
  rowwise() %>%
  mutate(
    resolved_raster_path = if (file.exists(path)) unpack_tcc_path(path, year) else NA_character_
  ) %>%
  ungroup()

write_csv(tcc_manifest, file.path(base_dir, "tree_canopy_download_manifest.csv"))

missing_files <- tcc_manifest %>%
  filter(!file.exists(resolved_raster_path))
if (nrow(missing_files) > 0) {
  stop("Missing tree canopy rasters after download/local resolution.")
}

ca_sf <- read_community_areas(ca_shp_path)

tree_canopy_year_ca <- lapply(seq_len(nrow(tcc_manifest)), function(i) {
  message("Aggregating tree canopy for ", tcc_manifest$year[i], ": ", tcc_manifest$resolved_raster_path[i])
  aggregate_tcc(tcc_manifest$resolved_raster_path[i], tcc_manifest$year[i], ca_sf)
}) %>%
  bind_rows() %>%
  arrange(community, year)

tree_canopy_year_ca_qa <- tree_canopy_year_ca %>%
  group_by(year) %>%
  summarise(
    n_communities = n_distinct(community),
    n_missing_tree_canopy = sum(is.na(tree_canopy_pct)),
    min_tree_canopy_pct = min(tree_canopy_pct, na.rm = TRUE),
    median_tree_canopy_pct = median(tree_canopy_pct, na.rm = TRUE),
    max_tree_canopy_pct = max(tree_canopy_pct, na.rm = TRUE),
    .groups = "drop"
  )

write_csv(tree_canopy_year_ca, file.path(base_dir, "tree_canopy_year_ca.csv"))
write_csv(tree_canopy_year_ca_qa, file.path(base_dir, "tree_canopy_year_ca_qa.csv"))

assign("tree_canopy_year_ca", tree_canopy_year_ca, envir = .GlobalEnv)
assign("tree_canopy_year_ca_qa", tree_canopy_year_ca_qa, envir = .GlobalEnv)
assign("tree_canopy_download_manifest", tcc_manifest, envir = .GlobalEnv)

message("Tree canopy aggregation complete: ", file.path(base_dir, "tree_canopy_year_ca.csv"))
