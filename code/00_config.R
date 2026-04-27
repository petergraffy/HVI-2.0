# ================================================================================================
# HVI 2.0 | Shared configuration
# Centralizes private, intermediate, manuscript, and dashboard paths.
# ================================================================================================

hvi_find_repo_root <- function(start = getwd()) {
  cur <- normalizePath(start, winslash = "/", mustWork = FALSE)
  for (i in seq_len(8)) {
    if (file.exists(file.path(cur, "README.md")) && dir.exists(file.path(cur, "code"))) {
      return(cur)
    }
    parent <- dirname(cur)
    if (identical(parent, cur)) break
    cur <- parent
  }
  normalizePath(start, winslash = "/", mustWork = FALSE)
}

hvi_env <- function(name, default = NULL) {
  val <- Sys.getenv(name, unset = NA_character_)
  if (is.na(val) || !nzchar(val)) default else val
}

HVI_REPO_DIR <- hvi_find_repo_root()
HVI_BOX_DIR <- normalizePath(
  hvi_env("HVI_BOX_DIR", "C:/Users/Peter Graffy/Box/HVI2.0"),
  winslash = "/",
  mustWork = FALSE
)
HVI_PRIVATE_DIR <- normalizePath(
  hvi_env("HVI_PRIVATE_DIR", HVI_BOX_DIR),
  winslash = "/",
  mustWork = FALSE
)
HVI_PUBLIC_EXPORT_DIR <- normalizePath(
  hvi_env("HVI_PUBLIC_EXPORT_DIR", file.path(HVI_REPO_DIR, "public_exports")),
  winslash = "/",
  mustWork = FALSE
)
HVI_SMALL_CELL_THRESHOLD <- as.integer(hvi_env("HVI_SUPPRESS_SMALL_CELLS", "11"))

HVI_PATHS <- list(
  repo = HVI_REPO_DIR,
  box = HVI_BOX_DIR,
  private = HVI_PRIVATE_DIR,
  public_exports = HVI_PUBLIC_EXPORT_DIR,
  raw = list(
    mortality = file.path(HVI_PRIVATE_DIR, "Mortality", "all_deaths.csv"),
    ed = file.path(HVI_PRIVATE_DIR, "ED", "ed_outcomes_complete.csv"),
    ems = file.path(HVI_PRIVATE_DIR, "EMS", "new_files_deduplicated_clean.csv"),
    community_areas = file.path(HVI_PRIVATE_DIR, "comm_areas.geojson"),
    climate_hist = file.path(HVI_PRIVATE_DIR, "Climate", "CA_temps_rh_90-23.csv"),
    climate_2024_nc = file.path(HVI_PRIVATE_DIR, "Climate", "DAYMET.004_1km_aid0001.nc"),
    daymet_2024 = file.path(HVI_PRIVATE_DIR, "Climate", "ca_daymet_2024_aggregated.csv"),
    climate_full = file.path(HVI_PRIVATE_DIR, "Climate", "CA_temps_rh_90-24.csv")
  ),
  private_outputs = list(
    derived = file.path(HVI_PRIVATE_DIR, "derived"),
    model_outputs = file.path(HVI_PRIVATE_DIR, "09_model_outputs"),
    variable_selection = file.path(HVI_PRIVATE_DIR, "variable_selection"),
    publication_outputs = file.path(HVI_PRIVATE_DIR, "publication_outputs"),
    baseline_burden = file.path(HVI_PRIVATE_DIR, "Baseline Burden")
  ),
  public = list(
    aggregate = file.path(HVI_PUBLIC_EXPORT_DIR, "aggregate"),
    dashboard = file.path(HVI_PUBLIC_EXPORT_DIR, "dashboard"),
    manuscript = file.path(HVI_PUBLIC_EXPORT_DIR, "manuscript")
  )
)

hvi_path <- function(...) file.path(...)

hvi_source <- function(path) {
  candidates <- c(
    file.path(HVI_REPO_DIR, path),
    file.path(HVI_REPO_DIR, "code", path),
    path
  )
  hit <- candidates[file.exists(candidates)][1]
  if (is.na(hit)) stop("Could not locate source file: ", path)
  source(hit, chdir = TRUE)
}

hvi_dir_create <- function(path) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  invisible(path)
}

hvi_public_count <- function(x, threshold = HVI_SMALL_CELL_THRESHOLD) {
  ifelse(is.na(x), NA, ifelse(x > 0 & x < threshold, NA, x))
}

hvi_assert_no_phi_columns <- function(df, context = "public export") {
  phi_patterns <- c(
    "name", "fname", "lname", "dob", "ssn", "address", "street", "zip",
    "latitude", "longitude", "gps", "narrative", "phone", "study_id", "enc_id"
  )
  hits <- names(df)[grepl(paste(phi_patterns, collapse = "|"), names(df), ignore.case = TRUE)]
  if (length(hits) > 0) {
    stop(context, " contains potential PHI columns: ", paste(hits, collapse = ", "))
  }
  invisible(TRUE)
}

hvi_write_public_csv <- function(df, path, threshold = HVI_SMALL_CELL_THRESHOLD) {
  hvi_assert_no_phi_columns(df, path)
  hvi_dir_create(dirname(path))
  utils::write.csv(df, path, row.names = FALSE, na = "")
  invisible(path)
}
