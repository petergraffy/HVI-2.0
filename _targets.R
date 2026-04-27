# ================================================================================================
# HVI 2.0 | Reproducible pipeline scaffold
# This wraps the current script sequence while the analysis is migrated away from global state.
# ================================================================================================

library(targets)

tar_option_set(
  packages = c(
    "data.table", "dplyr", "readr", "tidyr", "lubridate", "stringr", "janitor",
    "sf", "mgcv", "slider", "jsonlite"
  ),
  error = "stop"
)

source("code/00_config.R")

tar_script_step <- function(path) {
  force(path)
  function() {
    source(path, chdir = TRUE)
    path
  }
}

list(
  tar_target(config_loaded, HVI_PATHS, cue = tar_cue(mode = "always")),

  tar_target(
    standardized_panel,
    tar_script_step("code/01_standardize_time_space.R")(),
    cue = tar_cue(mode = "thorough")
  ),

  tar_target(
    climate_panel,
    tar_script_step("code/02_climate_qc_daymet.R")(),
    cue = tar_cue(mode = "thorough")
  ),

  tar_target(
    baseline_vulnerability,
    tar_script_step("code/06_baseline_vulnerability.R")(),
    cue = tar_cue(mode = "thorough")
  ),

  tar_target(
    public_exports,
    tar_script_step("code/11_build_public_exports.R")(),
    cue = tar_cue(mode = "always")
  )
)
