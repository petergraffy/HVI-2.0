# ================================================================================================
# HVI 2.0 | Reproducible pipeline scaffold
# This wraps the current script sequence while the analysis is migrated away from global state.
# ================================================================================================

library(targets)

tar_option_set(
  packages = c(
    "data.table", "dplyr", "readr", "tidyr", "lubridate", "stringr", "janitor",
    "sf", "terra", "mgcv", "slider", "jsonlite"
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
    tree_canopy_coverage,
    tar_script_step("code/02b_tree_canopy_coverage.R")(),
    cue = tar_cue(mode = "thorough")
  ),

  tar_target(
    baseline_vulnerability,
    {
      tree_canopy_coverage
      tar_script_step("code/06_baseline_vulnerability.R")()
    },
    cue = tar_cue(mode = "thorough")
  ),

  tar_target(
    scenario_exports,
    tar_script_step("code/12_build_scenario_exports.R")(),
    cue = tar_cue(mode = "always")
  ),

  tar_target(
    long_term_forecast_exports,
    tar_script_step("code/13_build_long_term_forecast_exports.R")(),
    cue = tar_cue(mode = "always")
  ),

  tar_target(
    public_exports,
    {
      scenario_exports
      long_term_forecast_exports
      tar_script_step("code/11_build_public_exports.R")()
    },
    cue = tar_cue(mode = "always")
  )
)
