# ================================================================================================
# HVI 2.0 | Build public dashboard and manuscript export contract
# Copies only aggregate, dashboard-safe artifacts from the private analysis area into public_exports/.
# ================================================================================================

source(file.path(if (dir.exists("code")) "code" else ".", "00_config.R"))

dashboard_private_dirs <- c(
  file.path(HVI_PATHS$private_outputs$model_outputs, "frontend_exports_v2"),
  file.path(HVI_PATHS$private_outputs$model_outputs, "frontend_exports"),
  file.path(HVI_PATHS$private, "results", "frontend_exports"),
  file.path(HVI_PATHS$repo, "code", "09_model_outputs", "frontend_exports_v2"),
  file.path(HVI_PATHS$repo, "code", "09_model_outputs", "frontend_exports"),
  file.path(HVI_PATHS$repo, "results", "frontend_exports")
)
dashboard_private_dirs <- dashboard_private_dirs[dir.exists(dashboard_private_dirs)]

if (length(dashboard_private_dirs) == 0) {
  stop("No private frontend export directory found. Run 09e_build_frontend_exports.R first.")
}

dashboard_src <- dashboard_private_dirs[1]
dashboard_dst <- HVI_PATHS$public$dashboard
aggregate_dst <- HVI_PATHS$public$aggregate
manuscript_dst <- HVI_PATHS$public$manuscript

hvi_dir_create(dashboard_dst)
hvi_dir_create(aggregate_dst)
hvi_dir_create(manuscript_dst)

copy_public_csv <- function(src, dst) {
  dat <- utils::read.csv(src, check.names = FALSE)

  count_like <- grep(
    "(^observed_count$|^expected_events$|count$|events$|excess_events$|baseline_events$)",
    names(dat),
    ignore.case = TRUE,
    value = TRUE
  )
  count_like <- setdiff(count_like, c("relative_risk", "attributable_fraction"))
  for (col in count_like) {
    if (is.numeric(dat[[col]])) dat[[col]] <- hvi_public_count(dat[[col]])
  }

  hvi_write_public_csv(dat, dst)
}

dashboard_csv <- list.files(dashboard_src, pattern = "\\.csv$", full.names = TRUE)
for (src in dashboard_csv) {
  copy_public_csv(src, file.path(dashboard_dst, basename(src)))
}

dashboard_geojson <- list.files(dashboard_src, pattern = "\\.geojson$", full.names = TRUE)
for (src in dashboard_geojson) {
  file.copy(src, file.path(dashboard_dst, basename(src)), overwrite = TRUE)
}

scenario_src_dirs <- c(
  file.path(HVI_PATHS$private_outputs$model_outputs, "scenario_exports"),
  file.path(HVI_PATHS$repo, "public_exports", "dashboard", "scenarios")
)
scenario_src_dirs <- scenario_src_dirs[dir.exists(scenario_src_dirs)]
if (length(scenario_src_dirs) > 0) {
  scenario_dst <- file.path(dashboard_dst, "scenarios")
  hvi_dir_create(scenario_dst)
  scenario_csv <- list.files(scenario_src_dirs[1], pattern = "\\.csv$", full.names = TRUE)
  for (src in scenario_csv) {
    copy_public_csv(src, file.path(scenario_dst, basename(src)))
  }
}

model_public <- list.files(
  HVI_PATHS$private_outputs$model_outputs,
  pattern = "^(09a_endpoint_model_performance|09b_endpoint_weights|09c_temperature_grid_overall_risk|09c_temperature_grid_endpoint_risk|09d_ca_day_overall_operational_hvi)\\.csv$",
  full.names = TRUE
)
if (length(model_public) == 0) {
  model_public <- list.files(
    file.path(HVI_PATHS$repo, "code", "09_model_outputs"),
    pattern = "^(09a_endpoint_model_performance|09b_endpoint_weights|09c_temperature_grid_overall_risk|09c_temperature_grid_endpoint_risk|09d_ca_day_overall_operational_hvi)\\.csv$",
    full.names = TRUE
  )
}
for (src in model_public) {
  copy_public_csv(src, file.path(aggregate_dst, basename(src)))
}

publication_src <- HVI_PATHS$private_outputs$publication_outputs
if (!dir.exists(publication_src)) {
  publication_src <- file.path(HVI_PATHS$repo, "code", "outputs", "publication")
}
if (dir.exists(publication_src)) {
  pub_csv <- list.files(publication_src, pattern = "\\.csv$", recursive = TRUE, full.names = TRUE)
  for (src in pub_csv) {
    rel <- substring(src, nchar(publication_src) + 2)
    copy_public_csv(src, file.path(manuscript_dst, rel))
  }
}

json_escape <- function(x) {
  x <- gsub("\\\\", "\\\\\\\\", x)
  x <- gsub('"', '\\"', x)
  paste0('"', x, '"')
}

json_array <- function(x) {
  paste0("[", paste(vapply(x, json_escape, character(1)), collapse = ", "), "]")
}

manifest_lines <- c(
  "{",
  paste0('  "generated_at": ', json_escape(format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")), ","),
  paste0('  "source_dashboard_dir": ', json_escape(dashboard_src), ","),
  paste0('  "small_cell_threshold": ', HVI_SMALL_CELL_THRESHOLD, ","),
  '  "public_dirs": {',
  paste0('    "dashboard": ', json_escape(dashboard_dst), ","),
  paste0('    "aggregate": ', json_escape(aggregate_dst), ","),
  paste0('    "manuscript": ', json_escape(manuscript_dst)),
  "  },",
  '  "files": {',
  paste0('    "dashboard": ', json_array(list.files(dashboard_dst, recursive = TRUE, full.names = FALSE)), ","),
  paste0('    "aggregate": ', json_array(basename(list.files(aggregate_dst, full.names = FALSE))), ","),
  paste0('    "manuscript": ', json_array(list.files(manuscript_dst, recursive = TRUE, full.names = FALSE))),
  "  }",
  "}"
)

writeLines(manifest_lines, con = file.path(HVI_PUBLIC_EXPORT_DIR, "manifest.json"))

message("Public exports written to: ", HVI_PUBLIC_EXPORT_DIR)
