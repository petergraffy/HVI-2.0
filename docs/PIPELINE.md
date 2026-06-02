# Pipeline

The repository now has a first `targets` scaffold in `_targets.R`. It is intentionally conservative: it wraps the existing scripts while the analysis is gradually migrated from global-session objects into explicit target return values.

## Current Stages

1. `standardized_panel`: creates private standardized records and aggregate community-day panel.
2. `climate_panel`: standardizes Daymet and joins climate to the community-day panel.
3. `tree_canopy_coverage`: downloads or reads annual USDA/MRLC tree canopy cover rasters and aggregates them to community-area-year coverage.
4. `baseline_vulnerability`: builds community-area-year structural vulnerability covariates, including tree canopy coverage when available. Missing vulnerability covariates are median-imputed before z-scoring so communities such as O'Hare stay in downstream scoring.
5. `model_matrix`: builds endpoint-specific heat-dose features. All-cause ED uses a 32 C default MRT override for operational excess calculations, with the original DLNM MRT retained in metadata for audit.
6. `risk_scoring`: scores structural, temperature-grid, and historical daily risk. Public 0-100 temporal scores are scaled from positive modeled heat-attributable excess; signed excess remains in aggregate outputs for QA.
7. `public_exports`: creates public dashboard/manuscript-safe outputs, including temporal risk maps, context-only baseline vulnerability tables, all-cause excess cost files, and condition-specific excess drilldown files.
8. `scenario_exports`: creates precomputed dashboard slider outputs for temperature, NDVI, and AC loss-of-cooling scenarios, with tree canopy baseline values carried for frontend context.

Canonical production scripts live directly under `code/`. Prototype, exploratory, or superseded scripts are kept in `code/archive/` and should not be called by `_targets.R` or dashboard export automation.

## Run

```r
targets::tar_make()
```

Set paths with `.Renviron` copied from `.Renviron.example`.

## Next Migration

The model-fitting scripts still depend on in-memory objects. The next step is to make each script expose a function such as `build_hvi_model_matrix()` or `fit_endpoint_models()` and return explicit data frames/RDS paths to `targets`.
