# Pipeline

The repository now has a first `targets` scaffold in `_targets.R`. It is intentionally conservative: it wraps the existing scripts while the analysis is gradually migrated from global-session objects into explicit target return values.

## Current Stages

1. `standardized_panel`: creates private standardized records and aggregate community-day panel.
2. `climate_panel`: standardizes Daymet and joins climate to the community-day panel.
3. `tree_canopy_coverage`: downloads or reads annual USDA/MRLC tree canopy cover rasters and aggregates them to community-area-year coverage.
4. `baseline_vulnerability`: builds community-area-year structural vulnerability covariates, including tree canopy coverage when available.
5. `public_exports`: creates public dashboard/manuscript-safe outputs.
6. `scenario_exports`: creates precomputed dashboard slider outputs for temperature and NDVI scenarios, with tree canopy baseline values carried for frontend context.

Canonical production scripts live directly under `code/`. Prototype, exploratory, or superseded scripts are kept in `code/archive/` and should not be called by `_targets.R` or dashboard export automation.

## Run

```r
targets::tar_make()
```

Set paths with `.Renviron` copied from `.Renviron.example`.

## Next Migration

The model-fitting scripts still depend on in-memory objects. The next step is to make each script expose a function such as `build_hvi_model_matrix()` or `fit_endpoint_models()` and return explicit data frames/RDS paths to `targets`.
