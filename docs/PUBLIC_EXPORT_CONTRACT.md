# Public Export Contract

The public dashboard should read from `public_exports/`, not from `data/`, `results/`, `code/results/`, or Box raw-data folders.

## Directory Layout

- `public_exports/dashboard/`: files intended for the Replit or production web app.
- `public_exports/aggregate/`: model-level aggregate outputs useful for QA, API layers, and documentation.
- `public_exports/manuscript/`: publication tables and supplements that are safe to share.
- `public_exports/manifest.json`: generated file inventory, generation time, and small-cell threshold.

## Privacy Rules

- No row-level health records.
- No names, DOB, addresses, SSNs, encounter IDs, study IDs, EMS narratives, GPS points, or ZIP-level identifiers.
- Community-area aggregates are the minimum public geography.
- Positive observed counts below `HVI_SUPPRESS_SMALL_CELLS` are suppressed to blank values in public CSVs.
- Model scores, relative risks, percentiles, and non-count metadata can remain visible if they cannot identify individuals.

## Dashboard Files

The current dashboard-facing exports are copied from the private `09_model_outputs/frontend_exports_v2/` directory when available, falling back to `frontend_exports/`.

Expected core files:

- `frontend_structural_summary.csv`
- `frontend_structural_endpoint_long.csv`
- `frontend_temperature_query_endpoint_long.csv`
- `frontend_temperature_query_overall.csv`
- `frontend_temperature_query_wide.csv`
- `frontend_daily_map.csv`
- `frontend_daily_endpoint_long.csv`
- `frontend_daily_wide.csv`

## Scenario Slider Files

The slider backend is exported under `public_exports/dashboard/scenarios/`.

- `scenario_grid_endpoint.csv`: endpoint-specific predictions for combinations of temperature and NDVI change.
- `scenario_grid_overall.csv`: overall HVI scenario scores and dominant endpoints.
- `scenario_variable_metadata.csv`: slider ranges, units, and whether each variable is active in the current fitted models.
- `scenario_baseline_values.csv`: baseline community-area values used to anchor delta sliders.

Current controls:

- `temperature_f`: active; scored through endpoint-specific MRT and lag windows.
- `ndvi_delta`: active when endpoint models include `z_ndvi`; interpreted as an association-based scenario, not a causal effect.
- `humidity_scenario`: included in metadata; only changes predictions after models are fit with humidity terms.

Temperature is displayed in Fahrenheit but converted to Celsius for model scoring because the fitted MRT/heat-dose artifacts are Celsius-based. Endpoint rows include `heat_dose_capped` and `prediction_capped_for_display`; the app should use the 0-100 score fields for visual comparison and show a model-boundary note when either cap flag is true.

Build scenario exports with:

```r
source("code/12_build_scenario_exports.R")
```

Run:

```r
source("code/11_build_public_exports.R")
```
