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

Run:

```r
source("code/11_build_public_exports.R")
```
