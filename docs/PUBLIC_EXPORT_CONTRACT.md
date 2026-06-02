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
- `frontend_baseline_risk_factors.csv`
- `frontend_tree_canopy_year_ca.csv`
- `frontend_temperature_query_endpoint_long.csv`
- `frontend_temperature_query_overall.csv`
- `frontend_temperature_query_wide.csv`
- `frontend_daily_map.csv`
- `frontend_daily_endpoint_long.csv`
- `frontend_daily_wide.csv`
- `frontend_historical_daily_summary.csv`
- `frontend_historical_daily_endpoint_ranked.csv`
- `frontend_excess_endpoint_daily.csv`
- `frontend_excess_endpoint_annual.csv`
- `frontend_excess_endpoint_citywide.csv`
- `frontend_excess_allcause_daily.csv`
- `frontend_excess_allcause_annual.csv`
- `frontend_excess_allcause_citywide.csv`
- `frontend_excess_cost_assumptions.csv`

## Historical Risk Files

The historical dashboard view should use the frontend exports for 2019-2022 warm-season community-area days. The primary map score is temporal heat-health risk, scaled from positive modeled heat-attributable excess; signed excess fields remain available for QA but should not be used directly as display scores.

- `frontend_daily_map.csv`: community-area-day map table. It includes `temperature_c`, `temperature_f`, `humidity`, `tree_canopy_pct`, and `tree_canopy_fraction` where available, daily 0-100 temporal risk score fields, and dominant endpoint labels. Count-like model fields are suppressed by the public export builder when they are below the small-cell threshold.
- `frontend_historical_daily_summary.csv`: public-safe community-area-day summary for tooltips, map cards, and time-series views. It includes temperature, tree canopy coverage, daily risk tier, and dominant health endpoint. It intentionally omits observed, predicted, reference, and excess count values.
- `frontend_historical_daily_endpoint_ranked.csv`: public-safe endpoint ranking by community-area-day. It includes endpoint labels, source family, rank, endpoint 0-100 risk score, qualitative impact tier, relative risk, and proportional contribution (`impact_share_pct`). It intentionally omits observed, predicted, reference, and excess count values.
- `frontend_baseline_risk_factors.csv`: context-only community-area-year baseline risk explanation table. It reports tree canopy coverage, the baseline vulnerability score/tier, dominant baseline endpoint, and top modeled structural or environmental contributors using human-readable labels. It is retained for explanatory panels and QA, not as the primary map score.
- `frontend_tree_canopy_year_ca.csv`: standalone community-area-year tree canopy coverage export with `tree_canopy_pct` and `tree_canopy_fraction` for map layers, filters, and frontend QA.
- `frontend_excess_endpoint_daily.csv`: modeled positive heat-attributable excess by community-area-day and condition-specific endpoint. This is the main file for disease/condition drilldowns. It includes endpoint family, domain, temperature, tree canopy coverage, population when available, modeled positive excess events, excess rate per 100,000 residents, endpoint risk score, and relative risk.
- `frontend_excess_endpoint_annual.csv`: community-area-year rollup of the condition-specific excess file.
- `frontend_excess_endpoint_citywide.csv`: citywide year-level and all-year rollup of condition-specific excess events and rates.
- `frontend_excess_allcause_daily.csv`: modeled heat-attributable excess for all-cause ED visits (`ed_visits`), EMS calls (`ems_calls`), and deaths (`deaths`) by community-area-day. It includes temperature, tree canopy coverage, population when available, modeled positive excess events, excess rate per 100,000 residents, unit cost, and estimated cost. Public export suppression applies to count-like fields.
- `frontend_excess_allcause_annual.csv`: community-area-year rollup of the all-cause excess file.
- `frontend_excess_allcause_citywide.csv`: citywide year-level and all-year rollup of all-cause excess events, rates, and estimated costs.
- `frontend_excess_cost_assumptions.csv`: unit-cost assumptions used to estimate health/medical/economic burden. Defaults can be overridden with `HVI_UNIT_COST_ED_VISIT_USD`, `HVI_UNIT_COST_EMS_CALL_USD`, `HVI_UNIT_COST_DEATH_USD`, `HVI_UNIT_COST_YEAR`, or a CSV path in `HVI_EXCESS_COST_ASSUMPTIONS`.

These files support statements such as: the temperature observed in a specific community area on a given historical day; which modeled health endpoint contributed most to that day's heat-health risk; which condition categories account for modeled excess burden; and which baseline community characteristics contribute to structural vulnerability. They should not be used to display raw health event counts.

Excess events are model-derived differences between heat-day predictions and reference no-heat predictions. Negative excess is clamped to zero for the cost files so estimated cost represents positive heat-attributable burden only. Mortality cost defaults are valuation assumptions, not direct medical bills; public displays should label them accordingly.

## Scenario Slider Files

The slider backend is exported under `public_exports/dashboard/scenarios/`.

- `scenario_grid_endpoint.csv`: endpoint-specific predictions for combinations of temperature, NDVI change, and AC loss-of-cooling scenarios, with humidity held at community baseline values unless explicitly enabled for research QA.
- `scenario_grid_overall.csv`: overall HVI scenario scores and dominant endpoints.
- `scenario_variable_metadata.csv`: slider ranges, units, and whether each variable is active in the current fitted models.
- `scenario_baseline_values.csv`: baseline community-area values used to anchor delta sliders, including tree canopy coverage when available.

Current controls:

- `temperature_f`: active; scored through endpoint-specific MRT and lag windows.
- `humidity_scenario`: inactive by default; humidity is retained as a model adjustment term but held at each community area's baseline value in public scenarios.
- `ndvi_delta`: active when endpoint models include `z_ndvi`; interpreted as an association-based scenario, not a causal effect.
- `ac_delta`: active when endpoint models include `z_ac_prob`; interpreted as an association-based loss-of-cooling scenario. Defaults include baseline and negative values only to simulate reduced effective AC access, including outage-like conditions.

Temperature is displayed in Fahrenheit but converted to Celsius for model scoring because the fitted MRT/heat-dose artifacts are Celsius-based. Humidity is not exposed as a public intervention slider because the fitted observational association is an adjustment term, not a causal response function. AC scenarios are applied as reductions from each community area's baseline `ac_baseline` and clamped to the 0-1 prevalence range before scoring. Endpoint rows include `heat_dose_capped` and `prediction_capped_for_display`; the app should use the 0-100 score fields for visual comparison and show a model-boundary note when either cap flag is true.

Build scenario exports with:

```r
source("code/12_build_scenario_exports.R")
```

Run:

```r
source("code/11_build_public_exports.R")
```
