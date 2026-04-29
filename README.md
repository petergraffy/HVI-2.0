# Chicago Health-Informed Heat Vulnerability Index 2.0

**PI:** Peter Graffy  
**Institution:** University of Chicago, Center for Computational Medicine and Clinical AI  
**Project type:** Spatiotemporal environmental health modeling and public decision-support tool  
**Primary geography:** Chicago community areas  
**Primary temporal unit:** Community area-day, warm season

## Overview

The Chicago Health-Informed Heat Vulnerability Index 2.0 (HVI 2.0) is a next-generation heat vulnerability framework for estimating neighborhood-level health risk during extreme heat. The project is designed to support both a scientific manuscript and a public-facing dashboard/map for Chicago.

Most heat vulnerability indices are built from static social, demographic, housing, and environmental indicators. HVI 2.0 takes a different approach: it calibrates vulnerability to observed health outcomes. The index integrates emergency medical services calls, emergency department visits, and mortality records to estimate where heat is most likely to translate into measurable health burden.

The core scientific target is expected excess health burden under heat conditions, not raw outcome rates. This distinction matters because high baseline utilization or mortality alone does not necessarily identify heat-specific vulnerability.

## Scientific Aims

1. Estimate historical heat-health relationships across EMS, ED, and mortality endpoints.
2. Identify community areas where health burden increases disproportionately during heat.
3. Construct a multi-endpoint, health-informed vulnerability index.
4. Produce endpoint-specific and composite risk scores for mapping and public communication.
5. Evaluate how environmental and structural covariates, such as greenness and air conditioning access, modify predicted heat-related health burden.
6. Provide dashboard-ready scenario outputs so users can explore how temperature and modifiable covariates may affect endpoint-specific risk.

## Data Streams

### Health Outcomes

| Source | Approximate Years | Role |
| --- | --- | --- |
| CAPriCORN emergency department records | 2011-2022 | Acute care utilization and cause-specific ED outcomes |
| Chicago mortality records | 1993-2022 | All-cause and cause-specific mortality outcomes |
| Chicago EMS records | 2018-2023 | Pre-hospital acute health events and syndromic response |
| CLIF ICU data | Optional extension | Severe illness validation or added endpoint layer |

All record-level health data are treated as private data. Public-facing exports are aggregated to community-area level and must not contain direct identifiers, addresses, narratives, coordinates, encounter IDs, or other record-level fields.

### Environmental Exposure

The primary historical heat exposure is daily maximum temperature from Daymet, harmonized to Chicago community areas. The current model matrix uses Celsius-scale temperature features internally. Public dashboard controls can display Fahrenheit, but backend scoring converts Fahrenheit to Celsius before applying model terms.

Other exposure or adaptive-capacity inputs include humidity where available, NDVI, PM2.5, NO2, and air-conditioning probability.

### Structural Vulnerability Covariates

The model matrix includes community-area-year covariates representing:

- age structure
- race and ethnicity
- income, education, employment, and poverty
- housing burden and household composition
- vehicle access, insurance, internet access, and social vulnerability indicators
- NDVI and other built or natural environment measures
- air conditioning probability
- air pollution burden

Candidate covariates are screened and selected using a combination of missingness checks, correlation pruning, random forest importance, and interaction-aware generalized additive models.

## Spatial and Temporal Design

The primary analytic unit is the Chicago community area by day. The main HVI modeling window currently emphasizes the 2019-2022 overlap period, when EMS, ED, mortality, climate, and baseline covariate data are jointly available. Warm-season analyses are prioritized because the target estimand is heat-related health burden.

The project is structured to support future sensitivity analyses at alternative spatial resolutions or under future climate scenarios.

## Endpoint Families

HVI 2.0 models three broad endpoint families:

- **Mortality:** all-cause, cardiovascular, respiratory, renal, neurologic, mental health, gastrointestinal, injury, and other heat-sensitive categories.
- **Emergency department visits:** all-cause and cause-specific acute care outcomes including cardiovascular, respiratory, renal, dehydration, heat illness, injury, syncope, neurologic, mental health, and gastrointestinal visits.
- **EMS calls:** all-cause and syndromic EMS outcomes including cardiovascular, respiratory, syncope, neurologic, mental health, injury, gastrointestinal, bleeding, and heat-related calls.

The endpoint dictionary is maintained in [docs/ENDPOINT_DICTIONARY.md](docs/ENDPOINT_DICTIONARY.md).

## Modeling Framework

### 1. Outcome Harmonization

Raw EMS, ED, and mortality records are standardized to a common event date and community-area geography. Cause-specific indicators are then derived from diagnosis codes, underlying cause of death codes, or EMS symptom/impression fields. Daily community-area panels are generated for all-cause and endpoint-specific counts.

Record-level standardized files are private working artifacts only. Aggregate panels and public exports are the safe downstream products.

### 2. Minimum-Risk Temperature and Heat Dose

Endpoint-specific temperature-response relationships are estimated to identify minimum-risk temperatures (MRTs) and relevant lag windows. Heat dose is then defined as cumulative excess temperature above the endpoint-specific MRT over the endpoint's lag window.

This design allows different health endpoints to have different heat thresholds and lag structures. For example, EMS and ED endpoints may respond over shorter windows, while mortality endpoints may require longer distributed lag structures.

### 3. Endpoint-Specific Health Models

Endpoint models are fit using generalized additive modeling machinery with count-family likelihoods. The production HVI model uses endpoint-specific outcome counts as the response, heat dose as the primary exposure term, selected vulnerability covariates, and heat-dose-by-covariate interaction terms.

Core adjustment variables include seasonal timing, day of week, year, and population offset where available. Model outputs include fitted endpoint models, coefficients, interaction terms, cross-validation predictions, and endpoint-level performance summaries.

### 4. Variable Selection and Effect Modification

The vulnerability modeling pipeline screens candidate structural and environmental covariates in stages:

1. Drop variables with excessive missingness or near-zero variance.
2. Use random forest screening to identify variables that improve endpoint prediction.
3. Fit interaction-aware GAMs to estimate whether covariates modify heat-dose response.
4. Rank variables by predictive signal, interaction strength, sign consistency, and redundancy.
5. Apply correlation pruning to keep a smaller, interpretable final covariate set.

This produces endpoint-specific and cross-endpoint covariate sets for structural vulnerability scoring.

### 5. Structural Vulnerability Scoring

For each community area and year, endpoint-specific vulnerability scores are computed from selected covariates and their heat interaction coefficients. These endpoint scores are combined into family-level and overall structural HVI scores.

The structural score answers:

> Which communities have baseline characteristics associated with greater heat-sensitive health risk?

### 6. Temperature and Operational Risk Scoring

The pipeline also scores risk across temperature scenarios and observed daily heat conditions. For each community area, endpoint, and temperature scenario, fitted models estimate predicted counts, reference counts, excess events, relative risk, endpoint risk, family risk, and overall risk.

The operational score answers:

> Given a specific heat scenario or observed day, which communities are expected to experience the largest heat-related health burden?

### 7. Scenario Slider Backend

The dashboard backend includes precomputed scenario exports for user controls:

- temperature in degrees Fahrenheit
- NDVI change relative to baseline
- humidity metadata, currently inactive unless endpoint models are refit with humidity terms

Scenario exports are generated by [code/12_build_scenario_exports.R](code/12_build_scenario_exports.R). The script converts Fahrenheit to Celsius for model scoring, caps heat dose at endpoint-specific observed support, and marks model-boundary rows using `heat_dose_capped` and `prediction_capped_for_display`.

For public visualization, the app should use the 0-100 score fields rather than interpreting scenario predicted counts as causal estimates. NDVI slider effects are association-based scenarios, not causal intervention estimates.

## Validation Strategy

The validation framework includes:

- spatial cross-validation by held-out community areas
- temporal validation by held-out years
- endpoint-specific predictive performance summaries
- comparison of high-HVI areas against observed heat-event burden
- evaluation of targeting efficiency, such as the share of health burden captured in the highest-risk communities
- comparison against traditional socioeconomic-only vulnerability indices

Performance artifacts are exported for manuscript tables and dashboard metadata.

## Public Dashboard Exports

Dashboard-safe files are generated under `public_exports/`. The public application should read from:

- `public_exports/dashboard/`
- `public_exports/dashboard/scenarios/`
- `public_exports/manifest.json`

The dashboard should not read from `data/`, `results/`, `code/results/`, `code/09_model_outputs/`, or private Box folders.

Public export details are documented in [docs/PUBLIC_EXPORT_CONTRACT.md](docs/PUBLIC_EXPORT_CONTRACT.md).

## Reproducibility and Configuration

Paths are centralized in [code/00_config.R](code/00_config.R). Copy `.Renviron.example` to `.Renviron` and adjust local paths as needed:

```r
HVI_BOX_DIR="C:/Users/Peter Graffy/Box/HVI2.0"
HVI_PRIVATE_DIR="C:/Users/Peter Graffy/Box/HVI2.0"
HVI_PUBLIC_EXPORT_DIR="public_exports"
HVI_SUPPRESS_SMALL_CELLS=11
```

Build public dashboard and manuscript-safe exports with:

```r
source("code/11_build_public_exports.R")
```

Build scenario slider exports with:

```r
source("code/12_build_scenario_exports.R")
```

A conservative `targets` scaffold is provided in [_targets.R](_targets.R). It wraps the current script sequence while the project is migrated away from global R session state and toward explicit target return values.

## Data Governance

This repository should contain code, documentation, manuscript figures/tables when appropriate, and public aggregate exports only. Raw health records and record-level standardized files must remain outside Git.

The `.gitignore` excludes private data, generated model outputs, large artifacts, and record-level health files. Additional governance notes are in [docs/DATA_GOVERNANCE.md](docs/DATA_GOVERNANCE.md).

## Key Documentation

- [docs/README.md](docs/README.md): documentation map for analysts, collaborators, and dashboard developers
- [docs/DATA_GOVERNANCE.md](docs/DATA_GOVERNANCE.md): private data handling, public-release checklist, and Git safety rules
- [docs/PUBLIC_EXPORT_CONTRACT.md](docs/PUBLIC_EXPORT_CONTRACT.md): dashboard file contract and scenario slider outputs
- [docs/PIPELINE.md](docs/PIPELINE.md): current pipeline stages and migration plan
- [docs/ENDPOINT_DICTIONARY.md](docs/ENDPOINT_DICTIONARY.md): endpoint names and families

## Scientific Contribution

HVI 2.0 advances heat vulnerability mapping by linking neighborhood risk directly to observed health outcomes across pre-hospital care, emergency care, and mortality. It provides both a manuscript-ready epidemiologic framework and a dashboard-ready public health tool for comparing structural vulnerability, temperature-driven risk, endpoint composition, and modifiable scenario inputs across Chicago communities.
