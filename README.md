
# Chicago Health-Informed Heat Vulnerability Index (HVI) 2.0

**PI:** Peter Graffy  
**Project Type:** Spatiotemporal Environmental Health Modeling  
**Scope:** Chicago, IL (Community Area–level primary geography)  

---

## Overview

The **Health-Informed Heat Vulnerability Index (HVI) 2.0** is a next-generation, multi-endpoint framework designed to quantify neighborhood-level health risk associated with extreme heat in Chicago.

Unlike traditional vulnerability indices based solely on socioeconomic characteristics, HVI 2.0 is explicitly **calibrated to observed health outcomes**, integrating:

- Emergency Department (ED) visits  
- Mortality  
- Emergency Medical Services (EMS) calls  
- (Optional) ICU admissions via CLIF  

The framework estimates **heat-attributable health burden**, models **spatial heterogeneity in vulnerability**, and projects **future risk under climate change scenarios** using downscaled temperature data.

---

## Core Objectives

1. **Quantify historical heat-health relationships** across multiple endpoints  
2. **Estimate neighborhood-specific vulnerability** to heat exposure  
3. **Construct a multi-endpoint composite vulnerability index**  
4. **Project future heat-related health burden** under climate scenarios  
5. **Validate index performance against observed heat events**

---

## Data Sources

### Health Outcomes

| Source | Years | Description |
|------|------|------------|
| CAPriCORN | 2011–2022 | ED visits across Chicago health systems |
| CDPH Mortality | 1993–2022 | Death records (all-cause and cause-specific) |
| EMS | ~2018–2023 | Emergency medical service calls |
| CLIF (optional) | TBD | ICU admissions and severe illness endpoints |

---

### Environmental Exposure

- **Daymet** (historical daily temperature)
- **Downscaled climate projections** (QDM / Earthmover dataset)

---

### Covariates

#### Structural Vulnerability
- Age distribution  
- Socioeconomic status  
- Race/ethnicity  
- Housing characteristics  
- Insurance coverage  

#### Baseline Health Burden
- Chronic disease prevalence  
- Prior healthcare utilization  
- Frailty proxies  

#### Built Environment / Adaptive Capacity
- Tree canopy  
- Impervious surface  
- Green space  
- Housing age  
- Cooling access proxies  

---

## Spatial and Temporal Resolution

- **Primary unit:** Community Area × Day  
- **Secondary analyses:** Census tract (sensitivity)  
- **Temporal scope:** Warm season (May–September primary)  

---

## Exposure Definition

Multiple complementary heat metrics will be used:

- Daily max temperature (Tmax)  
- Daily mean temperature  
- Heat index / apparent temperature  
- Heatwave indicators (duration, intensity)  
- Early vs late season heat  
- Consecutive hot days  
- Lagged exposure (0–10 days depending on endpoint)  

---

## Outcomes

### Mortality
- All-cause  
- Cardiovascular  
- Respiratory  
- Heat-related (if coded)

### Emergency Department
- All-cause visits  
- Heat illness  
- Dehydration  
- Renal injury  
- Cardiovascular  
- Respiratory  

### EMS Calls
- All calls  
- Heat-related dispatch  
- Cardiopulmonary  
- Syncope / collapse  
- Altered mental status  

### ICU (Optional)
- ICU admission  
- Mechanical ventilation  
- Vasopressor use  
- Acute respiratory failure  

---

## Analytical Framework

### Key Principle

The index is based on **heat-attributable burden**, not raw outcome rates.

---

## Modeling Strategy

Two complementary modeling approaches will be implemented and compared.

---

### Approach 1: DLNM + Second-Stage Pooling (Primary)

#### Stage 1: Endpoint-Specific Models

- Model: Generalized Additive Models (GAMs) with DLNM
- Family: Quasi-Poisson / Negative Binomial
- Exposure: Cross-basis temperature functions
- Lag structure:
  - ED / EMS: 0–5 days  
  - Mortality: 0–10 days  

#### Adjustments

- Long-term time trends  
- Seasonality  
- Day of week  
- Holidays  
- Humidity (if available)  
- Autocorrelation structures  

---

#### Stage 2: Hierarchical Pooling

- Estimate spatial heterogeneity across community areas  
- Model variation in heat-response curves  
- Regress heterogeneity on neighborhood covariates  
- Use mixed-effects or Bayesian meta-regression  

---

### Approach 2: Bayesian Spatiotemporal Hierarchical Model (Secondary)

A unified model incorporating:

- Poisson / Negative Binomial likelihood  
- Nonlinear temperature-response functions  
- Distributed lag structure  
- Spatial random effects  
- Temporal random effects  
- Area-level effect modification  

Potential implementations:

- **INLA (SPDE)**  
- **brms / Stan**  

---

### Comparison Between Approaches

We will evaluate:

- Predictive performance  
- Stability of estimates  
- Interpretability  
- Computational efficiency  
- Spatial smoothness vs overfitting  

---

## Heat-Attributable Burden Estimation

For each endpoint and geography:

- Baseline outcome rate  
- Exposure-response function  
- Lagged effects  
- Attributable number of events  
- Attributable rate per 100,000  

---

## Index Construction

### Definition

HVI 2.0 represents:

> Expected excess health burden under standardized heat conditions.

---

### Strategy A: Transparent Composite (Primary)

1. Standardize endpoint-specific attributable risks  
2. Apply weights based on:
   - Severity  
   - Reliability  
   - Public health relevance  
3. Combine into composite score  

---

### Strategy B: Latent Index (Secondary)

- Hierarchical Bayesian model combining endpoints  
- Derive latent “heat vulnerability” factor  

---

### Outputs

- Composite score  
- Percentile ranking  
- Domain-specific sub-scores:
  - Mortality burden  
  - Acute care burden  
  - EMS burden  
  - (Optional) ICU burden  

---

## Climate Projection Framework

### Pipeline

1. Fit historical exposure-response relationships  
2. Input projected daily temperatures  
3. Generate future heat metrics  
4. Estimate future attributable burden  

---

### Scenarios

- Climate-only change  
- Climate + demographic shifts  
- Climate + adaptation improvements  
- Climate + worsening inequities  

---

### Outputs

- Future attributable events  
- Change from baseline  
- Extreme event burden projections  

---

## Validation Strategy

### Internal Validation

- Train/test split across years  
- Leave-one-year-out validation  
- Leave-one-community-area-out validation  

---

### Endpoint Validation

Evaluate whether high-index areas show:

- Higher excess ED visits  
- Higher excess mortality  
- Higher EMS surge during heat events  

---

### Event-Based Validation

- Test performance during major Chicago heatwaves  
- Compare predicted vs observed burden  

---

### Comparative Validation

Compare HVI 2.0 against:

- HVI 1.0  
- Socioeconomic-only indices  
- Individual covariates (e.g., poverty, age)  

---

### Policy-Relevant Validation

- Proportion of total burden captured in top decile areas  
- Targeting efficiency vs baseline strategies  

---

## Outputs

### 1. Epidemiologic Results
- Exposure-response curves  
- Lag structures  
- Effect modification  

### 2. Attributable Burden Maps
- Events per year  
- Rates per population  

### 3. HVI 2.0 Maps
- Composite index  
- Subdomain indices  
- Uncertainty estimates  

### 4. Future Projection Maps
- Scenario-based risk  
- Temporal trends  

---

## Implementation Stack

- `targets` – pipeline orchestration  
- `data.table` / `arrow` – scalable data processing  
- `sf` / `terra` / `exactextractr` – spatial analysis  
- `mgcv` – GAM modeling  
- `dlnm` – distributed lag models  
- `glmmTMB` – mixed models  
- `INLA` / `brms` – Bayesian modeling  

---

## Reproducibility and Public Exports

This repository should be treated as code plus public aggregate outputs only. Raw and record-level EMS, ED, and mortality data should remain in the private Box-backed project directory configured by `HVI_PRIVATE_DIR`.

Copy `.Renviron.example` to `.Renviron` and adjust paths as needed:

```r
HVI_BOX_DIR="C:/Users/Peter Graffy/Box/HVI2.0"
HVI_PRIVATE_DIR="C:/Users/Peter Graffy/Box/HVI2.0"
HVI_PUBLIC_EXPORT_DIR="public_exports"
HVI_SUPPRESS_SMALL_CELLS=11
```

Build the dashboard/manuscript-safe export bundle with:

```r
source("code/11_build_public_exports.R")
```

The public-facing app should read from `public_exports/dashboard/` and `public_exports/manifest.json`, not from `data/`, `results/`, `code/results/`, or Box raw-data folders. See `docs/DATA_GOVERNANCE.md`, `docs/PUBLIC_EXPORT_CONTRACT.md`, and `docs/PIPELINE.md`.

---

## Project Phases

### Phase 1: Data Engineering
- Harmonize spatial units  
- Construct panel dataset  
- Assign exposures  
- Define outcomes  

### Phase 2: Endpoint Modeling
- Fit DLNM models  
- Estimate attributable burden  

### Phase 3: Spatial Vulnerability Modeling
- Estimate heterogeneity  
- Identify modifiers  

### Phase 4: Index Construction
- Build composite index  
- Compare methodologies  

### Phase 5: Climate Projections
- Apply future temperature scenarios  

### Phase 6: Validation and Deployment
- Validate against observed events  
- Produce maps and dashboards  

---

## Optional Extension: ICU / CLIF Integration

If feasible:

- Incorporate ICU admissions and severe illness  
- Evaluate as:
  - Additional endpoint  
  - Validation layer  
  - Separate analysis  

---

## Expected Contribution

HVI 2.0 moves beyond traditional indices by:

- Directly linking heat exposure to **observed health outcomes**  
- Integrating **multiple healthcare system signals**  
- Capturing **spatiotemporal heterogeneity**  
- Enabling **future climate risk projections**  

---

## Contact

**Peter Graffy**  
University of Chicago  
Center for Computational Medicine & Clinical AI  






