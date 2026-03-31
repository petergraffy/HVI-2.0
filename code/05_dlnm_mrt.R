# =====================================================================================
# HVI 2.0 | Minimum Risk Temperature (MRT) estimation + response curves
# Warm season only: April-October
# Endpoint-specific observation windows
# =====================================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(lubridate)
  library(ggplot2)
  library(dlnm)
  library(splines)
  library(data.table)
  library(purrr)
  library(stringr)
})

# -------------------------------------------------------------------------------------
# Output directories
# -------------------------------------------------------------------------------------

dir.create("results/mrt", recursive = TRUE, showWarnings = FALSE)
dir.create("figures/mrt", recursive = TRUE, showWarnings = FALSE)

# -------------------------------------------------------------------------------------
# Helper: save plots
# -------------------------------------------------------------------------------------

save_mrt_plot <- function(plot_obj, filename, width = 8, height = 6, dpi = 300) {
  ggsave(
    filename = file.path("figures/mrt", filename),
    plot = plot_obj,
    width = width,
    height = height,
    dpi = dpi
  )
}

# -------------------------------------------------------------------------------------
# 0) Prepare warm-season panel
# -------------------------------------------------------------------------------------

panel_mrt <- panel_causes %>%
  mutate(
    event_date = as.Date(event_date),
    year = year(event_date),
    month = month(event_date),
    dow = wday(event_date),
    warm_season = month %in% 5:9,
    has_deaths = year >= 1993 & year < 2023,
    has_ed     = year >= 2011 & year < 2024,
    has_ems    = year >= 2019
  ) %>%
  filter(warm_season)

# -------------------------------------------------------------------------------------
# 1) Build citywide daily totals for each endpoint window
# -------------------------------------------------------------------------------------

city_deaths <- panel_mrt %>%
  filter(has_deaths) %>%
  group_by(event_date, year, month, dow) %>%
  summarise(
    tmax = mean(tmax, na.rm = TRUE),
    tmin = mean(tmin, na.rm = TRUE),
    tmean = mean(tmean, na.rm = TRUE),
    humidity = mean(humidity, na.rm = TRUE),
    
    deaths = sum(deaths, na.rm = TRUE),
    death_cvd = sum(death_cvd, na.rm = TRUE),
    death_respiratory = sum(death_respiratory, na.rm = TRUE),
    death_heat = sum(death_heat, na.rm = TRUE),
    death_dehydration = sum(death_dehydration, na.rm = TRUE),
    .groups = "drop"
  )

city_ed <- panel_mrt %>%
  filter(has_ed) %>%
  group_by(event_date, year, month, dow) %>%
  summarise(
    tmax = mean(tmax, na.rm = TRUE),
    tmin = mean(tmin, na.rm = TRUE),
    tmean = mean(tmean, na.rm = TRUE),
    humidity = mean(humidity, na.rm = TRUE),
    
    ed_visits = sum(ed_visits, na.rm = TRUE),
    ed_cvd = sum(ed_cvd, na.rm = TRUE),
    ed_respiratory = sum(ed_respiratory, na.rm = TRUE),
    ed_heat = sum(ed_heat, na.rm = TRUE),
    ed_dehydration = sum(ed_dehydration, na.rm = TRUE),
    .groups = "drop"
  )

city_ems <- panel_mrt %>%
  filter(has_ems) %>%
  group_by(event_date, year, month, dow) %>%
  summarise(
    tmax = mean(tmax, na.rm = TRUE),
    tmin = mean(tmin, na.rm = TRUE),
    tmean = mean(tmean, na.rm = TRUE),
    humidity = mean(humidity, na.rm = TRUE),
    
    ems_calls = sum(ems_calls, na.rm = TRUE),
    ems_cvd = sum(ems_cvd, na.rm = TRUE),
    ems_respiratory = sum(ems_respiratory, na.rm = TRUE),
    ems_heat = sum(ems_heat, na.rm = TRUE),
    ems_syncope = sum(ems_syncope, na.rm = TRUE),
    ems_neuro = sum(ems_neuro, na.rm = TRUE),
    .groups = "drop"
  )

# -------------------------------------------------------------------------------------
# 2) Generic DLNM MRT function
# -------------------------------------------------------------------------------------
# Notes:
# - var = temperature variable to model, default tmax
# - lag = choose based on endpoint
# - family = quasipoisson
# - control for humidity, day of week, year, month, and long-term time
# -------------------------------------------------------------------------------------

estimate_mrt_dlnm <- function(data, outcome, var = "tmax", lag = 5,
                              outcome_label = outcome) {
  
  dat <- data %>%
    filter(
      !is.na(.data[[outcome]]),
      !is.na(.data[[var]]),
      !is.na(humidity)
    ) %>%
    arrange(event_date) %>%
    mutate(
      time_num = as.numeric(event_date - min(event_date)) + 1,
      dow = factor(lubridate::wday(event_date, label = TRUE, abbr = TRUE)),
      month = factor(lubridate::month(event_date)),
      year = factor(lubridate::year(event_date))
    )
  
  if (nrow(dat) < 50) {
    stop(paste("Not enough data for outcome:", outcome))
  }
  
  is_ems <- stringr::str_detect(outcome, "^ems")
  is_sparse <- mean(dat[[outcome]], na.rm = TRUE) < 0.2 ||
    sum(dat[[outcome]] > 0, na.rm = TRUE) < 100
  
  force_no_lag <- outcome %in% c(
    "death_dehydration",
    "ed_dehydration",
    "death_heat",
    "ed_heat",
    "death_injury",
    "ed_renal",
    "ed_neurologic",
    "ed_mental",
    "ed_gi",
    "ed_injury",
    "ed_syncope"
  )
  
  pred_grid <- seq(
    floor(min(dat[[var]], na.rm = TRUE)),
    ceiling(max(dat[[var]], na.rm = TRUE)),
    by = 0.1
  )
  
  # ------------------------------------------------------------------
  # Generic no-lag spline fallback
  # ------------------------------------------------------------------
  fit_nolag_fallback <- function(model_type = "no_lag_spline_fallback") {
    
    form <- as.formula(
      paste0(
        outcome,
        " ~ ns(", var, ", df = 3) + ns(humidity, 3) + factor(dow) + ns(time_num, df = 4)"
      )
    )
    
    mod <- glm(form, family = quasipoisson(link = "log"), data = dat)
    
    pred_dat <- tibble::tibble(
      temp = pred_grid
    ) %>%
      rename(!!var := temp) %>%
      mutate(
        humidity = median(dat$humidity, na.rm = TRUE),
        dow = factor(levels(dat$dow)[1], levels = levels(dat$dow)),
        time_num = median(dat$time_num, na.rm = TRUE)
      )
    
    pred <- predict(mod, newdata = pred_dat, type = "link", se.fit = TRUE)
    
    fit <- as.numeric(pred$fit)
    se <- as.numeric(pred$se.fit)
    
    mrt <- pred_grid[which.min(fit)]
    ref <- min(fit, na.rm = TRUE)
    
    curve_df <- tibble::tibble(
      outcome = outcome,
      outcome_label = outcome_label,
      temp = pred_grid,
      rr = exp(fit - ref),
      rr_low = exp((fit - 1.96 * se) - ref),
      rr_high = exp((fit + 1.96 * se) - ref),
      mrt = mrt,
      model_type = model_type
    )
    
    mrt_row <- tibble::tibble(
      outcome = outcome,
      outcome_label = outcome_label,
      mrt = mrt,
      n_days = nrow(dat),
      mean_daily_count = mean(dat[[outcome]], na.rm = TRUE),
      model_type = model_type
    )
    
    list(curve = curve_df, mrt = mrt_row, model = mod)
  }
  
  # ------------------------------------------------------------------
  # Forced no-lag branch for dehydration outcomes
  # ------------------------------------------------------------------
  if (force_no_lag) {
    return(fit_nolag_fallback(model_type = "no_lag_forced"))
  }
  
  # ------------------------------------------------------------------
  # DLNM fit helper
  # ------------------------------------------------------------------
  fit_dlnm_once <- function(simple = FALSE) {
    
    if (simple || is_sparse) {
      
      cb <- dlnm::crossbasis(
        dat[[var]],
        lag = min(lag, 3),
        argvar = list(fun = "ns", df = 2),
        arglag = list(fun = "ns", df = 2)
      )
      
      form <- as.formula(
        paste0(
          outcome,
          " ~ cb + ns(humidity, 3) + factor(dow) + ns(time_num, df = 4)"
        )
      )
      
      model_type <- "DLNM_simple"
      
    } else if (is_ems) {
      
      cb <- dlnm::crossbasis(
        dat[[var]],
        lag = min(lag, 3),
        argvar = list(fun = "ns", df = 3),
        arglag = list(fun = "ns", df = 2)
      )
      
      form <- as.formula(
        paste0(
          outcome,
          " ~ cb + ns(humidity, 3) + factor(dow) + ns(time_num, df = 4)"
        )
      )
      
      model_type <- "DLNM_ems"
      
    } else {
      
      temp_knots <- quantile(
        dat[[var]],
        probs = c(0.10, 0.50, 0.90),
        na.rm = TRUE
      )
      
      cb <- dlnm::crossbasis(
        dat[[var]],
        lag = lag,
        argvar = list(fun = "ns", knots = temp_knots),
        arglag = list(fun = "ns", knots = dlnm::logknots(lag, 3))
      )
      
      form <- as.formula(
        paste0(
          outcome,
          " ~ cb + ns(humidity, 3) + factor(dow) + factor(month) + factor(year) + ns(time_num, df = 6)"
        )
      )
      
      model_type <- "DLNM_full"
    }
    
    mod <- glm(form, family = quasipoisson(link = "log"), data = dat)
    
    coef_names <- names(coef(mod))
    if (!any(stringr::str_detect(coef_names, "cb"))) {
      stop("No crossbasis coefficients retained in fitted model.")
    }
    
    cp0 <- dlnm::crosspred(
      cb,
      mod,
      at = pred_grid,
      cen = median(dat[[var]], na.rm = TRUE)
    )
    
    rr_raw <- as.numeric(cp0$allRRfit)
    mrt <- pred_grid[which.min(rr_raw)]
    
    cp <- dlnm::crosspred(
      cb,
      mod,
      at = pred_grid,
      cen = mrt
    )
    
    curve_df <- tibble::tibble(
      outcome = outcome,
      outcome_label = outcome_label,
      temp = pred_grid,
      rr = as.numeric(cp$allRRfit),
      rr_low = as.numeric(cp$allRRlow),
      rr_high = as.numeric(cp$allRRhigh),
      mrt = mrt,
      model_type = model_type
    )
    
    mrt_row <- tibble::tibble(
      outcome = outcome,
      outcome_label = outcome_label,
      mrt = mrt,
      n_days = nrow(dat),
      mean_daily_count = mean(dat[[outcome]], na.rm = TRUE),
      model_type = model_type
    )
    
    list(curve = curve_df, mrt = mrt_row, model = mod)
  }
  
  # ------------------------------------------------------------------
  # Try primary DLNM, then simplified DLNM, then no-lag fallback
  # ------------------------------------------------------------------
  out <- tryCatch(
    fit_dlnm_once(simple = FALSE),
    error = function(e1) {
      message("Primary DLNM failed for ", outcome, ": ", e1$message)
      
      tryCatch(
        fit_dlnm_once(simple = TRUE),
        error = function(e2) {
          message("Simplified DLNM failed for ", outcome, ": ", e2$message)
          fit_nolag_fallback(model_type = "no_lag_fallback_after_dlnm_failure")
        }
      )
    }
  )
  
  out
}

# -------------------------------------------------------------------------------------
# 3) Outcomes to estimate
# Use endpoint-specific windows
# -------------------------------------------------------------------------------------

death_outcomes <- tribble(
  ~outcome,              ~label,                    ~lag,
  "deaths",              "All-Cause Mortality",     7,
  "death_cvd",           "CVD Mortality",           7,
  "death_respiratory",   "Respiratory Mortality",   7,
  "death_heat",          "Heat-Related Mortality",  3,
  "death_dehydration",   "Dehydration Mortality",   3
)

ed_outcomes <- tribble(
  ~outcome,             ~label,                  ~lag,
  "ed_visits",          "All-Cause ED Visits",   5,
  "ed_cvd",             "CVD ED Visits",         5,
  "ed_respiratory",     "Respiratory ED Visits", 5,
  "ed_heat",            "Heat-Related ED Visits",3,
  "ed_dehydration",     "Dehydration ED Visits", 3
)

ems_outcomes <- tribble(
  ~outcome,             ~label,                     ~lag,
  "ems_calls",          "All-Cause EMS Calls",      3,
  "ems_cvd",            "CVD EMS Calls",            3,
  "ems_respiratory",    "Respiratory EMS Calls",    3,
  "ems_heat",           "Heat-Related EMS Calls",   3,
  "ems_syncope",        "Syncope EMS Calls",        3,
  "ems_neuro",          "Neurologic EMS Calls",     3
)

# -------------------------------------------------------------------------------------
# 4) Run MRT estimation
# -------------------------------------------------------------------------------------

run_outcome_set <- function(df, outcomes_tbl, temp_var = "tmax") {
  purrr::pmap(
    list(outcomes_tbl$outcome, outcomes_tbl$label, outcomes_tbl$lag),
    function(outcome, label, lag) {
      tryCatch(
        estimate_mrt_dlnm(
          data = df,
          outcome = outcome,
          var = temp_var,
          lag = lag,
          outcome_label = label
        ),
        error = function(e) {
          message("Failed for outcome: ", outcome, " | ", e$message)
          list(
            model = NULL,
            crossbasis = NULL,
            curve = tibble(
              outcome = character(),
              outcome_label = character(),
              temp = numeric(),
              rr = numeric(),
              rr_low = numeric(),
              rr_high = numeric(),
              mrt = numeric(),
              lag = numeric(),
              temp_var = character()
            ),
            mrt = tibble(
              outcome = outcome,
              outcome_label = label,
              temp_var = temp_var,
              lag = lag,
              mrt = NA_real_,
              min_temp = NA_real_,
              max_temp = NA_real_,
              n_days = nrow(df),
              mean_daily_count = NA_real_
            )
          )
        }
      )
    }
  )
}

screen_outcome <- function(data, outcome, min_nonzero_days = 25, min_total_events = 100) {
  y <- data[[outcome]]
  
  tibble(
    outcome = outcome,
    n_days = sum(!is.na(y)),
    nonzero_days = sum(y > 0, na.rm = TRUE),
    total_events = sum(y, na.rm = TRUE),
    keep_for_dlnm = nonzero_days >= min_nonzero_days & total_events >= min_total_events
  )
}

death_screen <- bind_rows(lapply(death_outcomes$outcome, function(x) screen_outcome(city_deaths, x)))
ed_screen    <- bind_rows(lapply(ed_outcomes$outcome, function(x) screen_outcome(city_ed, x)))
ems_screen   <- bind_rows(lapply(ems_outcomes$outcome, function(x) screen_outcome(city_ems, x)))

death_screen
ed_screen
ems_screen

death_outcomes_keep <- death_outcomes %>%
  inner_join(death_screen %>% filter(keep_for_dlnm), by = c("outcome"))

ed_outcomes_keep <- ed_outcomes %>%
  inner_join(ed_screen %>% filter(keep_for_dlnm), by = c("outcome"))

ems_outcomes_keep <- ems_outcomes %>%
  inner_join(ems_screen %>% filter(keep_for_dlnm), by = c("outcome"))

death_results <- run_outcome_set(city_deaths, death_outcomes_keep, temp_var = "tmax")
ed_results    <- run_outcome_set(city_ed, ed_outcomes_keep, temp_var = "tmax")
ems_results   <- run_outcome_set(city_ems, ems_outcomes_keep, temp_var = "tmax")

all_results <- c(death_results, ed_results, ems_results)

# -------------------------------------------------------------------------------------
# 5) Collect tables
# -------------------------------------------------------------------------------------

mrt_table <- bind_rows(lapply(all_results, function(x) x$mrt))
curve_table <- bind_rows(lapply(all_results, function(x) x$curve))

fwrite(mrt_table, "results/mrt/mrt_table_tmax_warmseason.csv")
fwrite(curve_table, "results/mrt/mrt_response_curves_tmax_warmseason.csv")

# -------------------------------------------------------------------------------------
# 6) Plot response curves
# -------------------------------------------------------------------------------------

plot_one_curve <- function(curve_df, title_text) {
  mrt_val <- unique(curve_df$mrt)[1]
  
  ggplot(curve_df, aes(x = temp, y = rr)) +
    geom_ribbon(aes(ymin = rr_low, ymax = rr_high), alpha = 0.2, fill = "steelblue") +
    geom_line(linewidth = 1, color = "steelblue4") +
    geom_hline(yintercept = 1, linetype = "dashed", color = "gray40") +
    geom_vline(xintercept = mrt_val, linetype = "dotted", color = "firebrick", linewidth = 1) +
    annotate("text", x = mrt_val, y = max(curve_df$rr_high, na.rm = TRUE),
             label = paste0("MRT = ", round(mrt_val, 1), "°C"),
             vjust = -0.4, color = "firebrick", size = 4) +
    labs(
      title = title_text,
      subtitle = "Cumulative relative risk from warm-season DLNM",
      x = "Daily Maximum Temperature (°C)",
      y = "Relative Risk"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold")
    )
}

for (i in seq_along(all_results)) {
  obj <- all_results[[i]]
  curve_df <- obj$curve
  lab <- unique(curve_df$outcome_label)[1]
  
  p <- plot_one_curve(curve_df, lab)
  
  file_stub <- curve_df$outcome[1]
  
  save_mrt_plot(p, paste0(file_stub, "_mrt_curve.png"), width = 8, height = 6)
  save_mrt_plot(p, paste0(file_stub, "_mrt_curve.pdf"), width = 8, height = 6)
}

# -------------------------------------------------------------------------------------
# 7) Faceted plots by endpoint class
# -------------------------------------------------------------------------------------

curve_table <- curve_table %>%
  mutate(
    domain = case_when(
      str_detect(outcome, "^death") ~ "Mortality",
      str_detect(outcome, "^ed") ~ "ED",
      str_detect(outcome, "^ems") ~ "EMS",
      TRUE ~ "Other"
    )
  )

plot_curves_mort <- curve_table %>%
  filter(domain == "Mortality") %>%
  ggplot(aes(x = temp, y = rr)) +
  geom_ribbon(aes(ymin = rr_low, ymax = rr_high), alpha = 0.15, fill = "steelblue") +
  geom_line(color = "steelblue4", linewidth = 0.9) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray40") +
  geom_vline(aes(xintercept = mrt), linetype = "dotted", color = "firebrick") +
  facet_wrap(~ outcome_label, scales = "free_y") +
  labs(
    title = "Mortality Temperature-Response Curves",
    subtitle = "Warm season only (April-October)",
    x = "Daily Maximum Temperature (°C)",
    y = "Relative Risk"
  ) +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold"))

plot_curves_ed <- curve_table %>%
  filter(domain == "ED") %>%
  ggplot(aes(x = temp, y = rr)) +
  geom_ribbon(aes(ymin = rr_low, ymax = rr_high), alpha = 0.15, fill = "darkorange") +
  geom_line(color = "darkorange4", linewidth = 0.9) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray40") +
  geom_vline(aes(xintercept = mrt), linetype = "dotted", color = "firebrick") +
  facet_wrap(~ outcome_label, scales = "free_y") +
  labs(
    title = "ED Temperature-Response Curves",
    subtitle = "Warm season only (April-October)",
    x = "Daily Maximum Temperature (°C)",
    y = "Relative Risk"
  ) +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold"))

plot_curves_ems <- curve_table %>%
  filter(domain == "EMS") %>%
  ggplot(aes(x = temp, y = rr)) +
  geom_ribbon(aes(ymin = rr_low, ymax = rr_high), alpha = 0.15, fill = "darkgreen") +
  geom_line(color = "darkgreen", linewidth = 0.9) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray40") +
  geom_vline(aes(xintercept = mrt), linetype = "dotted", color = "firebrick") +
  facet_wrap(~ outcome_label, scales = "free_y") +
  labs(
    title = "EMS Temperature-Response Curves",
    subtitle = "Warm season only (April-October)",
    x = "Daily Maximum Temperature (°C)",
    y = "Relative Risk"
  ) +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold"))

save_mrt_plot(plot_curves_mort, "mortality_mrt_curves_faceted.png", width = 12, height = 8)
save_mrt_plot(plot_curves_mort, "mortality_mrt_curves_faceted.pdf", width = 12, height = 8)

save_mrt_plot(plot_curves_ed, "ed_mrt_curves_faceted.png", width = 12, height = 8)
save_mrt_plot(plot_curves_ed, "ed_mrt_curves_faceted.pdf", width = 12, height = 8)

save_mrt_plot(plot_curves_ems, "ems_mrt_curves_faceted.png", width = 12, height = 8)
save_mrt_plot(plot_curves_ems, "ems_mrt_curves_faceted.pdf", width = 12, height = 8)

# -------------------------------------------------------------------------------------
# 8) Suggested knot placements for future iterations
# Based on all-cause MRTs + exposure distribution
# -------------------------------------------------------------------------------------

knot_suggestions <- mrt_table %>%
  filter(outcome %in% c("deaths", "ed_visits", "ems_calls")) %>%
  select(outcome_label, mrt) %>%
  mutate(
    suggestion = paste0(
      "Center future DLNM at MRT ≈ ", round(mrt, 1),
      "°C; place exposure knots around empirical quantiles and verify stability."
    )
  )

fwrite(knot_suggestions, "results/mrt/mrt_knot_suggestions.csv")

# -------------------------------------------------------------------------------------
# 9) Console output
# -------------------------------------------------------------------------------------

cat("\nMinimum Risk Temperatures (warm season):\n")
print(mrt_table %>% arrange(domain = outcome_label))

cat("\nSaved MRT tables to results/mrt/\n")
cat("Saved response-curve figures to figures/mrt/\n")


# =====================================================================================
# 10) EXPANDED CAUSE-SPECIFIC MRT MODELING
# Add broader cause-specific endpoints and domain-level plots
# This block piggybacks off objects/functions already created above:
#   - panel_mrt
#   - estimate_mrt_dlnm()
#   - run_outcome_set()
#   - screen_outcome()
#   - save_mrt_plot()
# =====================================================================================

dir.create("results/mrt_expanded", recursive = TRUE, showWarnings = FALSE)
dir.create("figures/mrt_expanded", recursive = TRUE, showWarnings = FALSE)

# -------------------------------------------------------------------------------------
# Rebuild citywide daily totals with ALL available cause-specific outcomes
# -------------------------------------------------------------------------------------

city_deaths_full <- panel_mrt %>%
  filter(has_deaths) %>%
  group_by(event_date, year, month, dow) %>%
  summarise(
    tmax = mean(tmax, na.rm = TRUE),
    tmin = mean(tmin, na.rm = TRUE),
    tmean = mean(tmean, na.rm = TRUE),
    humidity = mean(humidity, na.rm = TRUE),
    
    deaths = sum(deaths, na.rm = TRUE),
    death_cvd = sum(death_cvd, na.rm = TRUE),
    death_respiratory = sum(death_respiratory, na.rm = TRUE),
    death_renal = sum(death_renal, na.rm = TRUE),
    death_heat = sum(death_heat, na.rm = TRUE),
    death_dehydration = sum(death_dehydration, na.rm = TRUE),
    death_syncope = sum(death_syncope, na.rm = TRUE),
    death_mental = sum(death_mental, na.rm = TRUE),
    death_neurologic = sum(death_neurologic, na.rm = TRUE),
    death_injury = sum(death_injury, na.rm = TRUE),
    death_gi = sum(death_gi, na.rm = TRUE),
    .groups = "drop"
  )

city_ed_full <- panel_mrt %>%
  filter(has_ed) %>%
  group_by(event_date, year, month, dow) %>%
  summarise(
    tmax = mean(tmax, na.rm = TRUE),
    tmin = mean(tmin, na.rm = TRUE),
    tmean = mean(tmean, na.rm = TRUE),
    humidity = mean(humidity, na.rm = TRUE),
    
    ed_visits = sum(ed_visits, na.rm = TRUE),
    ed_cvd = sum(ed_cvd, na.rm = TRUE),
    ed_respiratory = sum(ed_respiratory, na.rm = TRUE),
    ed_renal = sum(ed_renal, na.rm = TRUE),
    ed_heat = sum(ed_heat, na.rm = TRUE),
    ed_dehydration = sum(ed_dehydration, na.rm = TRUE),
    ed_syncope = sum(ed_syncope, na.rm = TRUE),
    ed_mental = sum(ed_mental, na.rm = TRUE),
    ed_neurologic = sum(ed_neurologic, na.rm = TRUE),
    ed_injury = sum(ed_injury, na.rm = TRUE),
    ed_gi = sum(ed_gi, na.rm = TRUE),
    .groups = "drop"
  )

city_ems_full <- panel_mrt %>%
  filter(has_ems) %>%
  group_by(event_date, year, month, dow) %>%
  summarise(
    tmax = mean(tmax, na.rm = TRUE),
    tmin = mean(tmin, na.rm = TRUE),
    tmean = mean(tmean, na.rm = TRUE),
    humidity = mean(humidity, na.rm = TRUE),
    
    ems_calls = sum(ems_calls, na.rm = TRUE),
    ems_cvd = sum(ems_cvd, na.rm = TRUE),
    ems_respiratory = sum(ems_respiratory, na.rm = TRUE),
    ems_heat = sum(ems_heat, na.rm = TRUE),
    ems_syncope = sum(ems_syncope, na.rm = TRUE),
    ems_neuro = sum(ems_neuro, na.rm = TRUE),
    ems_mental = sum(ems_mental, na.rm = TRUE),
    ems_injury = sum(ems_injury, na.rm = TRUE),
    ems_gi = sum(ems_gi, na.rm = TRUE),
    ems_bleeding = sum(ems_bleeding, na.rm = TRUE),
    .groups = "drop"
  )

# -------------------------------------------------------------------------------------
# Expanded outcome sets
# Note:
# - lag 7 for mortality domains
# - lag 5 for most ED domains
# - lag 3 for most EMS domains
# - lag 1 for acute syndromes like injury/syncope
# - dehydration/heat stay in the table but your function may force no-lag
# -------------------------------------------------------------------------------------

death_outcomes_expanded <- tibble::tribble(
  ~outcome,              ~label,                     ~lag,
  "deaths",              "All-Cause Mortality",      7,
  "death_cvd",           "CVD Mortality",            7,
  "death_respiratory",   "Respiratory Mortality",    7,
  "death_renal",         "Renal Mortality",          7,
  "death_neurologic",    "Neurologic Mortality",     7,
  "death_mental",        "Mental Health Mortality",  7,
  "death_gi",            "GI Mortality",             5,
  "death_injury",        "Injury Mortality",         3,
  "death_syncope",       "Syncope Mortality",        3,
  "death_dehydration",   "Dehydration Mortality",    3,
  "death_heat",          "Heat-Related Mortality",   3
)

ed_outcomes_expanded <- tibble::tribble(
  ~outcome,              ~label,                    ~lag,
  "ed_visits",           "All-Cause ED Visits",     5,
  "ed_cvd",              "CVD ED Visits",           5,
  "ed_respiratory",      "Respiratory ED Visits",   5,
  "ed_renal",            "Renal ED Visits",         3,
  "ed_neurologic",       "Neurologic ED Visits",    3,
  "ed_mental",           "Mental Health ED Visits", 3,
  "ed_gi",               "GI ED Visits",            3,
  "ed_injury",           "Injury ED Visits",        1,
  "ed_syncope",          "Syncope ED Visits",       1,
  "ed_dehydration",      "Dehydration ED Visits",   3,
  "ed_heat",             "Heat-Related ED Visits",  3
)

ems_outcomes_expanded <- tibble::tribble(
  ~outcome,              ~label,                     ~lag,
  "ems_calls",           "All-Cause EMS Calls",      3,
  "ems_cvd",             "CVD EMS Calls",            3,
  "ems_respiratory",     "Respiratory EMS Calls",    3,
  "ems_neuro",           "Neurologic EMS Calls",     3,
  "ems_mental",          "Mental Health EMS Calls",  3,
  "ems_gi",              "GI EMS Calls",             3,
  "ems_bleeding",        "Bleeding EMS Calls",       3,
  "ems_injury",          "Injury EMS Calls",         1,
  "ems_syncope",         "Syncope EMS Calls",        1,
  "ems_heat",            "Heat-Related EMS Calls",   3
)

# -------------------------------------------------------------------------------------
# Only keep outcomes that actually exist in the rebuilt citywide datasets
# -------------------------------------------------------------------------------------

death_outcomes_expanded <- death_outcomes_expanded %>%
  filter(outcome %in% names(city_deaths_full))

ed_outcomes_expanded <- ed_outcomes_expanded %>%
  filter(outcome %in% names(city_ed_full))

ems_outcomes_expanded <- ems_outcomes_expanded %>%
  filter(outcome %in% names(city_ems_full))

# -------------------------------------------------------------------------------------
# Screen outcomes before modeling
# -------------------------------------------------------------------------------------

death_screen_expanded <- bind_rows(lapply(death_outcomes_expanded$outcome, function(x) screen_outcome(city_deaths_full, x)))
ed_screen_expanded    <- bind_rows(lapply(ed_outcomes_expanded$outcome,    function(x) screen_outcome(city_ed_full, x)))
ems_screen_expanded   <- bind_rows(lapply(ems_outcomes_expanded$outcome,   function(x) screen_outcome(city_ems_full, x)))

fwrite(death_screen_expanded, "results/mrt_expanded/death_screen_expanded.csv")
fwrite(ed_screen_expanded,    "results/mrt_expanded/ed_screen_expanded.csv")
fwrite(ems_screen_expanded,   "results/mrt_expanded/ems_screen_expanded.csv")

death_outcomes_keep_expanded <- death_outcomes_expanded %>%
  inner_join(death_screen_expanded %>% filter(keep_for_dlnm), by = "outcome")

ed_outcomes_keep_expanded <- ed_outcomes_expanded %>%
  inner_join(ed_screen_expanded %>% filter(keep_for_dlnm), by = "outcome")

ems_outcomes_keep_expanded <- ems_outcomes_expanded %>%
  inner_join(ems_screen_expanded %>% filter(keep_for_dlnm), by = "outcome")

# -------------------------------------------------------------------------------------
# Run MRT estimation
# -------------------------------------------------------------------------------------

death_results_expanded <- run_outcome_set(city_deaths_full, death_outcomes_keep_expanded, temp_var = "tmax")
ed_results_expanded    <- run_outcome_set(city_ed_full,     ed_outcomes_keep_expanded,    temp_var = "tmax")
ems_results_expanded   <- run_outcome_set(city_ems_full,    ems_outcomes_keep_expanded,   temp_var = "tmax")

all_results_expanded <- c(death_results_expanded, ed_results_expanded, ems_results_expanded)

# -------------------------------------------------------------------------------------
# Collect tables
# -------------------------------------------------------------------------------------

mrt_table_expanded <- bind_rows(lapply(all_results_expanded, function(x) x$mrt))
curve_table_expanded <- bind_rows(lapply(all_results_expanded, function(x) x$curve))

# add endpoint and physiologic domain labels
curve_table_expanded <- curve_table_expanded %>%
  mutate(
    endpoint = case_when(
      str_detect(outcome, "^death") ~ "Mortality",
      str_detect(outcome, "^ed")    ~ "ED",
      str_detect(outcome, "^ems")   ~ "EMS",
      TRUE ~ "Other"
    ),
    domain = case_when(
      outcome %in% c("deaths", "ed_visits", "ems_calls") ~ "All-cause",
      str_detect(outcome, "cvd")           ~ "Cardiovascular",
      str_detect(outcome, "respiratory")   ~ "Respiratory",
      str_detect(outcome, "renal")         ~ "Renal",
      str_detect(outcome, "neuro")         ~ "Neurologic",
      str_detect(outcome, "mental")        ~ "Mental Health",
      str_detect(outcome, "injury")        ~ "Injury",
      str_detect(outcome, "gi")            ~ "Gastrointestinal",
      str_detect(outcome, "bleeding")      ~ "Bleeding",
      str_detect(outcome, "syncope")       ~ "Syncope",
      str_detect(outcome, "dehydration|heat") ~ "Heat-related",
      TRUE ~ "Other"
    )
  )

mrt_table_expanded <- mrt_table_expanded %>%
  mutate(
    endpoint = case_when(
      str_detect(outcome, "^death") ~ "Mortality",
      str_detect(outcome, "^ed")    ~ "ED",
      str_detect(outcome, "^ems")   ~ "EMS",
      TRUE ~ "Other"
    ),
    domain = case_when(
      outcome %in% c("deaths", "ed_visits", "ems_calls") ~ "All-cause",
      str_detect(outcome, "cvd")           ~ "Cardiovascular",
      str_detect(outcome, "respiratory")   ~ "Respiratory",
      str_detect(outcome, "renal")         ~ "Renal",
      str_detect(outcome, "neuro")         ~ "Neurologic",
      str_detect(outcome, "mental")        ~ "Mental Health",
      str_detect(outcome, "injury")        ~ "Injury",
      str_detect(outcome, "gi")            ~ "Gastrointestinal",
      str_detect(outcome, "bleeding")      ~ "Bleeding",
      str_detect(outcome, "syncope")       ~ "Syncope",
      str_detect(outcome, "dehydration|heat") ~ "Heat-related",
      TRUE ~ "Other"
    )
  )

fwrite(mrt_table_expanded, "results/mrt_expanded/mrt_table_tmax_warmseason_expanded.csv")
fwrite(curve_table_expanded, "results/mrt_expanded/mrt_response_curves_tmax_warmseason_expanded.csv")

# -------------------------------------------------------------------------------------
# Domain-level faceted figures across all endpoints
# -------------------------------------------------------------------------------------

plot_by_domain <- function(df, domain_name, fill_col, line_col) {
  ggplot(df, aes(x = temp, y = rr)) +
    geom_ribbon(aes(ymin = rr_low, ymax = rr_high), alpha = 0.18, fill = fill_col) +
    geom_line(color = line_col, linewidth = 1) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "gray40") +
    geom_vline(aes(xintercept = mrt), linetype = "dotted", color = "firebrick") +
    facet_wrap(~ outcome_label, scales = "free_y") +
    labs(
      title = paste0(domain_name, " Temperature-Response Curves"),
      subtitle = "Expanded cause-specific outcomes",
      x = "Daily Maximum Temperature (°C)",
      y = "Relative Risk"
    ) +
    theme_minimal(base_size = 11) +
    theme(plot.title = element_text(face = "bold"))
}

domain_names <- unique(curve_table_expanded$domain)

domain_plot_list <- lapply(domain_names, function(dom) {
  dom_df <- curve_table_expanded %>% filter(domain == dom)
  
  fill_col <- dplyr::case_when(
    dom == "All-cause" ~ "gray70",
    dom == "Cardiovascular" ~ "firebrick2",
    dom == "Respiratory" ~ "steelblue",
    dom == "Renal" ~ "purple",
    dom == "Neurologic" ~ "darkgreen",
    dom == "Mental Health" ~ "orchid4",
    dom == "Injury" ~ "tan4",
    dom == "Gastrointestinal" ~ "darkorange",
    dom == "Bleeding" ~ "brown3",
    dom == "Syncope" ~ "goldenrod",
    dom == "Heat-related" ~ "red3",
    TRUE ~ "gray70"
  )
  
  line_col <- dplyr::case_when(
    dom == "All-cause" ~ "gray25",
    dom == "Cardiovascular" ~ "firebrick4",
    dom == "Respiratory" ~ "steelblue4",
    dom == "Renal" ~ "purple4",
    dom == "Neurologic" ~ "darkgreen",
    dom == "Mental Health" ~ "orchid4",
    dom == "Injury" ~ "tan4",
    dom == "Gastrointestinal" ~ "darkorange4",
    dom == "Bleeding" ~ "brown4",
    dom == "Syncope" ~ "goldenrod4",
    dom == "Heat-related" ~ "red4",
    TRUE ~ "gray30"
  )
  
  p <- plot_by_domain(dom_df, dom, fill_col, line_col)
  
  file_stub <- str_replace_all(str_to_lower(dom), "[^a-z0-9]+", "_")
  
  save_mrt_plot(p, paste0("expanded_", file_stub, "_mrt_curves.png"), width = 12, height = 8)
  save_mrt_plot(p, paste0("expanded_", file_stub, "_mrt_curves.pdf"), width = 12, height = 8)
  
  p
})

# -------------------------------------------------------------------------------------
# Endpoint-wide expanded figures
# -------------------------------------------------------------------------------------

plot_curves_mort_expanded <- curve_table_expanded %>%
  filter(endpoint == "Mortality") %>%
  ggplot(aes(x = temp, y = rr)) +
  geom_ribbon(aes(ymin = rr_low, ymax = rr_high), alpha = 0.15, fill = "steelblue") +
  geom_line(color = "steelblue4", linewidth = 0.9) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray40") +
  geom_vline(aes(xintercept = mrt), linetype = "dotted", color = "firebrick") +
  facet_wrap(~ outcome_label, scales = "free_y") +
  labs(
    title = "Expanded Mortality Temperature-Response Curves",
    subtitle = "Warm season only",
    x = "Daily Maximum Temperature (°C)",
    y = "Relative Risk"
  ) +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold"))

plot_curves_ed_expanded <- curve_table_expanded %>%
  filter(endpoint == "ED") %>%
  ggplot(aes(x = temp, y = rr)) +
  geom_ribbon(aes(ymin = rr_low, ymax = rr_high), alpha = 0.15, fill = "darkorange") +
  geom_line(color = "darkorange4", linewidth = 0.9) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray40") +
  geom_vline(aes(xintercept = mrt), linetype = "dotted", color = "firebrick") +
  facet_wrap(~ outcome_label, scales = "free_y") +
  labs(
    title = "Expanded ED Temperature-Response Curves",
    subtitle = "Warm season only",
    x = "Daily Maximum Temperature (°C)",
    y = "Relative Risk"
  ) +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold"))

plot_curves_ems_expanded <- curve_table_expanded %>%
  filter(endpoint == "EMS") %>%
  ggplot(aes(x = temp, y = rr)) +
  geom_ribbon(aes(ymin = rr_low, ymax = rr_high), alpha = 0.15, fill = "darkgreen") +
  geom_line(color = "darkgreen", linewidth = 0.9) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray40") +
  geom_vline(aes(xintercept = mrt), linetype = "dotted", color = "firebrick") +
  facet_wrap(~ outcome_label, scales = "free_y") +
  labs(
    title = "Expanded EMS Temperature-Response Curves",
    subtitle = "Warm season only",
    x = "Daily Maximum Temperature (°C)",
    y = "Relative Risk"
  ) +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold"))

save_mrt_plot(plot_curves_mort_expanded, "expanded_mortality_mrt_curves_faceted.png", width = 14, height = 10)
save_mrt_plot(plot_curves_mort_expanded, "expanded_mortality_mrt_curves_faceted.pdf", width = 14, height = 10)

save_mrt_plot(plot_curves_ed_expanded, "expanded_ed_mrt_curves_faceted.png", width = 14, height = 10)
save_mrt_plot(plot_curves_ed_expanded, "expanded_ed_mrt_curves_faceted.pdf", width = 14, height = 10)

save_mrt_plot(plot_curves_ems_expanded, "expanded_ems_mrt_curves_faceted.png", width = 14, height = 10)
save_mrt_plot(plot_curves_ems_expanded, "expanded_ems_mrt_curves_faceted.pdf", width = 14, height = 10)

# -------------------------------------------------------------------------------------
# Optional: a single physiologic-domain summary table
# -------------------------------------------------------------------------------------

mrt_domain_summary <- mrt_table_expanded %>%
  arrange(endpoint, domain, outcome_label)

fwrite(mrt_domain_summary, "results/mrt_expanded/mrt_domain_summary.csv")

cat("\nExpanded MRT modeling complete.\n")
cat("Saved expanded tables to results/mrt_expanded/\n")
cat("Saved expanded figures to figures/mrt_expanded/\n")
print(mrt_table_expanded %>% select(endpoint, domain, outcome_label, mrt, mean_daily_count, model_type) %>% arrange(endpoint, domain, outcome_label))

