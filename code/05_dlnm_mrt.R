# =====================================================================================
# HVI 2.0 | 05_dlnm_mrt.R
# Minimum Risk Temperature (MRT) estimation for HVI 2.0
# Primary analysis: May-September, 2019-2022
# Sensitivity analysis: full warm-season history by source availability
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
  library(tibble)
})

# -------------------------------------------------------------------------------------
# Directories
# -------------------------------------------------------------------------------------

dir.create("results/mrt", recursive = TRUE, showWarnings = FALSE)
dir.create("figures/mrt", recursive = TRUE, showWarnings = FALSE)

# -------------------------------------------------------------------------------------
# Config
# -------------------------------------------------------------------------------------

analysis_year_start  <- 2019L
analysis_year_end    <- 2022L
analysis_month_start <- 5L
analysis_month_end   <- 9L

primary_label     <- "2019_2022_overlap"
sensitivity_label <- "full_history"

min_nonzero_days <- 25
min_total_events <- 100

# -------------------------------------------------------------------------------------
# Helpers
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

in_warm_season <- function(date) {
  mo <- lubridate::month(date)
  mo >= analysis_month_start & mo <= analysis_month_end
}

in_primary_window <- function(date) {
  yr <- lubridate::year(date)
  mo <- lubridate::month(date)
  yr >= analysis_year_start &
    yr <= analysis_year_end &
    mo >= analysis_month_start &
    mo <= analysis_month_end
}

classify_endpoint <- function(outcome) {
  case_when(
    str_detect(outcome, "^death") ~ "Mortality",
    str_detect(outcome, "^ed")    ~ "ED",
    str_detect(outcome, "^ems")   ~ "EMS",
    TRUE ~ "Other"
  )
}

classify_domain <- function(outcome) {
  case_when(
    outcome %in% c("deaths", "ed_visits", "ems_calls") ~ "All-cause",
    str_detect(outcome, "cvd") ~ "Cardiovascular",
    str_detect(outcome, "respiratory") ~ "Respiratory",
    str_detect(outcome, "renal") ~ "Renal",
    str_detect(outcome, "neuro") ~ "Neurologic",
    str_detect(outcome, "mental") ~ "Mental Health",
    str_detect(outcome, "injury") ~ "Injury",
    str_detect(outcome, "gi") ~ "Gastrointestinal",
    str_detect(outcome, "bleeding") ~ "Bleeding",
    str_detect(outcome, "syncope") ~ "Syncope",
    str_detect(outcome, "dehydration|heat") ~ "Heat-related",
    TRUE ~ "Other"
  )
}

screen_outcome <- function(data, outcome,
                           min_nonzero_days = 25,
                           min_total_events = 100) {
  y <- data[[outcome]]
  
  tibble(
    outcome = outcome,
    n_days = sum(!is.na(y)),
    nonzero_days = sum(y > 0, na.rm = TRUE),
    total_events = sum(y, na.rm = TRUE),
    keep_for_dlnm = nonzero_days >= min_nonzero_days &
      total_events >= min_total_events
  )
}

build_city_daily <- function(data, outcome_vars) {
  data %>%
    group_by(event_date, year, month) %>%
    summarise(
      tmax = mean(tmax, na.rm = TRUE),
      tmin = mean(tmin, na.rm = TRUE),
      tmean = mean(tmean, na.rm = TRUE),
      humidity = mean(humidity, na.rm = TRUE),
      across(all_of(outcome_vars), ~ sum(.x, na.rm = TRUE)),
      .groups = "drop"
    ) %>%
    mutate(
      dow = factor(
        wday(event_date),
        levels = levels(wday(seq.Date(as.Date("2023-01-01"), by = "day", length.out = 7)))
      )
    ) %>%
    arrange(event_date)
}

plot_one_curve <- function(curve_df, title_text, subtitle_text = NULL) {
  if (nrow(curve_df) == 0 || all(is.na(curve_df$mrt))) return(NULL)
  
  mrt_val <- unique(curve_df$mrt)[1]
  
  ggplot(curve_df, aes(x = temp, y = rr)) +
    geom_ribbon(aes(ymin = rr_low, ymax = rr_high), alpha = 0.2, fill = "steelblue") +
    geom_line(linewidth = 1, color = "steelblue4") +
    geom_hline(yintercept = 1, linetype = "dashed", color = "gray40") +
    geom_vline(xintercept = mrt_val, linetype = "dotted", color = "firebrick", linewidth = 1) +
    labs(
      title = title_text,
      subtitle = subtitle_text,
      x = "Daily Maximum Temperature (°C)",
      y = "Relative Risk"
    ) +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(face = "bold"))
}

# -------------------------------------------------------------------------------------
# Outcome specifications
# -------------------------------------------------------------------------------------

death_outcomes <- tribble(
  ~outcome,              ~outcome_label,               ~lag_input,
  "deaths",              "All-Cause Mortality",        7L,
  "death_cvd",           "CVD Mortality",              7L,
  "death_respiratory",   "Respiratory Mortality",      7L,
  "death_renal",         "Renal Mortality",            7L,
  "death_neurologic",    "Neurologic Mortality",       7L,
  "death_mental",        "Mental Health Mortality",    7L,
  "death_gi",            "GI Mortality",               5L,
  "death_injury",        "Injury Mortality",           3L,
  "death_syncope",       "Syncope Mortality",          3L,
  "death_dehydration",   "Dehydration Mortality",      3L,
  "death_heat",          "Heat-Related Mortality",     3L
)

ed_outcomes <- tribble(
  ~outcome,              ~outcome_label,               ~lag_input,
  "ed_visits",           "All-Cause ED Visits",        5L,
  "ed_cvd",              "CVD ED Visits",              5L,
  "ed_respiratory",      "Respiratory ED Visits",      5L,
  "ed_renal",            "Renal ED Visits",            3L,
  "ed_neurologic",       "Neurologic ED Visits",       3L,
  "ed_mental",           "Mental Health ED Visits",    3L,
  "ed_gi",               "GI ED Visits",               3L,
  "ed_injury",           "Injury ED Visits",           1L,
  "ed_syncope",          "Syncope ED Visits",          1L,
  "ed_dehydration",      "Dehydration ED Visits",      3L,
  "ed_heat",             "Heat-Related ED Visits",     3L
)

ems_outcomes <- tribble(
  ~outcome,              ~outcome_label,               ~lag_input,
  "ems_calls",           "All-Cause EMS Calls",        3L,
  "ems_cvd",             "CVD EMS Calls",              3L,
  "ems_respiratory",     "Respiratory EMS Calls",      3L,
  "ems_neuro",           "Neurologic EMS Calls",       3L,
  "ems_mental",          "Mental Health EMS Calls",    3L,
  "ems_gi",              "GI EMS Calls",               3L,
  "ems_bleeding",        "Bleeding EMS Calls",         3L,
  "ems_injury",          "Injury EMS Calls",           1L,
  "ems_syncope",         "Syncope EMS Calls",          1L,
  "ems_heat",            "Heat-Related EMS Calls",     3L
)

all_outcomes_spec <- bind_rows(
  death_outcomes %>% mutate(endpoint = "Mortality"),
  ed_outcomes    %>% mutate(endpoint = "ED"),
  ems_outcomes   %>% mutate(endpoint = "EMS")
) %>%
  mutate(domain = classify_domain(outcome))

# -------------------------------------------------------------------------------------
# Prepare panel
# Requires panel_causes in memory
# -------------------------------------------------------------------------------------

panel_mrt_all <- panel_causes %>%
  mutate(
    event_date = as.Date(event_date),
    year = year(event_date),
    month = month(event_date),
    has_deaths = year >= 1993 & year <= 2022,
    has_ed     = year >= 2011 & year <= 2023,
    has_ems    = year >= 2019
  ) %>%
  filter(in_warm_season(event_date))

panel_mrt_primary <- panel_mrt_all %>%
  filter(in_primary_window(event_date))

# -------------------------------------------------------------------------------------
# Build citywide daily panels
# -------------------------------------------------------------------------------------

death_vars <- death_outcomes$outcome[death_outcomes$outcome %in% names(panel_mrt_all)]
ed_vars    <- ed_outcomes$outcome[ed_outcomes$outcome %in% names(panel_mrt_all)]
ems_vars   <- ems_outcomes$outcome[ems_outcomes$outcome %in% names(panel_mrt_all)]

city_deaths_primary <- panel_mrt_primary %>%
  filter(has_deaths) %>%
  build_city_daily(outcome_vars = death_vars)

city_ed_primary <- panel_mrt_primary %>%
  filter(has_ed) %>%
  build_city_daily(outcome_vars = ed_vars)

city_ems_primary <- panel_mrt_primary %>%
  filter(has_ems) %>%
  build_city_daily(outcome_vars = ems_vars)

city_deaths_full <- panel_mrt_all %>%
  filter(has_deaths) %>%
  build_city_daily(outcome_vars = death_vars)

city_ed_full <- panel_mrt_all %>%
  filter(has_ed) %>%
  build_city_daily(outcome_vars = ed_vars)

city_ems_full <- panel_mrt_all %>%
  filter(has_ems) %>%
  build_city_daily(outcome_vars = ems_vars)

# -------------------------------------------------------------------------------------
# MRT estimation
# -------------------------------------------------------------------------------------

estimate_mrt_dlnm <- function(data,
                              outcome,
                              outcome_label,
                              endpoint,
                              domain,
                              lag_input,
                              study_window_label,
                              temp_var = "tmax") {
  
  dat <- data %>%
    filter(
      !is.na(.data[[outcome]]),
      !is.na(.data[[temp_var]]),
      !is.na(humidity)
    ) %>%
    arrange(event_date) %>%
    mutate(
      time_num = as.numeric(event_date - min(event_date)) + 1,
      dow = factor(lubridate::wday(event_date, label = TRUE, abbr = TRUE)),
      month_factor = factor(lubridate::month(event_date)),
      year_factor = factor(lubridate::year(event_date))
    )
  
  if (nrow(dat) < 50) {
    stop(paste("Not enough data for outcome:", outcome))
  }
  
  is_ems <- str_detect(outcome, "^ems")
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
    floor(min(dat[[temp_var]], na.rm = TRUE)),
    ceiling(max(dat[[temp_var]], na.rm = TRUE)),
    by = 0.1
  )
  
  fit_nolag_fallback <- function(model_type = "no_lag_fallback") {
    form <- as.formula(
      paste0(
        outcome,
        " ~ ns(", temp_var, ", df = 3) + ns(humidity, 3) + ",
        "factor(dow) + ns(time_num, df = 4)"
      )
    )
    
    mod <- glm(form, family = quasipoisson(link = "log"), data = dat)
    
    pred_dat <- tibble(temp = pred_grid) %>%
      rename(!!temp_var := temp) %>%
      mutate(
        humidity = median(dat$humidity, na.rm = TRUE),
        dow = factor(levels(dat$dow)[1], levels = levels(dat$dow)),
        time_num = median(dat$time_num, na.rm = TRUE)
      )
    
    pred <- predict(mod, newdata = pred_dat, type = "link", se.fit = TRUE)
    
    fit <- as.numeric(pred$fit)
    se  <- as.numeric(pred$se.fit)
    
    mrt <- pred_grid[which.min(fit)]
    ref <- min(fit, na.rm = TRUE)
    
    curve_df <- tibble(
      outcome = outcome,
      outcome_label = outcome_label,
      endpoint = endpoint,
      domain = domain,
      temp_var = temp_var,
      temp = pred_grid,
      rr = exp(fit - ref),
      rr_low = exp((fit - 1.96 * se) - ref),
      rr_high = exp((fit + 1.96 * se) - ref),
      mrt = mrt,
      model_type = model_type,
      max_lag = 0L,
      lag_input = lag_input,
      study_year_start = min(year(dat$event_date), na.rm = TRUE),
      study_year_end = max(year(dat$event_date), na.rm = TRUE),
      study_window_label = study_window_label
    )
    
    mrt_row <- tibble(
      outcome = outcome,
      outcome_label = outcome_label,
      endpoint = endpoint,
      domain = domain,
      temp_var = temp_var,
      mrt = mrt,
      n_days = nrow(dat),
      mean_daily_count = mean(dat[[outcome]], na.rm = TRUE),
      model_type = model_type,
      max_lag = 0L,
      lag_input = lag_input,
      study_year_start = min(year(dat$event_date), na.rm = TRUE),
      study_year_end = max(year(dat$event_date), na.rm = TRUE),
      study_window_label = study_window_label
    )
    
    list(curve = curve_df, mrt = mrt_row, model = mod)
  }
  
  if (force_no_lag) {
    return(fit_nolag_fallback(model_type = "no_lag_forced"))
  }
  
  fit_dlnm_once <- function(simple = FALSE) {
    if (simple || is_sparse) {
      max_lag <- min(lag_input, 3L)
      cb <- crossbasis(
        dat[[temp_var]],
        lag = max_lag,
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
      max_lag <- min(lag_input, 3L)
      cb <- crossbasis(
        dat[[temp_var]],
        lag = max_lag,
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
      max_lag <- lag_input
      temp_knots <- quantile(dat[[temp_var]], probs = c(0.10, 0.50, 0.90), na.rm = TRUE)
      
      cb <- crossbasis(
        dat[[temp_var]],
        lag = max_lag,
        argvar = list(fun = "ns", knots = temp_knots),
        arglag = list(fun = "ns", knots = logknots(max_lag, 3))
      )
      
      form <- as.formula(
        paste0(
          outcome,
          " ~ cb + ns(humidity, 3) + factor(dow) + factor(month_factor) + ",
          "factor(year_factor) + ns(time_num, df = 6)"
        )
      )
      
      model_type <- "DLNM_full"
    }
    
    mod <- glm(form, family = quasipoisson(link = "log"), data = dat)
    
    coef_names <- names(coef(mod))
    if (!any(str_detect(coef_names, "cb"))) {
      stop("No crossbasis coefficients retained in fitted model.")
    }
    
    cp0 <- crosspred(
      cb,
      mod,
      at = pred_grid,
      cen = median(dat[[temp_var]], na.rm = TRUE)
    )
    
    rr_raw <- as.numeric(cp0$allRRfit)
    mrt <- pred_grid[which.min(rr_raw)]
    
    cp <- crosspred(
      cb,
      mod,
      at = pred_grid,
      cen = mrt
    )
    
    curve_df <- tibble(
      outcome = outcome,
      outcome_label = outcome_label,
      endpoint = endpoint,
      domain = domain,
      temp_var = temp_var,
      temp = pred_grid,
      rr = as.numeric(cp$allRRfit),
      rr_low = as.numeric(cp$allRRlow),
      rr_high = as.numeric(cp$allRRhigh),
      mrt = mrt,
      model_type = model_type,
      max_lag = max_lag,
      lag_input = lag_input,
      study_year_start = min(year(dat$event_date), na.rm = TRUE),
      study_year_end = max(year(dat$event_date), na.rm = TRUE),
      study_window_label = study_window_label
    )
    
    mrt_row <- tibble(
      outcome = outcome,
      outcome_label = outcome_label,
      endpoint = endpoint,
      domain = domain,
      temp_var = temp_var,
      mrt = mrt,
      n_days = nrow(dat),
      mean_daily_count = mean(dat[[outcome]], na.rm = TRUE),
      model_type = model_type,
      max_lag = max_lag,
      lag_input = lag_input,
      study_year_start = min(year(dat$event_date), na.rm = TRUE),
      study_year_end = max(year(dat$event_date), na.rm = TRUE),
      study_window_label = study_window_label
    )
    
    list(curve = curve_df, mrt = mrt_row, model = mod)
  }
  
  tryCatch(
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
}

run_outcome_set <- function(df, outcomes_tbl, study_window_label, temp_var = "tmax") {
  purrr::pmap(
    list(
      outcomes_tbl$outcome,
      outcomes_tbl$outcome_label,
      outcomes_tbl$endpoint,
      outcomes_tbl$domain,
      outcomes_tbl$lag_input
    ),
    function(outcome, outcome_label, endpoint, domain, lag_input) {
      tryCatch(
        estimate_mrt_dlnm(
          data = df,
          outcome = outcome,
          outcome_label = outcome_label,
          endpoint = endpoint,
          domain = domain,
          lag_input = lag_input,
          study_window_label = study_window_label,
          temp_var = temp_var
        ),
        error = function(e) {
          message("Failed for outcome: ", outcome, " | ", e$message)
          
          list(
            model = NULL,
            curve = tibble(
              outcome = character(),
              outcome_label = character(),
              endpoint = character(),
              domain = character(),
              temp_var = character(),
              temp = numeric(),
              rr = numeric(),
              rr_low = numeric(),
              rr_high = numeric(),
              mrt = numeric(),
              model_type = character(),
              max_lag = integer(),
              lag_input = integer(),
              study_year_start = integer(),
              study_year_end = integer(),
              study_window_label = character()
            ),
            mrt = tibble(
              outcome = outcome,
              outcome_label = outcome_label,
              endpoint = endpoint,
              domain = domain,
              temp_var = temp_var,
              mrt = NA_real_,
              n_days = nrow(df),
              mean_daily_count = NA_real_,
              model_type = "failed",
              max_lag = NA_integer_,
              lag_input = lag_input,
              study_year_start = suppressWarnings(min(year(df$event_date), na.rm = TRUE)),
              study_year_end = suppressWarnings(max(year(df$event_date), na.rm = TRUE)),
              study_window_label = study_window_label
            )
          )
        }
      )
    }
  )
}

# -------------------------------------------------------------------------------------
# Screen outcomes
# -------------------------------------------------------------------------------------

death_screen_primary <- bind_rows(lapply(death_vars, function(x) screen_outcome(city_deaths_primary, x, min_nonzero_days, min_total_events)))
ed_screen_primary    <- bind_rows(lapply(ed_vars,    function(x) screen_outcome(city_ed_primary, x, min_nonzero_days, min_total_events)))
ems_screen_primary   <- bind_rows(lapply(ems_vars,   function(x) screen_outcome(city_ems_primary, x, min_nonzero_days, min_total_events)))

death_screen_full <- bind_rows(lapply(death_vars, function(x) screen_outcome(city_deaths_full, x, min_nonzero_days, min_total_events)))
ed_screen_full    <- bind_rows(lapply(ed_vars,    function(x) screen_outcome(city_ed_full, x, min_nonzero_days, min_total_events)))
ems_screen_full   <- bind_rows(lapply(ems_vars,   function(x) screen_outcome(city_ems_full, x, min_nonzero_days, min_total_events)))

death_outcomes_keep_primary <- death_outcomes %>%
  filter(outcome %in% death_vars) %>%
  mutate(endpoint = "Mortality", domain = classify_domain(outcome)) %>%
  inner_join(death_screen_primary %>% filter(keep_for_dlnm), by = "outcome")

ed_outcomes_keep_primary <- ed_outcomes %>%
  filter(outcome %in% ed_vars) %>%
  mutate(endpoint = "ED", domain = classify_domain(outcome)) %>%
  inner_join(ed_screen_primary %>% filter(keep_for_dlnm), by = "outcome")

ems_outcomes_keep_primary <- ems_outcomes %>%
  filter(outcome %in% ems_vars) %>%
  mutate(endpoint = "EMS", domain = classify_domain(outcome)) %>%
  inner_join(ems_screen_primary %>% filter(keep_for_dlnm), by = "outcome")

death_outcomes_keep_full <- death_outcomes %>%
  filter(outcome %in% death_vars) %>%
  mutate(endpoint = "Mortality", domain = classify_domain(outcome)) %>%
  inner_join(death_screen_full %>% filter(keep_for_dlnm), by = "outcome")

ed_outcomes_keep_full <- ed_outcomes %>%
  filter(outcome %in% ed_vars) %>%
  mutate(endpoint = "ED", domain = classify_domain(outcome)) %>%
  inner_join(ed_screen_full %>% filter(keep_for_dlnm), by = "outcome")

ems_outcomes_keep_full <- ems_outcomes %>%
  filter(outcome %in% ems_vars) %>%
  mutate(endpoint = "EMS", domain = classify_domain(outcome)) %>%
  inner_join(ems_screen_full %>% filter(keep_for_dlnm), by = "outcome")

fwrite(death_screen_primary, "results/mrt/death_screen_primary.csv")
fwrite(ed_screen_primary,    "results/mrt/ed_screen_primary.csv")
fwrite(ems_screen_primary,   "results/mrt/ems_screen_primary.csv")

fwrite(death_screen_full, "results/mrt/death_screen_full_history.csv")
fwrite(ed_screen_full,    "results/mrt/ed_screen_full_history.csv")
fwrite(ems_screen_full,   "results/mrt/ems_screen_full_history.csv")

# -------------------------------------------------------------------------------------
# Run primary analysis
# -------------------------------------------------------------------------------------

death_results_primary <- run_outcome_set(
  city_deaths_primary,
  death_outcomes_keep_primary,
  study_window_label = primary_label,
  temp_var = "tmax"
)

ed_results_primary <- run_outcome_set(
  city_ed_primary,
  ed_outcomes_keep_primary,
  study_window_label = primary_label,
  temp_var = "tmax"
)

ems_results_primary <- run_outcome_set(
  city_ems_primary,
  ems_outcomes_keep_primary,
  study_window_label = primary_label,
  temp_var = "tmax"
)

all_results_primary <- c(
  death_results_primary,
  ed_results_primary,
  ems_results_primary
)

mrt_table_expanded <- bind_rows(lapply(all_results_primary, function(x) x$mrt)) %>%
  arrange(endpoint, domain, outcome_label)

curve_table_expanded <- bind_rows(lapply(all_results_primary, function(x) x$curve)) %>%
  arrange(endpoint, domain, outcome_label, temp)

fwrite(mrt_table_expanded, "results/mrt/mrt_table_expanded_2019_2022.csv")
fwrite(curve_table_expanded, "results/mrt/mrt_response_curves_expanded_2019_2022.csv")

saveRDS(mrt_table_expanded, "results/mrt/mrt_table_expanded_2019_2022.rds")
saveRDS(curve_table_expanded, "results/mrt/mrt_response_curves_expanded_2019_2022.rds")

# -------------------------------------------------------------------------------------
# Run sensitivity analysis
# -------------------------------------------------------------------------------------

death_results_full <- run_outcome_set(
  city_deaths_full,
  death_outcomes_keep_full,
  study_window_label = sensitivity_label,
  temp_var = "tmax"
)

ed_results_full <- run_outcome_set(
  city_ed_full,
  ed_outcomes_keep_full,
  study_window_label = sensitivity_label,
  temp_var = "tmax"
)

ems_results_full <- run_outcome_set(
  city_ems_full,
  ems_outcomes_keep_full,
  study_window_label = sensitivity_label,
  temp_var = "tmax"
)

all_results_full <- c(
  death_results_full,
  ed_results_full,
  ems_results_full
)

mrt_table_expanded_full_history <- bind_rows(lapply(all_results_full, function(x) x$mrt)) %>%
  arrange(endpoint, domain, outcome_label)

curve_table_expanded_full_history <- bind_rows(lapply(all_results_full, function(x) x$curve)) %>%
  arrange(endpoint, domain, outcome_label, temp)

fwrite(mrt_table_expanded_full_history, "results/mrt/mrt_table_expanded_full_history.csv")
fwrite(curve_table_expanded_full_history, "results/mrt/mrt_response_curves_expanded_full_history.csv")

saveRDS(mrt_table_expanded_full_history, "results/mrt/mrt_table_expanded_full_history.rds")
saveRDS(curve_table_expanded_full_history, "results/mrt/mrt_response_curves_expanded_full_history.rds")

# -------------------------------------------------------------------------------------
# Optional primary plots
# -------------------------------------------------------------------------------------

for (i in seq_along(all_results_primary)) {
  curve_df <- all_results_primary[[i]]$curve
  if (nrow(curve_df) == 0) next
  
  lab <- unique(curve_df$outcome_label)[1]
  stub <- unique(curve_df$outcome)[1]
  
  p <- plot_one_curve(
    curve_df,
    title_text = lab,
    subtitle_text = "Primary analysis: May-September, 2019-2022"
  )
  
  if (!is.null(p)) {
    save_mrt_plot(p, paste0(stub, "_mrt_curve_2019_2022.png"))
  }
}

# -------------------------------------------------------------------------------------
# Console output
# -------------------------------------------------------------------------------------

cat("\nPrimary MRT modeling complete.\n")
cat("Saved: results/mrt/mrt_table_expanded_2019_2022.csv\n")
cat("Saved: results/mrt/mrt_response_curves_expanded_2019_2022.csv\n")

cat("\nSensitivity MRT modeling complete.\n")
cat("Saved: results/mrt/mrt_table_expanded_full_history.csv\n")
cat("Saved: results/mrt/mrt_response_curves_expanded_full_history.csv\n")

cat("\nPrimary MRT summary:\n")
print(
  mrt_table_expanded %>%
    select(endpoint, domain, outcome_label, mrt, max_lag, lag_input, mean_daily_count, model_type) %>%
    arrange(endpoint, domain, outcome_label))
  
  

# =====================================================================================
# CREATE HEAT INDEX IN panel_causes
# Assumptions:
#   - tmax is in degrees C
#   - humidity is relative humidity in percent (0-100)
# Output:
#   - panel_causes$heat_index in degrees C
# =====================================================================================

# -------------------------------------------------------------------------------------
# Helper: calculate NOAA/NWS heat index from temperature (C) and RH (%)
# Returns heat index in C
# -------------------------------------------------------------------------------------

compute_heat_index_c <- function(temp_c, rh) {
  # initialize output
  hi_c <- rep(NA_real_, length(temp_c))
  
  # valid observations
  ok <- !is.na(temp_c) & !is.na(rh)
  
  if (!any(ok)) {
    return(hi_c)
  }
  
  # clamp RH to plausible bounds
  rh_use <- pmin(pmax(rh[ok], 0), 100)
  
  # convert C -> F
  temp_f <- temp_c[ok] * 9 / 5 + 32
  
  # simple approximation used by NOAA as a screening step
  hi_simple_f <- 0.5 * (
    temp_f + 61.0 + ((temp_f - 68.0) * 1.2) + (rh_use * 0.094)
  )
  
  # full Rothfusz regression
  hi_full_f <- -42.379 +
    2.04901523 * temp_f +
    10.14333127 * rh_use -
    0.22475541 * temp_f * rh_use -
    0.00683783 * temp_f^2 -
    0.05481717 * rh_use^2 +
    0.00122874 * temp_f^2 * rh_use +
    0.00085282 * temp_f * rh_use^2 -
    0.00000199 * temp_f^2 * rh_use^2
  
  # use full HI only in standard regime; otherwise use air temperature
  use_full <- temp_f >= 80 & rh_use >= 40
  
  hi_f <- ifelse(use_full, hi_full_f, temp_f)
  
  # low-RH adjustment
  low_rh_adj <- ((13 - rh_use) / 4) * sqrt((17 - abs(temp_f - 95)) / 17)
  low_rh_idx <- use_full & rh_use < 13 & temp_f >= 80 & temp_f <= 112
  hi_f[low_rh_idx] <- hi_f[low_rh_idx] - low_rh_adj[low_rh_idx]
  
  # high-RH adjustment
  high_rh_adj <- ((rh_use - 85) / 10) * ((87 - temp_f) / 5)
  high_rh_idx <- use_full & rh_use > 85 & temp_f >= 80 & temp_f <= 87
  hi_f[high_rh_idx] <- hi_f[high_rh_idx] + high_rh_adj[high_rh_idx]
  
  # convert back to C
  hi_c[ok] <- (hi_f - 32) * 5 / 9
  
  hi_c
}

# -------------------------------------------------------------------------------------
# Create heat index in panel_causes
# -------------------------------------------------------------------------------------

panel_causes <- panel_causes %>%
  mutate(
    heat_index = compute_heat_index_c(
      temp_c = tmax,
      rh = humidity
    )
  )

# -------------------------------------------------------------------------------------
# Quick QC
# -------------------------------------------------------------------------------------

cat("\nHeat index created in panel_causes.\n")
cat("Non-missing heat_index rows:", sum(!is.na(panel_causes$heat_index)), "\n")
cat("Missing heat_index rows:", sum(is.na(panel_causes$heat_index)), "\n")

print(
  panel_causes %>%
    summarise(
      tmax_min = min(tmax, na.rm = TRUE),
      tmax_median = median(tmax, na.rm = TRUE),
      tmax_max = max(tmax, na.rm = TRUE),
      rh_min = min(humidity, na.rm = TRUE),
      rh_median = median(humidity, na.rm = TRUE),
      rh_max = max(humidity, na.rm = TRUE),
      hi_min = min(heat_index, na.rm = TRUE),
      hi_median = median(heat_index, na.rm = TRUE),
      hi_max = max(heat_index, na.rm = TRUE)
    )
)

# optional spot check
print(
  panel_causes %>%
    select(event_date, tmax, humidity, heat_index) %>%
    arrange(desc(tmax)) %>%
    head(10)
)


  # =====================================================================================
  # HEAT INDEX MRT MODELING BRANCH
  # Parallel exposure analysis using heat_index instead of tmax
  # Primary analysis: May-September, 2019-2022
  # Sensitivity analysis: full warm-season history by source availability
  # =====================================================================================
  
  # -------------------------------------------------------------------------------------
  # Sanity check
  # -------------------------------------------------------------------------------------
  
  if (!"heat_index" %in% names(panel_causes)) {
    stop(
      "panel_causes does not contain a 'heat_index' column. ",
      "Create heat_index upstream before running this block."
    )
  }
  
  # -------------------------------------------------------------------------------------
  # Build citywide daily panels with heat index
  # -------------------------------------------------------------------------------------
  
  build_city_daily_hi <- function(data, outcome_vars) {
    data %>%
      group_by(event_date, year, month) %>%
      summarise(
        heat_index = mean(heat_index, na.rm = TRUE),
        tmax = mean(tmax, na.rm = TRUE),
        tmin = mean(tmin, na.rm = TRUE),
        tmean = mean(tmean, na.rm = TRUE),
        humidity = mean(humidity, na.rm = TRUE),
        across(all_of(outcome_vars), ~ sum(.x, na.rm = TRUE)),
        .groups = "drop"
      ) %>%
      mutate(
        dow = factor(
          wday(event_date),
          levels = levels(
            wday(
              seq.Date(as.Date("2023-01-01"), by = "day", length.out = 7)
            )
          )
        )
      ) %>%
      arrange(event_date)
  }

# -------------------------------------------------------------------------------------
# Rebuild MRT panels after adding heat_index to panel_causes
# -------------------------------------------------------------------------------------

panel_mrt_all <- panel_causes %>%
  mutate(
    event_date = as.Date(event_date),
    year = lubridate::year(event_date),
    month = lubridate::month(event_date),
    has_deaths = year >= 1993 & year <= 2022,
    has_ed     = year >= 2011 & year <= 2023,
    has_ems    = year >= 2019
  ) %>%
  filter(in_warm_season(event_date))

panel_mrt_primary <- panel_mrt_all %>%
  filter(in_primary_window(event_date))

# quick check
stopifnot("heat_index" %in% names(panel_mrt_all))
stopifnot("heat_index" %in% names(panel_mrt_primary))
  
  city_deaths_primary_hi <- panel_mrt_primary %>%
    filter(has_deaths) %>%
    build_city_daily_hi(outcome_vars = death_vars)
  
  city_ed_primary_hi <- panel_mrt_primary %>%
    filter(has_ed) %>%
    build_city_daily_hi(outcome_vars = ed_vars)
  
  city_ems_primary_hi <- panel_mrt_primary %>%
    filter(has_ems) %>%
    build_city_daily_hi(outcome_vars = ems_vars)
  
  city_deaths_full_hi <- panel_mrt_all %>%
    filter(has_deaths) %>%
    build_city_daily_hi(outcome_vars = death_vars)
  
  city_ed_full_hi <- panel_mrt_all %>%
    filter(has_ed) %>%
    build_city_daily_hi(outcome_vars = ed_vars)
  
  city_ems_full_hi <- panel_mrt_all %>%
    filter(has_ems) %>%
    build_city_daily_hi(outcome_vars = ems_vars)
  
  # -------------------------------------------------------------------------------------
  # Heat-index-specific MRT estimation function
  # Note:
  # - exposure = heat_index
  # - no separate humidity adjustment because HI already incorporates humidity
  # -------------------------------------------------------------------------------------
  
  estimate_mrt_dlnm_hi <- function(data,
                                   outcome,
                                   outcome_label,
                                   endpoint,
                                   domain,
                                   lag_input,
                                   study_window_label,
                                   exposure_var = "heat_index") {
    
    dat <- data %>%
      filter(
        !is.na(.data[[outcome]]),
        !is.na(.data[[exposure_var]])
      ) %>%
      arrange(event_date) %>%
      mutate(
        time_num = as.numeric(event_date - min(event_date)) + 1,
        dow = factor(lubridate::wday(event_date)),
        month_factor = factor(lubridate::month(event_date)),
        year_factor = factor(lubridate::year(event_date))
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
      floor(min(dat[[exposure_var]], na.rm = TRUE)),
      ceiling(max(dat[[exposure_var]], na.rm = TRUE)),
      by = 0.1
    )
    
    fit_nolag_fallback <- function(model_type = "no_lag_fallback") {
      
      form <- as.formula(
        paste0(
          outcome,
          " ~ ns(", exposure_var, ", df = 3) + factor(dow) + ns(time_num, df = 4)"
        )
      )
      
      mod <- glm(form, family = quasipoisson(link = "log"), data = dat)
      
      pred_dat <- tibble::tibble(exposure_value = pred_grid) %>%
        rename(!!exposure_var := exposure_value) %>%
        mutate(
          dow = factor(levels(dat$dow)[1], levels = levels(dat$dow)),
          time_num = median(dat$time_num, na.rm = TRUE)
        )
      
      pred <- predict(mod, newdata = pred_dat, type = "link", se.fit = TRUE)
      
      fit <- as.numeric(pred$fit)
      se  <- as.numeric(pred$se.fit)
      
      mrt <- pred_grid[which.min(fit)]
      ref <- min(fit, na.rm = TRUE)
      
      curve_df <- tibble(
        outcome = outcome,
        outcome_label = outcome_label,
        endpoint = endpoint,
        domain = domain,
        temp_var = exposure_var,
        temp = pred_grid,
        rr = exp(fit - ref),
        rr_low = exp((fit - 1.96 * se) - ref),
        rr_high = exp((fit + 1.96 * se) - ref),
        mrt = mrt,
        model_type = model_type,
        max_lag = 0L,
        lag_input = lag_input,
        study_year_start = min(year(dat$event_date), na.rm = TRUE),
        study_year_end = max(year(dat$event_date), na.rm = TRUE),
        study_window_label = study_window_label
      )
      
      mrt_row <- tibble(
        outcome = outcome,
        outcome_label = outcome_label,
        endpoint = endpoint,
        domain = domain,
        temp_var = exposure_var,
        mrt = mrt,
        n_days = nrow(dat),
        mean_daily_count = mean(dat[[outcome]], na.rm = TRUE),
        model_type = model_type,
        max_lag = 0L,
        lag_input = lag_input,
        study_year_start = min(year(dat$event_date), na.rm = TRUE),
        study_year_end = max(year(dat$event_date), na.rm = TRUE),
        study_window_label = study_window_label
      )
      
      list(curve = curve_df, mrt = mrt_row, model = mod)
    }
    
    if (force_no_lag) {
      return(fit_nolag_fallback(model_type = "no_lag_forced"))
    }
    
    fit_dlnm_once <- function(simple = FALSE) {
      
      if (simple || is_sparse) {
        max_lag <- min(lag_input, 3L)
        
        cb <- dlnm::crossbasis(
          dat[[exposure_var]],
          lag = max_lag,
          argvar = list(fun = "ns", df = 2),
          arglag = list(fun = "ns", df = 2)
        )
        
        form <- as.formula(
          paste0(
            outcome,
            " ~ cb + factor(dow) + ns(time_num, df = 4)"
          )
        )
        
        model_type <- "DLNM_simple"
        
      } else if (is_ems) {
        max_lag <- min(lag_input, 3L)
        
        cb <- dlnm::crossbasis(
          dat[[exposure_var]],
          lag = max_lag,
          argvar = list(fun = "ns", df = 3),
          arglag = list(fun = "ns", df = 2)
        )
        
        form <- as.formula(
          paste0(
            outcome,
            " ~ cb + factor(dow) + ns(time_num, df = 4)"
          )
        )
        
        model_type <- "DLNM_ems"
        
      } else {
        max_lag <- lag_input
        
        exposure_knots <- quantile(
          dat[[exposure_var]],
          probs = c(0.10, 0.50, 0.90),
          na.rm = TRUE
        )
        
        cb <- dlnm::crossbasis(
          dat[[exposure_var]],
          lag = max_lag,
          argvar = list(fun = "ns", knots = exposure_knots),
          arglag = list(fun = "ns", knots = dlnm::logknots(max_lag, 3))
        )
        
        form <- as.formula(
          paste0(
            outcome,
            " ~ cb + factor(dow) + factor(month_factor) + ",
            "factor(year_factor) + ns(time_num, df = 6)"
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
        cen = median(dat[[exposure_var]], na.rm = TRUE)
      )
      
      rr_raw <- as.numeric(cp0$allRRfit)
      mrt <- pred_grid[which.min(rr_raw)]
      
      cp <- dlnm::crosspred(
        cb,
        mod,
        at = pred_grid,
        cen = mrt
      )
      
      curve_df <- tibble(
        outcome = outcome,
        outcome_label = outcome_label,
        endpoint = endpoint,
        domain = domain,
        temp_var = exposure_var,
        temp = pred_grid,
        rr = as.numeric(cp$allRRfit),
        rr_low = as.numeric(cp$allRRlow),
        rr_high = as.numeric(cp$allRRhigh),
        mrt = mrt,
        model_type = model_type,
        max_lag = max_lag,
        lag_input = lag_input,
        study_year_start = min(year(dat$event_date), na.rm = TRUE),
        study_year_end = max(year(dat$event_date), na.rm = TRUE),
        study_window_label = study_window_label
      )
      
      mrt_row <- tibble(
        outcome = outcome,
        outcome_label = outcome_label,
        endpoint = endpoint,
        domain = domain,
        temp_var = exposure_var,
        mrt = mrt,
        n_days = nrow(dat),
        mean_daily_count = mean(dat[[outcome]], na.rm = TRUE),
        model_type = model_type,
        max_lag = max_lag,
        lag_input = lag_input,
        study_year_start = min(year(dat$event_date), na.rm = TRUE),
        study_year_end = max(year(dat$event_date), na.rm = TRUE),
        study_window_label = study_window_label
      )
      
      list(curve = curve_df, mrt = mrt_row, model = mod)
    }
    
    tryCatch(
      fit_dlnm_once(simple = FALSE),
      error = function(e1) {
        message("Primary HI DLNM failed for ", outcome, ": ", e1$message)
        
        tryCatch(
          fit_dlnm_once(simple = TRUE),
          error = function(e2) {
            message("Simplified HI DLNM failed for ", outcome, ": ", e2$message)
            fit_nolag_fallback(model_type = "no_lag_fallback_after_dlnm_failure")
          }
        )
      }
    )
  }
  
  run_outcome_set_hi <- function(df, outcomes_tbl, study_window_label, exposure_var = "heat_index") {
    purrr::pmap(
      list(
        outcomes_tbl$outcome,
        outcomes_tbl$outcome_label,
        outcomes_tbl$endpoint,
        outcomes_tbl$domain,
        outcomes_tbl$lag_input
      ),
      function(outcome, outcome_label, endpoint, domain, lag_input) {
        tryCatch(
          estimate_mrt_dlnm_hi(
            data = df,
            outcome = outcome,
            outcome_label = outcome_label,
            endpoint = endpoint,
            domain = domain,
            lag_input = lag_input,
            study_window_label = study_window_label,
            exposure_var = exposure_var
          ),
          error = function(e) {
            message("Failed for HI outcome: ", outcome, " | ", e$message)
            
            list(
              model = NULL,
              curve = tibble(
                outcome = character(),
                outcome_label = character(),
                endpoint = character(),
                domain = character(),
                temp_var = character(),
                temp = numeric(),
                rr = numeric(),
                rr_low = numeric(),
                rr_high = numeric(),
                mrt = numeric(),
                model_type = character(),
                max_lag = integer(),
                lag_input = integer(),
                study_year_start = integer(),
                study_year_end = integer(),
                study_window_label = character()
              ),
              mrt = tibble(
                outcome = outcome,
                outcome_label = outcome_label,
                endpoint = endpoint,
                domain = domain,
                temp_var = exposure_var,
                mrt = NA_real_,
                n_days = nrow(df),
                mean_daily_count = NA_real_,
                model_type = "failed",
                max_lag = NA_integer_,
                lag_input = lag_input,
                study_year_start = suppressWarnings(min(year(df$event_date), na.rm = TRUE)),
                study_year_end = suppressWarnings(max(year(df$event_date), na.rm = TRUE)),
                study_window_label = study_window_label
              )
            )
          }
        )
      }
    )
  }
  
  # -------------------------------------------------------------------------------------
  # Screen heat-index panels
  # -------------------------------------------------------------------------------------
  
  death_screen_primary_hi <- bind_rows(lapply(death_vars, function(x) screen_outcome(
    city_deaths_primary_hi, x, min_nonzero_days, min_total_events
  )))
  ed_screen_primary_hi <- bind_rows(lapply(ed_vars, function(x) screen_outcome(
    city_ed_primary_hi, x, min_nonzero_days, min_total_events
  )))
  ems_screen_primary_hi <- bind_rows(lapply(ems_vars, function(x) screen_outcome(
    city_ems_primary_hi, x, min_nonzero_days, min_total_events
  )))
  
  death_screen_full_hi <- bind_rows(lapply(death_vars, function(x) screen_outcome(
    city_deaths_full_hi, x, min_nonzero_days, min_total_events
  )))
  ed_screen_full_hi <- bind_rows(lapply(ed_vars, function(x) screen_outcome(
    city_ed_full_hi, x, min_nonzero_days, min_total_events
  )))
  ems_screen_full_hi <- bind_rows(lapply(ems_vars, function(x) screen_outcome(
    city_ems_full_hi, x, min_nonzero_days, min_total_events
  )))
  
  death_outcomes_keep_primary_hi <- death_outcomes %>%
    filter(outcome %in% death_vars) %>%
    mutate(endpoint = "Mortality", domain = classify_domain(outcome)) %>%
    inner_join(death_screen_primary_hi %>% filter(keep_for_dlnm), by = "outcome")
  
  ed_outcomes_keep_primary_hi <- ed_outcomes %>%
    filter(outcome %in% ed_vars) %>%
    mutate(endpoint = "ED", domain = classify_domain(outcome)) %>%
    inner_join(ed_screen_primary_hi %>% filter(keep_for_dlnm), by = "outcome")
  
  ems_outcomes_keep_primary_hi <- ems_outcomes %>%
    filter(outcome %in% ems_vars) %>%
    mutate(endpoint = "EMS", domain = classify_domain(outcome)) %>%
    inner_join(ems_screen_primary_hi %>% filter(keep_for_dlnm), by = "outcome")
  
  death_outcomes_keep_full_hi <- death_outcomes %>%
    filter(outcome %in% death_vars) %>%
    mutate(endpoint = "Mortality", domain = classify_domain(outcome)) %>%
    inner_join(death_screen_full_hi %>% filter(keep_for_dlnm), by = "outcome")
  
  ed_outcomes_keep_full_hi <- ed_outcomes %>%
    filter(outcome %in% ed_vars) %>%
    mutate(endpoint = "ED", domain = classify_domain(outcome)) %>%
    inner_join(ed_screen_full_hi %>% filter(keep_for_dlnm), by = "outcome")
  
  ems_outcomes_keep_full_hi <- ems_outcomes %>%
    filter(outcome %in% ems_vars) %>%
    mutate(endpoint = "EMS", domain = classify_domain(outcome)) %>%
    inner_join(ems_screen_full_hi %>% filter(keep_for_dlnm), by = "outcome")
  
  fwrite(death_screen_primary_hi, "results/mrt/death_screen_primary_heat_index.csv")
  fwrite(ed_screen_primary_hi,    "results/mrt/ed_screen_primary_heat_index.csv")
  fwrite(ems_screen_primary_hi,   "results/mrt/ems_screen_primary_heat_index.csv")
  
  fwrite(death_screen_full_hi, "results/mrt/death_screen_full_history_heat_index.csv")
  fwrite(ed_screen_full_hi,    "results/mrt/ed_screen_full_history_heat_index.csv")
  fwrite(ems_screen_full_hi,   "results/mrt/ems_screen_full_history_heat_index.csv")
  
  # -------------------------------------------------------------------------------------
  # Run heat-index primary analysis
  # -------------------------------------------------------------------------------------
  
  death_results_primary_hi <- run_outcome_set_hi(
    city_deaths_primary_hi,
    death_outcomes_keep_primary_hi,
    study_window_label = paste0(primary_label, "_heat_index"),
    exposure_var = "heat_index"
  )
  
  ed_results_primary_hi <- run_outcome_set_hi(
    city_ed_primary_hi,
    ed_outcomes_keep_primary_hi,
    study_window_label = paste0(primary_label, "_heat_index"),
    exposure_var = "heat_index"
  )
  
  ems_results_primary_hi <- run_outcome_set_hi(
    city_ems_primary_hi,
    ems_outcomes_keep_primary_hi,
    study_window_label = paste0(primary_label, "_heat_index"),
    exposure_var = "heat_index"
  )
  
  all_results_primary_hi <- c(
    death_results_primary_hi,
    ed_results_primary_hi,
    ems_results_primary_hi
  )
  
  mrt_table_expanded_heat_index <- bind_rows(lapply(all_results_primary_hi, function(x) x$mrt)) %>%
    arrange(endpoint, domain, outcome_label)
  
  curve_table_expanded_heat_index <- bind_rows(lapply(all_results_primary_hi, function(x) x$curve)) %>%
    arrange(endpoint, domain, outcome_label, temp)
  
  fwrite(
    mrt_table_expanded_heat_index,
    "results/mrt/mrt_table_expanded_heat_index_2019_2022.csv"
  )
  fwrite(
    curve_table_expanded_heat_index,
    "results/mrt/mrt_response_curves_expanded_heat_index_2019_2022.csv"
  )
  
  saveRDS(
    mrt_table_expanded_heat_index,
    "results/mrt/mrt_table_expanded_heat_index_2019_2022.rds"
  )
  saveRDS(
    curve_table_expanded_heat_index,
    "results/mrt/mrt_response_curves_expanded_heat_index_2019_2022.rds"
  )
  
  # -------------------------------------------------------------------------------------
  # Run heat-index full-history sensitivity analysis
  # -------------------------------------------------------------------------------------
  
  death_results_full_hi <- run_outcome_set_hi(
    city_deaths_full_hi,
    death_outcomes_keep_full_hi,
    study_window_label = paste0(sensitivity_label, "_heat_index"),
    exposure_var = "heat_index"
  )
  
  ed_results_full_hi <- run_outcome_set_hi(
    city_ed_full_hi,
    ed_outcomes_keep_full_hi,
    study_window_label = paste0(sensitivity_label, "_heat_index"),
    exposure_var = "heat_index"
  )
  
  ems_results_full_hi <- run_outcome_set_hi(
    city_ems_full_hi,
    ems_outcomes_keep_full_hi,
    study_window_label = paste0(sensitivity_label, "_heat_index"),
    exposure_var = "heat_index"
  )
  
  all_results_full_hi <- c(
    death_results_full_hi,
    ed_results_full_hi,
    ems_results_full_hi
  )
  
  mrt_table_expanded_heat_index_full_history <- bind_rows(lapply(all_results_full_hi, function(x) x$mrt)) %>%
    arrange(endpoint, domain, outcome_label)
  
  curve_table_expanded_heat_index_full_history <- bind_rows(lapply(all_results_full_hi, function(x) x$curve)) %>%
    arrange(endpoint, domain, outcome_label, temp)
  
  fwrite(
    mrt_table_expanded_heat_index_full_history,
    "results/mrt/mrt_table_expanded_heat_index_full_history.csv"
  )
  fwrite(
    curve_table_expanded_heat_index_full_history,
    "results/mrt/mrt_response_curves_expanded_heat_index_full_history.csv"
  )
  
  saveRDS(
    mrt_table_expanded_heat_index_full_history,
    "results/mrt/mrt_table_expanded_heat_index_full_history.rds"
  )
  saveRDS(
    curve_table_expanded_heat_index_full_history,
    "results/mrt/mrt_response_curves_expanded_heat_index_full_history.rds"
  )
  
  # -------------------------------------------------------------------------------------
  # Optional heat-index plots for primary analysis
  # -------------------------------------------------------------------------------------
  
  for (i in seq_along(all_results_primary_hi)) {
    curve_df <- all_results_primary_hi[[i]]$curve
    if (nrow(curve_df) == 0) next
    
    lab <- unique(curve_df$outcome_label)[1]
    stub <- unique(curve_df$outcome)[1]
    
    p <- plot_one_curve(
      curve_df,
      title_text = paste0(lab, " (Heat Index)"),
      subtitle_text = "Primary analysis: May-September, 2019-2022"
    )
    
    if (!is.null(p)) {
      save_mrt_plot(p, paste0(stub, "_heat_index_mrt_curve_2019_2022.png"))
    }
  }
  
  # -------------------------------------------------------------------------------------
  # Console output
  # -------------------------------------------------------------------------------------
  
  cat("\nHeat index MRT modeling complete.\n")
  cat("Saved: results/mrt/mrt_table_expanded_heat_index_2019_2022.csv\n")
  cat("Saved: results/mrt/mrt_response_curves_expanded_heat_index_2019_2022.csv\n")
  cat("Saved: results/mrt/mrt_table_expanded_heat_index_full_history.csv\n")
  cat("Saved: results/mrt/mrt_response_curves_expanded_heat_index_full_history.csv\n")
  
  cat("\nPrimary heat-index MRT summary:\n")
  print(
    mrt_table_expanded_heat_index %>%
      select(endpoint, domain, outcome_label, mrt, max_lag, lag_input, mean_daily_count, model_type) %>%
      arrange(endpoint, domain, outcome_label)
  )
  
  
  
  
  # =====================================================================================
  # FACETED DLNM PLOTS (TMAX + HEAT INDEX)
  # =====================================================================================
  
  # =====================================================================================
  # FACETED DLNM PLOTS BY OUTCOME TYPE AND EXPOSURE
  # Produces:
  #   - mortality_tmax_faceted
  #   - mortality_heat_index_faceted
  #   - ed_tmax_faceted
  #   - ed_heat_index_faceted
  #   - ems_tmax_faceted
  #   - ems_heat_index_faceted
  # =====================================================================================
  
  plot_endpoint_facets <- function(curve_df,
                                   endpoint_name,
                                   exposure_var,
                                   file_stub,
                                   subtitle_text = NULL) {
    
    df_plot <- curve_df %>%
      filter(endpoint == endpoint_name) %>%
      mutate(outcome_label = factor(outcome_label))
    
    if (nrow(df_plot) == 0) {
      message("No rows found for endpoint = ", endpoint_name,
              " and exposure = ", exposure_var)
      return(NULL)
    }
    
    # -----------------------------
    # Color mapping
    # -----------------------------
    if (endpoint_name == "Mortality") {
      line_col <- ifelse(exposure_var == "heat_index", "navy", "steelblue4")
      fill_col <- ifelse(exposure_var == "heat_index", "navy", "steelblue")
    } else if (endpoint_name == "ED") {
      line_col <- ifelse(exposure_var == "heat_index", "orangered4", "darkorange4")
      fill_col <- ifelse(exposure_var == "heat_index", "orangered", "darkorange")
    } else if (endpoint_name == "EMS") {
      line_col <- ifelse(exposure_var == "heat_index", "darkgreen", "darkgreen")
      fill_col <- ifelse(exposure_var == "heat_index", "forestgreen", "darkgreen")
    } else {
      line_col <- "gray30"
      fill_col <- "gray70"
    }
    
    # -----------------------------
    # Axis + title
    # -----------------------------
    x_label <- ifelse(
      exposure_var == "heat_index",
      "Heat Index (°C)",
      "Daily Maximum Temperature (°C)"
    )
    
    exposure_title <- ifelse(
      exposure_var == "heat_index",
      "Heat Index",
      "Temperature"
    )
    
    # -----------------------------
    # Plot
    # -----------------------------
    p <- ggplot(df_plot, aes(x = temp, y = rr)) +
      geom_ribbon(
        aes(ymin = rr_low, ymax = rr_high),
        alpha = 0.18,
        fill = fill_col
      ) +
      geom_line(
        color = line_col,
        linewidth = 1
      ) +
      geom_hline(yintercept = 1, linetype = "dashed", color = "gray40") +
      geom_vline(aes(xintercept = mrt), linetype = "dotted", color = "firebrick") +
      
      facet_wrap(~ outcome_label, scales = "free_y") +
      
      labs(
        title = paste0(endpoint_name, " ", exposure_title, " Response Curves"),
        subtitle = subtitle_text,
        x = x_label,
        y = "Relative Risk"
      ) +
      
      theme_minimal(base_size = 11) +
      theme(
        plot.title = element_text(face = "bold"),
        strip.text = element_text(size = 10),
        panel.spacing = unit(0.8, "lines")
      )
    
    save_mrt_plot(p, paste0(file_stub, ".png"), width = 14, height = 10)
    save_mrt_plot(p, paste0(file_stub, ".pdf"), width = 14, height = 10)
    
    return(p)
  }
  
  # -------------------------------------------------------------------------------------
  # TMAX FACETED PLOTS
  # Uses curve_table_expanded
  # -------------------------------------------------------------------------------------
  
  plot_mortality_tmax_faceted <- plot_endpoint_facets(
    curve_df = curve_table_expanded,
    endpoint_name = "Mortality",
    exposure_var = "tmax",
    file_stub = "mortality_tmax_faceted",
    subtitle_text = "DLNM-based cumulative relative risk"
  )
  
  plot_ed_tmax_faceted <- plot_endpoint_facets(
    curve_df = curve_table_expanded,
    endpoint_name = "ED",
    exposure_var = "tmax",
    file_stub = "ed_tmax_faceted",
    subtitle_text = "DLNM-based cumulative relative risk"
  )
  
  plot_ems_tmax_faceted <- plot_endpoint_facets(
    curve_df = curve_table_expanded,
    endpoint_name = "EMS",
    exposure_var = "tmax",
    file_stub = "ems_tmax_faceted",
    subtitle_text = "DLNM-based cumulative relative risk"
  )
  
  # -------------------------------------------------------------------------------------
  # HEAT INDEX FACETED PLOTS
  # Uses curve_table_expanded_heat_index
  # -------------------------------------------------------------------------------------
  
  plot_mortality_hi_faceted <- plot_endpoint_facets(
    curve_df = curve_table_expanded_heat_index,
    endpoint_name = "Mortality",
    exposure_var = "heat_index",
    file_stub = "mortality_heat_index_faceted",
    subtitle_text = "DLNM-based cumulative relative risk"
  )
  
  plot_ed_hi_faceted <- plot_endpoint_facets(
    curve_df = curve_table_expanded_heat_index,
    endpoint_name = "ED",
    exposure_var = "heat_index",
    file_stub = "ed_heat_index_faceted",
    subtitle_text = "DLNM-based cumulative relative risk"
  )
  
  plot_ems_hi_faceted <- plot_endpoint_facets(
    curve_df = curve_table_expanded_heat_index,
    endpoint_name = "EMS",
    exposure_var = "heat_index",
    file_stub = "ems_heat_index_faceted",
    subtitle_text = "DLNM-based cumulative relative risk"
  )
  
  # -------------------------------------------------------------------------------------
  # Console output
  # -------------------------------------------------------------------------------------
  
  cat("\nSaved faceted endpoint-specific DLNM plots:\n")
  cat(" - mortality_tmax_faceted.png/.pdf\n")
  cat(" - mortality_heat_index_faceted.png/.pdf\n")
  cat(" - ed_tmax_faceted.png/.pdf\n")
  cat(" - ed_heat_index_faceted.png/.pdf\n")
  cat(" - ems_tmax_faceted.png/.pdf\n")
  cat(" - ems_heat_index_faceted.png/.pdf\n")
