# =====================================================================================
# HVI 2.0 | Data Audit Script
# Purpose: Inspect structure, completeness, and harmonization potential of datasets
# =====================================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(lubridate)
  library(skimr)
  library(janitor)
})

source(file.path(if (dir.exists("code")) "code" else ".", "00_config.R"))

# -------------------------------------------------------------------------------------
# File paths
# -------------------------------------------------------------------------------------

path_deaths <- HVI_PATHS$raw$mortality
path_ed     <- HVI_PATHS$raw$ed
path_ems    <- HVI_PATHS$raw$ems

# -------------------------------------------------------------------------------------
# Helper function: dataset audit
# -------------------------------------------------------------------------------------

audit_dataset <- function(df, name) {
  
  cat("\n=====================================================\n")
  cat("DATASET:", name, "\n")
  cat("=====================================================\n\n")
  
  # Basic structure
  cat("Dimensions:\n")
  print(dim(df))
  
  cat("\nColumn names:\n")
  print(names(df))
  
  cat("\nData types:\n")
  print(sapply(df, class))
  
  # Missingness
  cat("\nMissingness (%):\n")
  miss <- sapply(df, function(x) mean(is.na(x))) * 100
  print(sort(miss, decreasing = TRUE))
  
  # Unique counts (for categorical insight)
  cat("\nUnique values (first 20 vars):\n")
  uniq <- sapply(df[ , 1:min(20, ncol(df))], function(x) length(unique(x)))
  print(uniq)
  
  # Head
  cat("\nPreview:\n")
  print(head(df, 5))
  
  # skimr summary
  cat("\nDetailed summary (skimr):\n")
  print(skim(df))
}

# -------------------------------------------------------------------------------------
# Load datasets
# -------------------------------------------------------------------------------------

cat("\nLoading datasets...\n")

deaths <- fread(path_deaths)
ed     <- fread(path_ed)
ems    <- fread(path_ems)

# Clean column names for consistency
deaths <- clean_names(deaths)
ed     <- clean_names(ed)
ems    <- clean_names(ems)

# -------------------------------------------------------------------------------------
# Run audits
# -------------------------------------------------------------------------------------

audit_dataset(deaths, "Mortality")
audit_dataset(ed, "ED Visits")
audit_dataset(ems, "EMS Calls")

# -------------------------------------------------------------------------------------
# Identify key variables (semi-automated)
# -------------------------------------------------------------------------------------

detect_key_fields <- function(df, name) {
  
  cat("\n-----------------------------------------------------\n")
  cat("Key Field Detection:", name, "\n")
  cat("-----------------------------------------------------\n")
  
  vars <- names(df)
  
  cat("\nPotential DATE/TIME fields:\n")
  print(vars[grepl("date|time|dttm|dt", vars, ignore.case = TRUE)])
  
  cat("\nPotential LOCATION fields:\n")
  print(vars[grepl("community|zip|lat|lon|tract|geo|address", vars, ignore.case = TRUE)])
  
  cat("\nPotential OUTCOME fields:\n")
  print(vars[grepl("death|visit|call|result|dx|diag|chief|complaint", vars, ignore.case = TRUE)])
}

detect_key_fields(deaths, "Mortality")
detect_key_fields(ed, "ED Visits")
detect_key_fields(ems, "EMS Calls")

# -------------------------------------------------------------------------------------
# Date parsing check (you will need to customize column names after inspection)
# -------------------------------------------------------------------------------------

check_date <- function(df, col) {
  cat("\nChecking date parsing for:", col, "\n")
  print(head(df[[col]]))
  print(summary(as.Date(df[[col]])))
}

# Example placeholders (update after inspection)
check_date(deaths, "dod")
check_date(ed, "enc_admit_date")
#check_date(ems, "patient_assessment_date_time_eexam_03")

# -------------------------------------------------------------------------------------
# Output quick summary tables
# -------------------------------------------------------------------------------------

summary_table <- function(df, name) {
  data.frame(
    dataset = name,
    n_rows = nrow(df),
    n_cols = ncol(df)
  )
}

summary_df <- bind_rows(
  summary_table(deaths, "mortality"),
  summary_table(ed, "ed"),
  summary_table(ems, "ems")
)

print(summary_df)

# =====================================================================================
# End script
# =====================================================================================






