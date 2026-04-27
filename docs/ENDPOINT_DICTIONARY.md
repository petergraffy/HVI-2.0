# Endpoint Dictionary

HVI 2.0 combines mortality, ED, and EMS endpoints at the Chicago community-area by day level.

## Mortality

- `deaths`: all-cause mortality
- `death_cvd`: cardiovascular mortality
- `death_respiratory`: respiratory mortality
- `death_renal`: renal mortality
- `death_heat`: heat-coded mortality
- `death_dehydration`: dehydration mortality
- `death_syncope`: syncope mortality
- `death_mental`: mental health mortality
- `death_neurologic`: neurologic mortality
- `death_injury`: injury mortality
- `death_gi`: gastrointestinal mortality

## Emergency Department

- `ed_visits`: all-cause ED visits
- `ed_cvd`: cardiovascular ED visits
- `ed_respiratory`: respiratory ED visits
- `ed_renal`: renal ED visits
- `ed_heat`: heat illness ED visits
- `ed_dehydration`: dehydration ED visits
- `ed_syncope`: syncope ED visits
- `ed_mental`: mental health ED visits
- `ed_neurologic`: neurologic ED visits
- `ed_injury`: injury ED visits
- `ed_gi`: gastrointestinal ED visits

## EMS

- `ems_calls`: all-cause EMS calls
- `ems_cvd`: cardiovascular EMS calls
- `ems_respiratory`: respiratory EMS calls
- `ems_heat`: heat-related EMS calls
- `ems_syncope`: syncope/collapse EMS calls
- `ems_neuro`: neurologic EMS calls
- `ems_mental`: mental health EMS calls
- `ems_injury`: injury EMS calls
- `ems_gi`: gastrointestinal EMS calls
- `ems_bleeding`: bleeding EMS calls

Endpoint definitions should be kept synchronized with `code/03_panel_causes.R` and `hvi_endpoint_metadata.csv`.
