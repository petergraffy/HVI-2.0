# Archived Scripts

This folder keeps exploratory, prototype, or superseded scripts out of the production pipeline while preserving their history for reference.

## Archived in April 2026 Cleanup

- `00_hvi_prototype_output.R`: prototype frontend export bundle superseded by `09e_build_frontend_exports.R`, `11_build_public_exports.R`, and `12_build_scenario_exports.R`.
- `04_modeling_prototype.R`: early heat-response and prototype HVI workflow superseded by the MRT, model-matrix, variable-selection, and endpoint-modeling stages.
- `08_final_variables.R`: alternate variable-selection script. The canonical production script is `code/08_variable_selection_hvi.R`.
- `09_fit_hvi_model.R`: monolithic HVI model fit superseded by modular `09a` through `09e` scripts.
- `10_pub_outputs.R`: alternate publication-output script. The canonical production script is `code/10_publication_output.R`.

Archived scripts should not be called by `_targets.R` or public export automation. If a useful function remains in one of these files, move that function into a production script or shared helper before using it in the pipeline.
