# HVI 2.0 Documentation

This directory contains the operational documentation that supports the scientific overview in the main `README.md`. The intent is to keep the repository understandable to three audiences at once: analysts rerunning the pipeline, collaborators reviewing methods, and dashboard developers consuming public exports.

## Start Here

- `DATA_GOVERNANCE.md`: privacy boundaries, PHI handling rules, storage zones, and public-release checklist.
- `PIPELINE.md`: current reproducible workflow scaffold and the staged migration toward `targets`.
- `ENDPOINT_DICTIONARY.md`: mortality, emergency department, and EMS endpoint names used across the analysis.
- `PUBLIC_EXPORT_CONTRACT.md`: dashboard-safe export layout, small-cell suppression rules, and scenario slider files.

## Scientific Workflow

HVI 2.0 starts with private record-level EMS, ED, and mortality data. Those records are standardized in time and space, joined to environmental exposure data, and aggregated to a Chicago community-area by day panel. Public analysis products are generated only after aggregation and suppression.

The scientific pipeline has four major modeling layers:

1. Health endpoint harmonization across EMS calls, ED visits, and mortality records.
2. Environmental heat exposure scoring using endpoint-specific heat dose and MRT-based temperature summaries.
3. Structural vulnerability modeling using community-level demographic, built-environment, health, and socioeconomic covariates.
4. Public dashboard export generation, including precomputed temperature and NDVI scenario grids for interactive exploration.

## Public Tool Boundary

The public dashboard should consume only files under `public_exports/`. That directory is intentionally untracked by git because it is generated output, but its structure is described in `PUBLIC_EXPORT_CONTRACT.md` so the frontend can rely on stable filenames and fields.

Record-level files, private derived panels, model RDS files, and Box raw-data folders should never be read directly by the public app.

## Maintenance Notes

When adding a new endpoint, covariate family, model output, or dashboard field, update the relevant documentation in the same commit as the code change. At minimum, check:

- `ENDPOINT_DICTIONARY.md` for endpoint naming.
- `PIPELINE.md` for new pipeline stages or changed script responsibilities.
- `PUBLIC_EXPORT_CONTRACT.md` for any public CSV schema changes.
- `DATA_GOVERNANCE.md` if the change affects privacy boundaries or release criteria.
