# Data Governance

HVI 2.0 uses protected health data from EMS, ED, and mortality records. Treat the repository as code and public aggregate output only.

## Storage Zones

- Private raw data: `HVI_PRIVATE_DIR`, normally `C:/Users/Peter Graffy/Box/HVI2.0`.
- Private derived data: `HVI_PRIVATE_DIR/derived`.
- Private model outputs: `HVI_PRIVATE_DIR/09_model_outputs`.
- Public exports: `HVI_PUBLIC_EXPORT_DIR`, normally `public_exports`.

## Git Rules

Do not commit:

- `data/`, `code/data/`, `results/`, `code/results/`, `code/09_model_outputs/`
- record-level `*_standardized.csv` files
- raw ED, EMS, mortality, or narrative files
- `.rds`, `.nc`, parquet, or other large binary analysis artifacts

The `.gitignore` now enforces these defaults for future untracked files. Already tracked publication artifacts remain tracked until intentionally removed with `git rm --cached`.

## Public Release Checklist

Before sharing or deploying:

- Rebuild public outputs with `code/11_build_public_exports.R`.
- Check `public_exports/manifest.json`.
- Search exported headers for PHI-like fields.
- Confirm small-cell suppression threshold.
- Confirm dashboard reads from `public_exports/dashboard/`.
