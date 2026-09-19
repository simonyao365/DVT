# Code Audit Report (Stage 1: Static Review)

## Scope
Static review of existing scripts prior to module expansion:

- `/home/runner/work/DVT/DVT/Figure2_ABC.R`
- `/home/runner/work/DVT/DVT/Figure2_ABC_5cluster.R`
- `/home/runner/work/DVT/DVT/Figure2_ABC_spatial3cluster.R`
- `/home/runner/work/DVT/DVT/analyze_bli_ms_fbln7.R`

## Key Findings Summary

- **Undefined/missing objects**: Critical blockers identified (`metadata_d14`, `feature_d14`, missing `src/data_loader.R`).
- **Hardcoded local paths**: Multiple absolute Windows paths for inputs and databases; not reproducible.
- **Contrast governance**: Controlled and mixed contrasts are combined; several mixed contrasts need explicit exploratory labeling.
- **Identifier handling**: `str_to_title()` and `make.unique()` are used in ways that can distort symbols or conflate internal IDs with biological identifiers.
- **Statistical consistency risks**:
  - Potential confusion between p-value and FDR thresholds across steps.
  - BLI-MS currently relies on a single imputation + t-test pathway without broader QC evidence integration.
- **Section-title mismatch**:
  - `Figure2_ABC_spatial3cluster.R` labels PCA as D14 spatial while data matrix includes D2/D7 groups.

## Detailed Audit Artifacts

- `/home/runner/work/DVT/DVT/results/logs/Undefined_Objects.csv`
- `/home/runner/work/DVT/DVT/results/logs/Hardcoded_Paths.csv`
- `/home/runner/work/DVT/DVT/results/logs/Contrast_Audit.csv`
- `/home/runner/work/DVT/DVT/results/logs/Identifier_Audit.csv`
- `/home/runner/work/DVT/DVT/results/logs/Priority_Fixes.md`

## Immediate Recommendations for Next Step

1. Land Module 0 config-first framework and contrast manifest.
2. Refactor scripts to read all paths/thresholds from config.
3. Enforce controlled-vs-exploratory contrast labeling at data and figure generation layers.
4. Freeze biological IDs (`Gene_symbol_raw`, `Gene_symbol_clean`) and separate from internal keys.
5. Rebuild loader as explicit, validated import interface before Module 1 QC.
