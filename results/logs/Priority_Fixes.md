# Priority Fixes (Static Audit)

1. **Critical stop errors**
   - Define/replace `metadata_d14` and `feature_d14` in `Figure2_ABC_spatial3cluster.R`.
   - Add missing loader module currently referenced as `src/data_loader.R`.

2. **Scope-control violations**
   - Split controlled comparisons and mixed comparisons in outputs/titles.
   - Keep `D14_Collagen vs D2_Fibrin`, `D14_Wall vs D2_Fibrin`, `D14_Collagen vs D7_Fibrin` as exploratory only.

3. **Identifier integrity**
   - Stop applying `str_to_title()` to biological gene symbols.
   - Keep `UniqueSymbol` internal; never treat it as true gene symbol in biological interpretation.

4. **Configuration debt**
   - Move all hardcoded paths, thresholds, seed, and contrast definitions to `config/config.yml` + `config/contrasts.csv`.

5. **BLI-MS statistical robustness**
   - Current script performs fixed LOD substitution + Welch t-test only; add replicate-aware detection frequency, background controls, and sensitivity analyses in later modules.

6. **Reproducibility safeguards**
   - Replace package auto-install in analysis scripts with explicit dependency checks.
   - Ensure output file names are contrast-aware and cannot overwrite prior runs.
