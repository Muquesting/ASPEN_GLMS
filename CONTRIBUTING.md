# Contributing

Start with the synthetic demo and open an issue describing the behaviour you
want to improve. Small, focused pull requests are welcome.

1. Fork/clone the repository and create a branch.
2. Install a compatible R/Bioconductor environment and restore `renv.lock`.
3. Make the change; document public functions with roxygen2.
4. Run `roxygen2::roxygenise()`, `rcmdcheck::rcmdcheck(args = "--as-cran", error_on = "warning")`,
   then install the package and run `inst/scripts/demo_pipeline.R`.
5. Include the trigger, expected/actual result, relevant tests, and limitations in the PR.

Public CI uses synthetic data. Tests requiring a local `tests/data/` dataset
are optional and skipped when that data is absent; they are not evidence of
public reproducibility. Do not upload private data or change scientific defaults
to make a benchmark look better. Dependency updates require an updated snapshot
and a successful restoration check. Preserve existing licence attribution.

For scientific changes, explain the statistical assumptions and use a seeded
simulation with known truth. A small demo is not a substitute for calibration.
