# Synthetic-data walkthrough

This tutorial exercises software behaviour; it is not a validation study.

## 1. Run the reproducible example

After installing the package, run:

```r
source(system.file("scripts", "demo_pipeline.R", package = "ASPENGLMS"))
```

The demo uses seed 42 and generates beta-binomial counts with donor-level sex,
age, and random effects. It fits `~ sex + age + celltype_new` with `(1|sample)`.
The simulated inputs contain no human participant data.

## 2. Understand the input

```r
sce <- simulate_demo_sce() # defined by the example script
SummarizedExperiment::assayNames(sce)
head(as.data.frame(SummarizedExperiment::colData(sce)))
ASPENGLMS::check_sce(sce)
```

Rows are genes and columns are cells. `a1` counts successes for the chosen allele;
`tot` counts all informative trials. `a2` is optional and unused by fitting.
Gene and cell order must agree with the metadata. Sex and age are constant within
each synthetic donor. A custom formula can use different metadata names.

## 3. Fit and inspect diagnostics

```r
fit <- ASPENGLMS::fit_glmm_bb(sce, ~ sex + age + celltype_new,
  rand = "(1|sample)", min_trials = 5, min_cells = 30, ncores = 1)
table(fit$converged)
unique(fit[!fit$converged, c("gene", "error")])
```

Each successful gene contributes one row per coefficient. An unsuccessful gene
contributes a failure row with missing estimates and a reason. Filtering can
remove genes altogether; an empty result is reported with a warning. Investigate
failures and your design rather than treating a fit's existence as convergence.

## 4. Summarise a coefficient

```r
sex <- ASPENGLMS::shrink_with_ash(fit, "sexMale")
head(sex)
adjusted <- ASPENGLMS::tidy_contrasts(fit)
head(adjusted)
ASPENGLMS::plot_coefficient_diagnostics(fit, "sexMale")
```

`sexMale` is a conditional log odds difference relative to Female in this design.
`exp(estimate)` is an odds ratio. The `ashr` output provides posterior mean
coefficients and local false-sign rates, not ordinary adjusted p-values.
The `fdr` output is BH-adjusted within each term. These summaries depend on
model assumptions and do not establish biological discoveries in synthetic data.

## 5. Reproduce and report

Save `sessionInfo()`, the release tag, formula, coverage thresholds, and a count
of successful/failed genes. For issue reports, construct the smallest synthetic
example that exhibits the problem. Never attach restricted study data.
