# ASPENGLMS 0.1.0 (release candidate)

- Validate count domains and formula-specific covariates before fitting.
- Require optimizer success, positive-definite Hessian and finite statistics.
- Exclude failed or invalid rows from coefficient summaries.
- Preserve gene labels with positional indexing and restore the caller's future plan.
- Replace the broken apeglm call with a documented experimental interface that
  requires explicit beta-binomial concentration estimates.
- Add regression tests, generated function help, dependency capture and a seeded
  donor-aware synthetic tutorial, plus cross-platform package/demo checks.

This is early-stage research software. These engineering checks do not establish
scientific calibration, broad adoption, or superiority to other methods.
