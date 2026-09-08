# Roadmap

The first release focuses on installation, explicit model diagnostics, and a
reproducible synthetic workflow. These follow-up tasks are intentionally bounded:

1. **[Paired-donor tutorial](https://github.com/Muquesting/ASPEN_GLMS/issues/2):** demonstrate a repeated-measures design with a seeded
   synthetic dataset, explain the random effect and reference levels, and test
   that the documented commands run without private inputs.
2. **[Simulation calibration harness](https://github.com/Muquesting/ASPEN_GLMS/issues/3):** predefine null/effect scenarios and report
   convergence, type-I error and interval coverage without tuning to the results.
   This is needed before broader scientific-validity claims.
3. **[Sparse-input benchmark](https://github.com/Muquesting/ASPEN_GLMS/issues/4):** compare dense/sparse versions of identical synthetic
   counts, report session information and elapsed/memory observations, and test
   numerical equivalence within a declared tolerance.

No completion date, adoption level, or performance advantage is promised.
