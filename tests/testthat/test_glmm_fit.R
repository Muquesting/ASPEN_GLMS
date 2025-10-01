test_that("fit_glmm_bb returns coefficients", {
  skip_if_not_installed("glmmTMB")
  sce <- toy_sce()
  res <- fit_glmm_bb(sce, formula_fixed = ~ sex + age, min_cells = 10, min_trials = 1, ncores = 1)
  expect_s3_class(res, "tbl_df")
  expect_true(all(c("gene", "term", "estimate", "std_error") %in% colnames(res)))
  expect_true(any(res$converged))
})
