skip_on_cran()

real_data_path <- test_path("..", "data", "aspensce_sexupdated.rds")

test_that("fit_glmm_bb runs on ASPEN real dataset subset", {
  skip_if_not_installed("glmmTMB")
  skip_if_not(file.exists(real_data_path), "Real dataset not available")

  fit_info <- real_aspen_fit()
  res <- fit_info$fit

  expect_s3_class(res, "tbl_df")
  expect_true(nrow(res) > 0)
  expect_true(all(c("gene", "term", "estimate", "std_error", "converged") %in% names(res)))
  expect_true(any(res$converged))
})
