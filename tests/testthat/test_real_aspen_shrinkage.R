skip_on_cran()

real_data_path <- test_path("..", "data", "aspensce_sexupdated.rds")

test_that("shrink_with_ash runs on ASPEN real dataset subset", {
  skip_if_not_installed("glmmTMB")
  skip_if_not_installed("ashr")
  skip_if_not(file.exists(real_data_path), "Real dataset not available")

  fit_info <- real_aspen_fit()
  res <- fit_info$fit

  target_term <- "sexMale"
  skip_if(!(target_term %in% res$term), "Target term not present in fit")

  shrunk <- shrink_with_ash(res, term = target_term)

  expect_s3_class(shrunk, "tbl_df")
  expect_true(nrow(shrunk) > 0)
  expect_true(all(is.finite(shrunk$beta_shrunk)))
  expect_true(all(shrunk$lfsr >= 0 & shrunk$lfsr <= 1))
})
