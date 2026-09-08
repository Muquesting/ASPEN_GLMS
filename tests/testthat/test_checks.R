test_that("check_sce validates required structure", {
  sce <- toy_sce()
  expect_invisible(check_sce(sce))

  bad <- sce
  SummarizedExperiment::assayNames(bad) <- c("a1", "a2", "missing")
  expect_error(check_sce(bad), "Missing assays")

  bad2 <- sce
  coldata <- SummarizedExperiment::colData(bad2)
  coldata$sex <- NULL
  SummarizedExperiment::colData(bad2) <- coldata
  expect_invisible(check_sce(bad2))
  expect_error(fit_glmm_bb(bad2, ~ sex, rand = NULL), "Missing colData columns: sex")
})

test_that("invalid counts fail before model fitting", {
  for (value in c(-1, 0.5, NA_real_, NaN, Inf)) {
    sce <- toy_sce()
    SummarizedExperiment::assay(sce, "a1")[1, 1] <- value
    expect_error(check_sce(sce), "finite nonnegative integer")
  }
  sce <- toy_sce()
  SummarizedExperiment::assay(sce, "a1")[1, 1] <-
    SummarizedExperiment::assay(sce, "tot")[1, 1] + 1
  expect_error(check_sce(sce), "cannot exceed")
})

test_that("sparse counts and formula-specific covariates are accepted", {
  sce <- toy_sce()
  for (name in c("a1", "tot")) {
    SummarizedExperiment::assay(sce, name) <-
      Matrix::Matrix(SummarizedExperiment::assay(sce, name), sparse = TRUE)
  }
  expect_invisible(check_sce(sce))
  cd <- SummarizedExperiment::colData(sce)
  cd$condition <- cd$age
  cd$sex <- cd$age <- cd$celltype_new <- cd$sample <- NULL
  SummarizedExperiment::colData(sce) <- cd
  res <- suppressWarnings(fit_glmm_bb(sce, ~ condition, rand = NULL, min_cells = 10))
  expect_s3_class(res, "tbl_df")
  expect_true(any(res$converged))
  expect_error(fit_glmm_bb(sce, ~ condition, rand = "(1|donor)"), "donor")
})

test_that("invalid configuration and missing covariates fail clearly", {
  sce <- toy_sce()
  expect_error(fit_glmm_bb(sce, sex ~ age), "one-sided")
  expect_error(fit_glmm_bb(sce, ncores = 0), "positive integer")
  expect_error(filter_coverage(sce, min_trials = 0), "positive integer")
  expect_error(filter_coverage(sce, min_cells = NA), "positive integer")
  SummarizedExperiment::colData(sce)$age[1] <- NA
  expect_error(fit_glmm_bb(sce, ~ age, rand = NULL), "missing values")
})
