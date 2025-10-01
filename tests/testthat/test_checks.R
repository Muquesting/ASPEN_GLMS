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
  expect_error(check_sce(bad2), "Missing colData")
})
