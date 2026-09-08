test_that("optimizer and Hessian failures cannot produce valid estimates", {
  coef <- matrix(c(0.2, 0.1, 2, 0.0455), nrow = 1,
                 dimnames = list("effect", c("Estimate", "Std. Error", "z value", "Pr(>|z|)")))
  model <- list(fit = list(convergence = 0L), sdr = list(pdHess = TRUE))
  expect_null(ASPENGLMS:::fit_diagnostic_error(model, coef))
  model$fit$convergence <- 1L
  expect_match(ASPENGLMS:::fit_diagnostic_error(model, coef), "Optimizer")
  model$fit$convergence <- 0L
  model$sdr$pdHess <- FALSE
  expect_match(ASPENGLMS:::fit_diagnostic_error(model, coef), "Hessian")
  model$sdr$pdHess <- TRUE
  coef[1, "Std. Error"] <- Inf
  expect_match(ASPENGLMS:::fit_diagnostic_error(model, coef), "non-finite")
})

test_that("failed model fitting returns an explicit failure row", {
  sce <- toy_sce()
  SummarizedExperiment::colData(sce)$constant <- factor(rep("only", ncol(sce)))
  res <- suppressWarnings(fit_glmm_bb(sce, ~ constant, rand = NULL, min_cells = 10))
  expect_equal(nrow(res), nrow(sce))
  expect_true(all(!res$converged))
  expect_true(all(is.na(res$estimate)))
  expect_true(all(nzchar(res$error)))
})

test_that("zero coverage returns an empty result with explanation", {
  sce <- toy_sce()
  expect_warning(res <- fit_glmm_bb(sce, ~ sex, rand = NULL, min_cells = 1000), "No genes")
  expect_equal(nrow(res), 0L)
})

test_that("gene positions work without row names", {
  sce <- toy_sce()
  rownames(sce) <- NULL
  res <- suppressWarnings(fit_glmm_bb(sce, ~ sex, rand = NULL, min_cells = 10))
  expect_true(all(res$gene %in% as.character(seq_len(nrow(sce)))))
  expect_true(any(res$converged))
})

test_that("invalid estimates are excluded from downstream inference", {
  tbl <- tibble::tibble(gene = letters[1:7], term = rep("effect", 7),
    estimate = c(0.2, 0.5, Inf, 1, 0.2, 0.1, 0.1),
    std_error = c(0.1, 0.2, 0.1, 0, 0.1, 0.1, 0.1),
    p_value = c(0.02, 0.04, Inf, -1, 0.01, 1.5, NA),
    converged = c(TRUE, TRUE, TRUE, TRUE, FALSE, NA, TRUE))
  shrunk <- shrink_with_ash(tbl, "effect")
  expect_setequal(shrunk$gene, c("a", "b", "g"))
  expect_true(all(is.finite(shrunk$beta_shrunk)))
  expect_true(all(shrunk$lfsr >= 0 & shrunk$lfsr <= 1))
  contrasts <- tidy_contrasts(tbl)
  expect_equal(contrasts$gene, c("a", "b"))
  expect_equal(contrasts$fdr, c(0.04, 0.04))
})

test_that("BH adjustment is performed separately for each coefficient", {
  tbl <- tibble::tibble(gene = rep(c("a", "b"), 2),
    term = rep(c("first", "second"), each = 2),
    p_value = c(0.01, 0.2, 0.04, 0.06), converged = TRUE)
  res <- tidy_contrasts(tbl)
  expect_equal(res$fdr, c(0.02, 0.2, 0.06, 0.06))
})
