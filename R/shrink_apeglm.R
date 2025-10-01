#' Beta-binomial shrinkage using apeglm
#'
#' @param Y_succ Matrix of successes (e.g. reference allele counts).
#' @param Y_tot Matrix of trials (total counts).
#' @param design Model matrix with rows matching columns of `Y_succ`.
#' @param coef Index or name of the coefficient to shrink.
#' @param ... Additional arguments passed to `apeglm::apeglm`.
#' @return Result object returned by `apeglm` containing MAP estimates and FSR.
apeglm_bb <- function(Y_succ, Y_tot, design, coef, ...) {
  if (!requireNamespace("apeglm", quietly = TRUE)) {
    stop("Package `apeglm` is required for shrinkage.", call. = FALSE)
  }

  if (!all(dim(Y_succ) == dim(Y_tot))) {
    stop("`Y_succ` and `Y_tot` must have the same dimensions.", call. = FALSE)
  }

  if (nrow(design) != ncol(Y_succ)) {
    stop("Number of rows in `design` must match number of cells.", call. = FALSE)
  }

  apeglm::apeglm(
    Y = Y_succ,
    x = design,
    log.lik = apeglm::logLik_bb,
    param = list(size = Y_tot),
    coef = coef,
    method = "betabinCR",
    ...
  )
}
