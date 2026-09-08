#' Beta-binomial shrinkage using apeglm
#'
#' @param Y_succ Matrix of successes (e.g. reference allele counts).
#' @param Y_tot Matrix of trials (total counts).
#' @param design Model matrix with rows matching columns of `Y_succ`.
#' @param coef Index or name of the coefficient to shrink.
#' @param ... Additional arguments passed to `apeglm::apeglm`. Must include
#'   `param = cbind(theta, Y_tot)`, where `theta` contains externally estimated
#'   positive beta-binomial concentration parameters, one per gene.
#' @details Experimental fixed-effects interface, not a mixed-model shrinker.
#'   Dispersion is never guessed. Supply `mle` through `...` for an adaptive
#'   prior; otherwise apeglm uses its default prior. See the apeglm ASE tutorial.
#'   Uses apeglm's R beta-binomial optimizer and rejects unsuccessful diagnostics
#'   or invalid posterior estimates instead of returning them for inference.
#' @return Result object returned by `apeglm` containing MAP estimates and FSR.
#' @export
apeglm_bb <- function(Y_succ, Y_tot, design, coef, ...) {
  if (!requireNamespace("apeglm", quietly = TRUE)) {
    stop("Package `apeglm` is required for shrinkage.", call. = FALSE)
  }

  if (!is.matrix(Y_succ) || !is.matrix(Y_tot) ||
      !identical(dim(Y_succ), dim(Y_tot)) || !length(Y_succ)) {
    stop("`Y_succ` and `Y_tot` must have the same dimensions.", call. = FALSE)
  }

  if (!is.numeric(Y_succ) || !is.numeric(Y_tot) ||
      any(!is.finite(Y_succ)) || any(!is.finite(Y_tot)) ||
      any(Y_succ < 0 | Y_succ != floor(Y_succ)) ||
      any(Y_tot < 0 | Y_tot != floor(Y_tot)) || any(Y_succ > Y_tot)) {
    stop("Counts must be finite nonnegative integers with successes <= totals.", call. = FALSE)
  }
  if (!is.matrix(design) || !is.numeric(design) ||
      nrow(design) != ncol(Y_succ) || any(!is.finite(design))) {
    stop("Number of rows in `design` must match number of cells.", call. = FALSE)
  }

  if (is.character(coef)) coef <- match(coef, colnames(design))
  if (length(coef) != 1L || !is.numeric(coef) || is.na(coef) ||
      coef < 1 || coef > ncol(design) || coef != floor(coef)) {
    stop("`coef` must identify one design coefficient.", call. = FALSE)
  }
  dots <- list(...)
  param <- dots$param
  if (!is.matrix(param) || !is.numeric(param) ||
      !identical(dim(param), c(nrow(Y_tot), ncol(Y_tot) + 1L)) ||
      any(!is.finite(param)) || any(param[, 1L] <= 0) ||
      !isTRUE(all.equal(unname(param[, -1L, drop = FALSE]), unname(Y_tot)))) {
    stop("Supply `param = cbind(theta, Y_tot)` with positive externally estimated theta per gene.",
         call. = FALSE)
  }
  reserved <- intersect(names(dots), c("Y", "x", "coef", "log.lik", "method", "log.link"))
  if (length(reserved)) stop("Do not override: ", paste(reserved, collapse = ", "), call. = FALSE)
  result <- do.call(apeglm::apeglm, c(list(Y = Y_succ, x = design, log.lik = NULL,
                               coef = coef, method = "betabinR", log.link = FALSE), dots))
  if (any(!is.finite(result$diag[, "conv"])) || any(result$diag[, "conv"] != 0) ||
      any(!is.finite(result$map)) || any(!is.finite(result$sd)) || any(result$sd <= 0) ||
      any(!is.finite(result$fsr)) || any(result$fsr < 0 | result$fsr > 1)) {
    stop("apeglm optimization failed or returned invalid posterior estimates.", call. = FALSE)
  }
  result
}
