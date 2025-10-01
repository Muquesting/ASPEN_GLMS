#' Empirical Bayes shrinkage of GLMM coefficients using ashr
#'
#' @param fit_tbl Tibble returned by `fit_glmm_bb`.
#' @param term Coefficient name to shrink (e.g. "sexMale").
#' @param ... Additional arguments passed to `ashr::ash`.
#' @return Tibble with original and shrunk estimates plus local false sign rate.
shrink_with_ash <- function(fit_tbl, term, ...) {
  if (!requireNamespace("ashr", quietly = TRUE)) {
    stop("Package `ashr` is required for shrinkage.", call. = FALSE)
  }

  required_cols <- c("gene", "term", "estimate", "std_error", "converged")
  missing_cols <- setdiff(required_cols, colnames(fit_tbl))
  if (length(missing_cols) > 0) {
    stop(
      "fit_tbl is missing columns: ", paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }

  subset <- fit_tbl[fit_tbl$converged & fit_tbl$term == term, , drop = FALSE]
  subset <- subset[!is.na(subset$estimate) & !is.na(subset$std_error), , drop = FALSE]

  if (nrow(subset) == 0) {
    stop("No converged coefficients found for the requested term.", call. = FALSE)
  }

  shrink_fit <- ashr::ash(subset$estimate, subset$std_error, ...)

  tibble::tibble(
    gene = subset$gene,
    term = term,
    beta = subset$estimate,
    std_error = subset$std_error,
    beta_shrunk = ashr::get_pm(shrink_fit),
    lfsr = ashr::get_lfsr(shrink_fit)
  )
}
