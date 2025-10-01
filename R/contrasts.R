#' Tidy coefficient table with multiple-testing correction
#'
#' @param fit_tbl Tibble returned by `fit_glmm_bb`.
#' @param method P-value adjustment method, passed to `stats::p.adjust`.
#' @return Tibble with per-term adjusted p-values.
tidy_contrasts <- function(fit_tbl, method = "BH") {
  required_cols <- c("gene", "term", "p_value", "converged")
  missing_cols <- setdiff(required_cols, colnames(fit_tbl))
  if (length(missing_cols) > 0) {
    stop(
      "fit_tbl is missing columns: ", paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }

  usable <- fit_tbl[fit_tbl$converged & !is.na(fit_tbl$term), , drop = FALSE]
  if (!nrow(usable)) {
    warning("No converged coefficients available for adjustment.")
    return(tibble::tibble())
  }

  usable <- dplyr::group_by(usable, term)
  usable <- dplyr::mutate(usable, fdr = stats::p.adjust(p_value, method = method))
  dplyr::ungroup(usable)
}
