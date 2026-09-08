#' Tidy coefficient table with multiple-testing correction
#'
#' @param fit_tbl Tibble returned by `fit_glmm_bb`.
#' @param method P-value adjustment method, passed to `stats::p.adjust`.
#' @return Tibble with per-term adjusted p-values.
#' @export
tidy_contrasts <- function(fit_tbl, method = "BH") {
  required_cols <- c("gene", "term", "p_value", "converged")
  missing_cols <- setdiff(required_cols, colnames(fit_tbl))
  if (length(missing_cols) > 0) {
    stop(
      "fit_tbl is missing columns: ", paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }

  valid <- fit_tbl$converged %in% TRUE & !is.na(fit_tbl$term) &
    is.finite(fit_tbl$p_value) & fit_tbl$p_value >= 0 & fit_tbl$p_value <= 1
  usable <- fit_tbl[which(valid), , drop = FALSE]
  if (!nrow(usable)) {
    warning("No converged coefficients available for adjustment.")
    return(tibble::tibble())
  }

  usable$fdr <- stats::ave(usable$p_value, usable$term,
                           FUN = function(p) stats::p.adjust(p, method = method))
  usable
}
