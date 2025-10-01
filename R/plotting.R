#' Quick diagnostic plot for a single coefficient
#'
#' @param fit_tbl Tibble returned by `fit_glmm_bb`.
#' @param term Coefficient name to visualise.
#' @return A `ggplot` object showing estimates vs. standard errors.
plot_coefficient_diagnostics <- function(fit_tbl, term) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package `ggplot2` is required for plotting.", call. = FALSE)
  }

  subset <- fit_tbl[fit_tbl$converged & fit_tbl$term == term, , drop = FALSE]
  subset <- subset[!is.na(subset$estimate) & !is.na(subset$std_error), , drop = FALSE]
  if (!nrow(subset)) {
    stop("No converged coefficients found for the requested term.", call. = FALSE)
  }

  subset$fdr <- stats::p.adjust(subset$p_value, method = "BH")

  ggplot2::ggplot(subset, ggplot2::aes(x = estimate, y = std_error, colour = fdr)) +
    ggplot2::geom_point(alpha = 0.8) +
    ggplot2::scale_colour_viridis_c(option = "B", end = 0.9) +
    ggplot2::labs(
      title = paste0("Diagnostic for term: ", term),
      x = "Estimate (logit scale)",
      y = "Standard error",
      colour = "FDR"
    ) +
    ggplot2::theme_minimal()
}
