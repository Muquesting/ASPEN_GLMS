#' Quick diagnostic plot for a single coefficient
#'
#' @param fit_tbl Tibble returned by `fit_glmm_bb`.
#' @param term Coefficient name to visualise.
#' @return A `ggplot` object showing estimates vs. standard errors.
#' @export
plot_coefficient_diagnostics <- function(fit_tbl, term) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package `ggplot2` is required for plotting.", call. = FALSE)
  }

  usable <- fit_tbl$converged %in% TRUE & !is.na(fit_tbl$term) &
    fit_tbl$term == term & is.finite(fit_tbl$estimate) &
    is.finite(fit_tbl$std_error) & fit_tbl$std_error > 0 &
    is.finite(fit_tbl$p_value) & fit_tbl$p_value >= 0 & fit_tbl$p_value <= 1
  subset <- fit_tbl[which(usable), , drop = FALSE]
  if (!nrow(subset)) {
    stop("No converged coefficients found for the requested term.", call. = FALSE)
  }

  subset$fdr <- stats::p.adjust(subset$p_value, method = "BH")

  ggplot2::ggplot(subset, ggplot2::aes(x = .data$estimate, y = .data$std_error, colour = .data$fdr)) +
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
