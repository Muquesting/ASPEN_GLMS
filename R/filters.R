#' Filter genes by coverage across cells
#'
#' Retains genes meeting minimum coverage and cell thresholds.
#' @param sce A `SingleCellExperiment` with assays `a1` and `tot`.
#' @param min_trials Minimum number of total allelic trials per cell to count.
#' @param min_cells Minimum number of cells meeting `min_trials` required per gene.
#' @return A filtered `SingleCellExperiment`.
#' @export
filter_coverage <- function(sce, min_trials = 5, min_cells = 50) {
  check_sce(sce)
  check_positive_integer(min_trials, "min_trials")
  check_positive_integer(min_cells, "min_cells")
  totals <- SummarizedExperiment::assay(sce, "tot")
  keep_gene <- Matrix::rowSums(totals >= min_trials) >= min_cells
  sce[keep_gene, , drop = FALSE]
}

check_positive_integer <- function(x, name) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) ||
      x < 1 || x != floor(x)) {
    stop("`", name, "` must be a positive integer.", call. = FALSE)
  }
}
