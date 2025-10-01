#' Filter genes by coverage across cells
#'
#' Retains genes meeting minimum coverage and cell thresholds.
#' @param sce A `SingleCellExperiment` with assays `a1` and `tot`.
#' @param min_trials Minimum number of total allelic trials per cell to count.
#' @param min_cells Minimum number of cells meeting `min_trials` required per gene.
#' @return A filtered `SingleCellExperiment`.
filter_coverage <- function(sce, min_trials = 5, min_cells = 50) {
  check_sce(sce)
  totals <- SummarizedExperiment::assay(sce, "tot")
  keep_gene <- Matrix::rowSums(totals >= min_trials) >= min_cells
  sce[keep_gene, , drop = FALSE]
}
