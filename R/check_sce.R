#' Validate a SingleCellExperiment input for ASE modelling
#'
#' Ensures the expected assays and covariates exist before modelling.
#' @param sce A `SingleCellExperiment` object.
#' @return Invisibly returns `TRUE` when validation passes.
check_sce <- function(sce) {
  if (!inherits(sce, "SingleCellExperiment")) {
    stop("`sce` must inherit from SingleCellExperiment", call. = FALSE)
  }

  available_assays <- SummarizedExperiment::assayNames(sce)
  required_assays <- c("a1", "tot")
  missing_assays <- setdiff(required_assays, available_assays)
  if (length(missing_assays) > 0) {
    stop(
      "Missing assays: ", paste(missing_assays, collapse = ", "),
      call. = FALSE
    )
  }

  available_cols <- colnames(SummarizedExperiment::colData(sce))
  required_cols <- c("sex", "age", "celltype_new", "sample")
  missing_cols <- setdiff(required_cols, available_cols)
  if (length(missing_cols) > 0) {
    stop(
      "Missing colData columns: ", paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }

  invisible(TRUE)
}
