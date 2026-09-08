#' Validate a SingleCellExperiment input for ASE modelling
#'
#' Ensures the count assays are finite, nonnegative integers and that successes
#' do not exceed totals. Model-specific covariates are checked by `fit_glmm_bb`.
#' @param sce A `SingleCellExperiment` object.
#' @return Invisibly returns `TRUE` when validation passes.
#' @export
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

  for (assay in required_assays) {
    counts <- SummarizedExperiment::assay(sce, assay)
    # Sparse implicit zeros are valid; inspect stored values without densifying.
    values <- if (inherits(counts, "sparseMatrix")) counts@x else as.vector(counts)
    if (!is.numeric(values) || any(!is.finite(values)) ||
        any(values < 0 | values != floor(values))) {
      stop("Assay `", assay, "` must contain finite nonnegative integer counts.",
           call. = FALSE)
    }
  }
  if (any(SummarizedExperiment::assay(sce, "a1") >
          SummarizedExperiment::assay(sce, "tot"))) {
    stop("Assay `a1` cannot exceed `tot`.", call. = FALSE)
  }

  invisible(TRUE)
}
