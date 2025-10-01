#' Fit per-gene beta-binomial GLMMs to allelic counts
#'
#' @param sce A `SingleCellExperiment` with `a1` and `tot` assays.
#' @param formula_fixed Fixed-effects portion supplied as a formula (rhs only).
#' @param rand Random-effects term (string) appended to the model.
#' @param family Family passed to `glmmTMB` (defaults to beta-binomial with logit link).
#' @param min_trials Minimum trials per cell used during filtering.
#' @param min_cells Minimum cells meeting the trial threshold per gene.
#' @param ncores Number of workers for parallel fitting (`future.apply`).
#' @return A tibble with one row per gene/term combination containing estimates.
#' @details The function filters genes by coverage, fits `glmmTMB` models and
#' returns Wald statistics. Genes failing to converge are reported with
#' `converged = FALSE` and contain an error message in the `error` column.
fit_glmm_bb <- function(
  sce,
  formula_fixed = ~ sex + age + celltype_new + sex:age,
  rand = "(1|sample)",
  family = glmmTMB::betabinomial(link = "logit"),
  min_trials = 5,
  min_cells = 50,
  ncores = 1
) {
  check_sce(sce)

  if (!requireNamespace("glmmTMB", quietly = TRUE)) {
    stop("Package `glmmTMB` is required for model fitting.", call. = FALSE)
  }

  if (ncores > 1) {
    if (!requireNamespace("future", quietly = TRUE)) {
      stop("Package `future` is required for parallel fitting.", call. = FALSE)
    }
    if (!requireNamespace("future.apply", quietly = TRUE)) {
      stop("Package `future.apply` is required for parallel fitting.", call. = FALSE)
    }
  }

  sce_filt <- filter_coverage(sce, min_trials = min_trials, min_cells = min_cells)
  if (nrow(sce_filt) == 0) {
    warning("No genes passed the coverage filter; returning empty tibble.")
    return(tibble::tibble())
  }

  col_df <- as.data.frame(SummarizedExperiment::colData(sce_filt))
  col_df$cell_id <- colnames(sce_filt)

  fixed_rhs <- paste(deparse(formula_fixed), collapse = " ")
  fixed_rhs <- trimws(sub("^~", "", fixed_rhs))
  if (!nzchar(fixed_rhs)) {
    fixed_rhs <- "1"
  }
  random_part <- if (!is.null(rand) && nzchar(rand)) paste("+", rand) else ""
  model_formula <- stats::as.formula(
    paste0("cbind(a1, tot - a1) ~ ", fixed_rhs, " ", random_part)
  )

  genes <- rownames(sce_filt)
  if (is.null(genes)) {
    genes <- as.character(seq_len(nrow(sce_filt)))
  }

  assay_a1 <- SummarizedExperiment::assay(sce_filt, "a1")
  assay_tot <- SummarizedExperiment::assay(sce_filt, "tot")

  fit_one <- function(g) {
    successes <- as.numeric(assay_a1[g, ])
    totals <- as.numeric(assay_tot[g, ])
    keep <- totals >= min_trials
    if (sum(keep) < min_cells) {
      return(NULL)
    }

    df <- col_df[keep, , drop = FALSE]
    df$a1 <- successes[keep]
    df$tot <- totals[keep]

    tryCatch({
      model <- glmmTMB::glmmTMB(model_formula, family = family, data = df)
      cond <- summary(model)$coefficients$cond
      tibble::tibble(
        gene = g,
        term = rownames(cond),
        estimate = cond[, "Estimate"],
        std_error = cond[, "Std. Error"],
        z_value = cond[, "z value"],
        p_value = cond[, "Pr(>|z|)"],
        converged = TRUE,
        error = NA_character_
      )
    }, error = function(e) {
      tibble::tibble(
        gene = g,
        term = NA_character_,
        estimate = NA_real_,
        std_error = NA_real_,
        z_value = NA_real_,
        p_value = NA_real_,
        converged = FALSE,
        error = conditionMessage(e)
      )
    })
  }

  if (ncores > 1) {
    future::plan(future::multisession, workers = ncores)
    on.exit({
      future::plan(future::sequential)
    }, add = TRUE)
    fits <- future.apply::future_lapply(genes, fit_one)
  } else {
    fits <- lapply(genes, fit_one)
  }

  non_null <- Filter(Negate(is.null), fits)
  if (length(non_null) == 0) {
    warning("All genes were dropped before fitting; returning empty tibble.")
    return(tibble::tibble())
  }

  dplyr::bind_rows(non_null)
}
