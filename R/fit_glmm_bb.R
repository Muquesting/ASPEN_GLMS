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
#' returns Wald statistics. Covariates must exist in `colData(sce)` and be
#' complete. Optimizer success, a positive-definite Hessian and finite
#' coefficient statistics are required. Genes failing these checks are reported with
#' `converged = FALSE` and contain an error message in the `error` column.
#' @export
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
  check_positive_integer(ncores, "ncores")
  if (!inherits(formula_fixed, "formula") || length(formula_fixed) != 2L) {
    stop("`formula_fixed` must be a one-sided formula.", call. = FALSE)
  }
  if (!is.null(rand) && (!is.character(rand) || length(rand) != 1L || is.na(rand))) {
    stop("`rand` must be NULL or a single string.", call. = FALSE)
  }
  fixed_rhs <- paste(deparse(formula_fixed[[2L]]), collapse = " ")
  random_part <- if (!is.null(rand) && nzchar(rand)) paste("+", rand) else ""
  model_formula <- stats::as.formula(
    paste0("cbind(a1, tot - a1) ~ ", fixed_rhs, " ", random_part),
    env = environment(formula_fixed)
  )
  required_cols <- all.vars(stats::as.formula(
    paste0("~ ", fixed_rhs, " ", random_part), env = environment(formula_fixed)
  ))
  col_df <- as.data.frame(SummarizedExperiment::colData(sce))
  missing_cols <- setdiff(required_cols, names(col_df))
  if (length(missing_cols)) {
    stop("Missing colData columns: ", paste(missing_cols, collapse = ", "), call. = FALSE)
  }
  if (length(required_cols) && any(!stats::complete.cases(col_df[, required_cols, drop = FALSE]))) {
    stop("Model covariates must not contain missing values.", call. = FALSE)
  }
  for (column in required_cols) {
    if (is.numeric(col_df[[column]]) && any(!is.finite(col_df[[column]]))) {
      stop("Model covariates must be finite.", call. = FALSE)
    }
  }

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

  genes <- rownames(sce_filt)
  if (is.null(genes)) {
    genes <- as.character(seq_len(nrow(sce_filt)))
  }

  assay_a1 <- SummarizedExperiment::assay(sce_filt, "a1")
  assay_tot <- SummarizedExperiment::assay(sce_filt, "tot")

  fit_one <- function(i) {
    g <- genes[[i]]
    successes <- as.numeric(assay_a1[i, ])
    totals <- as.numeric(assay_tot[i, ])
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
      problem <- fit_diagnostic_error(model, cond)
      if (!is.null(problem)) stop(problem, call. = FALSE)
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
    old_plan <- future::plan()
    on.exit({
      future::plan(old_plan)
    }, add = TRUE)
    future::plan(future::multisession, workers = ncores)
    fits <- future.apply::future_lapply(seq_along(genes), fit_one, future.seed = TRUE)
  } else {
    fits <- lapply(seq_along(genes), fit_one)
  }

  non_null <- Filter(Negate(is.null), fits)
  if (length(non_null) == 0) {
    warning("All genes were dropped before fitting; returning empty tibble.")
    return(tibble::tibble())
  }

  dplyr::bind_rows(non_null)
}

fit_diagnostic_error <- function(model, coefficients) {
  if (!isTRUE(model$fit$convergence == 0L)) {
    return(paste("Optimizer did not converge:", paste(model$fit$message, collapse = "; ")))
  }
  if (!isTRUE(model$sdr$pdHess)) return("Hessian is not positive definite.")
  needed <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
  if (!is.matrix(coefficients) || !nrow(coefficients) ||
      !all(needed %in% colnames(coefficients))) return("Coefficient statistics are missing.")
  if (any(!is.finite(coefficients[, needed, drop = FALSE])) ||
      any(coefficients[, "Std. Error"] <= 0) ||
      any(coefficients[, "Pr(>|z|)"] < 0 | coefficients[, "Pr(>|z|)"] > 1)) {
    return("Coefficient statistics are invalid or non-finite.")
  }
  NULL
}
