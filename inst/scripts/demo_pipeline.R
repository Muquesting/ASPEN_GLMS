#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Matrix)
  library(SingleCellExperiment)
})

if (!requireNamespace("ASPENGLMS", quietly = TRUE)) {
  stop("Please install/load the ASPENGLMS package before running this script.", call. = FALSE)
}

# helper to simulate a tiny dataset when bundled demo data is absent
simulate_demo_sce <- function(n_genes = 50L, n_cells = 120L, n_samples = 6L) {
  set.seed(42)
  cell_sample <- sample(paste0("sample", seq_len(n_samples)), n_cells, replace = TRUE)
  sex <- sample(c("Female", "Male"), n_cells, replace = TRUE)
  age <- sample(c("Young", "Aged"), n_cells, replace = TRUE)
  celltype <- sample(c("T", "B", "Mono"), n_cells, replace = TRUE)

  # simulate total counts and reference successes
  tot <- matrix(rnbinom(n_genes * n_cells, size = 10, mu = 20), nrow = n_genes)
  baseline <- matrix(stats::plogis(rnorm(n_genes * n_cells, sd = 0.5)), nrow = n_genes)
  male_shift <- ifelse(sex == "Male", 0.2, -0.2)
  prob_ref <- baseline * stats::plogis(male_shift)
  prob_ref <- pmin(pmax(prob_ref, 1e-4), 1 - 1e-4)
  a1 <- matrix(rbinom(n_genes * n_cells, size = as.vector(tot), prob = as.vector(prob_ref)),
               nrow = n_genes)
  a2 <- tot - a1

  rownames(tot) <- rownames(a1) <- rownames(a2) <- paste0("gene", seq_len(n_genes))
  colnames(tot) <- colnames(a1) <- colnames(a2) <- paste0("cell", seq_len(n_cells))

  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(a1 = a1, a2 = a2, tot = tot),
    colData = data.frame(
      sex = factor(sex, levels = c("Female", "Male")),
      age = factor(age, levels = c("Young", "Aged")),
      celltype_new = celltype,
      sample = factor(cell_sample)
    )
  )

  sce
}

sce_path <- system.file("extdata", "mini_ase.rds", package = "ASPENGLMS", mustWork = FALSE)
if (nzchar(sce_path) && file.exists(sce_path)) {
  sce <- readRDS(sce_path)
} else {
  message("mini_ase.rds not found; simulating demo data instead.")
  sce <- simulate_demo_sce()
}

fit_tbl <- ASPENGLMS::fit_glmm_bb(
  sce = sce,
  formula_fixed = ~ sex + age + celltype_new,
  rand = "(1|sample)",
  min_trials = 5,
  min_cells = 30,
  ncores = 1
)

if (!nrow(fit_tbl)) {
  stop("No model fits were produced. Check coverage thresholds or input data.", call. = FALSE)
}

sex_shrunk <- tryCatch({
  ASPENGLMS::shrink_with_ash(fit_tbl, term = "sexMale")
}, error = function(e) {
  warning("Shrinkage step failed: ", conditionMessage(e))
  NULL
})

contrasts <- ASPENGLMS::tidy_contrasts(fit_tbl)

if (!dir.exists("results")) {
  dir.create("results", recursive = TRUE)
}

utils::write.table(fit_tbl, file = file.path("results", "glmm_results.tsv"),
                   sep = "\t", quote = FALSE, row.names = FALSE)

if (!is.null(sex_shrunk)) {
  utils::write.table(sex_shrunk, file = file.path("results", "sex_shrinkage.tsv"),
                     sep = "\t", quote = FALSE, row.names = FALSE)
}

if (nrow(contrasts)) {
  utils::write.table(contrasts, file = file.path("results", "contrasts.tsv"),
                     sep = "\t", quote = FALSE, row.names = FALSE)
}

message("Demo pipeline completed. Results written to ./results.")
