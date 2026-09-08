#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Matrix)
  library(SingleCellExperiment)
})

if (!requireNamespace("ASPENGLMS", quietly = TRUE)) {
  stop("Please install/load the ASPENGLMS package before running this script.", call. = FALSE)
}

# helper to simulate a tiny dataset when bundled demo data is absent
simulate_demo_sce <- function(n_genes = 24L, n_cells = 480L, n_samples = 24L) {
  set.seed(42)
  sample_index <- rep(seq_len(n_samples), length.out = n_cells)
  donor_sex <- rep(c("Female", "Male"), length.out = n_samples)
  donor_age <- rep(rep(c("Young", "Aged"), each = 2), length.out = n_samples)
  sex <- factor(donor_sex[sample_index], levels = c("Female", "Male"))
  age <- factor(donor_age[sample_index], levels = c("Young", "Aged"))
  celltype <- factor(rep(c("T", "B", "Mono"), length.out = n_cells))
  total <- matrix(rnbinom(n_genes * n_cells, size = 20, mu = 40), nrow = n_genes)
  a1 <- matrix(0, nrow = n_genes, ncol = n_cells)
  for (gene in seq_len(n_genes)) {
    donor_effect <- rnorm(n_samples, sd = 0.6)
    eta <- -0.2 + 0.5 * (sex == "Male") - 0.2 * (age == "Aged") +
      0.15 * (celltype == "B") + donor_effect[sample_index]
    mean_prob <- plogis(eta)
    prob <- rbeta(n_cells, mean_prob * 20, (1 - mean_prob) * 20)
    a1[gene, ] <- rbinom(n_cells, total[gene, ], prob)
  }
  dimnames(total) <- dimnames(a1) <- list(paste0("gene", seq_len(n_genes)),
                                         paste0("cell", seq_len(n_cells)))
  SingleCellExperiment::SingleCellExperiment(
    assays = list(a1 = a1, a2 = total - a1, tot = total),
    colData = S4Vectors::DataFrame(sex = sex, age = age, celltype_new = celltype,
      sample = factor(paste0("sample", sample_index)))
  )
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

if (!nrow(fit_tbl) || !any(fit_tbl$converged)) {
  stop("No model fits were produced. Check coverage thresholds or input data.", call. = FALSE)
}

sex_shrunk <- ASPENGLMS::shrink_with_ash(fit_tbl, term = "sexMale")

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

if (requireNamespace("ggplot2", quietly = TRUE)) {
  plot <- ASPENGLMS::plot_coefficient_diagnostics(fit_tbl, "sexMale")
  ggplot2::ggsave(file.path("results", "coefficient-diagnostics.png"), plot,
                 width = 7, height = 4.5, dpi = 140)
}
