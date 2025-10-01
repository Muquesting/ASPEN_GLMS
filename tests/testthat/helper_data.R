toy_sce <- function() {
  set.seed(1)
  n_genes <- 5L
  n_cells <- 40L
  tot <- matrix(rnbinom(n_genes * n_cells, mu = 15, size = 8), nrow = n_genes)
  prob <- matrix(runif(n_genes * n_cells, min = 0.35, max = 0.65), nrow = n_genes)
  a1 <- matrix(rbinom(n_genes * n_cells, size = as.vector(tot), prob = as.vector(prob)),
               nrow = n_genes)
  a2 <- tot - a1

  colData <- data.frame(
    sex = factor(rep(c("Female", "Male"), length.out = n_cells), levels = c("Female", "Male")),
    age = factor(rep(c("Young", "Aged"), each = n_cells / 2), levels = c("Young", "Aged")),
    celltype_new = sample(c("T", "B"), n_cells, replace = TRUE),
    sample = factor(rep(paste0("sample", 1:4), length.out = n_cells))
  )

  SingleCellExperiment::SingleCellExperiment(
    assays = list(a1 = a1, a2 = a2, tot = tot),
    colData = colData
  )
}
