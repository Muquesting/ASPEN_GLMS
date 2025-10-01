toy_sce <- function() {
  set.seed(1)
  n_genes <- 5L
  n_cells <- 40L
  genes <- paste0("gene", seq_len(n_genes))
  cells <- paste0("cell", seq_len(n_cells))
  tot <- matrix(rnbinom(n_genes * n_cells, mu = 15, size = 8), nrow = n_genes)
  prob <- matrix(runif(n_genes * n_cells, min = 0.35, max = 0.65), nrow = n_genes)
  a1 <- matrix(rbinom(n_genes * n_cells, size = as.vector(tot), prob = as.vector(prob)),
               nrow = n_genes)
  a2 <- tot - a1

  dimnames(tot) <- list(genes, cells)
  dimnames(a1) <- list(genes, cells)
  dimnames(a2) <- list(genes, cells)

  colData <- data.frame(
    sex = factor(rep(c("Female", "Male"), length.out = n_cells), levels = c("Female", "Male")),
    age = factor(rep(c("Young", "Aged"), each = n_cells / 2), levels = c("Young", "Aged")),
    celltype_new = sample(c("T", "B"), n_cells, replace = TRUE),
    sample = factor(rep(paste0("sample", 1:4), length.out = n_cells)),
    row.names = cells
  )

  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(a1 = a1, a2 = a2, tot = tot),
    colData = colData
  )

  SummarizedExperiment::rowData(sce)$gene_id <- genes
  rownames(sce) <- genes
  colnames(sce) <- cells

  sce
}
