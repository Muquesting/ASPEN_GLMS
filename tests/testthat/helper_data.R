toy_sce <- function() {
  set.seed(1)
  n_genes <- 5L
  n_cells <- 120L
  n_donors <- 12L
  donor <- rep(seq_len(n_donors), length.out = n_cells)
  sex <- factor(rep(c("Female", "Male"), length.out = n_donors)[donor],
                levels = c("Female", "Male"))
  age <- factor(rep(c("Young", "Young", "Aged", "Aged"), length.out = n_donors)[donor])
  tot <- matrix(rnbinom(n_genes * n_cells, mu = 30, size = 10), nrow = n_genes)
  a1 <- matrix(0, nrow = n_genes, ncol = n_cells)
  for (gene in seq_len(n_genes)) {
    eta <- 0.4 * (sex == "Male") + 0.2 * (age == "Aged") + rnorm(n_donors, sd = 0.7)[donor]
    mu <- plogis(eta)
    prob <- rbeta(n_cells, mu * 12, (1 - mu) * 12)
    a1[gene, ] <- rbinom(n_cells, tot[gene, ], prob)
  }
  dimnames(tot) <- dimnames(a1) <- list(paste0("gene", seq_len(n_genes)), paste0("cell", seq_len(n_cells)))
  SingleCellExperiment::SingleCellExperiment(
    assays = list(a1 = a1, a2 = tot - a1, tot = tot),
    colData = S4Vectors::DataFrame(sex = sex, age = age,
      celltype_new = sample(c("T", "B"), n_cells, replace = TRUE),
      sample = factor(paste0("sample", donor)))
  )
}

.real_aspen_cache <- new.env(parent = emptyenv())

real_sce_subset <- function(
  path = testthat::test_path("..", "data", "aspensce_sexupdated.rds"),
  min_trials = 10,
  min_cells = 200,
  max_genes = 12
) {
  stopifnot(file.exists(path))
  sce <- readRDS(path)

  cd <- as.data.frame(SummarizedExperiment::colData(sce))
  if (!"age" %in% names(cd)) {
    cd$age <- cd$condition_new
  }
  if (!"sample" %in% names(cd)) {
    fallback <- if ("orig.ident" %in% names(cd)) cd$orig.ident else seq_len(ncol(sce))
    cd$sample <- fallback
  }

  cd$sex <- factor(cd$sex)
  cd$condition_new <- factor(cd$condition_new)
  cd$celltype_new <- factor(cd$celltype_new)
  cd$age <- factor(cd$age)
  cd$sample <- factor(cd$sample)

  SummarizedExperiment::colData(sce) <- S4Vectors::DataFrame(cd)

  totals <- SummarizedExperiment::assay(sce, "tot")
  keep_counts <- Matrix::rowSums(totals >= min_trials)
  top_genes <- head(order(keep_counts, decreasing = TRUE), max_genes)

  sce[top_genes, , drop = FALSE]
}

real_aspen_fit <- function(
  formula_fixed = ~ sex + condition_new + celltype_new,
  min_trials = 10,
  min_cells = 200,
  max_genes = 12
) {
  key <- paste(deparse(formula_fixed), min_trials, min_cells, max_genes, sep = "|")
  if (!exists(key, envir = .real_aspen_cache, inherits = FALSE)) {
    sce <- real_sce_subset(
      min_trials = min_trials,
      min_cells = min_cells,
      max_genes = max_genes
    )
    fit <- suppressWarnings(fit_glmm_bb(
      sce,
      formula_fixed = formula_fixed,
      min_trials = min_trials,
      min_cells = min_cells,
      ncores = 1
    ))
    assign(key, list(sce = sce, fit = fit), envir = .real_aspen_cache)
  }
  get(key, envir = .real_aspen_cache, inherits = FALSE)
}
