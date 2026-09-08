# Use the same public Bioconductor mirror as setup-r-dependencies. The locked
# dependencies use only the software repository, not annotation or experiment data.
bioc_version <- renv::lockfile_read("renv.lock")$Bioconductor$Version
stopifnot(is.character(bioc_version), length(bioc_version) == 1L)
bioc_repository <- c(BioCsoft = paste0(
  "https://bioconductor.posit.co/packages/", bioc_version, "/bioc"))
options(BioC_mirror = "https://bioconductor.posit.co",
        renv.bioconductor.repos = bioc_repository)
cran_repositories <- getOption("repos")
cran_repositories <- cran_repositories[!grepl("^BioC", names(cran_repositories))]
options(repos = c(bioc_repository, cran_repositories))
