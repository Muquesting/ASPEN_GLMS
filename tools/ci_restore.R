# Restore into a separate library and verify in a fresh R process.
library_path <- file.path(Sys.getenv("RUNNER_TEMP", tempdir()), "aspen-restored-library")
dir.create(library_path, recursive = TRUE, showWarnings = FALSE)
lockfile <- "renv.lock"
renv::restore(lockfile = lockfile, library = library_path, prompt = FALSE)
status <- system2(file.path(R.home("bin"), "Rscript"),
  c("--vanilla", "tools/verify_restored.R", shQuote(library_path)))
if (status != 0L) stop("Restored-library verification failed.")
