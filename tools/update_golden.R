# Regenerate the golden files used by tests/testthat/test-golden-*.R.
#
# Run from the package root, with CONIPHER installed from this checkout and PyClone on
# PATH (e.g. inside the `conipher` conda env):
#
#   R CMD INSTALL . && Rscript tools/update_golden.R
#
# Only regenerate when a change in output is intended, and say why in the commit message:
# the golden files define "the algorithm works as expected".

suppressPackageStartupMessages(library(CONIPHER))
source("tests/testthat/helper-golden.R")

if (!pyclone_available()) stop("PyClone not found on PATH")

for (case_id in SYNTHETIC_CASES) {
  message("== ", case_id)
  work <- file.path(tempdir(), case_id)
  unlink(work, recursive = TRUE)

  run_clustering(case_id, file.path(work, "clustering"))
  dest <- golden_path(case_id, "clustering")
  dir.create(dest, recursive = TRUE, showWarnings = FALSE)
  for (f in CLUSTERING_GOLDEN_FILES) {
    stopifnot(file.copy(clustering_output_file(file.path(work, "clustering"), case_id, f),
                        file.path(dest, f), overwrite = TRUE))
  }

  # Trees are built from the golden clustering output, so tree tests run without PyClone.
  run_treebuilding(case_id, file.path(dest, "SCoutput.CLEAN.tsv"), file.path(work, "trees"))
  dest <- golden_path(case_id, "trees")
  dir.create(dest, recursive = TRUE, showWarnings = FALSE)
  for (f in TREE_GOLDEN_FILES) {
    stopifnot(file.copy(file.path(work, "trees", f), file.path(dest, f), overwrite = TRUE))
  }
}

writeLines(c(
  paste("CONIPHER", as.character(utils::packageVersion("CONIPHER"))),
  R.version.string,
  # PyClone prints its version to stderr.
  trimws(system2("PyClone", "--version", stdout = TRUE, stderr = TRUE)),
  paste("Commit", system2("git", c("rev-parse", "--short", "HEAD"), stdout = TRUE)),
  paste("Platform", R.version$platform)
), golden_path("VERSIONS.txt"))
message("Golden files written to ", golden_path())
