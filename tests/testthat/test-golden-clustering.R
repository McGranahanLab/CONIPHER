# Clustering (PyClone) against golden files. Needs PyClone on PATH (the `conipher` conda
# env) and takes ~1 min per case, so it is skipped on CRAN and when PyClone is missing.
# PyClone runs with a fixed seed, so results are reproducible on the same platform.

for (case_id in SYNTHETIC_CASES) {
  test_that(paste(case_id, "clustering matches golden files"), {
    skip_on_cran()
    skip_if_not(pyclone_available(), "PyClone not on PATH")

    out_dir <- file.path(tempdir(), "golden-clustering", case_id)
    unlink(out_dir, recursive = TRUE)
    run_clustering(case_id, out_dir)

    for (f in CLUSTERING_GOLDEN_FILES) {
      expect_matches_golden(clustering_output_file(out_dir, case_id, f),
                            golden_path(case_id, "clustering", f), f)
    }
  })
}
