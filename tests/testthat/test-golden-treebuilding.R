# Tree building from the golden clustering output. Pure R, no PyClone needed, so these
# run everywhere. Regenerate golden files with tools/update_golden.R.

for (case_id in SYNTHETIC_CASES) {
  out_dir <- file.path(tempdir(), "golden-trees", case_id)
  unlink(out_dir, recursive = TRUE)
  run_treebuilding(case_id, golden_path(case_id, "clustering", "SCoutput.CLEAN.tsv"), out_dir)

  test_that(paste(case_id, "tree building matches golden files"), {
    for (f in TREE_GOLDEN_FILES) {
      expect_matches_golden(file.path(out_dir, f), golden_path(case_id, "trees", f), f)
    }
  })

  test_that(paste(case_id, "default tree recovers the true topology"), {
    inferred <- inferred_edges_as_clones(out_dir, case_id)
    expect_setequal(inferred, truth_edges(fixture_path(case_id, "truth_tree.tsv")))
  })

  test_that(paste(case_id, "every true clone maps to exactly one tree cluster"), {
    mapping <- cluster_to_clone(file.path(out_dir, "treeTable.tsv"),
                                fixture_path(case_id, "truth_mutations.tsv"))
    truth <- read_golden_table(fixture_path(case_id, "truth_tree.tsv"))
    expect_setequal(unname(mapping), truth$CLONE)
    expect_false(anyDuplicated(unname(mapping)) > 0)
  })
}
