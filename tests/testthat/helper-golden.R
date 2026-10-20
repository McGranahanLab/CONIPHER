# Shared helpers for the golden-file tests. testthat sources helper-*.R before the
# tests; tools/update_golden.R sources this file too, so tests and golden files are
# always produced by the same calls.

SYNTHETIC_CASES <- c("SYN001", "SYN002", "SYN003")

# Fixed seeds: PyClone sampling (passed to `PyClone run_analysis --seed`) and the R RNG
# used by bootstrapped cluster CIs during tree building.
PYCLONE_SEED <- 1024
TREE_SEED    <- 1024

# Test-only PyClone settings: 2000 MCMC iterations (default 10000) keeps a case at ~1 min.
PYCLONE_BURN_IN <- 500

CLUSTERING_GOLDEN_FILES <- c("SCoutput.CLEAN.tsv", "SCoutput.FULL.tsv", "results.tsv")
TREE_GOLDEN_FILES <- c("allTrees.txt", "consensusBranches.txt", "consensusRelationships.txt",
                       "clusterInfo.txt", "treeTable.tsv", "alternativeTreeMetrics.txt",
                       "cloneProportionsMinErrorTrees.txt",
                       "subclonalExpansionScoreMinErrorTrees.txt")

# Inside a test run paths are relative to tests/testthat; tools/update_golden.R runs
# from the package root.
testthat_dir <- function() {
  if (testthat::is_testing()) testthat::test_path() else "tests/testthat"
}
fixture_path <- function(...) file.path(testthat_dir(), "fixtures", ...)
golden_path  <- function(...) file.path(testthat_dir(), "golden", ...)

pyclone_available <- function() nzchar(Sys.which("PyClone"))

run_clustering <- function(case_id, out_dir) {
  suppressMessages(utils::capture.output(
    conipher_clustering(case_id = case_id,
                        out_dir = out_dir,
                        input_tsv_loc = fixture_path(case_id, "input.tsv"),
                        pyclone_yaml_loc = normalizePath(fixture_path("pyclone_fast.yaml")),
                        burn_in = PYCLONE_BURN_IN,
                        seed = PYCLONE_SEED)
  ))
  invisible(out_dir)
}

run_treebuilding <- function(case_id, input_tsv_loc, out_dir) {
  set.seed(TREE_SEED)
  suppressMessages(utils::capture.output(
    conipher_treebuilding(input_tsv_loc = input_tsv_loc, out_dir = out_dir, prefix = "SYN")
  ))
  invisible(out_dir)
}

clustering_output_file <- function(out_dir, case_id, file) {
  file.path(out_dir, paste0(case_id, ".", file))
}

# Tables are compared column by column: numbers with a relative tolerance (PyClone and
# boot results can differ in the last digits across platforms/BLAS), everything else
# exactly. Files that are not plain tables (allTrees.txt, ...) are compared line by line.
read_golden_table <- function(path) {
  utils::read.delim(path, sep = "\t", stringsAsFactors = FALSE, quote = "",
                    check.names = FALSE, colClasses = NA)
}

is_plain_table <- function(file) {
  !file %in% c("allTrees.txt", "consensusBranches.txt", "consensusRelationships.txt")
}

expect_matches_golden <- function(actual, golden, file, tolerance = 1e-6) {
  testthat::expect_true(file.exists(actual), info = paste("missing output:", file))
  if (!file.exists(actual)) return(invisible())
  if (is_plain_table(file)) {
    a <- read_golden_table(actual)
    g <- read_golden_table(golden)
    testthat::expect_identical(colnames(a), colnames(g), info = file)
    testthat::expect_identical(nrow(a), nrow(g), info = file)
    for (col in intersect(colnames(a), colnames(g))) {
      if (is.numeric(g[[col]])) {
        testthat::expect_equal(a[[col]], g[[col]], tolerance = tolerance,
                               info = paste0(file, "$", col))
      } else {
        testthat::expect_identical(a[[col]], g[[col]], info = paste0(file, "$", col))
      }
    }
  } else {
    testthat::expect_identical(readLines(actual), readLines(golden), info = file)
  }
}

# --- Truth checks -------------------------------------------------------------------

# Edges of the default (first) tree in allTrees.txt, as "parent->child" strings.
default_tree_edges <- function(all_trees_file) {
  lines <- readLines(all_trees_file)
  starts <- grep("^# tree ", lines)
  end <- if (length(starts) > 1) starts[2] - 1 else length(lines)
  edges <- lines[(starts[1] + 1):end]
  edges <- edges[nzchar(edges)]
  vapply(strsplit(edges, "\t"), function(e) paste0(e[1], "->", e[2]), character(1))
}

# Map each inferred cluster to the true clone holding most of its mutations.
cluster_to_clone <- function(tree_table_file, truth_mutations_file) {
  tree  <- read_golden_table(tree_table_file)
  truth <- read_golden_table(truth_mutations_file)
  tree  <- unique(tree[, c("CHR", "POS", "treeCLUSTER")])
  m     <- merge(tree, truth[, c("CHR", "POS", "CLONE")], by = c("CHR", "POS"))
  votes <- table(m$treeCLUSTER, m$CLONE)
  stats::setNames(colnames(votes)[apply(votes, 1, which.max)], rownames(votes))
}

# True edges as "parent->child" strings over clone names.
truth_edges <- function(truth_tree_file) {
  t <- read_golden_table(truth_tree_file)
  t <- t[t$PARENT != "diploid", ]
  paste0(t$PARENT, "->", t$CLONE)
}

# Inferred default-tree edges translated into true clone names.
inferred_edges_as_clones <- function(tree_dir, case_id) {
  mapping <- cluster_to_clone(file.path(tree_dir, "treeTable.tsv"),
                              fixture_path(case_id, "truth_mutations.tsv"))
  edges <- strsplit(default_tree_edges(file.path(tree_dir, "allTrees.txt")), "->")
  vapply(edges, function(e) paste0(mapping[[e[1]]], "->", mapping[[e[2]]]), character(1))
}
