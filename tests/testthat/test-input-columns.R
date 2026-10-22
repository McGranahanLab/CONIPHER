# Input tables are read by column name: column order must not matter, and the legacy
# COPY_NUMBER_A / COPY_NUMBER_B names must still work in place of COPY_NUMBER_MAJOR /
# COPY_NUMBER_MINOR. Outputs must match the golden files either way.

write_tsv <- function(df, path) {
  utils::write.table(df, file = path, sep = "\t", quote = FALSE, row.names = FALSE)
  invisible(path)
}

to_legacy_names <- function(df) {
  colnames(df)[colnames(df) == "COPY_NUMBER_MAJOR"] <- "COPY_NUMBER_A"
  colnames(df)[colnames(df) == "COPY_NUMBER_MINOR"] <- "COPY_NUMBER_B"
  df
}

test_that("legacy copy number columns are renamed", {
  df <- data.frame(SAMPLE = "S1", COPY_NUMBER_A = 2, COPY_NUMBER_B = 1)
  out <- standardise_copy_number_columns(df)
  expect_identical(colnames(out), c("SAMPLE", "COPY_NUMBER_MAJOR", "COPY_NUMBER_MINOR"))
  expect_identical(out$COPY_NUMBER_MAJOR, 2)
  expect_identical(out$COPY_NUMBER_MINOR, 1)
})

test_that("current copy number columns are left unchanged", {
  df <- data.frame(COPY_NUMBER_MINOR = 1, SAMPLE = "S1", COPY_NUMBER_MAJOR = 2)
  expect_identical(standardise_copy_number_columns(df), df)
})

test_that("current copy number columns win over legacy ones", {
  df <- data.frame(COPY_NUMBER_MAJOR = 3, COPY_NUMBER_A = 2, COPY_NUMBER_MINOR = 1)
  expect_warning(out <- standardise_copy_number_columns(df), "COPY_NUMBER_A")
  expect_identical(colnames(out), c("COPY_NUMBER_MAJOR", "COPY_NUMBER_MINOR"))
  expect_identical(out$COPY_NUMBER_MAJOR, 3)
})

test_that("seg position columns are found by name", {
  seg <- data.frame(ENDPOS = 200, COPY_NUMBER_MAJOR = 2, CHR = 1, SAMPLE = "S1", STARTPOS = 100)
  out <- standardise_seg_position_columns(seg)
  expect_identical(out$Chromosome, 1)
  expect_identical(out$StartPosition, 100)
  expect_identical(out$EndPosition, 200)

  internal <- data.frame(SampleID = "S1", chr = 1, startpos = 100, endpos = 200)
  expect_identical(standardise_seg_position_columns(internal)$StartPosition, 100)

  expect_error(standardise_seg_position_columns(seg[, c("CHR", "STARTPOS")]), "ENDPOS")
})

test_that("tree building is unaffected by column order and legacy names", {
  case_id <- "SYN001"
  input <- read_golden_table(golden_path(case_id, "clustering", "SCoutput.CLEAN.tsv"))
  input <- to_legacy_names(input)
  input <- input[, rev(colnames(input))]

  in_dir <- file.path(tempdir(), "input-columns-trees-in")
  dir.create(in_dir, showWarnings = FALSE, recursive = TRUE)
  input_tsv_loc <- write_tsv(input, file.path(in_dir, paste0(case_id, ".tsv")))

  out_dir <- file.path(tempdir(), "input-columns-trees", case_id)
  unlink(out_dir, recursive = TRUE)
  run_treebuilding(case_id, input_tsv_loc, out_dir)

  # treeTable.tsv echoes the input table, so it keeps the input column order
  tree_table <- read_golden_table(file.path(out_dir, "treeTable.tsv"))
  golden_cols <- colnames(read_golden_table(golden_path(case_id, "trees", "treeTable.tsv")))
  expect_setequal(colnames(tree_table), golden_cols)
  write_tsv(tree_table[, golden_cols], file.path(out_dir, "treeTable.tsv"))

  for (f in TREE_GOLDEN_FILES) {
    expect_matches_golden(file.path(out_dir, f), golden_path(case_id, "trees", f), f)
  }
})

test_that("clustering is unaffected by column order, new names and seg column order", {
  skip_on_cran()
  skip_if_not(pyclone_available(), "PyClone not on PATH")

  case_id <- "SYN001"
  input <- read_golden_table(fixture_path(case_id, "input.tsv"))
  input <- standardise_copy_number_columns(input)
  input <- input[, rev(colnames(input))]

  # one segment per sample and chromosome, columns deliberately out of the usual order
  seg <- do.call(rbind, lapply(split(input, list(input$SAMPLE, input$CHR), drop = TRUE), function(x) {
    data.frame(COPY_NUMBER_MINOR = x$COPY_NUMBER_MINOR[1], ENDPOS = max(x$POS),
               SAMPLE = x$SAMPLE[1], COPY_NUMBER_MAJOR = x$COPY_NUMBER_MAJOR[1],
               CHR = x$CHR[1], STARTPOS = min(x$POS))
  }))

  in_dir <- file.path(tempdir(), "input-columns-clustering-in")
  dir.create(in_dir, showWarnings = FALSE, recursive = TRUE)
  input_tsv_loc <- write_tsv(input, file.path(in_dir, paste0(case_id, ".tsv")))
  seg_tsv_loc   <- write_tsv(seg, file.path(in_dir, paste0(case_id, ".seg.tsv")))

  out_dir <- file.path(tempdir(), "input-columns-clustering", case_id)
  unlink(out_dir, recursive = TRUE)
  suppressMessages(utils::capture.output(
    conipher_clustering(case_id = case_id,
                        out_dir = out_dir,
                        input_tsv_loc = input_tsv_loc,
                        input_seg_tsv_loc = seg_tsv_loc,
                        pyclone_yaml_loc = normalizePath(fixture_path("pyclone_fast.yaml")),
                        burn_in = PYCLONE_BURN_IN,
                        seed = PYCLONE_SEED)
  ))

  for (f in CLUSTERING_GOLDEN_FILES) {
    expect_matches_golden(clustering_output_file(out_dir, case_id, f),
                          golden_path(case_id, "clustering", f), f)
  }
  expect_true(file.exists(file.path(out_dir, paste0(case_id, "_pyclone_cluster_assignment_copynumber_clean.pdf"))))
})
