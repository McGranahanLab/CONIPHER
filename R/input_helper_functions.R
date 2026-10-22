################################################################################
##################          INPUT COLUMN HELPER FUNCTIONS          #############
################################################################################

# Legacy copy number column names, kept for backwards compatibility.
# Names are the current column names, values the legacy ones.
LEGACY_COPY_NUMBER_COLUMNS <- c(COPY_NUMBER_MAJOR = "COPY_NUMBER_A",
                                COPY_NUMBER_MINOR = "COPY_NUMBER_B")

#' Standardise copy number column names
#'
#' The raw copy number columns are called COPY_NUMBER_MAJOR and
#' COPY_NUMBER_MINOR. For backwards compatibility, tables using the legacy
#' names COPY_NUMBER_A and COPY_NUMBER_B are accepted and renamed. If both the
#' legacy and the current name are present, the current one is used.
#' @param input_table A dataframe, e.g. the input mutation table or a seg table.
#' @returns The dataframe with copy number columns using the current names.

standardise_copy_number_columns <- function(input_table) {
  for (new_col in names(LEGACY_COPY_NUMBER_COLUMNS)) {
    old_col <- LEGACY_COPY_NUMBER_COLUMNS[[new_col]]
    if (!old_col %in% colnames(input_table)) next
    if (new_col %in% colnames(input_table)) {
      warning(paste0("Both ", new_col, " and legacy ", old_col, " columns found, using ", new_col, " and ignoring ", old_col))
      input_table[[old_col]] <- NULL
    } else {
      colnames(input_table)[colnames(input_table) == old_col] <- new_col
    }
  }
  return(input_table)
}

#' Standardise segment position column names for plotting
#'
#' Copy number plotting expects Chromosome, StartPosition and EndPosition
#' columns. This selects them by name (from either the internal chr/startpos/
#' endpos or the seg file CHR/STARTPOS/ENDPOS), so column order does not matter.
#' @param seg.mat A dataframe of copy number segments.
#' @returns The dataframe with Chromosome, StartPosition and EndPosition columns.

standardise_seg_position_columns <- function(seg.mat) {
  position_cols <- list(Chromosome    = c("chr", "CHR"),
                        StartPosition = c("startpos", "STARTPOS"),
                        EndPosition   = c("endpos", "ENDPOS"))
  for (new_col in names(position_cols)) {
    if (new_col %in% colnames(seg.mat)) next
    old_col <- intersect(position_cols[[new_col]], colnames(seg.mat))
    if (length(old_col) == 0) {
      stop(paste0("Copy number segment table requires a ", paste(position_cols[[new_col]], collapse = " or "), " column"))
    }
    colnames(seg.mat)[colnames(seg.mat) == old_col[1]] <- new_col
  }
  return(seg.mat)
}
