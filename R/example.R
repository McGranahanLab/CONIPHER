#' Example input table
#'
#' Example input table to input to CONIPHER tree building. The input table should
#' have the following columns: CASE_ID, SAMPLE, CHR, POS, REF,
#' ALT, REF_COUNT, VAR_COUNT, DEPTH, CLUSTER, CCF_PHYLO, CCF_OBS,
#' MUT_COPY, COPY_NUMBER_MAJOR, COPY_NUMBER_MINOR, ACF, and PLOIDY. Full description of
#' the input format can be found in our companion manuscript. Columns are matched
#' by name, so their order does not matter. The legacy copy number column names
#' COPY_NUMBER_A and COPY_NUMBER_B are still accepted in place of COPY_NUMBER_MAJOR
#' and COPY_NUMBER_MINOR.
#'
#' @docType data
#'
#' @usage data(input_table)
#'
#' @format An object of class \code{"data.frame"}
#'
#' @keywords datasets
#'
#' @examples
#' data(input_table)
#' head(input_table)
"input_table"
