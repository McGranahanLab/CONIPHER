#' Example input table
#'
#' Example input table to input to CONIPHER tree building. The input table should
#' have the following columns: CASE_ID, SAMPLE, CHR, POS, REF,
#' ALT, REF_COUNT, VAR_COUNT, DEPTH, CLUSTER, CCF_PHYLO, CCF_OBS,
#' MUT_COPY, COPY_NUMBER_A, COPY_NUMBER_B, ACF, and PLOIDY. Full description of
#' the input format can be found in our companion manuscript.
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
