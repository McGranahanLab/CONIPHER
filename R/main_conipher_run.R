#' Full CONIPHER run
#'
#' This function takes all the input options and runs the three main steps:
#' preprocess, tree building run and postprocess
#' @param opt a list of options
#' @returns NULL
#' @export conipher_run

conipher_run <- function(case_id, prefix, out_dir, input_tsv_loc, ...) {
    conipher_clustering(case_id = case_id, 
                        out_dir = out_dir, 
                        input_tsv_loc = input_tsv_loc, 
                        ...)
    tree_input_tsv_loc <- paste0(out_dir, "/Clustering/", case_id, ".SCoutput.CLEAN.tsv")
    conipher_treebuilding(input_tsv_loc = tree_input_tsv_loc,
                          out_dir = out_dir,
                          prefix = prefix,
                          ...)
}
