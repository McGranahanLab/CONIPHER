library(tidyverse)
library(CONIPHER)
library(coin)

source("/nemo/project/proj-tracerx-wgs/working/tholk/CONIPHER/R/main_clustering_functions.R")
source("/nemo/project/proj-tracerx-wgs/working/tholk/CONIPHER/R/functionsForSimpleClustering.v13.R")

cmdArgs           <- commandArgs(trailingOnly = TRUE)
patient           <- cmdArgs[1]
input_tsv_loc     <- cmdArgs[2]
input_seg_tsv_loc <- cmdArgs[3]
new.dir           <- cmdArgs[4]

#patient           <- "CRUKTOY001"
#input_tsv_loc     <- "/nemo/project/proj-tracerx-wgs/working/tholk/CONIPHER/inst/extdata/input_table.tsv"
# input_seg_tsv_loc <- "/nemo/project/proj-tracerx-lung/tctProjects/mcgranahanLab/tholk/TRACERX-wgs/tree_test_input/LTX0038/LTX0038_tree_input_seg.tsv"
#new.dir           <- "/nemo/project/proj-tracerx-wgs/working/tholk/CONIPHER/testing/results/"

subclonal_copy_correction              <- TRUE
only_truncal_subclonal_copy_correction <- TRUE
multiple_test_correction               <- TRUE
fix_absentCCFs                         <- TRUE

min_cluster_size                       <- 5
clean_clusters                         <- TRUE
clonal_cutOff                          <- 0.9
propClonal_threshold                   <- 0.25
driver_cat                             <- unlist(strsplit("1A,1,2A", split = ","))


cat("\nCONIPHER clustering analysis of the following tumour case:\n")
print(patient)
cat("\n")

# if (!file.exists(new.dir)) {
#     if (!dir.create(new.dir, recursive = TRUE)) {
#         stop("Unable to create root directory.\n")
#     }
# }

# if(!file.exists(input_tsv_loc)) {
#     stop("Unable to find input_tsv.\n")
# }
# input_tsv     <- read.delim(input_tsv_loc, sep = "\t", stringsAsFactors = FALSE, header = TRUE, fill = TRUE, quote = "")
# if (nrow(input_tsv) == 0) {
#     stop('No mutations passed filtering, stopping PyClone phylo clustering')
# }

# cat("\nRunning clustering preprocess\n")
# input_list     <- clustering_preprocess(input_tsv, new.dir = new.dir, subclonal_copy_correction = subclonal_copy_correction, multiple_test_correction = multiple_test_correction, only_truncal_subclonal_copy_correction = only_truncal_subclonal_copy_correction, fix_absentCCFs = fix_absentCCFs)

cat("\nLoading previous preprocess\n")
load(paste0(new.dir, patient, ".PyClone.RData"))

cat("\nRunning clustering\n")
sample.results <- clustering_run(input_list, nProcs = nProcs, new.dir = new.dir, burn_in = burn_in, pyclone_seed = seed, template.config.yaml = template.config.yaml)

cat("\nRunning clustering postprocess\n")
input_tsv$CHR <- as.numeric(input_tsv$CHR)
input_tsv$POS <- as.numeric(input_tsv$POS)
clustering_postprocess(input_list, sample.results, new.dir = new.dir, input_tsv = input_tsv, input_seg_tsv_loc = input_seg_tsv_loc, min_cluster_size = min_cluster_size, driver_cat = driver_cat, clean_clusters = clean_clusters, min_ccf_present = 0.1, clonal_cutOff = clonal_cutOff, propClonal_threshold = propClonal_threshold)


cat("\nRunning tree building\n")
conipher_treebuilding(input_tsv_loc = paste0(new.dir, patient, ".SCoutput.CLEAN.tsv"), out_dir = new.dir, prefix = "CRUK", multi_trees = FALSE)
