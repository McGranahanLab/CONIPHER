% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/main_clustering_functions.R
\name{clustering_postprocess}
\alias{clustering_postprocess}
\title{Postprocessing of clustering function}
\usage{
clustering_postprocess(
  input_list,
  sample.results,
  new.dir,
  input_tsv,
  input_seg_tsv_loc = NULL,
  min_cluster_size = 5,
  driver_cat = "1",
  clean_clusters = TRUE,
  min_ccf_present = 0.1,
  clonal_cutOff = 0.9,
  propClonal_threshold = 0.25
)
}
\arguments{
\item{input_list}{A list created by the clustering preprocess function
including patient id, regions to use, phylo region list and others.}

\item{sample.results}{which is the location of the pyclone output table.}

\item{new.dir}{A character specifying the directory where the pyclone
output should be saved.}

\item{input_tsv}{the input mutation tsv.}

\item{input_seg_tsv_loc}{path to a copy number segment tsv file that is used for
across genome copy number plotting.
Default NULL}

\item{min_cluster_size}{Minimum number of mutations needed for a cluster to be considered.
Default 5}

\item{driver_cat}{Which categories to use as driver mutations
Default "1"}

\item{clean_clusters}{should clusters be cleaned and merged?
Default TRUE}

\item{min_ccf_present}{minimum CCF to consider a mutation as present. 
Default 0.1}

\item{clonal_cutOff}{lower threshold CCF to consider mutations as clonal.
Default 0.9}

\item{propClonal_threshold}{Proportion of mutations in cluster which needs to be
considered clonal to merge.
Default 0.25}
}
\description{
This function takes the input tsv and formats the data to be compatible with
the main CONIPHER clustering function.
}
