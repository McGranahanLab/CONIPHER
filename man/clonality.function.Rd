% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/treebuilding_helper_functions.R
\name{clonality.function}
\alias{clonality.function}
\title{Function to determine cluster clonality in each tumour region}
\usage{
clonality.function(
  pyclone,
  trunk,
  ccf_buffer = 10,
  prefix = "LTX",
  min_cluster_size = 5,
  pval_cutoff = 0.01,
  min_ccf = 0.05,
  use_boot = FALSE
)
}
\arguments{
\item{pyclone}{An R list object containing information about the PhyloCCF
of each mutation in each tumour region.}

\item{trunk}{truncal cluster name}

\item{ccf_buffer}{PhyloCCF buffer permitted when checking tree level issue}

\item{min_cluster_size}{Threshold for minimum number of mutations required in
a mutation cluster}

\item{pval_cutoff}{A p-value significance threshold for testing whether
clusters can be nested. (i.e. a p-value < pval_cutoff is significant)}

\item{min_ccf}{Minimum threshold for cluster PhyloCCF allowed to be classified
as present}

\item{use_boot}{Whether to use bootstrapping to determine confidence intervals
for each mutation cluster}
}
\value{
a matrix of dimensions (n_clusters x n_regions) classifying each
cluster as 'clonal', 'subclonal', or 'absent' in each tumour region.
}
\description{
This function takes as input an R list containing information about PhyloCCF
of each mutation (pyclone) and computes confidence intervals for each mutation
cluster, of the PhyloCCF distributions of the mutations in that cluster. If
use_boot==TRUE, then confidence intervals are computed using bootstrapping.
The function then performs a statistical test (Wilcoxon) for every pair of
clusters to determine whether one cluster can be nested within another.
}
