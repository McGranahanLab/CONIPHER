% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/example.R
\docType{data}
\name{input_table}
\alias{input_table}
\title{Example input table}
\format{
An object of class \code{"data.frame"}
}
\usage{
data(input_table)
}
\description{
Example input table to input to CONIPHER tree building. The input table should
have the following columns: CASE_ID, SAMPLE, CHR, POS, REF,
ALT, REF_COUNT, VAR_COUNT, DEPTH, CLUSTER, CCF_PHYLO, CCF_OBS,
MUT_COPY, COPY_NUMBER_A, COPY_NUMBER_B, ACF, and PLOIDY. Full description of
the input format can be found in our companion manuscript.
}
\examples{
data(input_table)
head(input_table)
}
\keyword{datasets}
