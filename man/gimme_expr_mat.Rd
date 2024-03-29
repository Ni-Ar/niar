% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/vst_utils.R
\name{gimme_expr_mat}
\alias{gimme_expr_mat}
\title{Extract the gene expression counts from an vast-tools expression table and return them as a matrix}
\usage{
gimme_expr_mat(vst_expression_tbl, min_mean_count = 5, verbose = FALSE)
}
\arguments{
\item{vst_expression_tbl}{Path to a vast-tools expression table (either cRPKM or TPM).}

\item{min_mean_count}{Filter out low expressed genes in the table read from \code{inclusion_tbl}. Defines the minimum row mean expression value across all samples that a gene must have to be selected.}

\item{verbose}{Print out extra info.}
}
\value{
A matrix
}
\description{
The matrix contains the samples as row names and the gene ID as column name. The gene ID used are the one present in the first column of the table (usually ENSEMBL gene ID).
}
\examples{
gimme_expr_mat(vst_expression_tbl = /path/to/expr/table/file.tab)
}
