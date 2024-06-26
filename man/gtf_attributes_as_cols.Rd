% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/genomics_files_utils.R
\name{gtf_attributes_as_cols}
\alias{gtf_attributes_as_cols}
\title{Helper function to convert the last column list of a GTF file to a different separate columns}
\usage{
gtf_attributes_as_cols(data)
}
\arguments{
\item{data}{A data frame of a GTF file with the 9th column labelled 'attribute'}
}
\value{
A tibble
}
\description{
Helper function to convert the last column list of a GTF file to a different separate columns
}
\details{
This function works only if the "attributes" column was first converted to a list with \code{gtf_attributes_as_list()}
}
\examples{
gtf <- read_delim(file = path, delim = '\t', comment = '#!',
                  col_names = c('seqname', 'source', 'feature', 'start', 
                  'end', 'score', 'strand', 'frame', 'attribute'), 
                  col_types = c('i', 'c', 'c', 'i', 'i', 'c', 'c', 'c', 'c'),
                  n_max = 100) |>
                  gtf_attributes_as_list() |>
                  gtf_attributes_as_cols()
}
