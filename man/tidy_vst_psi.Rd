% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/vst_utils.R
\name{tidy_vst_psi}
\alias{tidy_vst_psi}
\title{Reshape a wide vast-tools inclusion table into a tidyverse-friendly long table.}
\usage{
tidy_vst_psi(
  vst_psi_tbl,
  num_id_cols = 6,
  num_of_Score_Types = 5,
  quality_col_suffix = "-Q",
  return_quality_scores = TRUE,
  return_S1_only = TRUE,
  verbose = FALSE,
  add_ID_col = FALSE,
  col_ID_name = "banana"
)
}
\arguments{
\item{vst_psi_tbl}{A dataframe generated with \code{grep_psi()} from a vast-tools inclusion table.}

\item{num_id_cols}{How many first num_id_cols to consider as info/metadata/IDs of the AS event ID in the table.}

\item{num_of_Score_Types}{How many quality scores are present in the columns with headers ending with \code{quality_col_suffix}.}

\item{quality_col_suffix}{Suffix identifying the quality control columns. Default is "-Q".}

\item{return_quality_scores}{Logical. Do you want the individual scores to be returned in the output data.frame? Defaul \code{TRUE}.}

\item{return_S1_only}{Logical. Return only the first quality score, default \code{TRUE}, use \code{FALSE} for returning ALL \code{num_of_Score_Types} quality scores for each event.}

\item{verbose}{Print out information}

\item{add_ID_col}{Logical. Do you want to add an extra \code{col_ID_name} to the output data.frame? Default \code{FALSE}.}

\item{col_ID_name}{Extra column that could be used to add a new identifier to the data. Default is "banana".}
}
\value{
A reshaped \code{data.frame} in long format as \code{tibble}.
}
\description{
This function works with inclusion tables generated both by vast-tools combine or compare. In the latter case, if there's a 'dPSI' column it will be automatically detect and included in the reshaped table.
}
\details{
This function works well in conjunction with \code{grep_psi()} function.
}
\examples{
grep_psi(path_to_vst_PSI_tbl, vst_id = "HsaEX0000001") |>
    tidy_vst_psi() -> psi_tbl
}
