% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/epiprofile_utils.R
\name{read_EpiProfile_histone_ratio}
\alias{read_EpiProfile_histone_ratio}
\title{Import EpiProfile version2.1 output file \code{histone_ratio.xls} into R}
\usage{
read_EpiProfile_histone_ratio(histone_ratio_path)
}
\arguments{
\item{histone_ratio_path}{Path to where the xls file is}
}
\value{
A tidy tibble
}
\description{
Import EpiProfile version2.1 output file \code{histone_ratio.xls} into R
}
\examples{
read_EpiProfile_histone_ratio("path/to/histone_ratio.xls") |>
    tidy_hPTMs()   
}
