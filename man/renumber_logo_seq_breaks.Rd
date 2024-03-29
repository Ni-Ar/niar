% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pwm_utils.R
\name{renumber_logo_seq_breaks}
\alias{renumber_logo_seq_breaks}
\title{Given a sequence of upper and lower cases creates the best breaks and labels for a ggplot sequence logo}
\usage{
renumber_logo_seq_breaks(
  seq,
  lowercase_spacer = 5,
  uppercase_spacer = 10,
  uppercase_adjustment
)
}
\arguments{
\item{seq}{A character vector with UPPER and lower cases}

\item{lowercase_spacer}{An integer indicating the distance between the label numbers in the lowercase sequence.}

\item{uppercase_spacer}{An integer indicating the distance between the label numbers in the uppercase sequence.}

\item{uppercase_adjustment}{An integer defining the actual length of the region in the coord2seq.sh script.
If the upper-case region has been trimmed upstream (e.g. -u -39) this parameter
corrects the X-axis numbering to match the real length before trimming. Can be omitted.}
}
\value{
A list the best breaks and labels
}
\description{
This is a helper function used in \code{plot_bits_logo()} and \code{plot_JSD_logo()}.
}
\note{
If proving a non-integer number for \verb{*case_spacer} parameters it will be rounded up to the closer integer.
}
