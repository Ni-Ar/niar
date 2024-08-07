% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils_rmats.R
\name{multiply_csv}
\alias{multiply_csv}
\title{Multiply comma-separated numbers in a character vector by a number, handling NAs}
\usage{
multiply_csv(x, k = 100)
}
\arguments{
\item{x}{A character vector containing comma-separated numbers.}

\item{k}{The number to multiply by. Default \code{100}.}
}
\value{
A character vector with the same number of elements as \code{x},
where each element is the product of the original number and \code{k},
or NA if the original element could not be converted to numeric.
}
\description{
This function takes a character vector containing comma-separated numbers and
multiplies each number by 100. It handles missing values (NAs) by returning \code{NA}
for elements that cannot be converted to numeric.
}
\examples{
x <- c("0.1,0.2,NA,0.3", "0.4,0.5,NA,0.66")
multiply_csv(x)
# "10,20,NA,30" "40,50,NA,66"
y <- c("0.2,!!,0.9,NA", "*,NA,0.015")
multiply_csv(y)
# "20,NA,90,NA"  "NA,NA,1.5,NA"
}
