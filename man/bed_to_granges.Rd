% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bed_to_granges.R
\name{bed_to_granges}
\alias{bed_to_granges}
\title{BED to GRanges}
\usage{
bed_to_granges(file)
}
\arguments{
\item{file}{Location of your file}
}
\description{
This function loads a BED-like file and stores it as a GRanges object.
The tab-delimited file must be ordered as 'chr', 'start', 'end', 'id', 'score', 'strand'.
The minimal BED file must have the 'chr', 'start', 'end' columns.
Any columns after the strand column are ignored.
}
\keyword{BED}
\keyword{GRanges}
