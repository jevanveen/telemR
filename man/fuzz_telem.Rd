% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/telemR.R
\name{fuzz_telem}
\alias{fuzz_telem}
\title{average telemetry data over user defined windows. very useful when data appear very noisy/oversampled}
\usage{
fuzz_telem(tidy_telem, output_window, sum_counts = F)
}
\arguments{
\item{tidy_telem}{a tidy telemetry tibble, as produced by read_starr or read_oddi}

\item{output_window}{what should final output increment be?}

\item{sum_counts}{logical. should counts be summed rather than averaged?}
}
\value{
a time fuzzed tidy telemetry tibble
}
\description{
average telemetry data over user defined windows. very useful when data appear very noisy/oversampled
}
