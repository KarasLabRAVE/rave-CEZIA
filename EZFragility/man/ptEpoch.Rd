% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.r
\docType{data}
\name{ptEpoch}
\alias{ptEpoch}
\title{Pt01 seizure 1 (-10:10s) around seizure onset}
\format{
A Matrix with 20001 rows (time points) and 84 columns (electrodes)
}
\source{
Fragility Multi-Center Retrospective Study
(\href{https://openneuro.org/datasets/ds003029/versions/1.0.0}{OpenNeuro})
}
\usage{
data(ptEpoch)
}
\description{
This data corresponds to the first seizure of patient PT01 from the Fragility Data Set.
The data contains only the good channels.
It has been notch filtered and common average referenced in RAVE
It has been epoched -10:10s around the seizure onset
The acquisition frequency is 1000 Hz
EcoG recording gathered in collaboration with the National Institue of Health
}
\examples{
data(ptEpoch)
}
\keyword{datasets}
