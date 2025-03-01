% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fragility.r
\name{fragStat}
\alias{fragStat}
\title{Compute quantiles, mean and standard deviation for two electrodes group marked as soz non marked as soz}
\usage{
fragStat(frag, sozID)
}
\arguments{
\item{frag}{Matrix or Fragility object. Either a matrix with row as Electrode names and Column as fragility index, or a Fragility object from \code{calc_adj_frag}}

\item{sozID}{Integer.  Vector soz electrodes (for good electrodes)}
}
\value{
list of 5 items with quantile matrix, mean and sdv from both electrodes groups
}
\description{
Compute quantiles, mean and standard deviation for two electrodes group marked as soz non marked as soz
}
\examples{
data("pt01Frag")
data("pt01Epoch")
sozindex<-attr(pt01Epoch,"sozindex")
pt01fragstat<-fragStat(frag=pt01Frag, sozID=sozindex)
}
