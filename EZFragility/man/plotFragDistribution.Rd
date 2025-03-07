% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/visualization.R
\name{plotFragDistribution}
\alias{plotFragDistribution}
\title{Plot Fragility time distribution for two electrodes group marked and non-marked as soz}
\usage{
plotFragDistribution(
  stat = NULL,
  timeRange = NULL,
  title = "Average Fragility over time"
)
}
\arguments{
\item{stat}{FragStat object, a FragStat object from \code{fragStat}}

\item{timeRange}{Numeric Vector of length 2. The time window to display at the x-axis}

\item{title}{String. Figure title}
}
\value{
plot fragility distribution
}
\description{
Plot Fragility time distribution for two electrodes group marked and non-marked as soz
}
\examples{
data("pt01Epochm1sp2s")
sozIndex<-attr(pt01Epochm1sp2s,"sozIndex")
# Load the precomputed fragility object
timeRange <- c(-10,10)
data("pt01Fragm1sp2s")
# compute fragility statistics evolution with time (mean and standard deviation) for soz and
# non soz groups
pt01fragstat <- fragStat(frag=pt01Fragm1sp2s, sozID=sozIndex)
# plot the statistical results
pfragstat<-plotFragDistribution(stat=pt01fragstat,timeRange=timeRange)
pfragstat
}
