#' Visualization iEEGData
#'
#' @param ieegts Numeric. A matrix of iEEG time series x(t),
#' with time points as rows and electrodes names as columns
#' @param scaling Scaling factor
#' @param displayChannels Channels to display
#'
#' @return plot raw signal
#' @export
#'
#' @examples
#' data("PT01Epochm30sp30s")
#' data("ElectrodesDataPT01")
#' displayChannels=which(ElectrodesDataPT01$insoz==TRUE)
#' visuiEEGdata(ieegts=PT01Epochm30sp30s,1000000, displayChannels = displayChannels)
visuiEEGdata<-function( ieegts, scaling, displayChannels){


  plotData<-ieegts[,displayChannels]/scaling
  gaps<-2
  displayNames=colnames(ieegts)[displayChannels]

  for(i in seq_along(plotData)){
    plotData[, i] <- (plotData[, i]- mean(plotData[, i]))+
      (ncol(plotData)-i)*gaps
  }

  plot(plotData[, 1],type="l" ,cex=0.1,
       ylim = range(plotData), yaxt = "n")
  for(i in 2:ncol(plotData)){
    lines(plotData[, i])
  }
  axis(2, at = rev(seq_along(displayChannels) - 1)*gaps,
       labels = displayNames,las=1)


}
