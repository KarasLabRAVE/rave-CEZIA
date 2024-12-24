#' Visualization functions (raw signal, fragility matrix)
#' 
#' plot fragility heatmaps with electrodes marked as soz colored
#'
#' @param frag fragility matrix results
#' @param ElectrodesData electrodes data
#' @param ieegts Numeric. A matrix of iEEG time series x(t), 
#' with time points as rows and electrodes names as columns
#' @param time_window Fragility heatmap time window around seizure onset
#' @param subject_code patient name
#' @param j seizure number
#'
#' @return Heatmap plot of the fragility matrix with soz electrodes in blue in the bottom
#' @export
#'
#' @examples
#' data("fragm3sp5s")
#' data("fragrankm3sp5s")
#' data("pt01Epochm3sp5s")
#' data("ElectrodesDataPT01")
#' time_window=c(-3:5)
#' heatmap_frag(frag=fragm3sp5s,ElectrodesData=ElectrodesDataPT01,ieegts=pt01Epochm3sp5s,time_window=c(-3,5),subject_code='pt01',j=1)
heatmap_frag<-function(frag,ElectrodesData,ieegts,time_window,option=NULL,subject_code,j){
  
  titlepng=paste(subject_code,'Seizure',as.character(j),sep=" ")
  elecsoz=which(ElectrodesData$insoz==TRUE)
  elecsozc=which(ElectrodesData$insoz==FALSE)
  elecsozsozc=c(elecsoz,elecsozc)
  
  elecnum <- colnames(ieegts)[elecsozsozc]
  n_elec <- ncol(ieegts)
  nw<- ncol(frag)
  colorelec<-elecnum
  nsoz=length(elecsoz)
  colorelec[1:n_elec]="black"
  colorelec[1:nsoz]="blue"

  fragord<-frag[elecsozsozc,]
  fragdf<-data.frame(fragord)
  stimes=c(1:nw)*(time_window[2]-time_window[1])/nw+time_window[1]
  colnames(fragdf)<-stimes
  rownames(fragdf)<-elecnum
  
  fragmap_data <- expand.grid(Time = stimes, Electrode = elecnum)
  fragmap_data$Value <- c(t(fragord))
  
  ggplot2::ggplot(fragmap_data, ggplot2::aes(x = Time, y = Electrode, fill = Value)) +
    ggplot2::geom_tile() +
    ggplot2::ggtitle(titlepng)+
    ggplot2::labs(x = "Time (s)", y = "Electrode",size=2) +
    viridis::scale_fill_viridis(option = "turbo") +  #
    
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.y = ggplot2::element_text(size=4,colour=colorelec),     # Adjust depending on electrodes
    )
  
}

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