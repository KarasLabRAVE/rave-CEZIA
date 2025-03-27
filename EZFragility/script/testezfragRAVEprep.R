##############################
library(R.matlab)
library(readxl)
library(parallel)
library(doSNOW)
library(scales)
data <- readMat('data-raw/pt01epochdata.mat')
pt01EpochRaw <- data$a

## add channel names to the rows
goodChannels <- c(1:4,7:36,42:43,46:69,72:95)
sozChannels<-c(33:34,62:69)
channelNames <- read_excel('data-raw/Pt01ictalRun01EcoGChannels.xls')
rownames(pt01EpochRaw) <- channelNames$name[goodChannels]
sozIndex<-which(goodChannels%in%sozChannels==TRUE)
sozNames<-channelNames$name[sozChannels]

## Add time stamps to the columns
times <- seq(-10, 10, length.out=ncol(pt01EpochRaw))
times_with_sign <- ifelse(times >= 0, paste0("+", times), as.character(times))
colnames(pt01EpochRaw)<-times_with_sign

pt01EcoG<-pt01EpochRaw
attr(pt01EcoG, "sozIndex") <- sozIndex
attr(pt01EcoG, "sozNames") <- sozNames

cl <- parallel::makeCluster(4, type = "SOCK")
doSNOW::registerDoSNOW(cl)

#data("pt01EcoG")
epoch <- Epoch(pt01EcoG)
window <- 250
step <- 125
title <- "PT01 seizure 1"
fragtest<-calcAdjFrag(
  epoch = epoch, window = window,
  step = step, parallel = TRUE, progress = TRUE
)

## stop the parallel backend
parallel::stopCluster(cl)

## Result visualization

sozIndex <- attr(pt01EcoG, "sozIndex")
display <- c(sozNames, "MLT1", "MLT2", "MLT3", "MLT4")

plotheatmap<-plotFragHeatmap(fragtest, sozIndex)
plotheatmap<-plotheatmap+ggplot2::ggtitle("Fragility heatmap for patient PT01") 
# Add a vertical line at Seizure Onset
#plotheatmap <- plotheatmap + ggplot2::geom_vline(xintercept = as.numeric(0), color = "black", linetype = "dashed", size = 1)
plotheatmap
ggplot2::ggsave("~/pt01FragHeatmap.png")

plotDistribution<-plotFragDistribution(fragtest, sozIndex)
plotDistribution<-plotDistribution+ggplot2::ggtitle("Pooled fragility distribution for patient PT01") 
plotDistribution

plotquantile<-plotFragQuantile(fragtest, sozIndex)
plotquantile<-plotquantile+ggplot2::ggtitle("Pooled fragility quantiles for patient PT01") 
plotquantile


## Visualize a subject of electrodes
sozIndex <- attr(pt01EcoG, "sozIndex")
display <- c(sozIndex, 77:80)

epoch <- Epoch(pt01EcoG)
plotiEEG<-visuIEEGData(epoch = epoch[display, ])
plotiEEG<-plotiEEG+ggplot2::ggtitle("iEEG traces for patient PT01 around seizure Onset")
plotiEEG

## Compute quantiles, mean and standard deviation for two electrodes group marked as soz non marked as soz

sozIndex <- attr(pt01EcoG, "sozIndex")
pt01fragstat <- fragStat(fragtest, sozIndex)
