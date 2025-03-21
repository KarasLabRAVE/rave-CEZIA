## A more realistic example with parallel computing
## Not run:
## Register a SNOW backend with 4 workers
library(parallel)
library(doSNOW)


library(R.matlab)
library(readxl)
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

cl <- makeCluster(4, type = "SOCK")
registerDoSNOW(cl)

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
stopCluster(cl)

sozNames <- attr(pt01EcoG, "sozNames")
display <- c(sozNames, "MLT1", "MLT2", "MLT3", "MLT4")
plotFragHeatmap(frag = fragtest[display], sozIndex = sozNames)

plotFragDistribution(frag = fragtest[display], sozIndex = sozNames)

plotFragQuantile(frag = fragtest[display], sozIndex = sozNames)
