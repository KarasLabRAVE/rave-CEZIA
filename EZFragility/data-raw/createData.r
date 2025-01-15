## load pt01epochdata.mat
## Patient PT01 from the Fragility data set

library(R.matlab)
library(readxl)
library(gifti)
library(rgl)
library(magick)

data <- readMat('data-raw/pt01epochdata.mat')
pt01Epoch <- data$a

## add channel names to the rows
goodChannels <- c(1:4,7:36,42:43,46:69,72:95)
sozChannels<-c(33:34,62:69)
channelNames <- read_excel('data-raw/Pt01ictalRun01EcoGChannels.xls')
rownames(pt01Epoch) <- channelNames$name[goodChannels]
sozindex<-which(goodChannels%in%sozChannels==TRUE)
soznames<-channelNames$name[sozChannels]

## Add time stamps to the columns
times <- seq(-10, 10, length.out=ncol(pt01Epoch))
times_with_sign <- ifelse(times >= 0, paste0("+", times), as.character(times))
# why this?
#colnames(ptEpoch) <- paste0('t', times_with_sign)
colnames(pt01Epoch)<-times_with_sign

pt01Epoch <- t(pt01Epoch)
attr(pt01Epoch, "sozindex") <- sozindex
attr(pt01Epoch, "soznames") <- soznames
usethis::use_data(pt01Epoch, overwrite = TRUE)

pt01Epochm1sp2s<-pt01Epoch[9001:12000,]
attr(pt01Epochm1sp2s, "sozindex") <- sozindex
attr(pt01Epochm1sp2s, "soznames") <- soznames
usethis::use_data(pt01Epochm1sp2s, overwrite = TRUE)

pt01Epochm3sp5s<-pt01Epoch[7001:15000,]
attr(pt01Epochm3sp5s, "sozindex") <- sozindex
attr(pt01Epochm3sp5s, "soznames") <- soznames
usethis::use_data(pt01Epochm3sp5s, overwrite = TRUE)

# plot 3d brain in gif
lhpt01 <- readgii('data-raw/pt01lh.pial.gii')
lhptb<-lhpt01$data$pointset
rhpt01 <- readgii('data-raw/pt01rh.pial.gii')
rhptb<-rhpt01$data$pointset

size<-dim(lhptb)
np=size[1]
sample=seq(1,np,by=20)
lhpts=lhptb[sample,]
plot3d(x=lhpts[,1],y=lhpts[,2],z=lhpts[,3],xlab='X in mm',ylab='Y in mm',zlab='Z in mm')
rhpts=rhptb[sample,]
#rhpts[,1]=rhpts[,1]*0.1
plot3d(x=rhpts[,1],y=rhpts[,2],z=rhpts[,3],xlab='X in mm',ylab='Y in mm',zlab='Z in mm')
aspect3d(0.6,1,1)

movie3d(spin3d(axis=c(0,0,1),rpm=1),duration=30,dir="./data-raw",movie="brainvisupt01")