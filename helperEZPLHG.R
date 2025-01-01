data("PT01Epochm30sp30s")
ieegts=PT01Epochm30sp30s
sizeWindow <- 3000
sizeSkip <- 333
fs=1000
tBaseline=20
time_window_ictal=c(-10,10)
time_window=c(-30,30)

start <- Sys.time()
print(start)
resPLHG<-calc_PLHG(ieegts = PT01Epochm30sp30s, sizeWindow=sizeWindow, sizeSkip=sizeSkip,fs=fs,tBaseline=tBaseline,time_window=time_window,time_window_ictal=time_window_ictal)
end <- Sys.time()
print(end - start)

data("ElectrodesDataPT01")
head(ElectrodesDataPT01,n=10)
displayChannels=which(ElectrodesDataPT01$insoz==TRUE)
scaling=1000000
visuiEEGdata(ieegts=PT01Epochm30sp30s,scaling=scaling, displayChannels = displayChannels)

data("PT01Epochm30sp30s")
data("ElectrodesDataPT01")
data("resPLHG")
ElectrodesData=ElectrodesDataPT01
time_window_ictal=c(-10,10)
subject_code='PT01'
j=1
heatmap_PLHG(resPLHG=resPLHG, ElectrodesData=ElectrodesDataPT01,ieegts=PT01Epochm30sp30s,time_window_ictal,subject_code,j)
