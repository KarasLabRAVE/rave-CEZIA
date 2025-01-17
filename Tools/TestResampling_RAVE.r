# updated 011725
# Anne-Cecile Lesage
# UTMB
library(gsignal)

###########################################
# ---- Analysis script --------------------------------------------------------


pts <- dipsaus::parse_svec("1-2")

pipeline_xls <- readxl::read_xlsx("/Volumes/bigbrain/ACL/RAVE_Projects/PipelineScripts/KDataEEGDataset_pipeline_update_011625.xlsx")

ravepreproc='/Volumes/bigbrain/rave_data/data_dir/Retrostudy/'
raveraw='/Volumes/bigbrain/rave_data/raw_dir/'
#scaling=200000

time_window <- c(-10, 10)
reference_name <- "car"


#for(i in pts){
  
  i=2
  
  patname <- pipeline_xls$subject[i]
  
  project_name <- pipeline_xls$project[i]
  
  gelectrodes <- dipsaus::parse_svec(pipeline_xls$good_electrodes[i])
  soz<-dipsaus::parse_svec(pipeline_xls$soz[i])
  resect<-dipsaus::parse_svec(pipeline_xls$resect[i])
  fs<-pipeline_xls$sample_rate[i]
  if(class(fs)=='character'){
    fs=as.numeric(fs)
  }
  ictal_runs <- dipsaus::parse_svec(pipeline_xls$ictal_runs[i])
  load_electrodes<-gelectrodes
  
  
  insoz=which(load_electrodes%in%soz==TRUE)
  
  epoch_name<-paste(patname,"_seizure",sep="")
  
  
 
  subject <- raveio::RAVESubject$new(project_name = project_name, subject_code = patname)
  repository <- raveio::prepare_subject_voltage_with_epoch(
    subject = subject,
    epoch_name = epoch_name,
    electrodes = load_electrodes,
    time_windows = time_window,
    reference_name = reference_name
  )
  
  
  #for(j in ictal_runs){
    
   j=1 
    
    condition=paste('sz',j,sep="")
    print(condition)
    # Get initial matrix structure
    voltage_for_analysis <- repository$voltage$data_list[[sprintf("e_%s", load_electrodes[1])]]
    repository$epoch_table$Condition
    voltage_trace <- voltage_for_analysis[, 1, 1]
    selector <- repository$epoch_table$Condition %in% condition
    trial_list <- repository$epoch_table$Trial[selector]
    selected_trial_data <- subset(voltage_for_analysis, Trial ~ Trial %in% trial_list)
    collapsed_trial <- raveio::collapse2(selected_trial_data, keep = 1)

    
    electrodes_by_time <- matrix(nrow = length(load_electrodes), ncol = length(collapsed_trial))
    
    # Loop through each electrode
    for (ii in 1:length(load_electrodes)) {
      electrode <- load_electrodes[ii]
      
      voltage_for_analysis <- repository$voltage$data_list[[sprintf("e_%s", electrode)]]
      
      repository$epoch_table$Condition
      
      voltage_trace <- voltage_for_analysis[, 1, 1]
      
      selector <- repository$epoch_table$Condition %in% condition
      trial_list <- repository$epoch_table$Trial[selector]
      selected_trial_data <- subset(voltage_for_analysis, Trial ~ Trial %in% trial_list)
      collapsed_trial <- raveio::collapse2(selected_trial_data, keep = 1)
      
      electrodes_by_time[ii, ] <- collapsed_trial
      
    }
    
  dumm<-t(electrodes_by_time)  
  dumm2<-gsignal::resample(t(electrodes_by_time),1000,fs)
  elec_resamp<-t(dumm2)
  
  nt=ncol(electrodes_by_time)
  nts=ncol(elec_resamp)
  
  stimesb=seq(time_window[1],time_window[2],length.out=nt)
  stimess=seq(time_window[1],time_window[2],length.out=nts)
  
  y=electrodes_by_time[insoz[1],]
  ys=elec_resamp[insoz[1],]
  
  plot(stimesb,y,type='l',xlim=c(0,0.1))
  lines(stimess,ys,col="red",type='o')
  
#}    
#}


