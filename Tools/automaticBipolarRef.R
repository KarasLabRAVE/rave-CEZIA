# updated 022525
# Anne-Cecile Lesage
# UTMB
# automatic bipolar referencing

library(readxl)
library(writexl)
library(dplyr)

# Load the input data
channels <- read_excel("/Users/aclesage/Documents/RAVEProjects/BipolarReferencing/Kishan/sub-pt01_ses-presurgery_task-ictal_acq-ecog_run-02_channels.xlsx")

# count number of channels
nelec<-nrow(channels)

Type1 <- rep("Bipolar Reference",nelec)
# Extract Group name
Group <- gsub("[^a-zA-Z]", "", channels$name) # set this before the df is cleaned so that we can account for all grid electrodes
# Give a number id to each electrode group
GroupID <- as.numeric(dense_rank(Group))

Electrode<-c(1:nelec)
# set up the final df for bipolar referencing and all of the relevant columns 
Reference <- gsub(" ", "", paste("ref_",nelec)) # set this as default for now 
name<-channels$name
status<-channels$status
type<-channels$type

# Build the bipolar reference dataframe and its column names
bipolar_ref <- cbind.data.frame(Electrode,Group,Reference,Type1,GroupID,name,status,type)
# Filter out bad electrodes
bipolar_ref=bipolar_ref[bipolar_ref$status=="good",] 
bipolar_ref=bipolar_ref[bipolar_ref$type=="ECOG",] 
# Identify all unique group names as encountered in the electrode ordering
Group<-bipolar_ref$Group
groupnames<-unique(Group)

bipolar_ref$status <- NULL 
bipolar_ref$type <- NULL 
bipolar_ref$Label <- as.numeric(gsub("[^0-9]", "", bipolar_ref$name))
# Max contact number by group
ngroup=length(groupnames)
for (i in 2:ngroup) {

  namegroup<-groupnames[i]
  if(namegroup=='G'){
    
  }else{
    groupid<-which(bipolar_ref$Group==namegroup)
    maxcontact=length(groupid)
    je=maxcontact-1
    for (j in 1:je) {
      #j=3
      if(bipolar_ref$Label[groupid[j]]!=maxcontact){
        jp=bipolar_ref$Label[groupid[j]]+1
        nextc<-paste(namegroup,jp,sep="")
        id<-which(bipolar_ref$name==nextc)
        bipolar_ref$Reference[groupid[j]]<-paste("ref_",bipolar_ref$Electrode[id],sep="")
      }
    }  
    bipolar_ref$Reference[groupid[j+1]]<-"noref"
  }
}  
  
#}  



