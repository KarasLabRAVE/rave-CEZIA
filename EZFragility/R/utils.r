isWholeNumber <- function(x) {
  return(x %% 1 == 0)
}

#' validate seizure onset data
#'
#' @param ieeg_ts Numeric. A matrix of iEEG time series x(t), 
#' with time points as rows and electrodes names as columns
#' @param sozindex Integer. Vector soz electrodes 
#' @param soznames Vector string. soz electrodes names 
#'
#' @return boolean
#' @export
#'
#' @examples
#' data("pt01Epochm3sp5s")
#'sozindex<-attr(pt01Epochm3sp5s,"sozindex")
#'soznames<-attr(pt01Epochm3sp5s,"soznames")
#'valid<-valid_soz(ieeg_ts=pt01Epochm3sp5s,sozindex=sozindex,soznames=soznames)
#'
valid_soz <- function( ieeg_ts, sozindex, soznames){
  
  valid<-TRUE
  elecnames<-colnames(ieeg_ts)
  elecind=c(1:ncol(ieeg_ts))
  
  diffelecind<-setdiff(sozindex,elecind)
  
  if(length(diffelecind)!=0){
    listelecmissing<-paste(as.character(diffelecind),collapse=" ")
    message<-paste("ERROR in soz electrodes names. Numbers ",listelecmissing,"are out of electrode number limit")
    warning(message)
    valid<-FALSE
  }
  
  a<-character(0)
  #soznames<-c("LB1","LB2")
  diffelecname<-setdiff(soznames,elecnames)
 
  if(!identical(diffelecname,a)){
    listelecmissing<-paste(diffelecname,collapse=" ")
    message<-paste("ERROR in soz electrodes names.",listelecmissing,"are not in the electrode name set")
    warning(message)
    valid<-FALSE
  }
  
  return(valid)
  
}
