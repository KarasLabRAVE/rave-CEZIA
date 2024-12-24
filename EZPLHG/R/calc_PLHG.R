#' Computes the Phase Locked High Gamma matrix for ictal iEEG
#' Implementation of PLHG as described in:
# Ictal high frequency oscillations distinguish two types of seizure territories in humans
# Shennan A. Weiss, Garrett P. Banks, Guy M. McKhann, Jr, Robert R. Goodman, Ronald G. Emerson, Andrew J. Trevelyan, and Catherine A. Schevon
# Brain. 2013 Dec; 136(12): 3796?3808.
#'
#' @param ieegts Numeric. A matrix of iEEG time series x(t),
#' with time points as rows and electrodes names as columns
#' @param sizeWindow Integer. The number of time points to use in each window
#' @param sizeSkip Integer.The number of time points to move the window each time
#' @param fs Acquisition frequency
#' @param tBaseline baseline duration in seconds
#'
#' @return plhgMaster matrix of PLHG values
#' @export
#'
#' @examples
#' data("PT01Epochm30sp30s")
#' sizeWindow <- 3000
#' sizeSkip <- 333
#' fs=1000
#' tBaseline=20
#' time_window_ictal=c(-10,10)
#' time_window=c(-30,30)
#' resPLHG<-calc_PHG(ieegts = PT01Epochm30sp30s, sizeWindow=sizeWindow, sizeSkip=sizeSkip,fs=fs,tBaseline=tBaseline,time_window=time_window,time_window_ictal=time_window_ictal)
calc_PLHG<- function(ieegts,sizeWindow,sizeSkip,fs,tBaseline,time_window_ictal,time_window) {

  np<-reticulate::import('numpy')

  ## Number of electrodes and time points
  n_tps <- nrow(ieegts)
  n_elec <- ncol(ieegts)

  electrode_list <- colnames(ieegts)

  # number of windows for PLHG analysis
  nw<-floor(n_tps-sizeWindow)/sizeSkip

  #
  ntb=tBaseline*fs
  tsBaseline=matrix(0,ntb,n_elec)

  tsBaseline[1:ntb,1:n_elec]=ieegts[1:ntb,1:n_elec]

  #Epoch Ictal signal
  its=as.integer((time_window_ictal[1]-time_window[1])*fs)
  nt=as.integer((time_window_ictal[2]-time_window_ictal[1])*fs)
  ite=its+nt

  dt <- as.numeric(1/fs)

  ts <- 1:nt
  ts <- ts*dt+time_window_ictal[1]
  tsIctal<-ieegts[its:ite,1:n_elec]

  # number of windows for PLHG analysis
  nw<-floor((nt-sizeWindow)/sizeSkip)




}
