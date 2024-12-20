#' Computes the Phase Locked High Gamma matrix for ictal iEEG
#'
#' @param ieegts Numeric. A matrix of iEEG time series x(t),
#' with time points as rows and electrodes names as columns
#' @param sizeWindow Integer. The number of time points to use in each window
#' @param sizeSkip Integer.The number of time points to move the window each time
#' @param fs Acquisition frequency
#' @param tBaseline baseline start from seizure onset
#'
#' @return plhgMaster matrix of PLHG values
#' @export
#'
#' @examples
calc_PLHG<- function(ieegts,sizeWindow,sizeSkip,fs,tBaseline) {

  ## Number of electrodes and time points
  n_tps <- nrow(ieegts)
  n_elec <- ncol(ieegts)

  electrode_list <- colnames(ieegts)

  # number of windows for PLHG analysis
  nw<-floor(nt-sizeWindow)/sizeSkip

}
