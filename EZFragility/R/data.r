#' Pt01 seizure 1 (-10:10s) around seizure onset
#' 
#' This data corresponds to the first seizure of patient PT01 from 
#' the Fragility Data Set.
#' The data contains only the good channels. 
#' It has been notch filtered and common average referenced in RAVE. 
#' The full data has been epoched -1:2s around the seizure onset.
#' The acquisition frequency is 1000 Hz
#' EcoG recording gathered in collaboration with the National Institute of Health
#'
#' @docType data
#'
#' @usage 
#' ## EEG data
#' data(pt01Epochm1sp2s)
#' 
#' @format 
#' pt01Epochm1sp2s: A Matrix with 3000 rows (time points) and 84 columns (electrodes)
#' 
#' pt01Fragm1sp2s: A fragility object results of applying the main function \code{calcAdjFrag}
#' to pt01Epochm1sp2s
#'
#' @keywords datasets
#'
#' @source Fragility Multi-Center Retrospective Study
#' (\href{https://openneuro.org/datasets/ds003029/versions/1.0.0}{OpenNeuro})
#' 
#' @aliases pt01Epochm1sp2s pt01Fragm1sp2s
"pt01Epochm1sp2s"

