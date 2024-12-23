#' Pt01 seizure 1 (-10:10s) around seizure onset
#' 
#' This data corresponds to the first seizure of patient PT01 from the Fragility Data Set.
#' The data contains only the good channels. 
#' It has been notch filtered and common average referenced in RAVE 
#' It has been epoched -10:10s around the seizure onset
#' The acquisition frequency is 1000 Hz
#' EcoG recording gathered in collaboration with the National Institue of Health
#'
#' @docType data
#'
#' @usage data(ptEpoch)
#'
#' @format A Matrix with 20001 rows (time points) and 84 columns (electrodes)
#'
#' @keywords datasets
#'
#' @source Fragility Multi-Center Retrospective Study
#' (\href{https://openneuro.org/datasets/ds003029/versions/1.0.0}{OpenNeuro})
#'
#' @examples
#' data(ptEpoch)
"ptEpoch"

#' Pt01 seizure 1 (-3:5s) around seizure onset
#' 
#' This data corresponds to the first seizure of patient PT01 from the Fragility Data Set.
#' The data contains only the good channels. 
#' It has been notch filtered and common average referenced in RAVE 
#' It has been epoched -[3:5]s around the seizure onset
#' The acquisition frequency is 1000 Hz
#' EcoG recording gathered in collaboration with the National Institue of Helath
#'
#' @docType data
#'
#' @usage data(pt01Epochm3sp5s)
#'
#' @format A Matrix with 8001 rows (time points) and 84 columns (electrodes)
#'
#' @keywords datasets
#'
#' @source Fragility Multi-Center Retrospective Study
#' (\href{https://openneuro.org/datasets/ds003029/versions/1.0.0}{OpenNeuro})
#'
#' @examples
#' data(pt01Epochm3sp5s)
"pt01Epochm3sp5s"

#'
#' @docType data
#'
#' @usage data(fragm3sp5s)
#'
#' @format fragility matrix result of example 2 for calc_adj_frag function help
#'  with and 84 columns (electrodes)
#'  -[3:5]s around the seizure onset
#'
#' @keywords datasets
#'
#' @source Fragility Multi-Center Retrospective Study
#' (\href{https://openneuro.org/datasets/ds003029/versions/1.0.0}{OpenNeuro})
#'
#' @examples
#' data(fragm3sp5s)
"fragm3sp5s"

