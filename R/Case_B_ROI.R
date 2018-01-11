#' Sarcoma tumor region of interest in the pelvis. 
#' This scan was acquired at University of Washington in Seattle, USA.
#' Uptake values are expressed in SUV scale for this test dataset
#' (Clinical SUVmax = 7.8). The ellipsoidal ROI was hand-drawn by expert 
#' at the UW School of Medicine. 
#' This is frame 0, gate 0, i.e. one static FDG-PET scan.
#' @format An array with 5 columns as follows:
#' \describe{
#'    \item{Value}{Uptake value}
#'    \item{Weight}{Voxel weight (dummy values set to 1)}
#'    \item{X (mm)}{x-coordinate in the scanner referential}
#'    \item{Y (mm)}{y-coordinate in the scanner referential}
#'    \item{Z (mm)}{z-coordinate in the scanner referential}
#' }
"Case_B_ROI"
