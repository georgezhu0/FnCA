#' Compute Correlation matrix for ROI signals
#'
#' Use the mean time signal for each ROI and find the correlation between
#' each ROI signals
#'
#' @param roiMat ROI average signal matrix
#' @param labels A set of labels identifying anatomical regions of interest
#'
#' @return ROI signal correlation matrix
#'
#' @export corrMatROI

corrMatROI <- function(
  roiMat,labels){
  #check requirements
  if (!usePkg("ANTsR")) {
    print("Need ANTsR package")
    return(NULL)
  }

  nLabels = max(labels)

  missingROIs = which(colMeans(roiMat)==0)
  goodROIs = (1:nLabels)
  if ( length(missingROIs) > 0 ) {
    goodROIs = goodROIs[-missingROIs]
  }

  connMat = suppressWarnings(cor(roiMat))
  diag(connMat) = rep(0, length(diag(connMat)) )
  if ( length(missingROIs) > 0 ) {
    connMat[missingROIs,] = 0
    connMat[,missingROIs] = 0
  }
  return(connMat)
}
