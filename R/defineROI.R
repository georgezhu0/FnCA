#' Defining ROI of functional network
#'
#' A set of ROI is needed for creating network, define ROIs that contain many
#' voxels and represent regions for with priori knowledge regarding the
#' functional network to which each region belongs. Default radius of ROI
#' sphere is 5.0mm
#'
#' @param mask mask for BOLD image (3D)
#' @param pts labels identifying anatomical regions of interest
#' @param boldMat BOLD signal matrix
#' @param goodtimes list of interger index of good timepoints
#' @param radius radius of each ROI sphere, default is 5.0mm
#'
#' @return ROI labels image
#' list of output containing:
#' \itemize{
#' \item{labelImg:}{ Labeled ROI img}
#' \item{labels:}{ ROI labels}
#' \item{labeledBoldMat:}{ Labeled BOLD signal matrix}
#' }
#'
#' @export defineROI
#'
defineROI <- function(
  mask,pts,boldMat,goodtimes,radius = 5.0){
  #check requirements
  if (!usePkg("ANTsR")) {
    print("Need ANTsR package")
    return(NULL)
  }

  if(!is.antsImage(mask)){
    stop("ERROR: Input mask must be of classs 'antsImages'")
  }

  #create inputs
  labelImg <- mask*0
  nPts <- dim(pts)[1]
  rad <- radius
  n <- ceiling(rad / antsGetSpacing(mask))

  #Defining ROI from each voxel
  for ( r in 1:nPts) {
    pt = as.numeric(c(pts$X[r], pts$Y[r], pts$Z[r] ))
    idx = antsTransformPhysicalPointToIndex(mask,pt)
    for ( i in c(-n[1]:n[1]) ) {
      for (j in c(-n[2]:n[2])) {
        for (k in c(-n[3]:n[3])) {
          local = idx + c(i,j,k)
          localpt = antsTransformIndexToPhysicalPoint(mask,local)
          dist = sqrt( sum( (localpt-pt)*(localpt-pt) ))
          inImage = ( prod(idx <= dim(mask))==1) && ( length(which(idx<1)) == 0 )
          if ( (dist <= rad) && ( inImage == TRUE ) ) {
            labelImg[local[1], local[2], local[3]] = pts$ROI[r]
          }
        }
      }
    }
  }

  #create a labeled mask
  labelMask = labelImg*1
  labelMask[labelMask > 0] = 1
  labelMask[mask == 0] = 0
  labelVox = which(subset(labelMask, mask > 0)==1)

  #roi labels
  labels = labelImg[labelMask > 0]

  #label bold matrix
  labeledBoldMat = boldMat[goodtimes,labelVox]


  return(list( labelImg = labelImg, 
              labels = labels, 
              labeledBoldMat = labeledBoldMat))
}
