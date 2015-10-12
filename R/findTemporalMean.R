#' Calculate the Temporal Mean of BOLD image
#' 
#' Find the Mean of BOLD data Time Series, this is done by calculating 
#' the mean of each voxel at all timepoints
#' 
#' @param  bold input 4D image
#' @param  boolean for plotting, default is FALSE
#' @return antsImage of temporal mean
#' 
#' @examples
#' fmri <- antsImageRead(getANTsRData("rsbold"))
#' tempMean <- rsBOLDTemporalMean(fmri)
#' 
#' 
#' @export findTemporalMean


findTemporalMean <- function(img,providePlot = FALSE){
  
  #check requirements
  if (!usePkg("ANTsR")) {
    print("Need ANTsR package")
    return(NULL)
  }
  
  if (!usePkg("R.cache")) {
    print("Need R.cache package")
    return(NULL)
  }
  
  if (!is.antsImage(img)) {
    stop("ERROR: Input image must be of classs 'antsImages'")
  }  
  
  boldMean <- apply.antsImage(img,c(1,2,3),mean)
  if(providePlot){
    plot(boldMean, axis=3, slices=1:30, ncolumns=10)
  }
  
  
  return(boldMean)
}
  

