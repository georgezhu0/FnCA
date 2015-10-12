#' Identify Nuisance Parameters in Timeseries
#'
#' Collect the nuisance parameters in BOLD global signals to regress out.
#' Possible nuisance parameters:
#' \itemize{
#' \item{detrended motion parameters, their squares, and the derivatives of both}
#' \item{mean signal in white Matter and its derivative}
#' \item{mean signal in CSF and its derivative}
#' \item{physiologocial noise estimated via compcor()}
#' }
#'
#' @param boldMat BOLD Time Series matrix
#' @param seg Tissue segmentation image that identifies: CSF, gray and white matter
#' @param mask mask for BOLD image (3D)
#' @param goodtimes list of interger index of good timepoints
#' @param badtimes list of interger index of bad timepoints
#' @param providePlot Boolean determine whether descriptive statistics should be provided
#' if TRUE, the plot will be cached and can be called
#' @return list of output containing:
#'
#' @export collectNuisParams


collectNuisParams <- function(
  boldMat, seg, mask, goodtimes, badtimes, providePlot = TRUE){
  #check requirements
  if (!usePkg("ANTsR")) {
    print("Need ANTsR package")
    return(NULL)
  }


  if (!usePkg("R.cache")) {
    print("Need R.cache package")
    return(NULL)
  }

  if (providePlot){
    if (!usePkg("ggplot2")) {
      print("Need ggplot2 package")
      return(NULL)
    }
  }

  if (!is.antsImage(seg)) {
    stop("ERROR: Input segmentaion must be of classs 'antsImages'")
  }

  if(!is.antsImage(mask)){
    stop("ERROR: Input mask must be of classs 'antsImages'")
  }
  #required input
  nTimes <- length(goodtimes) + length(badtimes)

  ## white matter
  #labeled as 3
  wmMask = seg*1
  wmMask[ wmMask != 3] = 0
  wmMask[ wmMask == 3 ] = 1
  wmMask = iMath( wmMask, "ME", 1)
  wmVox = which(subset(wmMask, mask > 0 )==1)
  wmMean = rowMeans(boldMat[,wmVox])

  ## CSF
  #labeled as 1
  csfMask = seg*1
  csfMask[ csfMask != 1] = 0
  csfVox = which(subset(csfMask, mask > 0)==1)
  csfMean= rowMeans(boldMat[,csfVox])

  globalMean = rowMeans(boldMat)

  compcor = compcor(boldMat[goodtimes,], ncompcor = 4)
  compcorNuis = matrix(0, nTimes, 4)
  compcorNuis[goodtimes,] = compcor

  #if there is bad timepoints
  if(!(is.integer(badtimes) && length(badtimes) == 0)){
    compcorNuis[badtimes,] = NA
  }

  ##tissue
  tissueNuis = cbind(globalMean, wmMean, csfMean)
  if ( length(badtimes) > 0 ) {
    for ( v in c(1:dim(tissueNuis)[2]) ) {
      tissueInterp = spline( c(1:nTimes)[goodtimes], tissueNuis[goodtimes,v], xout=badtimes )$y
      tissueNuis[badtimes,v]=tissueInterp
    }
  }
  tissueDeriv = rbind( rep(0,dim(tissueNuis)[2]), diff(tissueNuis,1) )

  ##mean cortex signal
  ctxMask = seg*1
  ctxMask[ ctxMask != 2] = 0
  ctxMask[ ctxMask == 2 ] = 1
  ctxVox = which(subset(ctxMask, mask > 0)==1)
  ctxMean = rowMeans(boldMat[,ctxVox])

  #plotting
  if(providePlot){
    tissueType = c( rep("Global", nTimes), rep("White matter",nTimes), rep("CSF",nTimes) )
    tissueType = c(tissueType, rep("CompCor1",nTimes), rep("CompCor2",nTimes))
    tissueType = c(tissueType, rep("CompCor3",nTimes), rep("CompCor4",nTimes) )

    tissueCategory = c(rep("Tissue", nTimes*3), rep("CompCor", nTimes*4))

    signal = c(globalMean, wmMean, csfMean, compcorNuis[,1], compcorNuis[,2])
    signal = c(signal, compcorNuis[,3], compcorNuis[,4])

    tissue.dat = data.frame( Time=rep(1:nTimes,7) )
    tissue.dat$Signal = signal
    tissue.dat$Type = tissueType
    tissue.dat$Category = tissueCategory

    tissuePlot <- ggplot(tissue.dat, aes(x=Time, y=Signal, group=Type, colour=Type) )
    tissuePlot <- tissuePlot + geom_line(size=0.5)
    tissuePlot <- tissuePlot + theme(text=element_text(size=10), legend.position="top")
    tissuePlot <- tissuePlot + facet_grid(Category ~ ., scales="free" )
    tissuePlot <- tissuePlot + ggtitle("Nuisance parameters")

    saveCache(tissuePlot,key = list("tp"))

  }
  return(list(
    ctxVox = ctxVox,
    ctxMean = ctxMean,
    wmVox = wmVox,
    wmMean = wmMean,
    csfVox = csfVox,
    csfMean = csfMean,
    globalMean = globalMean,
    compcorNuis = compcorNuis,
    tissueNuis = tissueNuis,
    tissueDeriv = tissueDeriv))
}







