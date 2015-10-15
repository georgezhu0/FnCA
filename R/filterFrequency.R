#' Filter and Smooth BOLD Timeseries Frequency
#'
#' Frequency filtering get rid of non-evenly sampled time-series data. After
#' filtering out bad timepoints, spatial smoothing is applied.The Default
#' frequency range for filter is 0.009 Hz - 0.08 Hz 0.009 Hz - 0.08 Hz, default
#' smooth factor is Gaussian kernel with FWHM=6.0mm
#'
#' @param img BOLD fMRI time series image
#' @param mask mask for BOLD image (3D)
#' @param boldMat BOLD time series matrix
#' @param ctxVox list of integer index of mean cortex signal vortex
#' @param goodtimes list of interger index of good timepoints
#' @param badtimes list of interger index of bad timepoints
#' @param lowFreq lower bound of frequency for filtering, default is 0.009 Hz
#' @param highFreq upper bound of frequency for filtering, default is 0.08 Hz
#' @param smoothFactor frequency smooth factor ,default FWHM is 6.0mm
#' @param providePlot Boolean determine whether descriptive statistics should be provided
#' if TRUE, the plot will be cached and can be called
#'
#' @return list of output containing:
#' \itemize{
#' \item{img:}{ Processed 4D fMRI image}
#' \item{boldMat:}{ Processed BOLD time series matrix}
#' }
#'
#' @export filterFrequency

filterFrequency <- function(
  img, mask, boldMat, ctxVox, goodtimes, badtimes,
  lowFreq = 0.009, highFreq = 0.08, smoothFactor = 6.0, providePlot = TRUE){
  #check requirements
  if (!usePkg("ANTsR")) {
    print("Need ANTsR package")
    return(NULL)
  }

  if (!usePkg("R.cache")) {
    print("Need R.cache package")
    return(NULL)
  }

  if(!usePkg("pracma")){
    print("Need pracma package")
    return(NULL)
  }

  if(!usePkg("mFilter")){
    print("Need pracma Package")
    return(NULL)
  }

  if (providePlot){
    if (!usePkg("ggplot2")) {
      print("Need ggplot2 package")
      return(NULL)
    }
  }


  if (!is.antsImage(img)) {
    stop("ERROR: Input image must be of classs 'antsImages'")
  }

  if(!is.antsImage(mask)){
    stop("ERROR: Input mask must be of classs 'antsImages'")
  }

  if (img@dimension != 4) {
    stop("ERROR: Dimension of input image is not 4D")
  }

  #required inputs
  nTimes <- length(goodtimes) + length(badtimes)
  nVox = length(which(ANTsR::as.array(mask) == 1))
  
  #test if cache for image spacing exists
  if (is.null(loadCache(key = list("tr"))) ){
    tr <- antsGetSpacing(img)[4]
    saveCache(tr,list("tr"))
  } else {
    tr <- loadCache(key = list("tr"))
  }

  ##FILTER##
  #if there is bad timepoints
  if(!(is.integer(badtimes) && length(badtimes) == 0)){
    if ( length(badtimes) > 0 ) {
      for ( v in c(1:nVox) ) {
        boldMat[badtimes,v]=spline( c(1:nTimes)[goodtimes], boldMat[goodtimes,v], xout=badtimes )$y
      }
    }
  }

  # save interpolated values for plotting
  if(providePlot){
    ctxMeanSpline = rowMeans(boldMat[,ctxVox])
  }

  boldMat <- frequencyFilterfMRI( boldMat, tr=tr, freqLo=lowFreq, freqHi=highFreq, opt="trig" )

  # save filtered values for plotting
  if(providePlot){
    ctxMeanFiltered = rowMeans(boldMat[,ctxVox])
    #if there is bad timepoints
    if(!(is.integer(badtimes) && length(badtimes) == 0)){
      ctxMeanFiltered[badtimes] = NA
    }
  }

  ##SMOOTH##
  #dont smooth time dimension
  sf <- c(smoothFactor,smoothFactor,smoothFactor,0)

  img <- matrix2timeseries( img, mask, boldMat )
  img <- smoothImage(img, sf, FWHM=TRUE )
  boldMat <- timeseries2matrix(img, mask)

  if(providePlot){
    ctxMeanSmoothed = rowMeans(boldMat[,ctxVox])
    if(!(is.integer(badtimes) && length(badtimes) == 0)){
      ctxMeanSmoothed[badtimes] = NA
    }
  }

  #plot
  if(providePlot){
    freq.dat =  data.frame( Time=rep(1:nTimes,3) )
    freq.dat$Values = c(ctxMeanSpline, ctxMeanFiltered, ctxMeanSmoothed)
    freq.dat$Type = c(rep("Original",nTimes), rep("Filtered",nTimes), rep("Smoothed",nTimes))
    freq.dat$Data = freq.dat$Type
    freq.dat$Data[badtimes] = "Interpolated"
    freq.dat$Type = factor(freq.dat$Type, levels=c("Original", "Filtered", "Smoothed"))
    freq.dat$Data = factor(freq.dat$Data, levels=c("Original", "Interpolated", "Filtered", "Smoothed"))
    freqPlot = ggplot(freq.dat, aes(x=Time, y=Values, group=Type, colour=Data))
    freqPlot = freqPlot + geom_line(size=0.5)
    freqPlot = freqPlot + theme(text=element_text(size=10), legend.position="top")
    freqPlot = freqPlot + ggtitle("Effect of bandpass filtering & spatial smoothing")
    freqPlot = freqPlot + facet_grid(Type ~ ., scales="free" )

    saveCache(freqPlot,key = list("fp"))
  }

  return(list(img=img, boldMat = boldMat))
}




