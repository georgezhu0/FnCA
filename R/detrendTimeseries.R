#' Demand and Detrend the Time-series Data
#'
#' Demean and detrend BOLD time series while excluding bad timepoints in BOLD
#' images
#'
#'
#' @param img BOLD fMRI time series image
#' @param mask for image (3D)
#' @param goodtimes list of interger index of good timepoints
#' @param badtimes list of interger index of bad timepoints
#' @param providePlot Boolean determine whether descriptive statistics should be provided
#' if TRUE, the plot will be cached and can be called
#'
#' @return BOLD signal matrix
#'
#' @export detrendTimeseries

detrendTimeseries <- function(
  img, mask , goodtimes, badtimes, providePlot = TRUE){

  #check requirements
  if (!usePkg("ANTsR")) {
    print("Need ANTsR package")
    return(NULL)
  }

  if (!usePkg("pracma")) {
    print("Need pracma package")
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

  if (!is.antsImage(img)) {
    stop("ERROR: Input image must be of classs 'antsImages'")
  }

  if(!is.antsImage(mask)){
    stop("ERROR: Input mask must be of classs 'antsImages'")
  }

  if (img@dimension != 4) {
    stop("ERROR: Dimension of input image is not 4D")
  }

  #creating bold matrix
  global_moco <- rowMeans(timeseries2matrix(img, mask))
  boldMat <- timeseries2matrix(img, mask)

  #detrend
  boldMat[goodtimes,] = detrend(boldMat[goodtimes,])

  #if there is bad timepoints
  if(!(is.integer(badtimes) && length(badtimes) == 0)){
    boldMat[badtimes,] = NA
  }


  #plot
  if(providePlot){
    #required plot inputs
    nTimes <- length(goodtimes)+length(badtimes)

    #global MOCO detrend 
    global_moco_detrend <- rowMeans(boldMat)
    if(!(is.integer(badtimes) && length(badtimes) == 0)){
      global_moco[badtimes] = NA
    }

    trend.dat = data.frame( Time=rep(1:nTimes,2) )
    trendType = c( rep("Motion-corrected",nTimes) )
    trendType = c(trendType, rep("Moco & Detrended",nTimes) )
    trendNames = c(rep("Original",nTimes), rep("Detrended", nTimes))
    trendCategory = factor(trendNames, levels=c("Original", "Detrended"))
    trend.dat$Signal = c(global_moco, global_moco_detrend)
    trend.dat$Type = trendType
    trend.dat$Category = trendCategory
    trendPlot <- ggplot(trend.dat, aes(x=Time, y=Signal, group=Type, colour=Type) )
    trendPlot <- trendPlot + geom_line(size=0.5)
    trendPlot <- trendPlot + theme(text=element_text(size=10), legend.position="top")
    trendPlot <- trendPlot + facet_grid(Category ~ ., scales="free" )
    trendPlot <- trendPlot + ggtitle("Detrending the time-series")
    saveCache(trendPlot, key = list("dp"))
  }

  return(boldMat)
}
