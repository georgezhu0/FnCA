#' Identify Bad Timpoints in BOLD Image
#'
#' Use Motion Registration Parameters to identify bad timepoints in time series.
#' The bad timepoints are deteremined by mean framewise displacement.
#'
#' @param img BOLD fMRI time series image
#' @param meanDisplacement mean framewise displacement of BOLD image
#' @param threshold threshold of framewise displacement, default is 0.1mm
#' @param providePlotBoolean determine whether descriptive statistics should be provided
#' if TRUE, the plot will be cached and can be called
#' @return list of outputs containing:
#' \itemize{
#' \item{good:}{ Index of good timepoints}
#' \item{bad:}{ Index of bad timepoints}
#' }
#'
#'
#' @export findBadTimepoints

findBadTimepoints <- function(
  img,meanDisplacement,threshold=0.1,providePlot = TRUE){

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

  if (!is.antsImage(img)) {
    stop("ERROR: Input image must be of classs 'antsImages'")
  }

  if (img@dimension != 4) {
    stop("ERROR: Dimension of input image is not 4D")
  }

  #required inputs

  if (is.null(loadCache(key = list("tr"))) ){
    tr <- antsGetSpacing(img)[4]
    saveCache(tr,list("tr"))
  } else {
    tr <- loadCache(key = list("tr"))
  }
  nTimes = length(meanDisplacement)

  #compare meandisplacement with threshhold
  flag <- TRUE
  badtimes = which(meanDisplacement > threshold)


  if(is.integer(badtimes) && length(badtimes) == 0){
    flag <- FALSE
  }

  if(flag){
    badtimes = sort(c(badtimes, badtimes+1))
    goodtimes = (1:nTimes)[-badtimes]

    badstarts = which(meanDisplacement > threshold)

    if(providePlot){
      bad.data = data.frame(Time=(1:nTimes)*tr)
      bad.data$FD = meanDisplacement

      bad.data.rect = data.frame(Start=badstarts*tr)
      bad.data.rect$Stop = (badstarts+1)*tr
      rect.aes = aes(xmin=Start,xmax=Stop,ymin=-Inf,ymax=Inf,fill="pink",alpha=0.2)

      badPlot <- ggplot(bad.data) + geom_line(aes(x=Time, y=FD))
      badPlot <- badPlot + geom_hline( yintercept=0.1, linetype="dashed", alpha=0.5 )
      badPlot <- badPlot + theme(text=element_text(size=10), legend.position="none")
      badPlot <- badPlot + ggtitle("Bad timepoints")
      badPlot <- badPlot + geom_rect(data=bad.data.rect, rect.aes)
      saveCache(badPlot, key = list("bp"))
    }
  }


  #if no ba timpoints found
  if(!flag){
    goodtimes = (1:nTimes)

    if(providePlot){
      warning("No bad time points found, no plot generated")
    }
  }

  return(list(good = goodtimes, bad = badtimes))
}
