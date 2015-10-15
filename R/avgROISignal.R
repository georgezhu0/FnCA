#' Get ROI Average Signals
#'
#' Compute the mean time singal for each ROI using ROI labels.
#'
#'
#' @param labels ROI labels
#' @param labeledBoldMat Labeled BOLD signal matrix
#' @param providePlot Boolean determine whether descriptive statistics should be provided
#' if TRUE, the plot will be cached and can be called
#'
#' @return  Average ROI signal matrix
#'
#' @export avgROISignal


avgROISignal <- function(
  labels,labeledBoldMat,providePlot=TRUE){
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

    #test if cahced tr exists
    if (is.null(loadCache(key = list("tr"))) ){
      stop("ERROR: Image spacing cache not found, run rsBOLDSteadyState First")
    }
  }

  #creating inputs
  nLabels = max(labels)
  roiMat = matrix(0, nrow=dim(labeledBoldMat)[1], ncol=nLabels)
  for ( i in c(1:nLabels) ) {
    if (length(which(labels==i)) > 1 ) {
      roiMat[,i] = rowMeans(labeledBoldMat[,(labels==i)])
    }
  }
  nActualTimes = dim(roiMat)[1]

  #plot
  if(providePlot){
    #required plot input
    tr <- loadCache(key = list("tr"))

    plotMat = roiMat + min(roiMat)
    plotMat = plotMat / max(plotMat)

    means.dat = data.frame(Time=rep( (1:nActualTimes)*tr, nLabels))
    yoffset = (rep( 1:nLabels, each=nActualTimes)-1)*0.5
    means.dat$Signal = (as.vector(plotMat)/2)+yoffset
    means.dat$ID = factor(rep( 1:nLabels, each=nActualTimes))

    meanPlot = ggplot(means.dat, aes(x=Time, y=Signal, group=ID, colour=ID))
    meanPlot = meanPlot + geom_line(size=0.5)
    meanPlot = meanPlot + theme(text=element_text(size=10), legend.position="none", axis.text.y=element_blank())
    meanPlot = meanPlot + ggtitle("Mean BOLD signal in network ROIs")

    saveCache(meanPlot,list("mp"))
  }
  return(roiMat)

}
