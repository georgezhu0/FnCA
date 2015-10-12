#' Comparing ROI signal between systems
#'
#' Given ROI system labels, provide the mean and the standard deviation for
#' BOLD signal within system. Provides insight on how the systems are working
#' together
#'
#' @param labels labels
#' @param labeledBoldMat labeled BOLD time series matrix
#' @param roiMat ROI matrix
#' @param pts labels identifying anatomical regions of interest
#' @param providePlotprovidePlot Boolean determine whether descriptive statistics should be provided
#' if TRUE, the plot will be cached and can be called
#'
#' @return list of output containing
#'
#' @export systemROI


systemROI <- function(
  labels, labeledBoldMat, roiMat, pts, providePlot = TRUE){
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
      stop("ERROR: Image spacing cache not found, run findSteadyState First")
    }
  }

  if(!is.antsImage(mask)){
    stop("ERROR: Input mask must be of classs 'antsImages'")
  }


  systemNames = levels(pts$System)

  nSystems = length(systemNames)
  sysMatMean = matrix(0, nrow=dim(labeledBoldMat)[1], ncol=nSystems)
  sysMatSD = matrix(0, nrow=dim(labeledBoldMat)[1], ncol=nSystems)

  systems = pts$System[labels]

  for ( i in 1:nSystems ) {
    sys = systemNames[i]
    sysIdx = which(systems==sys)
    if ( length(sysIdx) > 0)
    {
      sysMatMean[,i] = rowMeans(labeledBoldMat[,sysIdx])
      sysMatSD[,i] = apply(labeledBoldMat[,sysIdx], 1, sd)
    }
  }

  if(providePlot){
    #required inputs
    tr <- loadCache(key = list("tr"))
    nActualTimes <- dim(roiMat)[1]

    #name and color
    systemNickNames = c("Motor/Hand", "Motor/Mouth", "CO-Task", "Auditory", "Default", "Memory", "Visual", "FP-Task", "Salience", "Subcortical", "V Attention", "D Attention", "Cerebellar", "Uncertain" )
    lut = list("Motor/Hand"="cyan3", "Motor/Mouth"="orange", "CO-Task"="purple", "Auditory"="pink2", "Default"="red", "Memory"="gray50", "Visual"="blue", "FP-Task"="yellow2", "Salience"="black", "Subcortical"="chocolate4", "V Attention"="aquamarine4", "D Attention"="green", "Cerebellar"="cadetblue1", "Uncertain"="peachpuff2" )

    sys.dat = data.frame(Time=rep( (1:nActualTimes)*tr, nSystems))
    sys.dat$Signal = as.vector(sysMatMean)
    sys.dat$System = factor( rep( systemNickNames, foreach=nActualTimes), levels=systemNickNames)
    sys.dat$Lower = as.vector(sysMatMean) - as.vector(sysMatSD)
    sys.dat$Upper = as.vector(sysMatMean) + as.vector(sysMatSD)

    sysPlot = ggplot(sys.dat)
    sysPlot = sysPlot + geom_line(aes(x=Time, y=Signal, group=System), size=0.5)
    sysPlot = sysPlot + geom_ribbon(aes(x=Time, ymin=Lower, ymax=Upper, alpha=0.05, fill=System))
    sysPlot = sysPlot + scale_fill_manual( values = lut, na.value="gray80", name="System", breaks=systemNickNames, drop=FALSE)
    sysPlot = sysPlot + theme(text=element_text(size=10), legend.position="none")
    sysPlot = sysPlot + ggtitle("Mean BOLD signal in systems")
    sysPlot = sysPlot + facet_grid(System ~ ., scales="free" )

    saveCache(sysPlot,key = list("sp"))

  }
  return(list( sysMatMean =  sysMatMean, sysMatSD = sysMatSD))
}
