#' Calculate Framewise Displacement of BOLD Image
#'
#' Measures the average displacement of voxels betweeen consecutive time points
#'
#' @param  img motion corrected BOLD fMRI time series image
#' @param  mask mask for image(3D)
#' @param  param MOCO registraion parameters data frame
#' @param  providePlot Boolean determine whether descriptive statistics should be provided
#' if TRUE, the plot will be cached and can be called
#' @return a dataframe consisting mean displacement and max displacement
#'
#'
#'
#'
#' @export getDisplacement


getDisplacement <- function(
  img,mask,param,providePlot = TRUE){

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

  if (!is.antsImage(img)) {
    stop("ERROR: Input image must be of classs 'antsImages'")
  }

  if(!is.antsImage(mask)){
    stop("ERROR: Input mask must be of classs 'antsImages'")
  }

  if (img@dimension != 4) {
    stop("ERROR: Dimension of input image is not 4D")
  }

  tsimg <- antsImageClone(img,"float")
  mocostats<-.Call("antsMotionCorrStats",tsimg, mask, as.matrix(param),PACKAGE = "ANTsR")
  displacement <- as.data.frame(mocostats$Displacements)
  names(displacement) <- c("MeanDisplacement","MaxDisplacement")

  #plot
  if(providePlot){
    #required plot inputs
    meanDisplacement <- displacement$MeanDisplacement
    reg_params <- as.matrix(param[,3:8])
    nTimes <- dim(reg_params)[1]
    tr <- loadCache(key = list("tr"))

    orderedBreaks = c("Framewise", "X", "Y", "Z", "Pitch", "Roll", "Yaw" )
    moco.dat <- data.frame(Time=rep(1:nTimes, 7)*tr)
    moco.dat$Values = c( as.vector(reg_params), meanDisplacement )
    moco.dat$Category = c( rep("Angle", 3*nTimes), rep("Displacement", 4*nTimes) )
    moco.dat$Type = rep(c("Pitch", "Roll", "Yaw","X", "Y", "Z", "Framewise"), each=nTimes)

    regPlot <- ggplot(moco.dat, aes(x=Time, y=Values, group=Type, colour=Type) )
    regPlot <- regPlot + geom_line(size=0.5)
    regPlot <- regPlot + theme(text=element_text(size=10), legend.position="top")
    regPlot <- regPlot + ggtitle("Motion correction parameters")
    regPlot <- regPlot + facet_grid(Category ~ ., scales="free" )
    regPlot <- regPlot + scale_color_discrete(breaks=orderedBreaks)
    saveCache(regPlot, key = list("rp"))
  }

  return(displacement)
}
