#' Remove pre Steady State Timepoints
#'
#' Identify initial time points that occur before magnetization steady state
#' in BOLD fMRI and remove pre steady state timepoints. Default is set to
#' excluding the first 10 seconds of the data.
#'
#' @param img BOLD fMRI time series image
#' @param mask mask for BOLD image (3D)
#' @param mocoParams if moco parameters data includes non steady state timepoints, the data
#' should be provided. Otherwise, no input is required
#' @param preSteady time points to exclude, default is set to 10 seconds
#' @param providePlot Boolean determine whether descriptive statistics should be provided
#' if TRUE, the plot will be cached and can be called
#' @return list of outputs containing:
#' \itemize{
#' \item{ss_img}{Steady State BOLD image}
#' \item{OPTIONAL: ss_MOCOParams}{Steady State MocoParameters}
#' }
#'
#' @examples
#' fmri <- antsImageRead(getANTsRData("rsbold"))
#' msk <- antsImageRead(getANTsRData("rsboldmask"))
#' ssResult <- rsBOLDSteadyState(fmri, msk)
#'
#' @export findSteadyState

findSteadyState <- function(
  img,mask,mocoParams = NULL,preSteady = 10,providePlot = TRUE){
  #check requirements
  if (!usePkg("ANTsR")) {
    print("Need ANTsR package")
    return(NULL)
  }

  if (providePlot){
    if (!usePkg("ggplot2")) {
      print("Need ggplot2 package")
      return(NULL)
   }
  }

  if (!usePkg("R.cache")) {
    print("Need R.cache package")
    return(NULL)
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

  if (!is.numeric(preSteady)){
    print("'preSteady' is not a numeric type")
    return(NULL)
  }

  if(preSteady <= 0){
    print("'preSteady' needs to be greater than 0")
    return(NULL)
  }

# Find first steady state timepoint

##CACHE FOR GLOBAL USAGE#######
  tr = antsGetSpacing(img)[4]
  saveCache(tr,list("tr"))
###############################

  steady = floor(preSteady / tr) + 1

  fullmean = rowMeans(timeseries2matrix(img, mask))
  allTimes = dim(img)[4]

  # Eliminate non steady-state timepoints in img and MOCO params
  img = cropIndices(img, c(1,1,1,steady), dim(img) )
  if(!is.null(mocoParams)){
    mocoCrop <- allTimes - dim(img)[4]
    if (mocoCrop != 0){
      mocoParams <- mocoParams[-1:-mocoCrop,]
    }
  }


  #creating dataframe for plot
  if(providePlot){
    noss.data <- data.frame(Start=0)
    noss.data$Stop <- (steady-1)*tr
    noss.rect.aes <- aes(xmin=Start,xmax=Stop,ymin=-Inf,ymax=Inf,fill="pink",alpha=0.2)

    ss.dat <- data.frame(Time=rep(1:allTimes)*tr)
    ss.dat$Values <- fullmean

    #generating plot
    ssPlot <- ggplot(ss.dat)
    ssPlot <- ssPlot + geom_line(aes(x=Time, y=Values), size=0.5)
    ssPlot <- ssPlot + geom_rect(data=noss.data, noss.rect.aes)
    ssPlot <- ssPlot + theme(text=element_text(size=10), legend.position="none")
    ssPlot <- ssPlot + ggtitle("Points Prior to Steady State")

    saveCache(ssPlot,key = list("ss"))
  }
  #if params exist
  if(!is.null(mocoParams)){
    return(list(ss_img = img,ss_MOCOParams = mocoParams))
  }

  return(list(ss_img = img))
}
