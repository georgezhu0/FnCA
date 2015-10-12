#' Output Descriptive Statistics from BOLD processing
#'
#' Provide all the plots associated with each steps of BOLD processing.
#' The plots are cached from function call
#'
#'
#' @param outdir string indicating the output directory
#' @param prefix string indicating the prefix of output filename
#'
#' @export exportStats

exportStats <- function(
  outdir,prefix){
  #check requirements
  if (!usePkg("R.cache")) {
    print("Need R.cache package")
    return(NULL)
  }

  if (!usePkg("ggplot2")) {
    print("Need ggplot2 package")
    return(NULL)
  }

  #hard code in the keys
  cacheKey <- c("ss","rp","bp","dp","tp","cp","fp","mp","sp","ap","np")

  #steadystate
  if(!is.null(loadCache(key = list("ss")))){
    temp <- loadCache(key = list("ss"))
    ggsave(plot = temp, filename = paste(outdir,prefix,"_Steady_State_Plot.png",sep = ""),width = 5, height = 3,dpi = 144)
  }

  #MOCO parameters
  if(!is.null(loadCache(key = list("rp")))){
    temp <- loadCache(key = list("rp"))
    ggsave(plot = temp ,filename = paste(outdir,prefix,"_MOCO_Params_Plot.png",sep = ""),width = 6, height = 3,dpi = 144)
  }

  #bad time points
  if(!is.null(loadCache(key = list("bp")))){
    temp <- loadCache(key = list("bp"))
    ggsave(plot = temp ,filename = paste(outdir,prefix,"_Bad_Timepoints_Plot.png",sep = ""),width = 6, height = 3,dpi = 144)
  }

  #detrend plot
  if(!is.null(loadCache(key = list("dp")))){
    temp <- loadCache(key = list("dp"))
    ggsave(plot = temp ,filename = paste(outdir,prefix,"_Time_Series_Detrended_Plot.png",sep = ""),width = 6, height = 3,dpi = 144)
  }

  #Nuisance Params
  if(!is.null(loadCache(key = list("tp")))){
    temp <- loadCache(key = list("tp"))
    ggsave(plot = temp ,filename = paste(outdir,prefix,"_Nuisance_Params_Plot.png",sep = ""),width = 6, height = 3,dpi = 144)
  }

  #regressed nuisance params
  if(!is.null(loadCache(key = list("cp")))){
    temp <- loadCache(key = list("cp"))
    ggsave(plot = temp ,filename = paste(outdir,prefix, "__Nuisance_Params_Regression_Plot.png",sep = ""),width = 6, height = 3,dpi = 144)
  }

  #filtered
  if(!is.null(loadCache(key = list("fp")))){
    temp <- loadCache(key = list("fp"))
    ggsave(plot = temp ,filename = paste(outdir,prefix,"_Filter_Effect_Plot.png",sep = ""),width = 6, height = 4,dpi = 144)
  }

  #Mean BOLD signal
  if(!is.null(loadCache(key = list("mp")))){
    temp <- loadCache(key = list("mp"))
    ggsave(plot = temp ,filename = paste(outdir,prefix,"_BOLD_Signal_ROI_Plot.png",sep = ""),width = 6, height = 6,dpi = 144)
  }

  if(!is.null(loadCache(key = list("sp")))){
    temp <- loadCache(key = list("sp"))
    ggsave(plot = temp ,filename = paste(outdir,prefix,"_BOLD_Signal_System_Plot.png",sep = ""),width = 6, height = 6,dpi = 144)
  }

  if(!is.null(loadCache(key = list("ap")))){
    temp <- loadCache(key = list("ap"))
    ggsave(plot = temp ,filename = paste(outdir,prefix,"_System_Conn_Mat_Plot.png",sep = ""),width = 6, height = 4,dpi = 144)
  }

  if(!is.null(loadCache(key = list("np")))){
    temp <- loadCache(key = list("np"))
    ggsave(plot = temp ,filename = paste(outdir,prefix,"_Node_Metrics_Plot.png",sep = ""),width = 6, height = 4,dpi = 144)
  }

  return(NULL)
}
