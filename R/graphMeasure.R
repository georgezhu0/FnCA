#' Examine Graph-metrics from Node and/or Global Measures
#'
#' Obtain graph-metrics from node measures to examine various properties of
#' functional network
#'
#' @param adjMat Network adjacency matrix
#' @param globalMeasure Boolean determining whether the entire network in a single measure
#' should bed summarized
#' @param providePlot Boolean determining whether descriptive statistics should be provided
#' if TRUE, the plot will be cached and can be called
#' @return list of node measures, containing:
#' \itemize{
#' \item{Degree：}{ Number of links connected to a node}
#' \item{localClustCoeff：}{ Local neighborhood connectivity}
#' \item{pathLen:}{ Mean shortest distance between one node and all others}
#' \item{localEff:}{ Closeness of nodes in a neighborhood}
#' \item{pageRank: }{ Google page rank measure}
#' \item{(OPTIONAL) globalClustCoeff： }{ Closeness of nodes to all other nodes}
#' \item{(OPTIONAL) globalEff： }{ Global neighborhood connectivity}
#'}
#'
#' @export graphMeasure

graphMeasure <- function(
  adjMat,globalMeasure = TRUE, providePlot = TRUE){
  #check requirements
  if (!usePkg("ANTsR")) {
    print("Need ANTsR package")
    return(NULL)
  }
  if (!usePkg("R.cache")) {
    print("Need R.cache package")
    return(NULL)
  }
  if (!usePkg("igraph")) {
    print("Need igraph package")
    return(NULL)
  }

  if (!usePkg("dplyr")){
    print("Need dplyr package")
    return(NULL)
  }

  if (providePlot){
    if (!usePkg("ggplot2")) {
      print("Need ggplot2 package")
      return(NULL)
    }
  }

  adj <- adjMat
  graph = graph.adjacency( adj, mode="undirected", weighted=NULL )

  #degree
  deg = degree(graph)
  deg[deg==0] = NA

  #path length
  pathsmat =  shortest.paths(graph, weights=NA)
  pathsmat[!is.finite(pathsmat)] = NA
  paths = rowMeans(pathsmat, na.rm=TRUE)
  paths[paths==0] = NA

  #clustering coefficient
  clust = transitivity(graph, type="local")
  clust[deg < 2] = NA

  #pagerank
  pager = page.rank(graph)$vector
  pager[deg < 2] = NA

  #local efficiency
  leff <- numeric(length(deg))
  goodnodes <- which(deg > 1)
  leff[goodnodes] <- sapply(goodnodes, function(x) {
    neighbs <- neighbors(graph, v=x)
    g.sub <- induced.subgraph(graph, neighbs)
    Nv <- vcount(g.sub)

    lpaths <- shortest.paths(g.sub, weights=NA)
    lpaths <- paths[upper.tri(lpaths)]

    pathsup <- lpaths[upper.tri(lpaths)]
    2 / Nv / (Nv - 1) * sum(1 / lpaths[which(is.na(lpaths)==FALSE)])
  }
  )
  leff[deg < 2] = NA
  leff[which(is.na(deg)==TRUE)] = NA

  #plot
  if (providePlot){
    nNodes = length(deg)

    cnode.dat = data.frame(Node=rep(1:nNodes,5))
    cnode.dat$Value = c( deg, paths, leff, clust, pager )
    cnode.dat$Metric = c( rep("Degree", nNodes), rep("Shortest Path", nNodes), rep("Local Efficiency", nNodes), rep("Clustering Coefficient", nNodes), rep("Page-Rank", nNodes) )

    cnodePlot = ggplot(cnode.dat, aes(x=Node, y=Value, group=Metric, fill=Metric, colour=Metric))
    cnodePlot = cnodePlot + geom_point()
    cnodePlot = cnodePlot + theme(text=element_text(size=10), legend.position="none")
    cnodePlot = cnodePlot + ggtitle("Node metrics")
    cnodePlot = cnodePlot + facet_grid(Metric~ ., scales="free")
    saveCache(cnodePlot, key = list("np"))
  }

  #global measures
  if(globalMeasure){
    #global efficiency
    geff<-1/(shortest.paths(graph))
    geff[!is.finite(geff)]<-NA
    geff<-mean(geff,na.rm=TRUE)
    #global cluster coeff
    gcc = transitivity(graph)

    return(list(Degree = deg,
                localClustCoeff = clust,
                pathLen = paths,
                localEff = leff,
                pageRank = pager,
                globalClustCoef = gcc,
                globalEff = geff))
  }
  #if global is optional
  return(list(Degree = deg,
              localClustCoeff = clust,
              pathLen = paths,
              localEff = leff,
              pageRank = pager))
}

