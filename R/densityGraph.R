#' Create Network Density Graph
#'
#' Binarize ROI correlation matrix by keeping certain edgeds and create an
#' adjacency matrix. The adjacency matrix is then used to create network
#' density graph from graph operations. Default density is 0.1
#'
#' @param corMat ROI correlation matrix
#' @param pts labels identifying anatomical regions of interest
#' @param % of edges to retain in the correlation matrix, default is 0.1 (10%)
#' @param Boolean determine whether descriptive statistics should be provided
#' if TRUE, the plot will be cached and can be called
#' @return adjacency matrix
#'
#' @export densityGraph

densityGraph <- function(
  corMat,pts,density = 0.1, providePlot = TRUE){
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

  #creating adjacency matrix
  connMat <- corMat
  nEdges = length(upper.tri(connMat))*density
  thresh = sort( connMat[upper.tri(connMat)], decreasing=T)[nEdges]
  adj = 1*(connMat >= thresh)

  bingraph = graph.adjacency(adj, mode="undirected", weighted=NULL)
  components = clusters(bingraph)
  maxID = which(components$csize == max(components$csize))[1]

  adj[components$membership!=maxID,] = 0
  adj[,components$membership!=maxID] = 0
  bingraph = graph.adjacency(adj, mode="undirected", weighted=NULL)

  if (providePlot){
    pts$System = factor(pts$System, levels=c("Sensory/Somatomotor Hand", "Sensory/Somatomotor Mouth", "Cingulo-opercular Task Control", "Auditory", "Default Mode", "Memory Retrieval", "Visual", "Fronto-parietal Task Control", "Salience", "Subcortical", "Ventral Attention", "Dorsal Attention", "Cerebellar", "Uncertain"))

    graph = graph.adjacency( adj, mode="directed", weighted=NULL )
    V(graph)$name = pts$ROI
    V(graph)$comm = pts$System
    V(graph)$degree = degree(graph)

    systems = levels(pts$System)
    systemNames = as.character(systems)

    node_list <- get.data.frame(graph, what = "vertices")
    edge_list <- get.data.frame(graph, what = "edges") %>%
      inner_join(node_list %>% select(name, comm), by = c("from" = "name")) %>%
      inner_join(node_list %>% select(name, comm), by = c("to" = "name")) %>%
      mutate(group = ifelse(comm.x == comm.y, comm.x, NA) %>% factor())

    all_nodes <- sort(node_list$name)
    plot_data <- edge_list %>% mutate(
      to = factor(to, levels = all_nodes),
      from = factor(from, levels = all_nodes))

    name_order <- (node_list %>% arrange(comm))$name
    plot_data <- edge_list %>% mutate(
      to = factor(to, levels = name_order),
      from = factor(from, levels = name_order))

    plot_data$group = as.integer(plot_data$group)
    for ( i in 1:length(systems) ) {
      plot_data$group[ which( plot_data$group == i) ] = as.character( systems[i] )
    }

    #legends
    lut = c("Sensory/Somatomotor Hand"="cyan3", "Sensory/Somatomotor Mouth"="orange", "Cingulo-opercular Task Control"="purple", "Auditory" = "pink2", "Default Mode"="red", "Memory Retrieval"="gray50", "Visual"="blue", "Fronto-parietal Task Control"="yellow2", "Salience"="black", "Subcortical"="chocolate4", "Ventral Attention"="aquamarine4", "Dorsal Attention"="green", "Cerebellar"="cadetblue1", "Uncertain"="peachpuff2" )
    adjplot = ggplot(plot_data, aes(x = from, y = to, fill = group)) + geom_raster() + theme_bw() + scale_x_discrete(drop = FALSE) + scale_y_discrete(drop = FALSE) + theme( axis.title=element_blank(), axis.ticks=element_blank(), axis.text = element_blank(),  aspect.ratio = 1 ) + scale_fill_manual( values = lut, na.value="gray80", name="System",  breaks=systemNames, drop=FALSE )

    saveCache(adjplot, key = list("ap"))
  }
  return(adj)
}
