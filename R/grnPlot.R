#' Gene Regulatory Network Plot
#'
#' This function plot the a cell/tissue-specific gene regulatory network with genes in the five categories,
#' alterntively reducing the network size by plotting given specific genes as nodes.
#' @usage grnPlot(grn.data, cate.gene, filter = TRUE, nodes = NULL, centrality.score,
#' col = NULL, main= NULL, vertex.label =NULL, vertex.label.dist = 0, vertex.label.font = 1,
#' vertex.label.cex = 0.5, vertex.label.color="black", edge.arrow.size = 0.4,
#' edge.color = "grey", layout ="layout_with_lgl", legend.labels = NULL, ...)
#' @param grn.data a data frame with two columns named "TF" and "TG" to specify the genes as transcription regulators (TF)
#' and target genes being regulated (TG).
#' @param cate.gene a list of the five gene categories as nodes in the network, alternatively output by
#' \code{\link{categorizeGene}}.
#' @param filter logical to specify if the network is reduced to less nodes by a filter, logical to \code{TRUE}
#' for a clear visualization.
#' @param nodes a character vector of genes to kept in the network, only applied when \code{filter} is \code{TRUE}.
#' @param centrality.score a vector or data frame of centrality scores for genes in the network, alternatively calculated
#'  by \code{\link{networkAnalyze}} function.
#' @param col colors of gene vertex in each \code{cate.gene}.
#' @param main title of the plot.
#' @param vertex.label,vertex.label.dist,vertex.label.font,vertex.label.cex,vertex.label.color parameters for
#' vertex labels, to specify the labels of vertex, the position of labels on the vertex, font size, cex and color
#' for label, see details in \code{\link[igraph]{igraph.plotting}}.
#' @param edge.arrow.size,edge.color parameters for edge to specify the arrow size and color of edge, see details
#' in \code{\link[igraph]{igraph.plotting}}.
#' @param layout the layout of network plot, see details in \code{\link[igraph]{layout}}.
#' @param legend.labels vector of label names for each \code{cate.gene}.
#' @param ... other parameters used in \code{\link[igraph]{igraph.plotting}}.
#' @return a igraph plot for gene regulatory network.
#' @export
#' @examples
#' \dontrun{
#'   # select genes to shown their regulation with others
#'   node.genes = c("ZNF641", "BCL6")
#'   # enlarge the centrality
#'   centrality.score = degree$centrality*100
#'   names(centrality.score) = degree$Gene
#'   par(mar = c(2,2,3,2))
#'   grnPlot(grn.data = human.grn[[tissue]], cate.gene = cate.gene, filter = TRUE,
#'          nodes = node.genes, centrality.score = centrality.score,
#'          main = "Gene regulatory network")
#' }

grnPlot = function(grn.data, cate.gene, filter = TRUE, nodes = NULL, centrality.score,
                   col = NULL, main= NULL, vertex.label =NULL, vertex.label.dist = 0,
                   vertex.label.font = 1, vertex.label.cex = 0.5, vertex.label.color="black",
                   edge.arrow.size = 0.4, edge.color = "grey", layout ="layout_with_lgl",
                   legend.labels = NULL, ...){
  grn = graph.data.frame(grn.data[,c("TF", "TG")],directed =TRUE)

  if(filter){
    vertex = nodes
    grn.data = grn.data[grn.data[,"TF"] %in% vertex | grn.data[,"TG"] %in% vertex, ]
    grn = graph.data.frame(grn.data[,c("TF","TG")],directed =TRUE)
  }
  vertex= V(grn)$name
  if(is.data.frame(centrality.score)){
    centrality.score = centrality.score[match(vertex, centrality.score$Gene),"centrality.score"]
  }
  else if (is.vector(centrality.score)){
    centrality.score = centrality.score[match(vertex, names(centrality.score))]
  }
  V(grn)$size  = centrality.score

  if(is.null(col)){
    col =  c("#abd9e9", "#2c7bb6", "#fee090", "#d7191c", "#fdae61")
  }

  V(grn)$color = ifelse(V(grn)$name %in% cate.gene[[1]],col[1],
                        ifelse(V(grn)$name %in% cate.gene[[2]],col[2],
                               ifelse(V(grn)$name %in% cate.gene[[3]],col[3],
                                      ifelse(V(grn)$name %in% cate.gene[[4]],col[4],
                                             ifelse(V(grn)$name %in% cate.gene[[5]],col[5],"grey")))))
  V(grn)$frame.color=V(grn)$color

  # plot the GRN
  l = do.call(layout, list(grn))
  if(is.null(vertex.label)){
    vertex.label = V(grn)$name
  }
  plot(grn,layout=l,
       vertex.label.dist = vertex.label.dist,			#puts the name labels slightly off the dots
       vertex.label.color = vertex.label.color,		#the color of the name labels
       vertex.label.font = vertex.label.font,			#the font of the name labels
       vertex.label = vertex.label,		#specifies the lables of the vertices. in this case the 'name' attribute is used
       vertex.label.cex= vertex.label.cex, #specifies the size of the font of the labels
       edge.arrow.size = edge.arrow.size,
       edge.color = edge.color,...)
  if(is.null(main)){
    main = layout
  }
  title(main = main, cex.main =1.2, font=2)

  if(is.null(legend.labels)){
    legend.labels =c("Reversed","Inactive","Insufficient","Successful","Over")
  }
  legend("topleft",legend = legend.labels, pch =20, col= col,
         bty="n", x.intersp=0.2)
}
