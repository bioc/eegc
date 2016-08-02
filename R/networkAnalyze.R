#' Network Topological Analysis
#'
#' This function analyzes the topological structure of gene regulation network (GRN) by calculating the "degree",
#' "betweenness", "closeness" and "stress" parameters, and output the centrality values for given genes in each
#' gene categories.
#' @usage networkAnalyze(grn.data, cate.gene, centrality = c("degree", "betweenness",
#' "stress", "closeness"), mode = c("all","in", "out", "total"))
#' @param grn.data a data frame with two columns named "TF" and "TG" to specify the genes as transcription regulators (TF)
#' and target genes (TG) being regulated.
#' @param cate.gene a list of the five gene categories as nodes in the network, alternatively output by
#' \code{\link{categorizeGene}}.
#' @param centrality charactor string of "degree", "betweenness", "closeness" and "stress" to calculate the
#' centrality of network built from input \code{grn.data}, see \code{\link[igraph]{degree}}.
#' \code{\link[igraph]{betweenness}}, \code{\link[igraph]{closeness}}, and \code{\link[sna]{stresscent}} for details.
#' @param mode character string of "all", "in", "out" and "total", only used when \code{centrality} is "degree"
#' or "closeness", see \code{\link[igraph]{degree}}, \code{\link[igraph]{closeness}} for details.
#' @importFrom sna stresscent
#' @return data frame with genes and centrality scores.
#' @export
#' @examples
#' # load the CellNet GRN and gene categories
#' data(human.grn)
#' data(cate.gene)
#'
#' # specify a tissue-specifc network
#' tissue = "Hspc"
#' degree = networkAnalyze(human.grn[[tissue]], cate.gene = cate.gene,
#'                        centrality = "degree", mode ="all")

networkAnalyze = function(grn.data, cate.gene, centrality = c("degree", "betweenness", "stress", "closeness"),
                          mode = c("all", "in", "out", "total")){
  #classify the genes into TF, TG, or TF&TG
  tf = unique(grn.data[,"TF"])
  tg = unique(grn.data[,"TG"])

  # the gene is both TF and TG
  tfg = intersect(tf,tg)
  diff.tf = setdiff(tf, tfg)
  diff.tg = setdiff(tg, tfg)
  if(length(tfg) == 0){
    gene = list(diff.tf, diff.tg)
  }
  else{
    gene = list(tfg, diff.tf, diff.tg)
  }

  frame = function(x,type ="TF",name="Type"){
    data = data.frame(Gene =x,stringsAsFactors =FALSE)
    data[,ncol(data)+1] = type
    colnames(data)[ncol(data)] = name
    data
  }

  if(length(tfg) == 0){
    type = c("TF","TG")
  }
  else{
    type = c("TF;TG","TF","TG")
  }
  for(i in 1:length(type)){
    if(length(gene[[i]])> 0){
      gene[[i]] = frame(gene[[i]],type = type[i],"Type")
    }
  }
  gene = do.call("rbind",gene)

  #classify the genes into five categories (Successful, Inactive et al.)

  cate.names = names(cate.gene)
  if(is.null(cate.names)){
    category = c("Reversed","Inactive","Insufficient","Successful","Over")
  }
  else{
    category = cate.names
  }
  for(i in 1:length(cate.gene)){
    cate.gene[[i]] = frame(cate.gene[[i]],type = category[i],"Category")
  }
  cate.gene = do.call("rbind",cate.gene)
  gene = merge(gene,cate.gene,by="Gene",all=TRUE)
  gene = gene[!is.na(gene$Type),]

  #calculate the degree, betweenness, closeness, stress
  all.ppi = graph.data.frame(grn.data[,c("TF","TG")],directed =TRUE)

  mode = match.arg(mode,  c("all", "in", "out","total"))
  if(centrality == "degree"){
    centrality.score  = igraph::degree(all.ppi, mode = mode,normalized = TRUE)
  }
  if(centrality == "closeness"){
    centrality.score = igraph::closeness(all.ppi,mode= mode, normalized= TRUE )
  }
  if(centrality == "betweenness"){
    centrality.score = igraph::betweenness(all.ppi,normalized =TRUE)
  }
  if(centrality == "stress"){
    mat= as_adjacency_matrix(all.ppi,type= "both",attr=NULL, edges=FALSE,
                             names=TRUE,sparse= FALSE)
    centrality.score = 2 * stresscent(mat, g=1, gmode= "graph",
                                      cmode= "undirected", ignore.eval=TRUE)
    names(centrality.score) = rownames(mat)
  }

  centrality.score = data.frame(centrality.score,stringsAsFactors = FALSE)
  centrality.score$Gene = rownames(centrality.score)
  centrality.score = merge(centrality.score,gene,by="Gene")
  centrality.score = centrality.score[order(centrality.score$centrality.score, decreasing = TRUE),]

  return(centrality = centrality.score)
}



