#' Heatmap Plot of Enriched Terms
#'
#' This function plot the significantly enriched terms in a heatmap by calling \pkg{pheatmap} package.
#' @usage heatmapPlot(enrichresult, GO = FALSE, terms = NULL, padjust =TRUE, pvalue = 0.05,
#' top= NA, filter = FALSE, main = NA, annotation = NULL, annotation_col =NULL,
#' annotated_row = FALSE, annotation_row =NULL,annotation_colors=NULL,
#' display_numbers=FALSE, annotation_legend=FALSE,...)
#' @param enrichresult a list of data frames with enrichment results, alternatively output by \code{\link{functionEnrich}}
#' or \code{\link{enrichment}}.
#' @param GO logical to determine whether the terms are Gene Ontology(GO) terms enriched by \code{\link{functionEnrich}}.
#' @param terms a character vector to specify the terms chosen to be listed in the heatmap, selected from the \code{enrichment} result,
#' or a data frame with terms and corresponding term annotations used for \code{annotation_row}.
#' @param padjust logical to determine whether the significantly enriched terms were selected by adjusted p.value
#' rather than the p.value, default to \code{TRUE}.
#' @param pvalue a cutoff value to select the significantly enriched terms.
#' @param top a number to specify the most significantly enriched terms to be drawn in each category, default to \code{NA}
#' without specifying.
#' @param filter logical to specify whether the enriched terms need to be filered with the ones which are significant among
#' the first four categories.
#' @param main a character of main title on the heatmap plot.
#' @param annotated_row logical to determine whether the the rows are annotated by \code{annotation_row}, default to \code{FALSE}.
#' When it's \code{TRUE}, \code{annotation_row} should be specified or using annotations in a data frame of \code{terms}.
#' @param annotation,annotation_row,annotation_col,annotation_colors,annotation_legend,... see details in \code{\link[pheatmap]{pheatmap}}.
#' @param display_numbers logical to determine whether the gene counts number is shown on the heatmap.
#' @importFrom gplots colorpanel
#' @return heatmap plot and the terms with p.values for the heatmap
#' @export
#' @examples
#' # plot the enrichment results by the five gene categories
#' data(goenrich)
#' heatmaptable = heatmapPlot(goenrich, GO = TRUE, top = 5, filter = FALSE,
#'                            main = "Gene ontology enrichment",
#'                            display_numbers =  FALSE)


heatmapPlot = function(enrichresult, GO =FALSE, terms = NULL, padjust =TRUE, pvalue = 0.05, top= NA,
                       filter = FALSE, main = NA, annotation = NULL, annotation_col =NULL,
                       annotated_row = FALSE, annotation_row =NULL,annotation_colors=NULL,
                       display_numbers=FALSE, annotation_legend=FALSE,...){

  if(GO){
    enrichresult = lapply(enrichresult,function(x){
      if("Description" %in% colnames(x)){
        x$Term = apply(x[,c("ID","Description")],1,paste,collapse="~")
      }
      x})}
  else{
    enrichresult = lapply(enrichresult,function(x){
      x$Term = x$Row.names
      x})}
  if(!is.null(terms)){
    if(is.vector(terms)){
      enrichresult = lapply(enrichresult, function(x){
        x = x[x[,1] %in% terms,]
        x})}
    else{
    terms = terms[order(terms[,2]),]
      enrichresult = lapply(enrichresult, function(x){
        x = x[x[,1] %in% terms[,1],]
        x})}
  }
  if(padjust){
    term = lapply(enrichresult,function(x){
      x = x[order(x$p.adjust),]
      if(!is.na(top) ){
        x = x[1:top,]
        }
      else{
        x = x[x$p.adjust <= pvalue,]
      }
      return(x)})
    allterm = unique(unlist(lapply(term,"[","Term")))
    mapterm = lapply(enrichresult, function(x){
      x[x$Term %in% allterm,c("Term","p.adjust")]})
  }
  else{
    term = lapply(enrichresult,function(x){
      x = x[order(x$P.value),]
      if(!is.na(top)){
        x= x[1:top,]
      }
      else{
        x = x[x$P.value <= pvalue,]
      }
      return(x)})
    allterm = unique(unlist(lapply(term,"[","Term")))
    mapterm = lapply(enrichresult, function(x){
      x[x$Term %in% allterm,c("Term","P.value")]})
  }
  for(i in 1:length(term)){
    colnames(mapterm[[i]])[2] = names(term)[i]
  }
  tablebuild = function(data,filter,replace =1){
    table = Reduce(function(...) merge(..., by="Term", all=TRUE), data)
    rownames(table) = table$Term
    table = table[,-1]
    table[is.na(table)] = replace
    if(filter){
      check = apply(table, 1, function(x){all(x<= 0.05)})
      table = table[!check,]
      sumcount = rowSums(table[,2:4] < 0.05)
      table = table[sumcount < 3,]
      table = table[!grepl("UNKNOWN",rownames(table)),]
    }
    table = data.matrix(table)
    table
  }
  table = tablebuild(data= mapterm, filter = filter)
  log.table = abs(log(table,base=10))

  if(display_numbers){
    countterm = lapply(enrichresult, function(x){
      x[,c("Term","Count")]})
    for(i in 1:length(term)){
      colnames(countterm[[i]])[2] = names(term)[i]
    }
    counttable = tablebuild(data = countterm, filter =FALSE,replace =0)
    counttable = counttable[match(rownames(table),rownames(counttable)),]
  }
  if(is.null(annotation)){
    annotation = data.frame(Group = colnames(log.table))
    rownames(annotation) = colnames(log.table)
  }
  if(annotated_row){
    if(is.null(annotation_row)){
      anno.data = table(terms[,"Group"])
      len = length(anno.data)
      annotation_row = data.frame(Tissue = factor(rep(names(anno.data), anno.data)))
      rownames(annotation_row) = terms[,1]
    }
  }
  if(is.null(annotation_colors)){
    if(annotated_row){
      annotation_colors = list(Group = c("#abd9e9", "#2c7bb6", "#fee090", "#d7191c", "#fdae61"),
                               Tissue = c("#3B528BFF","#440154FF","#21908CFF","#FDE725FF","#5DC863FF"))
      names(annotation_colors[[1]]) = colnames(log.table)
      names(annotation_colors[[2]]) = names(anno.data)
    }
    else{
      annotation_colors = list(Group = c("#abd9e9", "#2c7bb6", "#fee090", "#d7191c", "#fdae61"))
      names(annotation_colors[[1]]) = colnames(log.table)
    }
  }
  hm.colors = colorpanel(100,"darkblue","white","darkred")

  if(display_numbers){
    display = counttable
  }
  else{
    display = display_numbers
  }
  log.table[log.table > 4] = 4
  pheatmap(log.table,color = hm.colors, annotation= annotation,
           annotation_colors =annotation_colors,main = main,
           display_numbers=display,annotation_legend=annotation_legend,
           cluster_cols =FALSE, annotation_col = annotation_col,
           annotation_row = annotation_row,...)
  return(table)
}








