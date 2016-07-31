#' Funtional Enrichment Analysis
#'
#' This function performs Gene Ontology (GO) and KEGG pathways functional enrichment for the five gene categories by calling
#' \pkg{clusterProfiler} package.
#' @usage functionEnrich(cate.gene, organism = "human", convert = TRUE, from = "SYMBOL", ont = "BP",
#' pAdjustMethod = "bonferroni", GO = TRUE, KEGG = FALSE, raw = FALSE)
#' @param cate.gene a list of the five gene categories, alternatively output by \code{\link{categorizeGene}}.
#' @param organism a character of organism "human" or "mouse" to indicate the species of background genes.
#' @param convert logical to determine whether the gene ID should be converted to "ENTREZID", default to \code{TRUE}.
#' @param from the gene id type of input data, see
#' @param ont One of "MF", "BP", and "CC" subontologies, see \code{\link[clusterProfiler]{enrichGO}}.
#' @param pAdjustMethod correction method for p-value, one of "holm", "hochberg", "hommel", "bonferroni", "BH",
#' "BY", "fdr", "none", default to "fdr", see details in \code{\link[stats]{p.adjust}}.
#' @param GO logical to determine whether the functional enrichment is performed on Gene Ontology, default to \code{TRUE}.
#' @param KEGG logical to determine whether the functional enrichment is performed on KEGG pathways, default to \code{FALSE}.
#' @param raw logical to determine if the raw "enrichResult" is output, default to \code{FALSE} to output a summary of the "enrichResult".
#' @return Function enrichment analysis results.
#' @export

functionEnrich = function(cate.gene,organism = "human", convert=TRUE,from="SYMBOL", ont = "BP",
                          pAdjustMethod ="bonferroni",GO=TRUE,KEGG= FALSE,raw = FALSE){

  if(organism == "human"){
    OrgDb = "org.Hs.eg.db"
  }
  else if(organism == "mouse"){
    OrgDb = "org.Mm.eg.db"
  }

  if(convert){
    ids.list = lapply(cate.gene,function(x){
      bitr(x, fromType=from, toType=c("ENTREZID"), OrgDb = OrgDb)})
  }
  else
    ids.list = cate.gene

  if(GO){
    gobp=lapply(ids.list,function(x){
      if(convert){
        x = x$ENTREZID
      }
      else{
        x=x
      }
      enrichGO(x, OrgDb = OrgDb,ont =ont,
               pAdjustMethod= pAdjustMethod,
               pvalueCutoff = 1,qvalueCutoff=1,readable = TRUE)
    })
    gobp.result =lapply(gobp,summary)
    gobp.raw = gobp
  }

  if(KEGG){
    org2org <- c( 'human'= 'hsa','mouse'='mmu' )
    if(organism %in% names(org2org)){
      org = org2org[organism]
    }
    kegg = lapply(ids.list,function(x){
      if(convert){
        x= x$ENTREZID
      }
      else{
        x = x
      }
      enrichKEGG(x, organism = org, pAdjustMethod=pAdjustMethod,
                 pvalueCutoff = 1,qvalueCutoff=1,use_internal_data =F)
    })
    kegg.result = lapply(kegg,summary)
    kegg.raw = kegg
  }

  if(GO && KEGG){
    if(raw){
      return(list(GO=gobp, KEGG = kegg))
    }
    return (list(GO = gobp.result,KEGG = kegg.result))
  }
  else if(GO) {
    if(raw){
      return(gobp)
    }
    else{
      return (gobp.result)
    }
  }
  else if (KEGG){
    if(raw){
      return(kegg)
    }
    else{
      return (kegg.result)
    }
  }
  else{
    stop ("No GO or KEGG enrichment is selected.")
  }
}


