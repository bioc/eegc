#' Enrichment Analysis by Hypergeometric Distribution Test
#'
#' @usage enrichment(cate.gene, annotated.gene, background.gene, padjust.method = "fdr" )
#' @param cate.gene a list of the five gene categories, alternatively output by \code{\link{categorizeGene}}.
#' @param annotated.gene a list of annotated gene sets which the \code{cate.gene} are enriched in.
#' @param background.gene vector of background genes, e.g. all genes screened by microarray or RNA-sequencing.
#' @param padjust.method correction method for enrichment p-values, one of "holm", "hochberg", "hommel", "bonferroni",
#' "BH", "BY", "fdr", "none", default to "fdr", see details in \code{\link[stats]{p.adjust}}.
#' @return A list of enrichment results for the five gene categories.
#' @export
#' @examples
#' # load the cell/tissue-specific genes
#' data(tissueGenes)
#' # load the mapping file of cells/tissues to grouped cells/tissues
#' data(tissueGroup)
#'
#' # get the background genes
#' data(expr.filter)
#' genes = rownames(expr.filter)
#' # enrichment analysis for the five gene categories
#' data(cate.gene)
#' tissueenrich = enrichment(cate.gene = cate.gene, annotated.gene = tissueGenes,
#'                          background.gene = genes, padjust.method = "fdr")
#' # select a group of cells/tissues
#' tissueGroup.selec = c("stem cells","B cells","T cells","Myeloid","Endothelial CD105+")
#' tissues.selec = tissueGroup[tissueGroup[,"Group"] %in% tissueGroup.selec,c(2,3)]
#' # tissuetable = heatmapPlot(tissueenrich, terms = tissues.selec, GO=FALSE,
#' #                           annotated_row = TRUE,annotation_legend = TRUE,
#' #                           main = "Tissue-specific enrichment")

enrichment = function(cate.gene,annotated.gene,background.gene, padjust.method = "fdr"){
  if(!is.list(cate.gene)){
    stop("'cate.gene' is not a list.")
  }

  # define the gene count function
  gene.count = function(gene, annotated.gene,index){
    if(is.data.frame(gene)){
      if(!is.null(rownames(gene))){
        genes = rownames(gene)
      }
      else{
        genes = gene[,1]
      }
    }
    else{
      genes =gene
    }
    gene.num = length(genes)
    
    enrich.num = length(annotated.gene)
    enrich.names = names(annotated.gene)
    enrich.new = matrix(NA,enrich.num+1,2)
    rownames(enrich.new) = c(enrich.names,"Others")

    for(i in seq_len(enrich.num)){
      tis = annotated.gene[[i]]
      enrich.genes = tis[tis %in% genes]
      enrich.new[i,1] = paste(enrich.genes,collapse=";")
      enrich.new[i,2] = length(enrich.genes)
    }
    enrich.new[enrich.num+1,2] = length(genes[!genes %in% unlist(annotated.gene)])
    colnames = paste(index,c("Genes","Count"),sep="")
    colnames(enrich.new) = colnames
    return(list(count = data.frame(enrich.new),number = gene.num))
  }

  #count the genes included in each annotated.gene
  enrichall = gene.count(background.gene,annotated.gene,index="All")
  enrichall.count = enrichall[["count"]]
  enrichall.num = enrichall[["number"]]

  name = names(cate.gene)
  if(is.null(name))
    stop("'cate.gene' shoule be named.")

  enrichgene.count = list()
  enrichgene.num = list()
  for(i in seq_along(cate.gene)){
    enrichgene = gene.count(cate.gene[[i]],annotated.gene,index=name[i])
    enrichgene.count[[i]] = enrichgene[["count"]]
    enrichgene.num[[i]] = enrichgene[["number"]]
  }
  names(enrichgene.count) = names(enrichgene.num) = name

  enrichcount = list()
  enrichcount = lapply(enrichgene.count,function(x){
    merge(x,enrichall.count,by="row.names")}
  )
  names(enrichcount) = name

  # define the function of enrichment analysis
  enrichment = function(data,gene.num,diff.num,index,method){
    factor = sapply(data, is.factor)
    data[factor] = lapply(data[factor], as.character)

    cols = paste(index,"Count",sep="")
    data$P.value = apply(data[,cols], 1, function(x)
      phyper(as.numeric(x[cols[2]])-1,as.numeric(x[cols[1]]),
             gene.num -as.numeric(x[cols[1]]),diff.num,lower.tail=FALSE))
    data$p.adjust = p.adjust(data$P.value, method = method)
    data = data[with(data, order(P.value)), ]
    return(data)
  }

  # perform the enrichment analysis
  enrich.result = list()
  for(i in seq_along(enrichcount)){
    enrich.result[[i]] = enrichment(enrichcount[[i]],gene.num=enrichall.num, diff.num = enrichgene.num[[i]],
                                    index=c("All",name[i]),method = padjust.method)
  }
  names(enrich.result) = name
  return(enrich.result)
}





