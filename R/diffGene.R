#' Gene Expression Differential Analysis
#'
#' Identify the differentially expressed genes for each pair-wise comparison of given three types of samples.
#' @usage diffGene(expr, array = TRUE, fpkm = FALSE, counts =FALSE, method =c("limma","DESeq2"),
#' from.sample, to.sample, target.sample, filter = FALSE, filter.perc = 0.4,
#' padjust ="fdr", signif = TRUE, pvalue = 0.05)
#' @param expr a data frame with gene expression data.
#' @param array,fpkm,counts logical, specifying the type of input gene expression data.
#' @param method differential analysis method, alternatively to "limma" and "DESeq2", default to "limma".
#' "DESeq2" can be chosen only when \code{counts} is \code{TRUE}.
#' @param from.sample,to.sample,target.sample character to specify the name of initiating sample, derived sample and primary
#' sample during a cellular engineering.
#' @param filter logical to indicate whether the genes need to be filtered when match the parameter \code{filter.perc}
#' , only applied to fpkm and counts data.
#' @param filter.perc a 0 to 1 number to specify the gene filter criteria by the percentage of samples with non-zero expression.
#' Only used to fpkm and counts data when \code{filter} is \code{TRUE}, and filter the genes with non-zero expression in less than
#' filter.perc samples.
#' @param padjust indicate the method to do p.value correction, default to "fdr". See \code{\link[stats]{p.adjust}}.
#' @param signif logical to indicate whether only the significantly differential genes are output, default to FALSE.
#' @param pvalue a cutoff p.value for the significant genes, default to 0.05, only used when \code{signif} is TRUE.
#' @details This function can be applied on both microarray and RNA-seq data for differential analysis when one of the "array",
#' "fpkm", or "counts" is specified. It does differential analysis to each pair-wise sample comparison among the from.sample,
#' to.sample and target.sample.
#' @return A list with components :
#' a list with differential analysis result for each pair-wise comparison;
#' a list with differential gene names for each pair-wise comparison;
#' a data frame with filtered/unfiltered gene expression.
#' @export
#' @examples
#' data(SandlerFPKM)
#'
#' # differential expression analysis:
#' diffgene = diffGene(expr = SandlerFPKM, array=FALSE,  fpkm=TRUE,  counts=FALSE,
#'                    from.sample="DMEC", to.sample="rEChMPP", target.sample="CB",
#'                    filter=TRUE, filter.perc =0.4, pvalue = 0.05 )
#'
#' # differential analysis results
#' diffgene.result = diffgene[[1]]
#' # differential genes
#' diffgene.genes = diffgene[[2]]
#' # filtered expression data
#' expr.filter = diffgene[[3]]


diffGene = function(expr, array= TRUE, fpkm = FALSE, counts = FALSE, method = c("limma","DESeq2"),
                    from.sample, to.sample, target.sample, filter = FALSE, filter.perc = 0.4, padjust = "fdr",
                    signif = TRUE, pvalue = 0.05){
  if(!is.data.frame(expr)){
    expr = data.frame(expr, stringsAsFactors = FALSE)
  }
  #filter the genes with low expression
  n.sample = ncol(expr)
  if(!array){
    if(filter){
      expr = expr[rowSums(expr >= 1) >= n.sample*filter.perc,]
    }
  }
  expr.all = expr

  expr.comb = list()
  expr.comb[[1]] = expr[,gsub("[^[:alpha:]]","",colnames(expr)) %in% from.sample]
  expr.comb[[2]] = expr[,gsub("[^[:alpha:]]","",colnames(expr)) %in% to.sample]
  expr.comb[[3]] = expr[,gsub("[^[:alpha:]]","",colnames(expr)) %in% target.sample]
  n = lapply(expr.comb, ncol)

  expr.comb2 = list()
  n.comb = list()
  expr.comb2[[1]] = cbind(expr.comb[[1]],expr.comb[[2]])
  n.comb[[1]] = c(n[[1]],n[[2]])
  expr.comb2[[2]] = cbind(expr.comb[[1]],expr.comb[[3]])
  n.comb[[2]] = c(n[[1]],n[[3]])
  expr.comb2[[3]] = cbind(expr.comb[[3]],expr.comb[[2]])
  n.comb[[3]] = c(n[[3]],n[[2]])

  #differential analysis
  method = match.arg(method)
  diff.result = list()
  for(i in 1:3){
    expr = expr.comb2[[i]]
    n = n.comb[[i]]
    if(method =="limma"){
      if(array){
        expr = log(expr+1,base =2)
      }
      if(fpkm){
        #log transform expr+2 to avoid the negative logFPKM
        expr =log(expr+2,base=2)
      }
      if(counts){
        expr.count = DGEList(counts = expr,genes = rownames(expr))
        norm.expr = calcNormFactors(expr.count)
        v = voom(norm.expr, design, plot=FALSE)
        expr = v
      }

      comp = factor(c(rep("control",n[1]),rep("treat",n[2])))
      design = model.matrix(~0+comp)
      colnames(design)=c("control","treat")
      fit <- lmFit(expr, design)
      contrast.matrix = makeContrasts(treat-control,levels=design)
      fit2= contrasts.fit(fit,contrast.matrix)
      fit2 <- eBayes(fit2,trend=TRUE)
      if(signif){
        tablesign = topTable(fit2,p.value=pvalue,coef=1,adjust.method=padjust,number=Inf)
        sign = merge(tablesign,expr,by="row.names")
        sign = sign[sign$logFC>=1 | sign$logFC <=-1,]
        sign = sign[order(sign$P.Value),]
        diff.result[[i]] = sign
      }
      else{
        table = topTable(fit2,coef=1,adjust.method=padjust,number=Inf)
        all = merge(table,expr,by="row.names")
        all = all[order(all$P.Value),]
        diff.result[[i]] = all
      }
    }
    if(method == "DESeq2"){
      sample1 = "control"
      sample2 = "treat"
      comp = factor(c(rep(sample1,n[1]),rep(sample2,n[2])))
      coldata <- DataFrame(comp)
      exprFilt = expr[rowSums(expr) > 10,]
      dds <- DESeqDataSetFromMatrix(countData=exprFilt, colData=coldata, design=~comp)
      dds <- DESeq(dds)
      res = results(dds)
      res = as.data.frame(res)
      res[,`padjust`] = p.adjust(res$pval, padjust)
      if(sign){
        reSign = res[res[,`padjust`] < pvalue,]
        reSign = reSign[reSign$log2FoldChange >= 1 | reSign$log2FoldChange <= -1,]
        diff.result[[i]] = all
      }
      else{
        diff.result[[i]] = all
      }
    }
  }
  names(diff.result) = c(paste(to.sample,from.sample,sep="_"),paste(target.sample, from.sample, sep="_"),
                         paste(to.sample, target.sample,sep="_"))

  diffgene = lapply(diff.result, function(x) x$Row.names)
  return(list(difflist = diff.result, diffgene = diffgene, expr = expr.all))
}



