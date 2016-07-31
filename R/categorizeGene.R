#' Gene Categorization
#'
#' This function categorizes differential genes of each pair-wise comparison among e.g. initiating A, derived B and primary C
#' samples during a cellular engineering, into five categories named from \emph{Reversed},\emph{Inactive},\emph{Insufficient},
#' \emph{Successful} and \emph{Over} representing the gene reprogrammed states, and calculates the ratio of expression difference (ED)
#'  between B and A to the ED between C and A.
#' @usage categorizeGene(expr, diffGene, from.sample,to.sample, target.sample)
#' @param expr a data frame with expression for all genes in \code{diffGene}, alternatively output from \code{\link{diffGene}}
#' function.
#' @param diffGene a list of differential genes in three comparisons, alternatively output by \code{\link{diffGene}} function.
#' @param from.sample,to.sample,target.sample character to specify the name of initiating sample, derived sample and primary
#' sample during a cellular engineering, must be consistent with sample names in the \code{expr} data frame.
#' @details Gene (g) categorization is achieved by considering the pair-wise comparisons (Expression Difference, ED eq. 1) among
#' the three types of samples and the ratio of such differences (EDg ratio, eq. 2).
#' \eqn{EDg(B, A) = average expression of g in B minus average expression of g in A}    (1)
#' \eqn{EDg ratio = EDg(B, A)/EDg(C, A)}                                            (2)
#'
#' \emph{Reversed}, EDg(B,A) and EDg(B,C) are significantly differential, while EDg(C,A) is not limited by differential or not,
#'  EDg ratio is smaller than 0;
#' \emph{Inactive}, EDg(C,A) and EDg(B,C) are significantly differential, while EDg(B,A) is not differential;
#' \emph{Insufficient}, EDg(B,A) ,EDg(C,A) and EDg(B,C) are all significantly differential, EDg ratio is between 0 and 1;
#' \emph{Successful}, EDg(B,A) and EDg(C,A) are significantly differential, while EDg(B,C) is not differential;
#' \emph{Over}, EDg(B,A) and EDg(B,C) are significantly differential, while EDg(C,A) is not limited by differential or not,
#' EDg ratio is greater than 1.
#' For the \emph{Inactive} and \emph{Successful} categories, genes with bottom and top 5 percentage ED ratios are
#' removed to avoid the ambiguous categorization with \emph{Reversed}, \emph{Insufficient} or \emph{Over} categories.
#' @return A list with components:
#' a list of the five gene categories
#' a list of the ED ratios for the five gene categories.
#' @export


categorizeGene = function(expr, diffGene, from.sample, to.sample, target.sample){

  samples = c(from.sample, to.sample, target.sample)
  expr = expr[,gsub("[^[:alpha:]]","",colnames(expr)) %in% samples]
  expr = expr[rowSums(expr) >0,]

  group= c(to.sample, target.sample)
  comp = c(paste(group,from.sample,sep="_"), paste(group,collapse = "_"))
  from.mean = rowMeans(expr[,gsub("[^[:alpha:]]","",colnames(expr)) %in% from.sample])
  for(i in 1:2){
    group.mean = rowMeans(expr[,gsub("[^[:alpha:]]","",colnames(expr)) %in% group[i]])
    expr[[comp[i]]] = group.mean- from.mean
  }
  # calculate the ratio of expression difference based two comparisons
  expr$ratio = expr[,comp[1]]/expr[,comp[2]]
  diff.ratio = expr$ratio

  # considering the differential genes between to.sample and target.sample
  if(!all(names(diffGene) %in% comp) ){
    stop ("Please correct the names of 'diffGene' list.")
  }

  #get the five gene categories
  diff.to = diffGene[[comp[1]]]
  diff.target = diffGene[[comp[2]]]
  diff.totarget = diffGene[[comp[3]]]

  suc.group = expr[rownames(expr) %in% diff.target & rownames(expr) %in% diff.to & !rownames(expr) %in% diff.totarget,]
  inact.group = expr[rownames(expr) %in% diff.target & !rownames(expr) %in% diff.to & rownames(expr) %in% diff.totarget,]
  insuff.group = expr[rownames(expr) %in% diff.target & rownames(expr) %in% diff.to & rownames(expr) %in% diff.totarget,  ]
  insuff.group = insuff.group[insuff.group$ratio <1 & insuff.group$ratio >0,]

  amb.group = expr[rownames(expr) %in% diff.totarget & rownames(expr) %in% diff.to,]
  rev.group =  amb.group[amb.group$ratio<0, ]
  over.group = amb.group[amb.group$ratio >1,]

  fivelist = list(Reverse = rev.group,Inactive = inact.group, Insufficient = insuff.group,
                  Successful= suc.group,  Over = over.group)
  vali.gene = lapply(fivelist, function(x){rownames(x)})
  vali.gene = Reduce(union, vali.gene)
  diff.ratio = expr[rownames(expr) %in% vali.gene, ]
  ratio = diff.ratio$ratio
  ratio = data.frame(ratio)
  colnames(ratio) = "ED_ratio"
  rownames(ratio) = rownames(diff.ratio)

  fivelist = lapply(fivelist, function(x){rownames(x)})
  ratio.list = list()
  for(i in 1:length(fivelist)){
    ratio.list[[i]] = ratio[rownames(ratio) %in% fivelist[[i]],,drop=FALSE]
  }
  names(ratio.list) = names(fivelist)

  #remove the bottom quantile 5% and top quantile 5% genes based on the ED ratio for the Successful and Inactive genes
  ratio.list.filter = lapply(ratio.list[c(2,4)], function(x){
    y = x[,1]
    z= x[x[,1] >= quantile(y,0.05) & x[,1] <= quantile(y,0.95),,drop=F]
    z
  })
  ratio.listold = ratio.list
  ratio.list[2] = ratio.list.filter[1]
  ratio.list[4] = ratio.list.filter[2]
  fivelist = lapply(ratio.list, rownames)
  return(list(fivelist = fivelist, ratiolist = ratio.list, ratio.listold = ratio.listold))
}


