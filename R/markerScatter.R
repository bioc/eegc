#' Scatter Plot for Gene Expression
#'
#' Generates an expression profile of each gene catetory in one sample against another, alternatively plot the regression
#' line from linear modeling fitting.
#' @usage markerScatter(expr, log = FALSE, samples, cate.gene, markers, pch = 19, cex = 0.5,
#' col = NULL, xlab = NULL, ylab = NULL, main = NULL, add.line = TRUE, legend.labels = NULL,
#' ...  )
#' @param expr a data frame with gene expression.
#' @param log logical to determine if the gene expression data is log converted (add a small constant 2), default to FALSE.
#' @param samples a vector of samples to compare on the x axis and y axis.
#' @param cate.gene a list of the gene categories, alternatively output by \code{\link{categorizeGene}}.
#' @param markers vector of marker genes to be highlighted in the plot. No gene is highlighted when it's \code{NULL}.
#' @param pch,cex,col,xlab,ylab,main plot parameters, see details in \code{\link[graphics]{par}}.
#' @param add.line logical to determine if the linear model fitting line is added on the figure.
#' @param legend.labels vector of labels for the marker legend.
#' @param ... other parameters in \code{\link[graphics]{plot}}.
#' @importFrom wordcloud textplot
#' @details Visualization of gene expression in the five categories under each pair-wised comparison.
#' @export

markerScatter = function(expr, log =FALSE, samples, cate.gene, markers = NULL, pch = 19,
                         cex = 0.5, col = NULL, xlab = NULL, ylab = NULL, main = NULL,
                         add.line =TRUE, legend.labels = NULL,...){
  if(log){
    expr = log(expr+2,base = 2)
  }
  if(!is.list(cate.gene)){
    stop("'cate.gene' should be a list.")
  }

  expr = expr[,grepl(samples[1],colnames(expr)) | grepl(samples[2],colnames(expr))]
  expr[[samples[1]]] = apply(expr[,grepl(samples[1],colnames(expr))],1,mean)
  expr[[samples[2]]] = apply(expr[,grepl(samples[2],colnames(expr))],1,mean)

  expr.selec = list()
  n = length(cate.gene)
  for(i in 1:n){
    expr.selec[[i]] = expr[rownames(expr) %in% cate.gene[[i]],samples]
  }

  if(is.null(xlab)){
    if(log){
      xlab= bquote(paste('log'[2], " expression in ", .(samples[1]), sep=""))
    }
    else{
      xlab = paste("Expression in ",samples[1],sep="")
    }
  }
  if(is.null(ylab)){
    if(log){
      ylab = bquote(paste('log'[2], " expression in ", .(samples[2]), sep=""))
    }
    else{
      ylab = paste("Expression in ",samples[2],sep="")
    }
  }
  xlim.max = max(sapply(expr.selec, function(x){max(x[,samples[1]])}))
  ylim.max = max(sapply(expr.selec, function(x){max(x[,samples[2]])}))

  if(is.null(main)){
    main = "Scatter plot of gene category"
  }

  if(is.null(col)){
    col = c("#abd9e9", "#2c7bb6", "#fee090", "#d7191c", "#fdae61")
  }

  plot(expr.selec[[1]], pch = pch, cex =cex, col=col[1],
       xlim= c(0,ceiling(xlim.max)), ylim=c(0,ceiling(ylim.max)),
       xlab = xlab, ylab = ylab, main= main, font.lab =2,
       font.main =2,...)

  if(add.line){
    line1 = lm(expr.selec[[1]][,2] ~ expr.selec[[1]][,1])
    abline(line1, col = col[1], lty = 2, lwd = 2)

    for(i in 2:n){
      points(expr.selec[[i]],col = col[i], pch =pch, cex = cex)
      line = lm(expr.selec[[i]][,2] ~ expr.selec[[i]][,1])
      abline(line, col = col[i], lty = 2, lwd = 2)
    }
  }

  if(!is.null(markers)){
    expr.selec.all = do.call("rbind",expr.selec)
    expr.marker = expr.selec.all[rownames(expr.selec.all) %in% markers,]
    if(nrow(expr.marker) != 0){
      textplot(expr.marker[,1],expr.marker[,2],
               rownames(expr.marker), new=FALSE)
    }
  }
  if(is.null(legend.labels)){
    legend.labels = names(cate.gene)
  }
  legend("topleft",legend = legend.labels, col = col, pch =pch,
         cex= 1 , bty ="n")
  #return(expr.selec)
}



