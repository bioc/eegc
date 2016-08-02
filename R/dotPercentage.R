#' Percentage Calculation and Visualization
#'
#' This function calculate the percentage of genes in each category over given annotated gene sets and plot
#' the percentages.
#' @usage dotPercentage(cate.gene, annotated.gene, order.by = NULL, type = "l", lty = 1,
#' pch = NULL, col = NULL, srt = 50, font = 1, adj = c(1,1), cex = 1, add.line = TRUE,
#' legend = TRUE, legend.label = NULL, ...)
#' @param cate.gene a list of the five gene categories, alternatively output by \code{\link{categorizeGene}}.
#' @param annotated.gene a list of the annotated gene sets which the \code{cate.gene} are proportioned in.
#' @param order.by one character out of of "Reversed","Inactive","Insufficient","Successful" and "Over" to
#' specify a gene category the percentage is ordered by.
#' @param type,lty,pch,col parameters for the plotting, specifying the type of plotting; the line type when
#' \code{type} is "l"; the symbol of points on the line; and the color of lines, see graphic parameters
#' in \code{\link[graphics]{par}}().
#' @param srt,font,cex,adj parameters for the text labeled on x-axis, specifying the string rotation in degrees;
#' the font of text; the text size, see graphic parameters in \code{\link[graphics]{par}}().
#' @param add.line logical to determine if to add lines on the dots, logical to \code{TRUE}.
#' @param legend logical to determine whether the legend is added on the figure, default to \code{TRUE}.
#' @param legend.label labels of the legend, applied only when \code{legend} is \code{TRUE}.
#' @param ... other parameters see \code{\link[graphics]{plot}}.
#' @return a data frame with the percentage of \code{cate.gene} in the \code{annotated.gene}.
#' @importFrom R.utils capitalize
#' @export
#' @examples
#' # load the C/T-specific genes in 16 cells/tissues
#' data(human.gene)
#' data(cate.gene)
#' # perc = dotPercentage(cate.gene = cate.gene, annotated.gene = human.gene,
#' #                     order.by = "Successful")

dotPercentage = function(cate.gene, annotated.gene, order.by =NULL, type = "l", lty = 1,
                         pch= NULL, col = NULL, srt =50, font = 1, adj = c(1,1), cex = 1,
                         add.line = TRUE, legend = TRUE, legend.label = NULL, ...){
    net.len = length(annotated.gene)
    inter = function(x,y){
        int = intersect(x,y)
        perc = length(int)/length(y)*100
        perc
    }
    data = matrix(NA,length(cate.gene),net.len)

    gene.len = length(cate.gene)
    for(i in 1:gene.len){
        data[i,] = unlist(lapply(annotated.gene,function(x){inter(cate.gene[[i]],x)}))
    }
    rownames(data) = names(cate.gene)
    colnames(data) = names(annotated.gene)

    if(!is.null(order.by)){
      data= t(data)
      order.by  = match.arg(order.by, c("Reversed","Inactive","Insufficient","Successful","Over"))
    }
    else{
      order.by = "Inactive"
    }
      data = data[order(data[,order.by],decreasing =TRUE),]
      data = t(data)
      cate.name = rownames(data)
      annotated.name = colnames(data)

    plot(data[1,], ylim = c(0,ceiling(max(data))), xaxt ="n",yaxt = "n",las = 0,
         xlab="", ylab = "Gene percentage(%)",type ="n",...)

    if(is.null(pch)){
      pch = c(15:19)
    }

    if(is.null(col)){
    	col = c("#abd9e9", "#2c7bb6", "#fee090", "#d7191c", "#fdae61")
    }
    if(is.null(type)){
      type = "b"
    }
    if(add.line){
      lty = 1
    }
    else{
      lty = 0
    }
    for(i in 1:gene.len){
      lines(data[i,],type= "b",lwd = 1.5,lty =lty,col = col[i],pch = pch[i])
    }

    axis(1, at = seq(1,net.len, 1), labels =rep("",net.len))
    text(seq(1,net.len, 1), par("usr")[3] - 0.4, labels = annotated.name,
            srt =srt, font = font, adj =adj, cex = cex, xpd = TRUE)

    axis(2, at = seq(0,ceiling(max(data)),1), labels = seq(0,ceiling(max(data)),1),
         cex.axis = 1)

    if(legend){
      if(is.null(legend.label)){
        legend.label = cate.name
      }
      legend("topright", legend = legend.label, bty = "n",
             col = col, pch = pch,lty= lty, cex = cex)
    }
    return(data)
}
