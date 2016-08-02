#' Quantify Genes and Corresponding ED ratios in Each Category
#'
#' Quantify genes in each gene category and their expression difference (ED) ratios in a density plot.
#' @usage densityPlot(ratio, color = NULL, main = NA, xlab = NA, ylab = "Density",
#' legend.labels = NULL, shade = TRUE, transparency = TRUE, proportion = TRUE,
#' out.file = NULL, ...)
#' @param ratio a list of ED ratios for five gene categories, alternatively output by \code{\link{categorizeGene}}.
#' @param color vector of colors for the five gene categories.
#' @param main,xlab,ylab the overall title, tile for x axis, and title for y axis, see \code{\link[graphics]{title}}.
#' @param legend.labels vector of labels for the legend.
#' @param shade logical to determine if the five categories are filled with shades.
#' @param transparency logical to determine if the density plot is transparent.
#' @param proportion logical to determine whether the proportion of each category genes over the all genes is drawn on the
#'  density plot.
#' @param out.file a character string naming the output file with density plot.
#' @param ... parameters in \code{\link[graphics]{plot}}.
#' @return a density plot
#' @export
#' @examples
#' data(cate.ratio)
#' names(cate.ratio)
#' # make the extreme ED ratios in Reversed and Over categories to the median values
#' reverse = cate.ratio[[1]]
#' over = cate.ratio[[5]]
#' reverse[reverse[,1] <= median(reverse[,1]), 1]  = median(reverse[,1])
#' over[over[,1] >= median(over[,1]),1] = median(over[,1])
#' cate.ratio[[1]] = reverse
#' cate.ratio[[5]] = over
#'
#' # densityPlot(cate.ratio, xlab = "ED ratio", ylab = "Density", proportion = TRUE)



densityPlot = function(ratio, color = NULL, main = NA, xlab = NA,ylab = "Density",
                       legend.labels = NULL, shade = TRUE,transparency = TRUE,
                       proportion = TRUE, out.file = NULL,...){

  n = length(ratio)
  ratio.all = do.call("rbind", ratio)
  den = list()
  for(i in 1:n){
    data = list(ratio[[i]][,1])
    den[[i]] = lapply(data,density,n =512)
  }

  density.y = unlist(lapply(den, function(x){
    max(sapply(x, "[[","y"))
  }))
  density.max = max(density.y)

  x.min = min(unlist(lapply(den, function(x){
    min(sapply(x, "[[","x"))
  })))
  x.max = max(unlist(lapply(den, function(x){
    max(sapply(x, "[[","x"))
  })))
  xlim = c(x.min, x.max)

  if(is.null(color)){
    color = c("#abd9e9", "#2c7bb6", "#fee090", "#d7191c", "#fdae61")
  }
  plot(den[[1]][[1]],col = color[1],lwd=2,font.lab=2, xlim = xlim,
       ylim=c(0,density.max+ 0.2*density.max),xaxt="n",
       xlab =xlab,ylab=ylab,main =main,...)
  for(i in 2:n){
    lines(den[[i]][[1]],col = color[i],lwd=2)
  }
  axis(side=1,at = seq(floor(xlim[1]), ceiling(xlim[2]), by =1))

  if(shade){
    if(transparency & !is.null(out.file)){
      color = adjustcolor(color, alpha.f = 0.8)
    }
    for(i in 1:n){
      shade.at = c(min(den[[i]][[1]]$x),max(den[[i]][[1]]$x))
      shade.y = sapply(shade.at,function(z){which.min(abs(den[[i]][[1]]$x-z))})

      with(den[[i]][[1]], polygon(x=c(x[c(shade.y[1],shade.y[1]:shade.y[2],shade.y[2])]),
                                  y= c(0,y[shade.y[1]:shade.y[2]],0), col= color[i],border = NA))
    }
    if(is.null(legend.labels))
      legend.labels= names(ratio)
    legend("topleft",legend=legend.labels,fill= color,border = color,
           density=rep(NA,5),bty= "n", x.intersp=0.3,horiz=TRUE,...)
  }

  #add proportion of each gene category over all genes
  if(proportion){
    gene = nrow(ratio.all)
    gene.pro = lapply(ratio,function(x){
      gene.numb = nrow(x)
      pro = gene.numb/gene
      pro
    })
    gene.pro = unlist(gene.pro)
    labels = paste(round(100*gene.pro,2),"%",sep="")

    pro.label = list()
    pro.label[1] = lapply(den[1], function(x){
      quantile(sapply(x, "[[","x"), 0.25)
    })
    pro.label[c(2,3,4)] = lapply(den[c(2,3,4)], function(x){
      mean(sapply(x, "[[","x"))
    })
    pro.label[5] = lapply(den[5], function(x){
      quantile(sapply(x, "[[","x"), 0.75)
    })

    y= density.y+0.05*density.max
    delta.y = abs(y[3] - y[4])
    if(delta.y < 0.05*density.max){
      y[4] = y[4]+0.05*density.max
    }
  text(x = pro.label, y = y ,labels, adj=0.5,col="black",font=1,...)
  }
  if(!is.null(out.file))
    dev.copy2pdf(file = out.file)
}



