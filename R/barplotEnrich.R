#' Barplot the enrichResult
#'
#' This function is revised from the \code{\link[DOSE]{barplot.enrichResult}} function DOSE package, and used to
#' perform a barplot of the enrichResult object.
#' @usage barplotEnrich(height, x = "Count", colorBy = "p.adjust", top = 5,
#' font.size = 12, title = "", color = NULL,...)
#' @param height enrichResult object, alternatively output from \code{\link{functionEnrich}} function.
#' @param x one of "Count" and "GeneRatio" to specify the x axis.
#' @param colorBy one of 'pvalue', 'p.adjust', 'qvalue'.
#' @param top number of top categories to show.
#' @param font.size,title font size and title.
#' @param color the color of the bar.
#' @param ... see parameters in \code{\link[ggplot2]{fortify}} function.
#' @export

barplotEnrich = function (height, x = "Count", colorBy = "p.adjust", top = 5,
          font.size = 12, title = "", color = NULL,...)
{
  object <- height
  colorBy <- match.arg(colorBy, c("pvalue", "p.adjust", "qvalue"))
  if (x == "geneRatio" || x == "GeneRatio") {
    x <- "GeneRatio"
  }
  else if (x == "count" || x == "Count") {
    x <- "Count"
  }
  Description <- Count <- NULL
  df <- fortify(object, top = top, by = x,
                ...)
  p <- ggplot(df, aes_string(x = "Description", y = x))
  p <- p + geom_bar(stat = "identity") + coord_flip() + theme_dose(font.size)
  if ("pvalue" %in% colnames(p$data)) {
    pvalue <- NULL
    p <- p + aes_string(fill = colorBy) + scale_fill_continuous(low = color[1],
                                                                high = "grey" )
  }
  else {
    p <- p + aes(fill = Description) + theme(legend.position = "none")
  }
  p <- p + ggtitle(title) + xlab("") + ylab("")
  return(p)
}

