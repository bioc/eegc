#' eegc
#' @name eegc
#' @docType package
#' @import igraph org.Hs.eg.db org.Mm.eg.db DOSE
#' @importFrom R.utils capitalize
#' @importFrom S4Vectors DataFrame
#' @importFrom gplots colorpanel
#' @importFrom sna stresscent
#' @importFrom wordcloud textplot
#' @importFrom ggplot2 ggplot fortify geom_bar aes ggtitle xlab
#' ylab theme scale_fill_continuous coord_flip aes_string
#' @importFrom limma makeContrasts contrasts.fit lmFit topTable eBayes voom
#' @importFrom pheatmap pheatmap
#' @importFrom edgeR DGEList calcNormFactors
#' @importFrom DESeq2 DESeqDataSetFromMatrix DESeq
#' @importFrom clusterProfiler bitr enrichKEGG enrichGO plotGOgraph plot
NULL

#' RNA-seq data in FPKM
#' A data frame containing 16692 gene in 7 samples from paper published by Sandler et.al. in 2014.,
#' samples are named from somatic dermal microvascular endothelial cell (DMEC) to induced
#' multipotent haematopoietic progenitor (rEChMPP) cells and primary cord blood (CB) cells.
#'
#' @docType data
#' @keywords datasets
#' @name SandlerFPKM
#' @usage data(SandlerFPKM)
#' @format A data frame with 16692 rows and 7 columns specifying for genes and samples.
#' @return A data frame with gene expression data.
NULL

#' Marker genes
#'
#' A vector containing 65 genes specific in endothelial and haematopoietic cells as
#'  listed in Sandler's paper.
#'
#' @docType data
#' @keywords datasets
#' @name markers
#' @usage data(markers)
#' @format A vector
#' @return A vector of marker genes
NULL

#' Differential genes for three comparisons
#'
#' A list of three vectors with significantly differential genes among three comparisons.
#'
#' @docType data
#' @keywords datasets
#' @name diffgene.genes
#' @usage data(diffgene.genes)
#' @format A list
#' @return A list with differential genes in three comparisons
NULL

#' Filtered expression data
#'
#' A data frame with filtered RNASeq FPKM data.
#'
#' @docType data
#' @keywords datasets
#' @name expr.filter
#' @usage data(expr.filter)
#' @format A data frame
#' @return A data frame with filtered gene expression
NULL

#' Categorized genes
#'
#' A list with five gene categories.
#'
#' @docType data
#' @keywords datasets
#' @name cate.gene
#' @usage data(cate.gene)
#' @format A list
#' @return A list with five gene catogories.
NULL

#' Expression Difference ratios of categorized genes
#'
#' A list with five data frames of gene ED ratios
#'
#' @docType data
#' @keywords datasets
#' @name cate.ratio
#' @usage data(cate.ratio)
#' @format A list
#' @return A list with ED ratios for five gene catogories.
NULL

#' Cell/Tissue specific gene sets
#'
#' A list of 126 cells/tissues-specific gene sets identified by SpeCond algorithm
#'  from 126 cells/tissues in Gene Enrichment Profiler database.
#'
#' @docType data
#' @keywords datasets
#' @name tissueGenes
#' @usage data(tissueGenes)
#' @format A list
#' @return A list with 126 cells/tissues-specific gene sets
NULL

#' Tissue mapping to groups
#'
#' A data frame with 126 cells/tissues and 30 C/T groups they belong to.
#'
#' @docType data
#' @keywords datasets
#' @name tissueGroup
#' @usage data(tissueGroup)
#' @format A data frame with 126 rows (tissues) and 3 columns (Tissue, Tissue_abbr, Group)
#' @return A data frame with 126 cells/tissues and 30 C/T groups they mapped to
NULL

#' Gene regulatory network based human cell/tissue-specific gene sets
#'
#' A list of 16 human cells/tissues-specific gene sets summarized from the gene regulatory network
#' downloaded from the CellNet website.
#'
#' @docType data
#' @keywords datasets
#' @name human.gene
#' @usage data(human.gene)
#' @format A list
#' @return A list with 16 human cells/tissues-specific gene sets from CellNet.
NULL

#' Gene regulatory network based human cell/tissue-specific transcription factor (TF)
#' regulated gene sets
#'
#' A list of 1455 human cells/tissues-specific TF regulated gene sets summarized from the gene
#' regulatory network downloaded from the CellNet website.
#'
#' @docType data
#' @keywords datasets
#' @name human.tf
#' @usage data(human.tf)
#' @format A list
#' @return A list with 1455 human cells/tissues-specific TF regulated gene sets
NULL

#' Human cell/tissue-specific gene-gene regulation
#'
#' A list of 16 data frames (cells/tissues) with transcription factors (TF) to target genes (TG)
#' regulation information.
#'
#' @docType data
#' @keywords datasets
#' @name human.grn
#' @usage data(human.grn)
#' @format A list with 16 data frames
#' @return A list of human gene regulatory information.
NULL

#' Gene regulatory network based mouse cell/tissue-specific gene sets
#'
#' A list of 20 mouse cells/tissues-specific gene sets summarized from the gene regulatory network
#' downloaded from the CellNet website.
#'
#' @docType data
#' @keywords datasets
#' @name mouse.gene
#' @usage data(mouse.gene)
#' @format A list
#' @return A list of 20 mouse cells/tissues-specific gene sets
NULL

#' Gene regulatory network based mouse cell/tissue-specific transcription factor (TF) regulated
#' gene sets
#'
#' A list of 1744 mouse cells/tissues-specific TF regulated gene sets summarized from the gene
#' regulatory network downloaded from the CellNet website.
#'
#' @docType data
#' @keywords datasets
#' @name mouse.tf
#' @usage data(mouse.tf)
#' @format A list
#' @return A list of 1744 mouse cells/tissues-specific TF regulated gene sets.
NULL

#' Mouse cell/tissue-specific gene-gene regulation
#'
#' A list of 20 data frames (cells/tissues) with transcription factors (TF) to target genes (TG)
#' regulation information.
#'
#' @docType data
#' @keywords datasets
#' @name mouse.grn
#' @usage data(mouse.grn)
#' @format A list with 20 data frames
#' @return A list of mouse gene regulatory information.
NULL

#' GO enrichment results
#'
#' A list of 5 data frames with the Gene ontology enrichment results for the
#' five gene categories.
#'
#' @docType data
#' @keywords datasets
#' @name goenrich
#' @usage data(goenrich)
#' @format A list with 5 data frames
#' @return A list of GO enrichment results
NULL


