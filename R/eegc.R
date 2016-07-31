#' eegc
#' @name eegc
#' @docType package
#' @import igraph org.Hs.eg.db org.Mm.eg.db
#' @importFrom R.utils capitalize
#' @importFrom S4Vectors DataFrame
#' @importFrom gplots colorpanel
#' @importFrom sna stresscent
#' @importFrom wordcloud textplot
#' @importFrom ggplot2 ggplot fortify geom_bar aes ggtitle xlab
#' ylab theme scale_fill_continuous coord_flip aes_string
#' @importFrom DOSE theme_dose
#' @importFrom limma makeContrasts contrasts.fit lmFit topTable eBayes voom
#' @importFrom pheatmap pheatmap
#' @importFrom edgeR DGEList calcNormFactors
#' @importFrom DESeq2 DESeqDataSetFromMatrix DESeq
#' @importFrom clusterProfiler bitr enrichKEGG enrichGO plotGOgraph
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
NULL


