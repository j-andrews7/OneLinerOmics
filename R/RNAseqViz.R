#' Generate gene count heatmaps
#'
#'
#'
#'
#'
#'
#'
#'
PlotHeatmaps <- function(dds, res, rld, level, g1, g2, plot.annos = NULL) {
  # plot.annos should be a vector of at least two columns from the sample sheet to use for labeling plots. 
  # g1 and g2 should be strings for the groups to be compared for the level (eg. "Progression", "Response")
  # dds must be a DESeqExperiment object.
  # rld is the log-transformed dds object.
  
  ## Should really make another generic function for the actual heatmap plots and feed it in the pval/lfc values.
  
  # Set color breaks and palette.
  breaks <- c(seq(-3, -1.251, length=250), seq(-1.25, -0.1001, length=250), 
    seq(-0.1, 0.1, length=1), seq(0.1001, 1.25, length=250), 
    seq(1.251, 3, length=250))
  colors <- colorRampPalette(c("#053061","#2166ac", "#f5f5f5", 
    "#b2182b", "#67001f"))(n = 1000)
  
  ### Now let's throw in a LFC magnitude threshold. 
  resSig <- subset(res, padj <= 0.1)
  resSig <- subset(resSig, log2FoldChange >= 1 | log2FoldChange <= -1)
  if (nrow(resSig) > 5) {
    resDEG <- row.names(resSig)
    x <- assay(rld)
    matrix <- x[resDEG,]

    pheatmap(matrix, annotation_col=annotation_data, col=colors, scale="row", 
      breaks=breaks, show_rownames = F, 
      main="ALL DE Genes - padj <= 0.1 & LFC >|< 1 - RLD", fontsize_row=3,
      fontsize_col=5)
    pheatmap(matrix, annotation_col=annotation_data, col=colors, scale="row", 
      breaks=breaks, show_rownames = F, 
      main="ALL DE Genes - padj <= 0.1 & LFC >|< 1 - RLD", cluster_cols=F, 
      fontsize_row=3, fontsize_col=5)

    # Now the same thing but only with the wanted samples.
    x <- assay(rld)
    x.sub <- x[, colData(rld)[,level] %in% c(g1, g2)]
    if (!(ncol(x) == ncol(x.sub))) {

      # Set which columns we want to use for annotating samples.
      annotation_data <- as.data.frame(colData(rld)[plot.annos])

      matrix <- x.sub[res100,]
      # Plot the top 100.
      pheatmap(matrix, annotation_col=annotation_data, col=colors, scale="row", 
        breaks=breaks, show_rownames = T, 
        main="Top 100 DE Genes - padj <= 0.1 & LFC >|< 2 - RLD", fontsize_row=4, 
        fontsize_col=5)
      pheatmap(matrix, annotation_col=annotation_data, col=colors, scale="row", 
        breaks=breaks, show_rownames = T, 
        main="Top 100 DE Genes - padj <= 0.1 & LFC >|< 2 - RLD", cluster_cols=F, 
        fontsize_row=4, fontsize_col=5)

      # Now heatmaps for all DEGs with the variance-stabilized counts.
      resSig <- subset(res, padj <= 0.1)
      resSig <- subset(resSig, log2FoldChange >= 2 | log2FoldChange <= -2)
      if (nrow(resSig) > 5) {
        resDEG <- row.names(resSig)
        x <- assay(rld)
        x.sub <- x[, colData(rld)[,level] %in% c(g1, g2)]
        matrix <- x.sub[resDEG,]
        pheatmap(matrix, annotation_col=annotation_data, col=colors, 
          scale="row", breaks=breaks, show_rownames = F, 
          main="ALL DE Genes - padj <= 0.1 & LFC >|< 2 - RLD", fontsize_row=3, 
          fontsize_col=5)
        pheatmap(matrix, annotation_col=annotation_data, col=colors, 
          scale="row", breaks=breaks, show_rownames = F, 
          main="ALL DE Genes - padj <= 0.1 & LFC >|< 2 - RLD", cluster_cols=F, 
          fontsize_row=3, fontsize_col=5)
      }
    }
  }
}