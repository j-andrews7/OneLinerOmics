#' Generate gene count heatmaps
#'
#' @param outpath Path to directory to be used for output. 
#' @param level String defining variable of interest.
#' @param padj.thresh Number or numeric scalar indicating the adjusted p-value 
#'   cutoff(s) to be used for determining "significant" differential expression.
#'   If multiple are given, multiple tables/plots will be generated using all 
#'   combinations of \code{padj.thresh} and \code{fc.thresh}.
#' @param fc.thresh Number or numeric scalar indicating the log2 fold-change 
#'   cutoff(s) to be used for determining "significant" differential expression.
#'   If multiple are given, multiple tables/plots will be generated using all 
#'   combinations of \code{padj.thresh} and \code{fc.thresh}.
#' @param outpath Path to directory to be used for output. 
#' @param plot.annos String or character vector defining the column(s) to use to 
#'   annotate figures.
#'
#' @importFrom pheatmap pheatmap
#' @importFrom grDevices colorRampPalette pdf dev.off
#' @importFrom SummarizedExperiment assay
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
#'
#' @author Jared Andrews
#'
PlotHeatmaps <- function(dds, res, rld, vsd, level, g1, g2, padj.thresh, 
  fc.thresh, outpath, plot.annos = NULL) {
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

    pheatmap(matrix, annotation_col=annotation_data, color=colors, scale="row", 
      breaks=breaks, show_rownames = F, 
      main="ALL DE Genes - padj <= 0.1 & LFC >|< 1 - RLD", fontsize_row=3,
      fontsize_col=5)
    pheatmap(matrix, annotation_col=annotation_data, color=colors, scale="row", 
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
      pheatmap(matrix, annotation_col=annotation_data, color=colors, 
        scale="row", breaks=breaks, show_rownames = T, 
        main="Top 100 DE Genes - padj <= 0.1 & LFC >|< 2 - RLD", fontsize_row=4, 
        fontsize_col=5)
      pheatmap(matrix, annotation_col=annotation_data, color=colors, 
        scale="row", breaks=breaks, show_rownames = T, 
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
        pheatmap(matrix, annotation_col=annotation_data, color=colors, 
          scale="row", breaks=breaks, show_rownames = F, 
          main="ALL DE Genes - padj <= 0.1 & LFC >|< 2 - RLD", fontsize_row=3, 
          fontsize_col=5)
        pheatmap(matrix, annotation_col=annotation_data, color=colors, 
          scale="row", breaks=breaks, show_rownames = F, 
          main="ALL DE Genes - padj <= 0.1 & LFC >|< 2 - RLD", cluster_cols=F, 
          fontsize_row=3, fontsize_col=5)
      }
    }
  }
}

#' Plot gene counts as box plots
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'
#' @import ggplot2
#' @importFrom DESeq2 counts plotCounts
#'
#' @author Jared Andrews
#'
PlotBoxplots <- function(res.list, dds, outpath, padj.thresh, fc.thresh) {

  resdata <- merge(as.data.frame(res), 
    as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
  names(resdata)[1] <- "Gene"

  # Make table for just DEGs.
  resSig <- subset(resdata, padj <= 0.1)

  message("Creating boxplots for all genes with padj <= 0.1")
  # If you have >7 levels for your contrast, you need to add colors here.
  fill = c("#A6CEE3", "#B2DF8A", "#FDBF6F", "#CAB2D6", 
    "#f4cae4", "#f1e2cc", "#b3e2cd")
  line = c("#1F78B4", "#33A02C", "#FF7F00", "#6A3D9A", 
    "#e7298a", "#a6761d", "#1b9e77")
  if (nrow(resSig) > 1) {
    for (i in 1:nrow(resSig)) {
      if (!file.exists(paste0(base, "/GeneBoxPlots/", 
        gsub('/','-',resSig$Gene[i]),".BoxPlot.pdf"))) {
        pdf(paste0(base, "/GeneBoxPlots/", gsub('/','-',resSig$Gene[i]), 
          ".BoxPlot.pdf"))
        d <- plotCounts(dds, gene = resSig$Gene[i], intgroup = level, 
          returnData = T)
        p <- ggplot(d, aes(x = d[,level], y = count)) + 
          geom_boxplot(fill = fill[1:length(levels(colData(rld)[,level]))], 
            colour = line[1:length(levels(colData(rld)[,level]))]) + 
          ggtitle(resSig$Gene[i]) + coord_trans(y = "log10")
        print(p)

        d <- plotCounts(dds, gene = resSig$Gene[i], intgroup = level, 
          returnData = T)
        e <- subset(d, (get(level) == g1 | get(level) == g2))
        p <- ggplot(e, aes(x = e[,level], y = count)) + 
          geom_boxplot(fill=fill[1:2],colour=line[1:2]) + 
          ggtitle(resSig$Gene[i]) + coord_trans(y = "log10")
        print(p)
        dev.off()
      }
    }
  } else {
    message("No DEGs, skipping boxplots.")
  }
}