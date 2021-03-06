#' Generate per comparison DEG heatmaps for a list of DESeq2 results
#'
#' @param res.list Named List containing \linkS4class{DESeqResults} objects for 
#'   all comparisons generated by \code{\link{ProcessDEGs}}.
#' @param rld A \linkS4class{RangedSummarizedExperiment} object of 
#'   \code{\link[DESeq2]{rlog}} transformed counts as returned by
#'   \link{RunDESeq2}.
#' @param vsd A \linkS4class{RangedSummarizedExperiment} object of 
#'   \code{\link[DESeq2]{vst}} transformed counts as returned by
#'   \link{RunDESeq2}.
#' @param level String defining variable of interest.
#' @param outpath Path to directory to be used for output. 
#' @param padj.thresh Number indicating the adjusted p-value 
#'   cutoff to be used for determining "significant" differential expression.
#' @param fc.thresh Number indicating the log2 fold-change 
#'   cutoff to be used for determining "significant" differential expression.
#' @param plot.annos String or character vector defining the column(s) to use to 
#'   annotate figures.
#'
#' @importFrom pheatmap pheatmap
#' @importFrom grDevices colorRampPalette pdf dev.off
#' @importFrom SummarizedExperiment assay colData
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
#'
#' @export
#'
#' @author Jared Andrews
#'
PlotRNAHeatmaps <- function(res.list, rld, vsd, level, outpath, 
  padj.thresh, fc.thresh, plot.annos) {
  
  # Set color breaks and palette.
  breaks <- c(seq(-3, -1.251, length = 250), seq(-1.25, -0.1001, length = 250), 
    seq(-0.1, 0.1, length = 1), seq(0.1001, 1.25, length = 250), 
    seq(1.251, 3, length = 250))
  colors <- colorRampPalette(c("#053061", "#2166ac", "#f5f5f5", 
    "#b2182b", "#67001f"))(n = 1000)
  
  for (i in seq_along(res.list)) {
    res <- res.list[[i]]
    comp <- names(res.list[i])
    message("Comparison set: ", comp)
    message(paste0("Creating heatmaps for genes with padj <= ", 
      padj.thresh, " and abs(log2FoldChange) >= ", fc.thresh))
    g1 <- unlist(strsplit(comp, "-v-"))[1]
    g2 <- unlist(strsplit(comp, "-v-"))[2]
    
    ressig <- res[(res[, 'padj'] <= padj.thresh) %in% TRUE & 
      abs(res[, 'log2FoldChange']) >= fc.thresh, ]

    # Assay index holder.
    ind <- 1
    labs <- c("rlog", "vst")

    pdf(paste0(outpath, "/Heatmaps/", comp, ".padj.", padj.thresh, ".log2fc.", 
      fc.thresh, ".Heatmaps.pdf"), height = 7, width = 5)

    if (nrow(ressig) > 1) {
      for (x in list(rld, vsd)) {
        # Set which columns we want to use for annotating samples.
        annotation_data <- as.data.frame(colData(x)[,plot.annos])

        resdeg <- row.names(ressig)
        counts <- assay(x)
        matrix <- counts[resdeg,]

        pheatmap(matrix, annotation_col = annotation_data, color = colors, 
          scale = "row", breaks = breaks, show_rownames = FALSE, 
          main = paste0(comp, " - ", labs[ind]), 
          fontsize_col = 5)
        pheatmap(matrix, annotation_col = annotation_data, color = colors, 
          scale = "row", breaks = breaks, show_rownames = FALSE, 
          main = paste0(comp, " - ", labs[ind]), 
          cluster_cols = FALSE, fontsize_col = 5)

        counts.sub <- counts[, colData(x)[,level] %in% c(g1, g2)]
        matrix <- counts.sub[resdeg,]

        if (!(ncol(counts) == ncol(counts.sub))) {
          pheatmap(matrix, annotation_col = annotation_data, color = colors, 
            scale = "row", breaks = breaks, show_rownames = FALSE, 
            main = paste0(comp, " - ", labs[ind]), 
            fontsize_col = 5)
          pheatmap(matrix, annotation_col = annotation_data, color = colors, 
            scale = "row", breaks = breaks, show_rownames = FALSE, 
            main = paste0(comp, " - ", labs[ind]),  
            cluster_cols = FALSE, fontsize_col = 5)
        }
        ind <- ind + 1
      }
    } else {
      message("Skipping, no DEGs.")
    }
    dev.off()
  }
}


#' Generate gene heatmaps for all DEGs from a list of DESeq2 results
#'
#' @param res.list Named List containing \linkS4class{DESeqResults} objects for 
#'   all comparisons generated by \code{\link{ProcessDEGs}}.
#' @param rld A \linkS4class{RangedSummarizedExperiment} object of 
#'   \code{\link[DESeq2]{rlog}} transformed counts as returned by
#'   \link{RunDESeq2}.
#' @param vsd A \linkS4class{RangedSummarizedExperiment} object of 
#'   \code{\link[DESeq2]{vst}} transformed counts as returned by
#'   \link{RunDESeq2}.
#' @param outpath Path to directory to be used for output. 
#' @param padj.thresh Number indicating the adjusted p-value 
#'   cutoff to be used for determining "significant" differential expression.
#' @param fc.thresh Number indicating the log2 fold-change 
#'   cutoff to be used for determining "significant" differential expression.
#' @param plot.annos String or character vector defining the column(s) to use to 
#'   annotate figures.
#'
#' @importFrom pheatmap pheatmap
#' @importFrom grDevices colorRampPalette pdf dev.off
#' @importFrom SummarizedExperiment assay colData
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
#'
#' @export
#'
#' @author Jared Andrews
#'
PlotRNACombinedHeatmaps <- function(res.list, rld, vsd, outpath, 
  padj.thresh, fc.thresh, plot.annos) {

  message("Generating combined DEG heatmap using genes with p.adj <= ", 
      padj.thresh, " and abs(log2FoldChange) >= ", fc.thresh)

  # Set color breaks and palette.
  breaks <- c(seq(-2.5, -1.251, length=250), seq(-1.25, -0.1001, length=250), 
    seq(-0.1, 0.1, length=1), seq(0.1001, 1.25, length=250), 
    seq(1.251, 2.5, length=250))
  colors <- colorRampPalette(c("#053061","#2166ac", "#f5f5f5", 
    "#b2182b", "#67001f"))(n = 1000)

  # Get all DEGs from all comparisons.
  siggy <- unique(unlist(sapply(res.list, 
    function(x, qval, logfc) {
      rownames(x[(x[, 'padj'] <= qval) %in% TRUE & 
        abs(x[, 'log2FoldChange']) >= logfc, ])
    }, qval = padj.thresh, logfc = fc.thresh)))

  # Assay index holder.
  ind <- 1
  labs <- c("rlog", "vst")

  pdf(paste0(outpath, "/Heatmaps/AllComparisons.padj.", padj.thresh, ".log2fc.", 
      fc.thresh, ".Heatmaps.pdf"), height = 7, width = 5)

  # Actual plotting.
  for (x in list(rld, vsd)) {
    heat <- assay(x)[siggy,]
    df <- as.data.frame(colData(rld)[, plot.annos])
    p <- pheatmap(heat, cluster_rows = TRUE, show_rownames = FALSE,
      cluster_cols = FALSE, annotation_col = df, fontsize_col = 6, 
      fontsize = 6, scale = "row", color = colors, breaks = breaks,
      main = paste0("All Comparisons - ", labs[ind]))
    print(p)
    
    p <- pheatmap(heat, cluster_rows = TRUE, show_rownames = FALSE,
      cluster_cols = TRUE, annotation_col = df, fontsize_col = 6, 
      fontsize = 6, scale = "row", color = colors, breaks = breaks,
      main = paste0("All Comparisons - ", labs[ind]))
    print(p)

    # Save counts tables.
    write.table(heat, file = paste0(outpath, "/Heatmaps/AllComparisons.padj.", 
      padj.thresh, ".log2fc.", fc.thresh, ".", labs[ind], ".Heatmaps.txt"), 
      quote = FALSE, sep = "\t")

    ind <- ind + 1
  }

  dev.off()
}


#' Plot gene counts as box plots
#'
#'
#' @param res.list Named List containing \linkS4class{DESeqResults} objects for 
#'   all comparisons generated by \code{\link{ProcessDEGs}}.
#' @param dds A \linkS4class{DESeqDataSet} object as returned by 
#'   \code{\link[DESeq2]{DESeq}} or \link{RunDESeq2}.
#' @param rld A \linkS4class{RangedSummarizedExperiment} object of 
#'   \code{\link[DESeq2]{rlog}} transformed counts as returned by
#'   \link{RunDESeq2}.
#' @param outpath Path to directory to be used for output. 
#' @param padj.thresh Number indicating the adjusted p-value 
#'   cutoff to be used for determining "significant" differential expression.
#' @param fc.thresh Number indicating the log2 fold-change 
#'   cutoff to be used for determining "significant" differential expression.
#' @param top.n Number of differentially expressed genes to create boxplots for, 
#'   ranked by adj. p-value after applying \code{padj.thresh} and 
#'   \code{fc.thresh} thresholds.
#' @param level String defining variable of interest.
#'
#' @import ggplot2
#' @importFrom DESeq2 counts plotCounts
#' @importFrom SummarizedExperiment colData
#'
#' @export
#'
#' @author Jared Andrews
#'
PlotRNABoxplots <- function(res.list, dds, rld, outpath, padj.thresh, fc.thresh, 
  top.n, level) {

  for (i in seq_along(res.list)) {
    res <- res.list[[i]]
    comp <- names(res.list[i])
    message("Comparison set: ", comp)
    g1 <- unlist(strsplit(comp, "-v-"))[1]
    g2 <- unlist(strsplit(comp, "-v-"))[2]

    resdata <- merge(as.data.frame(res), 
      as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
    names(resdata)[1] <- "Gene"

    # Make table for just DEGs, order by padj and subset by top.n.
    ressig <- resdata[(resdata[, 'padj'] <= padj.thresh) %in% TRUE & 
      abs(resdata[, 'log2FoldChange']) >= fc.thresh, ]

    ressig <- ressig[order(ressig$padj),]
    if (nrow(ressig) >= top.n) {
      ressig <- ressig[1:top.n,]
    }

    message(paste0("Creating boxplots for top ", top.n, " genes with padj <= ", 
      padj.thresh, " and abs(log2FoldChange) >= ", fc.thresh))

    # If you have >7 levels for your contrast, you need to add colors here.
    fill = c("#A6CEE3", "#B2DF8A", "#FDBF6F", "#CAB2D6", 
      "#f4cae4", "#f1e2cc", "#b3e2cd")

    if (nrow(ressig) > 1) {
      for (i in 1:nrow(ressig)) {
        pdf(paste0(outpath, "/GeneBoxPlots/", gsub('/','-',ressig$Gene[i]), ".",
          comp, ".padj.", padj.thresh, ".log2fc.", fc.thresh, ".BoxPlot.pdf"), 
          height = 5)
        d <- plotCounts(dds, gene = ressig$Gene[i], intgroup = level, 
          returnData = T, pc = 1)
        p <- ggplot(d, aes(x = d[,level], y = count)) + 
          geom_boxplot(fill = fill[1:length(levels(colData(rld)[,level]))], 
          colour = "#000000") + ggtitle(ressig$Gene[i]) + 
          coord_trans(y = "log10") + theme_classic() +
          xlab("Group") + ylab("Normalized Counts") + theme(aspect.ratio=1)
        print(p)

        e <- d[d[, level] == g1 | d[, level] == g2, ]
        p <- ggplot(e, aes(x = e[,level], y = count)) + 
          geom_boxplot(fill = fill[1:2], colour = "#000000") + 
          ggtitle(ressig$Gene[i]) + coord_trans(y = "log10") + theme_classic() +
          xlab("Group") + ylab("Normalized Counts") + theme(aspect.ratio=1)
        print(p)

        # Make them for log2 counts as well.
        d$count <- log2(d$count)
        p <- ggplot(d, aes(x = d[,level], y = count)) + 
          geom_boxplot(fill = fill[1:length(levels(colData(rld)[,level]))], 
            colour = "#000000") + ggtitle(ressig$Gene[i]) + theme_classic() +
          xlab("Group") + ylab("log2(Normalized Counts)") + 
          theme(aspect.ratio=1)
        print(p)

        e <- d[d[, level] == g1 | d[, level] == g2, ]
        p <- ggplot(e, aes(x = e[,level], y = count)) + 
          geom_boxplot(fill = fill[1:2], colour = "#000000") + 
          ggtitle(ressig$Gene[i]) + theme_classic() +
          xlab("Group") + ylab("log2(Normalized Counts)") + 
          theme(aspect.ratio=1)
        print(p)
        dev.off()
      }
    } else {
      message("No DEGs, skipping boxplots.")
    }
  }
}


#' Plot PCAs from variance stabilized counts for differentially expressed genes
#'
#' @param res.list Named List containing \linkS4class{DESeqResults} objects for 
#'   all comparisons generated by \code{\link{ProcessDEGs}}.
#' @param rld A \linkS4class{RangedSummarizedExperiment} object of 
#'   \code{\link[DESeq2]{rlog}} transformed counts as returned by
#'   \link{RunDESeq2}.
#' @param vsd A \linkS4class{RangedSummarizedExperiment} object of 
#'   \code{\link[DESeq2]{vst}} transformed counts as returned by
#'   \link{RunDESeq2}.
#' @param outpath Path to directory to be used for output.
#' @param level String defining variable of interest.
#' @param plot.annos String or character vector defining the column(s) in 
#'   \code{samplesheet} to use to annotate figures.
#' @param padj.thresh Number indicating the adjusted p-value 
#'   cutoff to be used for determining "significant" differential expression.
#' @param fc.thresh Number indicating the log2 fold-change 
#'   cutoff to be used for determining "significant" differential expression.
#'
#' @importFrom grDevices pdf dev.off
#' @importFrom ggplot2 ggtitle theme_classic theme
#' @importFrom utils combn
#' @importFrom SummarizedExperiment colData
#'
#' @export
#'
#' @author Jared Andrews
#'
PlotRNADEGPCAs <- function(res.list, rld, vsd, outpath, level, plot.annos, 
  padj.thresh, fc.thresh) {

  for (i in seq_along(res.list)) {
    res <- res.list[[i]]
    comp <- names(res.list[i])
    message("Comparison set: ", comp)
    message("Generating PCAs using genes with p.adj <= ", 
      padj.thresh, " and abs(log2FoldChange) >= ", fc.thresh)
    g1 <- unlist(strsplit(comp, "-v-"))[1]
    g2 <- unlist(strsplit(comp, "-v-"))[2]

    pdf(paste0(outpath, "/DEGFigures/", comp, ".padj.", padj.thresh, ".log2fc.", 
      fc.thresh, ".PCA.pdf"), height = 5, width = 5)
    ind <- 1
    labs <- c("rlog", "vst")

    for (x in list(rld, vsd)) {
      resdata <- merge(as.data.frame(res), 
        as.data.frame(assay(x)), by="row.names", sort=FALSE)
      names(resdata)[1] <- "Gene"

      # Get DEGs.
      ressig <- resdata[(resdata[, 'padj'] <= padj.thresh) %in% TRUE & 
        abs(resdata[, 'log2FoldChange']) >= fc.thresh, ]

      x.sub <- x[ressig$Gene, colData(x)[,level] %in% c(g1, g2)]

      p <- DESeq2::plotPCA(x.sub, intgroup = level) +
        ggtitle(paste0("DEGs - ", labs[ind], " - ", comp)) + theme_classic() +
        theme(aspect.ratio=1)
      print(p)

      if (!identical(plot.annos, level)) {
        p <- DESeq2::plotPCA(x.sub, intgroup = plot.annos) +
          ggtitle(paste0("DEGs - ", labs[ind], " - ", comp)) + theme_classic() +
          theme(aspect.ratio=1)
        print(p)
      }

      ind <- ind + 1
    }
    dev.off()
  }
}


#' Create volcano plots from a list of DESeq2 Results objects
#'
#' @param res.list Named List containing \linkS4class{DESeqResults} objects for 
#'   all comparisons generated by \code{\link{ProcessDEGs}}.
#' @param dds A \linkS4class{DESeqDataSet} object as returned by 
#'   \code{\link[DESeq2]{DESeq}}, \link{RunDESeq2}, or \link{ProcessDEGs}.
#' @param outpath Path to directory to be used for output.
#' @param padj.thresh Number indicating the adjusted p-value 
#'   cutoff to be used for determining "significant" differential expression.
#' @param fc.thresh Number indicating the log2 fold-change 
#'   cutoff to be used for determining "significant" differential expression.
#' @param n.labels Number of genes to label on the volcano plot, ranked by 
#'   adjusted p-value.
#' @param use.labels Boolean indicating whether labels should be added to plot.
#'
#' @importFrom grDevices pdf dev.off
#' @importFrom plotly plot_ly as_widget layout
#' @importFrom htmlwidgets saveWidget
#' @importFrom EnhancedVolcano EnhancedVolcano
#'
#' @export
#'
#' @author Jared Andrews
#'
PlotRNAVolcanoes <- function(res.list, dds, outpath, padj.thresh, fc.thresh, 
  n.labels = 10, use.labels = TRUE) {

  for (i in seq_along(res.list)) {
    res <- res.list[[i]]
    comp <- names(res.list[i])
    g1 <- unlist(strsplit(comp, "-v-"))[1]
    g2 <- unlist(strsplit(comp, "-v-"))[2]
    message("Comparison set: ", comp)
    message("Generating volcano plots with p.adj <= ", 
      padj.thresh, " and abs(log2FoldChange) >= ", fc.thresh)

    # Set up results table.
    res <- res[order(res$padj),]
    resdata <- merge(as.data.frame(res), 
      as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
    names(resdata)[1] <- "Gene"
    
    diff_df <- resdata[c("Gene","log2FoldChange","padj")]

    # Add a grouping column; default value is "not significant"
    diff_df["group"] <- "Not Significant"

    # Change groupings as necessary.
    diff_df[which(diff_df['padj'] <= padj.thresh & 
      abs(diff_df['log2FoldChange']) <= fc.thresh), "group"] <- "Significant"
    diff_df[which(diff_df['padj'] >= padj.thresh & 
      abs(diff_df['log2FoldChange']) >= fc.thresh), "group"] <- "FoldChange"
    diff_df[which(diff_df['padj'] <= padj.thresh & 
      abs(diff_df['log2FoldChange']) >= fc.thresh), 
      "group"] <- "Significant & FoldChange"

    p <- plot_ly(
      diff_df, 
      x = ~log2FoldChange, 
      y = -log10(diff_df$padj), 
      text = ~Gene, 
      color = ~group,
      mode = "markers", 
      type = "scatter") %>%
      layout(title = paste0(g1, " versus ", g2), 
        xaxis = list(title = "log2 Fold Change"), 
        yaxis = list(title = "-log10(adjusted p-value)")) 

    # Save plot - has to use absolute path for some reason.
    saveWidget(as_widget(p), paste0(getwd(), "/", 
      substr(outpath, 2, 500), "/DEGFigures/", comp, 
      ".Volcano.padj.", padj.thresh, ".log2fc.", fc.thresh, ".html"))

    # Set axis limits based on max magnitude.
    min.lim <- ceiling(max(abs(diff_df$log2FoldChange)))

    pdf(paste0(outpath, "/DEGFigures/", comp, ".Volcano.padj", padj.thresh, 
      ".log2fc.", fc.thresh, ".pdf"), height = 6, width = 5)

    if (use.labels) {
      labels <- rownames(res[1:n.labels,])
    } else {
      labels <- NULL
    }

    p <- EnhancedVolcano(res,
      lab = rownames(res),
      x = "log2FoldChange",
      y = "padj",
      xlim = c(-min.lim, min.lim),
      title = paste0(g1, " versus ", g2),
      subtitle = "",
      xlab = bquote(~Log[2]~ 'fold change'),
      ylab = bquote(~-Log[10]~adjusted~italic(P)),
      pCutoff = padj.thresh,
      FCcutoff = fc.thresh,
      selectLab = labels,
      transcriptLabSize = 2.0,
      transcriptPointSize = 0.5,
      col = c("#8C8C8C", "#D55E00", "#0072B2", "#009E73"),
      legendLabSize = 8,
      legendIconSize = 2,
      legendLabels = c("NS", "Log2 FC", "adj. P", 
        "Log2 FC & adj. P"), 
      colAlpha = 1,
      legendPosition = "bottom",
      gridlines.major = FALSE,
      gridlines.minor = FALSE)
    print(p)
    dev.off()
  }
}