#' Generate gene heatmaps for a list of DESeq2 results
#'
#' @param level String defining variable of interest.
#' @param outpath Path to directory to be used for output. 
#' @param padj.thresh Number or numeric scalar indicating the adjusted p-value 
#'   cutoff(s) to be used for determining "significant" differential expression.
#' @param fc.thresh Number or numeric scalar indicating the log2 fold-change 
#'   cutoff(s) to be used for determining "significant" differential expression.
#' @param plot.annos String or character vector defining the column(s) to use to 
#'   annotate figures.
#'
#' @importFrom pheatmap pheatmap
#' @importFrom grDevices colorRampPalette pdf dev.off
#' @importFrom SummarizedExperiment assay colData
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
#'
#' @author Jared Andrews
#'
PlotHeatmaps <- function(res.list, rld, vsd, level, outpath, 
  padj.thresh, fc.thresh, plot.annos) {
  
  # Set color breaks and palette.
  breaks <- c(seq(-3, -1.251, length=250), seq(-1.25, -0.1001, length=250), 
    seq(-0.1, 0.1, length=1), seq(0.1001, 1.25, length=250), 
    seq(1.251, 3, length=250))
  colors <- colorRampPalette(c("#053061","#2166ac", "#f5f5f5", 
    "#b2182b", "#67001f"))(n = 1000)
  
  for (i in seq_along(res.list)) {
    res <- res.list[i]
    comp <- names(res.list[i])
    message("Comparison set: ", comp)
    message(paste0("Creating heatmaps for genes with padj <= ", 
      padj.thresh, " and abs(log2FoldChange) >= ", fc.thresh))
    g1 <- unlist(strsplit(comp, "-v-"))[1]
    g2 <- unlist(strsplit(comp, "-v-"))[2]
  
    ressig <- subset(res, padj <= padj.thresh & 
      abs(log2FoldChange) >= fc.thresh)

    # Assay index holder.
    ind <- 1
    labs <- c("rlog", "vst")

    pdf(paste0(outpath, "/Heatmaps/", comp, ".padj.", padj.thresh, ".log2fc.", 
      fc.thresh, ".Heatmaps.pdf"), height = 7, width = 5)

    if (nrow(ressig) > 1) {
      for (x in c(rld, vsd)) {
        # Set which columns we want to use for annotating samples.
        annotation_data <- as.data.frame(colData(x)[,plot.annos])

        resdeg <- row.names(ressig)
        counts <- assay(x)
        matrix <- counts[resdeg,]

        pheatmap(matrix, annotation_col = annotation_data, color = colors, 
          scale = "row", breaks = breaks, show_rownames = FALSE, 
          main = paste0(comp, " - DE Genes - p.adj <= ", padj.thresh, 
            " & abs(log2FC) >= ", fc.thresh, " - ", labs[ind]), 
          fontsize_col=5)
        pheatmap(matrix, annotation_col = annotation_data, color = colors, 
          scale = "row", breaks = breaks, show_rownames = FALSE, 
          main = paste0(comp, " - DE Genes - p.adj <= ", padj.thresh, 
            " & abs(log2FC) >= ", fc.thresh, " - ", labs[ind]), 
          cluster_cols = FALSE, fontsize_col = 5)

        counts.sub <- counts[, colData(x)[,level] %in% c(g1, g2)]
        matrix <- counts.sub[resdeg,]

        if (!(ncol(counts) == ncol(counts.sub))) {
          pheatmap(matrix, annotation_col = annotation_data, color = colors, 
            scale = "row", breaks = breaks, show_rownames = FALSE, 
            main = paste0(comp, " - DE Genes - p.adj <= ", padj.thresh, 
            " & abs(log2FC) >= ", fc.thresh, " - ", labs[ind]), 
            fontsize_col = 5)
          pheatmap(matrix, annotation_col = annotation_data, color = colors, 
            scale = "row", breaks = breaks, show_rownames = FALSE, 
            main = paste0(comp, " - DE Genes - p.adj <= ", padj.thresh, 
            " & abs(log2FC) >= ", fc.thresh, " - ", labs[ind]),  
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
#' @param level String defining variable of interest.
#' @param outpath Path to directory to be used for output. 
#' @param padj.thresh Number or numeric scalar indicating the adjusted p-value 
#'   cutoff(s) to be used for determining "significant" differential expression.
#' @param fc.thresh Number or numeric scalar indicating the log2 fold-change 
#'   cutoff(s) to be used for determining "significant" differential expression.
#' @param plot.annos String or character vector defining the column(s) to use to 
#'   annotate figures.
#'
#' @importFrom pheatmap pheatmap
#' @importFrom grDevices colorRampPalette pdf dev.off
#' @importFrom SummarizedExperiment assay colData
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
#'
#' @author Jared Andrews
#'
PlotCombinedHeatmaps <- function(res.list, rld, vsd, outpath, 
  padj.thresh, fc.thresh, plot.annos) {

  message("Generating combined DEG heatmap using genes with p.adj <= ", 
      padj.thresh, " and abs(log2FoldChange) >= ", fc.thresh)

  # Set color breaks and palette.
  mycol <- colorRampPalette(c("darkblue", "snow", "darkred"))(1000)
  colors <- c(seq(-3, -.11, length=500), seq(-.1, .1, length=1),
    seq(.11, 3, length=500))

  # Get all DEGs from all comparisons.
  siggy <- unique(unlist(sapply(res.list, 
    function(x, qval, logfc) {
      rownames(subset(x, padj <= qval & abs(log2FoldChange) >= logfc))
    }, qval = padj.thresh, logfc = fc.thresh)))

  # Assay index holder.
  ind <- 1
  labs <- c("rlog", "vst")

  pdf(paste0("/Heatmaps/AllComparisons.padj.", padj.thresh, ".log2fc.", 
      fc.thresh, ".Heatmaps.pdf"), height = 7, width = 5)

  # Actual plotting.
  for (x in c(rld, vsd)) {
    heat <- assay(x)[siggy,]
    df <- as.data.frame(colData(rld)[, plot.annos])
    p <- pheatmap(heat, cluster_rows = TRUE, show_rownames=  FALSE,
      cluster_cols = FALSE, annotation_col = df, fontsize_col = 6, 
      fontsize = 6, scale = "row", color = mycol, breaks = colors,
      main = paste0("All Comparisons - DE Genes - p.adj <= ", padj.thresh, 
            " & abs(log2FC) >= ", fc.thresh, " - ", labs[ind]))
    print(p)
    
    p <- pheatmap(heat, cluster_rows = TRUE, show_rownames=  FALSE,
      cluster_cols = TRUE, annotation_col = df, fontsize_col = 6, 
      fontsize = 6, scale = "row", color = mycol, breaks = colors,
      main = paste0("All Comparisons - DE Genes - p.adj <= ", padj.thresh, 
            " & abs(log2FC) >= ", fc.thresh, " - ", labs[ind]))
    print(p)

    ind <- ind + 1
  }
  dev.off()
}


#' Plot gene counts as box plots
#'
#' @import ggplot2
#' @importFrom DESeq2 counts plotCounts
#' @importFrom SummarizedExperiment colData
#'
#' @author Jared Andrews
#'
PlotBoxplots <- function(res.list, dds, rld, outpath, padj.thresh, fc.thresh, 
  top.n, level) {

  for (i in seq_along(res.list)) {
    res <- res.list[i]
    comp <- names(res.list[i])
    message("Comparison set: ", comp)
    g1 <- unlist(strsplit(comp, "-v-"))[1]
    g2 <- unlist(strsplit(comp, "-v-"))[2]

    resdata <- merge(as.data.frame(res), 
      as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
    names(resdata)[1] <- "Gene"

    # Make table for just DEGs, order by padj and subset by top.n.
    ressig <- subset(resdata, padj <= padj.thresh & abs(log2FoldChange) >= 
      fc.thresh)
    ressig <- ressig[order(ressig$padj),]
    if (nrow(ressig) >= top.n) {
      ressig <- ressig[1:top.n,]
    }

    message(paste0("Creating boxplots for top ", top.n, " genes with padj <= ", 
      padj.thresh, " and abs(log2FoldChange) >= ", fc.thresh))

    # If you have >7 levels for your contrast, you need to add colors here.
    fill = c("#A6CEE3", "#B2DF8A", "#FDBF6F", "#CAB2D6", 
      "#f4cae4", "#f1e2cc", "#b3e2cd")
    line = c("#1F78B4", "#33A02C", "#FF7F00", "#6A3D9A", 
      "#e7298a", "#a6761d", "#1b9e77")

    if (nrow(ressig) > 1) {
      for (i in 1:nrow(ressig)) {
        pdf(paste0(outpath, "/GeneBoxPlots/", gsub('/','-',ressig$Gene[i]), 
          comp, ".padj.", padj.thresh, ".log2fc.", fc.thresh, ".BoxPlot.pdf"), 
          height = 5)
        d <- plotCounts(dds, gene = ressig$Gene[i], intgroup = level, 
          returnData = T, pc = 1)
        p <- ggplot(d, aes(x = d[,level], y = count)) + 
          geom_boxplot(fill = fill[1:length(levels(colData(rld)[,level]))], 
            colour = line[1:length(levels(colData(rld)[,level]))]) + 
          ggtitle(ressig$Gene[i]) + coord_trans(y = "log10") + theme_classic() +
          xlab("Group") + ylab("Normalized Counts")
        print(p)

        e <- subset(d, (get(level) == g1 | get(level) == g2))
        p <- ggplot(e, aes(x = e[,level], y = count)) + 
          geom_boxplot(fill=fill[1:2],colour=line[1:2]) + 
          ggtitle(ressig$Gene[i]) + coord_trans(y = "log10") + theme_classic() +
          xlab("Group") + ylab("Normalized Counts")
        print(p)

        # Make them for log2 counts as well.
        d$count <- log2(d$count)
        p <- ggplot(d, aes(x = d[,level], y = count)) + 
          geom_boxplot(fill = fill[1:length(levels(colData(rld)[,level]))], 
            colour = line[1:length(levels(colData(rld)[,level]))]) + 
          ggtitle(ressig$Gene[i]) + theme_classic() +
          xlab("Group") + ylab("log2(Normalized Counts)")
        print(p)

        e <- subset(d, (get(level) == g1 | get(level) == g2))
        p <- ggplot(e, aes(x = e[,level], y = count)) + 
          geom_boxplot(fill=fill[1:2],colour=line[1:2]) + 
          ggtitle(ressig$Gene[i]) + theme_classic() +
          xlab("Group") + ylab("log2(Normalized Counts)")
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
#' @importFrom grDevices pdf dev.off
#' @importFrom ggplot2 ggtitle
#' @importFrom utils combn
#' @importFrom SummarizedExperiment colData
#'
#' @author Jared Andrews
#'
PlotDEGPCAs <- function(res.list, rld, vsd, outpath, level, plot.annos, 
  padj.thresh, fc.thresh) {

  for (i in seq_along(res.list)) {
    res <- res.list[i]
    comp <- names(res.list[i])
    message("Comparison set: ", comp)
    message("Generating PCAs using genes with p.adj <= ", 
      padj.thresh, " and abs(log2FoldChange) >= ", fc.thresh)
    g1 <- unlist(strsplit(comp, "-v-"))[1]
    g2 <- unlist(strsplit(comp, "-v-"))[2]

    pdf(paste0(outpath, "/DEGFigures/", comp, ".padj.", padj.thresh, ".log2fc.", 
      fc.thresh, ".BoxPlot.pdf"), height = 5, width = 5)
    ind <- 1
    labs <- c("rlog", "vst")

    for (x in c(rld, vsd)) {
      resdata <- merge(as.data.frame(res), 
        as.data.frame(assay(x)), by="row.names", sort=FALSE)
      names(resdata)[1] <- "Gene"

      # Get DEGs.
      ressig <- subset(resdata, padj <= padj.thresh & abs(log2FoldChange) >= 
        fc.thresh)


      x.sub <- x[ressig$Gene, colData(x)[,level] %in% c(g1, g2)]
      p <- DESeq2::plotPCA(x, intgroup = level) +
        ggtitle(paste0("DEGs - p.adj <= ", padj.thresh, " & abs(log2FC) > ", 
          fc.thresh, " - ", labs[ind], " - ", comp))
      print(p)

      p <- DESeq2::plotPCA(x.sub, intgroup = level) +
        ggtitle(paste0("DEGs - p.adj <= ", padj.thresh, " & abs(log2FC) > ", 
          fc.thresh, " - ", labs[ind], " - ", comp))
      print(p)

      if (plot.annos != level) {
        p <- DESeq2::plotPCA(x, intgroup = plot.annos) +
          ggtitle(paste0("DEGs - p.adj <= ", padj.thresh, " & abs(log2FC) > ", 
          fc.thresh, " - ", labs[ind], " - ", comp))
        print(p)

        p <- DESeq2::plotPCA(x.sub, intgroup = plot.annos) +
          ggtitle(paste0("DEGs - p.adj <= ", padj.thresh, " & abs(log2FC) > ", 
          fc.thresh, " - ", labs[ind], " - ", comp))
        print(p)
      }

      ind <- ind + 1
    }
    dev.off()
  }
}


#' Create volcano plots from results lists for all genes
#'
#' @importFrom grDevices pdf dev.off
#' @importFrom plotly plot_ly as_widget layout
#' @importFrom htmlwidgets saveWidget
#' @importFrom EnhancedVolcano EnhancedVolcano
#'
#' @author Jared Andrews
#'
PlotVolcanoes <- function(res.list, dds, outpath, padj.thresh, fc.thresh) {

  for (i in seq_along(res.list)) {
    res <- res.list[i]
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

    min.lim <- ceiling(max(abs(diff_df$log2FoldChange)))

    pdf(paste0(outpath, "/DEGFigures/", comp, ".Volcano.padj", padj.thresh, 
      ".log2fc.", fc.thresh, ".pdf"), height = 6, width = 5)
    p <- EnhancedVolcano(res,
      lab = rownames(res),
      x = "log2FoldChange",
      y = "padj",
      xlim = c(-min.lim, min.lim),
      title = paste0(g1, " versus ", g2),
      xlab = bquote(~Log[2]~ 'fold change'),
      ylab = bquote(~-Log[10]~adjusted~italic(P)),
      pCutoff = padj.thresh,
      FCcutoff = fc.thresh,
      transcriptLabSize = 2.0,
      col=c("#8C8C8C", "#D55E00", "#0072B2", "#009E73"),
      colAlpha = 1,
      legendPosition = "bottom")
    print(p)
  }
  dev.off()
}