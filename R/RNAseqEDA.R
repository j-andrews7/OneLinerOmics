#' Perform and visualize different variance stabilization transformations
#'
#' @param dds A \linkS4class{DESeqDataSet} object.
#' @param outpath Path to output file.
#'
#' @importFrom grDevices pdf dev.off
#' @importFrom DESeq2 rlog vst estimateSizeFactors normTransform
#' @importFrom ggplot2 ggplot aes_string geom_hex facet_grid coord_fixed ggtitle
#' @importFrom dplyr bind_rows as_data_frame mutate
#' @importFrom vsn meanSdPlot
#' @importFrom SummarizedExperiment assay
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
#' @importFrom magrittr %>%
#'
#' @export
#'
#' @author Jared Andrews
#'
PlotRNAVarianceTransformations <- function (dds, outpath) {
  message("This may take a while if you have many samples.")
  rld <- suppressMessages(rlog(dds, blind = FALSE))
  vsd <- suppressMessages(vst(dds, blind = FALSE))

  pdf(outpath)

  dds <- estimateSizeFactors(dds)

  df <- bind_rows(
    as_data_frame(log2(counts(dds, normalized=TRUE)[, 1:2])) %>%
           mutate(transformation = "log2(x + 1)"),
    as_data_frame(assay(rld)[, 1:2]) %>% mutate(transformation = "rlog"),
    as_data_frame(assay(vsd)[, 1:2]) %>% mutate(transformation = "vst"))

  colnames(df)[1:2] <- c("x", "y")

  p <- ggplot(df, aes_string(x = "x", y = "y")) + geom_hex(bins = 80) +
    coord_fixed() + facet_grid( . ~ transformation)
  print(p)

  # Plots to compare variance and counts.
  ntd <- normTransform(dds)

  ntd.p <- meanSdPlot(assay(ntd), plot = FALSE)
  rld.p <- meanSdPlot(assay(rld), plot = FALSE)
  vsd.p <- meanSdPlot(assay(vsd), plot = FALSE)

  ntd.p$gg <- ntd.p$gg + ggtitle("Transformation: log2(x + 1)")
  rld.p$gg <- rld.p$gg + ggtitle("Transformation: regularized log (rlog)")
  vsd.p$gg <- vsd.p$gg +
    ggtitle("Transformation: variance stabilizing transformation")

  print(ntd.p$gg)
  print(rld.p$gg)
  print(vsd.p$gg)

  dev.off()

  return(list(rld = rld, vsd = vsd))
}


#' Visualize sample distances from variance stabilized counts
#'
#' @param rld A \linkS4class{RangedSummarizedExperiment} object of 
#'   \code{\link[DESeq2]{rlog}} transformed counts as returned by
#'   \link{PlotRNAVarianceTransformations}.
#' @param vsd A \linkS4class{RangedSummarizedExperiment} object of 
#'   \code{\link[DESeq2]{vst}} transformed counts as returned by
#'   \link{PlotRNAVarianceTransformations}.
#' @param outpath Path to output file.
#' @param level String defining variable of interest.
#' @param plot.annos String or character vector defining the column(s) in 
#'   \code{samplesheet} to use to annotate figures.
#'
#' @importFrom grDevices pdf dev.off
#' @importFrom stats dist
#' @importFrom pheatmap pheatmap
#' @importFrom SummarizedExperiment assay colData
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
#' @importFrom grDevices colorRampPalette
#' @importFrom RColorBrewer brewer.pal
#'
#' @export
#'
#' @author Jared Andrews
#'
PlotRNASampleDistances <- function(rld, vsd, outpath, level, plot.annos) {

  pdf(outpath)
  i <- 1

  labs <- c("Sample Distances (rlog)",
    "Sample Distances (vst)")

  for (x in list(rld, vsd)) {
    sampleDists <- dist(t(assay(x)))
    sampleDistMatrix <- as.matrix(sampleDists)
    rownames(sampleDistMatrix) <- x$name
    colnames(sampleDistMatrix) <- NULL
    colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)

    # Get annotation data.
    annotation.data <- as.data.frame(colData(x)[ ,plot.annos])

    p <- pheatmap(sampleDistMatrix,
      clustering_distance_rows = sampleDists,
      clustering_distance_cols = sampleDists,
      color = colors, main = labs[i], annotation_row = annotation.data)
    print(p)

    i <- i + 1
  }

  dev.off()
}


#' Plot PCAs from variance stabilized counts for EDA
#'
#' @param rld A \linkS4class{RangedSummarizedExperiment} object of 
#'   \code{\link[DESeq2]{rlog}} transformed counts as returned by
#'   \link{PlotRNAVarianceTransformations}.
#' @param vsd A \linkS4class{RangedSummarizedExperiment} object of 
#'   \code{\link[DESeq2]{vst}} transformed counts as returned by
#'   \link{PlotRNAVarianceTransformations}.
#' @param outpath Path to output file.
#' @param level String defining variable of interest.
#' @param plot.annos String or character vector defining the column(s) in 
#'   \code{samplesheet} to use to annotate figures.
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
PlotRNAEDAPCAs <- function(rld, vsd, outpath, level, plot.annos) {

  pdf(outpath, height = 5, width = 5)
  i <- 1

  labs <- c("All Genes (rlog)",
    "All Genes (vst)")

  # Get all possible comparisons.
  combs <- combn(levels(colData(rld)[,level]), 2)
  combs.seq <- seq(1, length(combs), by = 2)

  for (x in list(rld, vsd)) {
    p <- DESeq2::plotPCA(x, intgroup = level) + ggtitle(labs[i]) +
      theme_classic() + theme(aspect.ratio = 1)
    print(p)

    if (plot.annos != level) {
      p <- DESeq2::plotPCA(x, intgroup = plot.annos) + ggtitle(labs[i]) +
        theme_classic() + theme(aspect.ratio = 1)
      print(p)
    }

    # PCA for all possible sample comparisons.
    for (samp in combs.seq) {
      x.sub <- x[, colData(x)[, level] %in% c(combs[samp], combs[samp + 1])]

      p <- DESeq2::plotPCA(x.sub, intgroup = level) +
        ggtitle(paste0(labs[i], " - ", combs[samp], " v ", combs[samp + 1])) +
        theme_classic() + theme(aspect.ratio = 1)
      print(p)

      if (plot.annos != level) {
        p <- DESeq2::plotPCA(x.sub, intgroup = plot.annos) +
          ggtitle(paste0(labs[i], " - ", combs[samp], " v ", combs[samp + 1])) +
          theme_classic() + theme(aspect.ratio = 1)
        print(p)
      }
    }

    i <- i + 1
  }
  dev.off()
}
