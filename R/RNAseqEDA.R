#' Perform and visualize different variance stabilization transformations
#'
#' @importFrom grDevices pdf dev.off
#' @importFrom DESeq2 rlog vst estimateSizeFactors normTransform
#' @importFrom ggplot2 ggplot aes geom_hex facet_grid coord_fixed ggtitle
#' @importFrom dplyr bind_rows as_data_frame mutate
#' @importFrom vsn meanSdPlot
#'
ApplyVarianceTransformations <- function (dds, outpath) {
  message("This may take a while if you have many samples.")
  rld <- rlog(dds, blind = FALSE)
  vsd <- vst(dds, blind = FALSE)

  pdf(outpath)

  dds <- estimateSizeFactors(dds)

  df <- bind_rows(
    as_data_frame(log2(counts(dds, normalized=TRUE)[, 1:2])) %>%
           mutate(transformation = "log2(x + 1)"),
    as_data_frame(assay(rld)[, 1:2]) %>% mutate(transformation = "rlog"),
    as_data_frame(assay(vsd)[, 1:2]) %>% mutate(transformation = "vst"))
  
  colnames(df)[1:2] <- c("x", "y")  

  p <- ggplot(df, aes(x = x, y = y)) + geom_hex(bins = 80) +
    coord_fixed() + facet_grid( . ~ transformation) 
  print(p)

  # Plots to compare variance and counts.
  ntd <- normTransform(dds)

  ntd.p <- meanSdPlot(assay(ntd))
  rld.p <- meanSdPlot(assay(rld))
  vsd.p <- meanSdPlot(assay(vsd))

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
#' @importFrom grDevices pdf dev.off
#' @importFrom stats dist
#' @importFrom dichromat colorRampPalette
#' @importFrom pheatmap pheatmap
#'
VizSampleDistances <- function(rld, vsd, outpath, level) {
  pdf(outpath)
  i <- 1

  for (x in c(rld, vsd)) {
    sampleDists <- dist(t(assay(x)))

    sampleDistMatrix <- as.matrix(sampleDists)
    rownames(sampleDistMatrix) <- paste(colData(x)[,level], x$name, 
      sep = " - ")
    colnames(sampleDistMatrix) <- NULL
    colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)
    p <- pheatmap(sampleDistMatrix,
             clustering_distance_rows = sampleDists,
             clustering_distance_cols = sampleDists,
             col = colors)
    print(p)
  }
  dev.off()
}