#' Perform and visualize different variance stabilization transformations
#'
#' @importFrom grDevices pdf dev.off
#' @importFrom DESeq2 rlog vst estimateSizeFactors
#' @importFrom ggplot2 ggplot aes geom_hex facet_grid coord_fixed
#' @importFrom dplyr bind_rows as_data_frame mutate
#'
RunVarianceTransformations <- function (dds, base) {
  rld <- rlog(dds, blind = FALSE)
  vsd <- vst(dds, blind = FALSE)

  pdf(paste0(base,"/GenericFigures/transformation.pdf"))

  dds <- estimateSizeFactors(dds)

  df <- bind_rows(
    as_data_frame(log2(counts(dds, normalized=TRUE)[, 1:2] + 1)) %>%
           mutate(transformation = "log2(x + 1)"),
    as_data_frame(assay(rld)[, 1:2]) %>% mutate(transformation = "rlog"),
    as_data_frame(assay(vsd)[, 1:2]) %>% mutate(transformation = "vst"))
  
  colnames(df)[1:2] <- c("x", "y")  

  p <- ggplot(df, aes(x = x, y = y)) + geom_hex(bins = 80) +
    coord_fixed() + facet_grid( . ~ transformation) 
  print(p)

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
VizSampleDistances <- function(rld, base, level, g1 = NULL, g2 = NULL) {
  dists.out <- paste0(base, "/GenericFigures/samp_dist.pdf")
  message(dists.out)
  pdf(dists.out)
  sampleDists <- dist(t(assay(rld)))

  sampleDistMatrix <- as.matrix(sampleDists)
  rownames(sampleDistMatrix) <- paste(colData(rld)[,level], rld$name, 
    sep = " - ")
  colnames(sampleDistMatrix) <- NULL
  colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
  p <- pheatmap(sampleDistMatrix,
           clustering_distance_rows = sampleDists,
           clustering_distance_cols = sampleDists,
           col = colors)
  print(p)

  if (!is.null(g1) & !is.null(g2)) {
    # Now with only the groups we want to compare.
    rld.sub = rld[, colData(rld)[, level] %in% c(g1, g2)]
    sampleDists <- dist(t(assay(rld.sub)))
    sampleDistMatrix <- as.matrix(sampleDists)
    rownames(sampleDistMatrix) <- paste(colData(rld.sub)[, level], rld.sub$name, 
      sep = " - ")
    p <- pheatmap(sampleDistMatrix,
             clustering_distance_rows = sampleDists,
             clustering_distance_cols = sampleDists,
             col = colors)
    print(p)
  }
  dev.off()
}