#' Perform GO/Pathway enrichment for DEGs via enrichR
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
#'
#' 
#' @importFrom EZscRNA RunEnrichr VizEnrichments
#'
#' @author Jared Andrews
#'
PlotEnrichments <- function(res.list, outpath, padj.thresh, 
  fc.thresh) {

  for (r in seq_along(res.list)) {
    res <- res.list[r]
    comp <- names(res.list[r])
    g1 <- unlist(strsplit(comp, "-v-"))[1]
    g2 <- unlist(strsplit(comp, "-v-"))[2]

    # Create directories.
    base <- paste0(oathpath, "/Enrichments/", comp)
    dir.create(file.path(base, "AllGenes"), showWarnings = FALSE, 
      recursive = TRUE)
    dir.create(file.path(base, paste0(g1, "up")), showWarnings = FALSE, 
      recursive = TRUE)
    dir.create(file.path(base, paste0(g2, "up")), showWarnings = FALSE, 
      recursive = TRUE)

    one.two <- subset(res, 
      padj <= padj.thresh & abs(log2FoldChange) >= fc.thresh)
    one.up <- subset(res, 
      padj <= padj.thresh & log2FoldChange >= fc.thresh)
    two.up <- subset(res, 
      padj <= padj.thresh & log2FoldChange <= -fc.thresh)
    
    one.two.terms <- RunEnrichr(rownames(one.two), 
      outdir = paste0(base, "AllGenes"))
    VizEnrichments(one.two.terms, outdir = paste0(base, "AllGenes"))

    one.up.terms <- RunEnrichr(rownames(one.up), 
      outdir = paste0(base, g1, "up"))
    VizEnrichments(one.up.terms, outdir = paste0(base, g1, "up"))

    two.up.terms <- RunEnrichr(rownames(two.up), 
      outdir = paste0(base, g2, "up"))
    VizEnrichments(two.up.terms, outdir = paste0(base, g2, "up"))
  }
}