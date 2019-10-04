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
DoEnrichments <- function(samp1, samp2, results, baseout, padj.thresh = 0.01, 
  lfc.thresh = 2, top.n = NULL) {

  base <- paste0(baseout, samp1, "v", samp2, "/padj.01.lfc2.read100/")
  dir.create(file.path(base, "AllGenes"), showWarnings = FALSE, 
    recursive = TRUE)
  dir.create(file.path(base, paste0(samp1, "up")), showWarnings = FALSE, 
    recursive = TRUE)
  dir.create(file.path(base, paste0(samp2, "up")), showWarnings = FALSE, 
    recursive = TRUE)

  One.Two <- subset(results, 
    padj <= padj.thresh & abs(log2FoldChange) >= lfc.thresh)
  One.up <- subset(results, 
    padj <= padj.thresh & log2FoldChange >= lfc.thresh)
  Two.up <- subset(results, 
    padj <= padj.thresh & log2FoldChange <= -lfc.thresh)
  
  One.Two.terms <- RunEnrichr(rownames(One.Two), 
    outdir = paste0(base, "AllGenes"))
  VizEnrichments(One.Two.terms, outdir = paste0(base, "AllGenes"))

  One.up.terms <- RunEnrichr(rownames(One.up), 
    outdir = paste0(base, samp1, "up"))
  VizEnrichments(One.up.terms, outdir = paste0(base, samp1, "up"))

  Two.up.terms <- RunEnrichr(rownames(Two.up), 
    outdir = paste0(base, samp2, "up"))
  VizEnrichments(Two.up.terms, outdir = paste0(base, samp2, "up"))
  
  # If top hits set, do it for them as well.
  if (!is.null(top.n)) {
    dir.create(file.path(base, "AllGenes", paste0("Top", top.n)), 
      showWarnings = FALSE, recursive = TRUE)
    dir.create(file.path(base, paste0(samp1, "up/", "Top", top.n)), 
      showWarnings = FALSE, recursive = TRUE)
    dir.create(file.path(base, paste0(samp2, "up/", "Top", top.n)), 
      showWarnings = FALSE, recursive = TRUE)
      
    # Top 100 for each category.
    One.Two.top.n <- One.Two[order(-(abs(One.Two$log2FoldChange))),]
    One.Two.top.n <- One.Two.top.n[1:top.n,]
    
    One.Two.top.n.terms <- RunEnrichr(rownames(One.Two.top.n), 
      outdir = paste0(base, "AllGenes/Top", top.n))
    VizEnrichments(One.Two.top.n.terms, outdir = paste0(base, "AllGenes/Top", 
      top.n))

    One.up.top.n <- One.Two[order(-One.Two$log2FoldChange),]
    One.up.top.n <- One.up.top.n[1:top.n,]
    
    One.up.top.n.terms <- RunEnrichr(rownames(One.up.top.n), 
      outdir = paste0(base, samp1, "up/Top", top.n))
    VizEnrichments(One.up.top.n.terms, outdir = paste0(base, samp1, "up/Top", 
      top.n))

    Two.up.top.n <- One.Two[order(One.Two$log2FoldChange),]
    Two.up.top.n <- Two.up.top.n[1:top.n,]
    
    Two.up.top.n.terms <- RunEnrichr(rownames(Two.up.top.n), 
      outdir = paste0(base, samp2, "up/Top", top.n))
    VizEnrichments(Two.up.top.n.terms, outdir = paste0(base, samp2, "up/Top", 
      top.n))
  }
}