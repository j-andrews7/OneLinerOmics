#' Run DiffBind analysis on a samplesheet
#'
#' \code{RunDiffBind} performs a high-level differential binding analysis 
#' with \code{DiffBind}. It, along with \link{ProcessDBRs}, form the crux of the
#' ChIP-seq portion of this package. 
#'
#' The default parameters should be an adequate starting place for most users,
#' but lazy folks can provide multiple thresholds to \code{fdr.thresh}
#' and/or \code{fc.thresh} if they aren't sure how stringent or lenient they 
#' need to be with their data.
#'
#' It's generally best to provide an empty directory as the output path, as
#' several directories will be generated. 
#'
#' @param outpath Path to directory to be used for output. Additional 
#'   directories will be generated within this folder.
#' @param samplesheet Path to samplesheet containing sample metadata.
#' @param txdb \code{TxDb} object to use for annotation.
#' @param level String defining variable of interest from \code{samplesheet}. 
#'   Must be one of: "Treatment", "Condition", "Tissue", or "Factor".
#' @param se Path to file containing consensus SEs, which will be used to
#'   to annotate whether individual peaks fall within an SE or not.
#' @param fdr.thresh Number or numeric scalar indicating the false discovery 
#'   rate (FDR) cutoff(s) to be used for determining "significant" differential 
#'   binding. If multiple are given, multiple tables/plots will be generated 
#'   using all combinations of \code{fdr.thresh} and \code{fc.thresh}.
#' @param fc.thresh Number or numeric scalar indicating the log2 fold-change 
#'   cutoff(s) to be used for determining "significant" differential binding.
#'   If multiple are given, multiple tables/plots will be generated using all 
#'   combinations of \code{padj.thresh} and \code{fc.thresh}.
#' @param block String or character vector defining the column(s) in 
#'   \code{samplesheet} to use to block for unwanted variance, e.g. batch or 
#'   technical effects. Must be one of: "Treatment", "Condition", "Tissue", 
#'   or "Factor".
#' @param n.consensus Number of samples in which peaks must overlap for the 
#'   peaks to be merged and included in the consensus peak set. 
#' @param breaks Vector of sequences to be used as breaks for signal heatmaps.
#' @param heatmap.colors Character vector containing custom colors to use for 
#'   heatmaps in hex (e.g. \code{c("#053061", "#f5f5f5", "#67001f")}).
#' @param heatmap.preset String indicating which of the color presets to use in
#'   heatmaps.
#'
#'   Available presets (low to high) are:
#'   \itemize{
#'      \item "BuRd" Blue to red.
#'      \item "OrPu" Orange to purple. 
#'      \item "BrTe" Brown to teal.
#'      \item "PuGr" Purple to green.
#'      \item "BuOr" Sea blue to orange.
#'   }
#' @param reverse Boolean indicating whether to flip heatmap color scheme 
#'   (high color will become low, etc).
#' @param plot.enrich Boolean indicating whether enrichment analyses for DBRs 
#'   should be run and plotted for each comparison. 
#' @param enrich.libs Vector of valid \code{enrichR} libraries to test the 
#'   genes against.
#'   
#'
#'   Available libraries can be viewed with 
#'   \code{\link[enrichR]{listEnrichrDbs}} from the \code{enrichR} package.
#' @param promoters Scalar vector containing how many basepairs up and 
#'   downstream of the TSS should be used to define gene promoters.
#' @param method String indicating method to be used for differential expression 
#'   analysis. Can be "DESeq2" or "edgeR".
#' @return A \code{DBA} object from \code{\link[DiffBind]{dba.analyze}}.
#'
#' @importFrom utils read.table
#' @import DiffBind
#'
#' @export
#'
#' @author Jared Andrews
#'
#' @seealso
#' \code{\link[DiffBind]{dba}}, \code{\link[DiffBind]{dba.count}}, 
#' \code{\link[DiffBind]{dba.contrast}}, \code{\link[DiffBind]{dba.analyze}},
#' \code{\link[DiffBind]{dba.report}} for more about ChIP-seq differential
#' binding analysis.
#'
#' \code{\link{ProcessDBRs}}, for analyzing and visualizing the results.
#'
RunDiffBind <- function(outpath, samplesheet, txdb,
  level = c("Treatment", "Condition", "Tissue", "Factor"), 
  se = NULL, 
  fdr.thresh = 0.05, 
  fc.thresh = 2, 
  block = NULL, 
  heatmap.colors = NULL, 
  heatmap.preset = NULL,
  reverse = FALSE,
  n.consensus = 2, 
  breaks = c(seq(-3, -1, length = 250), seq(-1, -0.1, length = 250),
    seq(-0.1, 0.1,length = 1), seq(0.1, 1, length = 250), 
    seq(1, 3, length = 250)),
  plot.enrich = TRUE, 
  enrich.libs = c("GO_Molecular_Function_2018", 
  "GO_Cellular_Component_2018", "GO_Biological_Process_2018", "KEGG_2019_Human",
  "Reactome_2016", "BioCarta_2016", "Panther_2016"),
  promoters = c(-2000, 2000),
  method = c("DESeq2", "edgeR")) {

  message("### SETUP ###\n")
  rblock = NULL

  # CHECKING ARGUMENTS #
  message("\n# CHECKING ARGUMENTS #\n")
  level <- match.arg(level)
  if (level == "Treatment") {
    rlevel <- DBA_TREATMENT
  } else if (level == "Condition") {
    rlevel <- DBA_CONDITION
  } else if (level == "Tissue") {
    rlevel <- DBA_TISSUE
  } else if (level == "Factor") {
    rlevel <- DBA_FACTOR
  }

  if (!is.null(block)) {
    block <- match.arg(block, choices = c("Treatment", "Condition", "Tissue", 
      "Factor"))

    if (block == "Treatment") {
      rblock <- DBA_TREATMENT
    } else if (block == "Condition") {
      rblock <- DBA_CONDITION
    } else if (block == "Tissue") {
      rblock <- DBA_TISSUE
    } else if (block == "Factor") {
      rblock <- DBA_FACTOR
    }
  }

  method <- match.arg(method)
  if (!is.null(rblock)) {
    if (method == "DESeq2") {
      method <- DBA_DESEQ2_BLOCK
    } else {
      method <- DBA_EDGER_BLOCK
    }
  } else {
    if (method == "DESeq2") {
      method <- DBA_DESEQ2
    } else {
      method <- DBA_EDGER
    }
  }

  # Set up colors/breaks for heatmaps later.
  hmap.colors <- .heatmap_colors(breaks = breaks, preset = heatmap.preset, 
    custom.colors = heatmap.colors, reverse = reverse)
  
  # CREATE DIRECTORY STRUCTURE #
  message("# CREATING DIRECTORY STRUCTURE #\n")
  setup <- CreateOutputStructure(block, level, base = outpath, chip = TRUE)
  base <- setup$base
  
  # LOAD SUPER ENHANCERS #
  if (!is.null(se)) {
    message("# LOADING SUPER ENHANCERS #\n")
    se <- read.table(se, sep = '\t', header = TRUE)
  }  

  # FINDING DIFFERENTIALLY BOUND REGIONS #
  message("### FINDING DIFFERENTIALLY BOUND REGIONS ###\n\n")
  samps <- dba(sampleSheet = samplesheet, minOverlap = n.consensus)
  count <- dba.count(samps, minOverlap = n.consensus)
  cont <- dba.contrast(count, categories = rlevel, block = rblock)
  results <- dba.analyze(cont)
  
  # CONSENSUS PEAKS #
  message("# DEALING WITH CONSENSUS PEAKS #\n")
  
  report <- dba.report(results, th = 1, bCalled = TRUE, 
    bCounts = TRUE, method = method)

  # ANNOTATION #
  message("# ANNOTATING & GENERATING CONSENSUS PLOTS #\n")
  peak.anno <- annotatePeak(report, tssRegion = promoters, TxDb = txdb, 
    annoDb = "org.Hs.eg.db")

  PlotChIPAnnos(peak.anno, outpath = base)
  PlotChIPPCAs(results, outpath = base, method = method)
  PlotChIPHeatmaps(results, outpath = base, method = method, breaks = breaks,
    colors = hmap.colors)
  
  # DB PEAKS #
  message("# ANNOTATING & GENERATING DBR PLOTS #\n")

  ProcessDBRs(results = results, outpath = base, txdb = txdb,
    fdr.thresh = fdr.thresh, fc.thresh = fc.thresh, method = method, 
    breaks = breaks, heatmap.colors = heatmap.colors, 
    heatmap.preset = heatmap.preset, reverse = reverse, 
    plot.enrich = plot.enrich, enrich.libs = enrich.libs)
  
  message("# SAVING RESULTS TABLES #\n")
  SaveResults(results, outpath = base, chip = TRUE, method = method, 
    promoters = promoters, se = se, txdb = txdb)

  message("# SAVING ROBJECTS #")
  saveRDS(results, file=paste0(base, "/Robjects/results.rds"))
  
  return(results)
}


#' Extract, visualize, and save DBRs
#'
#' \code{ProcessDBRs} wraps several plotting functions to generate figures 
#' specifically for differentially bound regions in a given comparison.
#' 
#' \code{ProcessDBRs} is called by \link{RunDiffBind} but can also be re-run 
#' with the \link{RunDiffBind} output if you want to save time and don't need to 
#' generate all of the EDA/consensus figures again.
#'
#' This function will generate many figures in addition to saving the results as
#' tables.
#'
#' @param results \code{DBA} object as returned by 
#'   \code{\link[DiffBind]{dba.analyze}}.
#' @param outpath Path to directory to be used for output.
#' @param txdb \code{TxDb} object to use for annotation.
#' @param fdr.thresh Number or numeric scalar indicating the false discovery 
#'   rate (FDR) cutoff(s) to be used for determining "significant" differential 
#'   binding. If multiple are given, multiple tables/plots will be generated 
#'   using all combinations of \code{fdr.thresh} and \code{fc.thresh}.
#' @param fc.thresh Number or numeric scalar indicating the log2 fold-change 
#'   cutoff(s) to be used for determining "significant" differential binding.
#'   If multiple are given, multiple tables/plots will be generated using all 
#'   combinations of \code{padj.thresh} and \code{fc.thresh}.
#' @param method Method used for differential binding e.g. \code{DBA_DESEQ2}.
#'
#'   Do not quote this parameter, it is a global variable from \code{DiffBind}.
#'   If a block was used, be sure to include that 
#'   (e.g. \code{DBA_DESEQ2_BLOCK}).
#' @param promoters Scalar vector containing how many basepairs up and 
#'   downstream of the TSS should be used to define gene promoters.
#' @param breaks Vector of sequences to be used as breaks for signal heatmaps.
#' @param heatmap.colors Character vector containing custom colors to use for 
#'   heatmaps in hex (e.g. \code{c("#053061", "#f5f5f5", "#67001f")}).
#' @param heatmap.preset String indicating which of the color presets to use in
#'   heatmaps.
#'
#'   Available presets (low to high) are:
#'   \itemize{
#'      \item "BuRd" Blue to red.
#'      \item "OrPu" Orange to purple. 
#'      \item "BrTe" Brown to teal.
#'      \item "PuGr" Purple to green.
#'      \item "BuOr" Sea blue to orange.
#'   }
#' @param reverse Boolean indicating whether to flip heatmap color scheme 
#'   (high color will become low, etc).
#' @param plot.enrich Boolean indicating whether enrichment analyses for DBRs 
#'   should be run and plotted for each comparison.
#' @param enrich.libs A vector of valid \code{enrichR} libraries to test the 
#'   genes against.
#'
#'   Available libraries can be viewed with 
#'   \code{\link[enrichR]{listEnrichrDbs}} from the \code{enrichR} package.
#'
#' @import DiffBind
#' @import parallel
#' @import GenomeInfoDb
#' @importFrom ChIPseeker annotatePeak
#' @importFrom GenomicRanges GRangesList
#'
#' @export
#'
#' @author Jared Andrews
#'
#' @seealso
#' \code{\link{RunDiffBind}}, for generating input for this function.
#'
ProcessDBRs <- function(results, outpath, txdb,
  fdr.thresh = 0.05,
  fc.thresh = 1,
  method = NULL,
  promoters = c(-2000, 2000),
  breaks = c(seq(-3, -1, length = 250), seq(-1, -0.1, length = 250),
    seq(-0.1, 0.1,length = 1), seq(0.1, 1, length = 250), 
    seq(1, 3, length = 250)), 
  heatmap.colors = NULL, 
  heatmap.preset = NULL,
  reverse = FALSE,
  plot.enrich = TRUE, 
  enrich.libs = c("GO_Molecular_Function_2018", 
  "GO_Cellular_Component_2018", "GO_Biological_Process_2018", "KEGG_2019_Human",
  "Reactome_2016", "BioCarta_2016", "Panther_2016")) {

  hmap.colors <- .heatmap_colors(breaks = breaks, preset = heatmap.preset, 
    custom.colors = heatmap.colors, reverse = reverse)

  for (fc in fc.thresh) {
    for (fdr in fdr.thresh) {
      for (i in seq_along(results$contrasts)) {
        reportdb = dba.report(results, th = fdr, fold = fc, bCalled = TRUE, 
          bCounts = TRUE, method = method, contrast = i)
      
        g1up = reportdb[(reportdb$Fold > 0)]
        g2up = reportdb[(reportdb$Fold < 0)]

        g1 <- results$contrasts[[i]]$name1
        g2 <- results$contrasts[[i]]$name2

        # Create a named list of our subsets.
        files <- GRangesList(reportdb, g1up, g2up)
        names(files) <- c(paste0(g1, "-v-", g2), paste0(g1, ".up"), 
          paste0(g2, ".up"))

        peak.anno.list <- lapply(files, annotatePeak, TxDb = txdb,
          tssRegion = promoters, verbose = FALSE, annoDb = "org.Hs.eg.db")
        
        PlotChIPAnnos(peak.anno.list, outpath, consensus = FALSE, comp = , 
          fc = fc, fdr = fdr)

        pdf(paste0(outpath, "/DBRFigures/MAPlots/", g1, "-v-", g2, ".fdr.", 
          fdr, ".log2fc.", fc, "MAplots.pdf"), width = 5, height = 5)
        dba.plotMA(results, report = reportdb, contrast = i, method = method)
        dba.plotMA(results, report = reportdb, contrast = i, method = method,
          bXY = TRUE)
        dev.off()

        pdf(paste0(outpath, "/DBRFigures/SignalBoxPlots/", g1, "-v-", g2, 
          ".fdr.", fdr, ".log2fc.", fc, "MAplots.pdf"), width = 7, height = 5)
        dba.plotBox(results, report = reportdb, contrast = i, method = method,
          bAll = TRUE)
        dev.off()

        pdf(paste0(outpath, "/DBRFigures/VolcanoPlots/", g1, "-v-", g2, 
          ".fdr.", fdr, ".log2fc.", fc, "MAplots.pdf"), width = 5, height = 5)
        dba.plotVolcano(results, th = fdr, fold = fc, contrast = i, 
          method = method)
        dev.off()

        if (plot.enrich) {
          PlotEnrichments(peak.anno.list, outpath, padj.thresh = fdr, 
            fc.thresh = fc, libraries = enrich.libs, chip = TRUE)
        }
      }

      PlotChIPPCAs(results, outpath, method = method, fdr.thresh = fdr, 
        fc.thresh = fc, consensus = FALSE)

      PlotChIPHeatmaps(results, outpath, method = method, fdr.thresh = fdr,
        fc.thresh = fc, breaks = breaks, consensus = FALSE)
    }
  }
  
}


#' @importFrom GenomicRanges makeGRangesFromDataFrame findOverlaps
#' @importFrom S4Vectors queryHits
.categorize_peaks <- function(dfs) {
  # dfs = Named list of two dataframes - first for peaks, second for SEs.
  
  # Make them GRanges objects for easy overlap.
  peaks <- makeGRangesFromDataFrame(dfs$peaks, keep.extra.columns = TRUE)

  peaks$cat <- "PROMOTER"
  peaks$cat[((startsWith(peaks$annotation, "Intron") | 
      startsWith(peaks$annotation, "Downstream") | 
      startsWith(peaks$annotation, "Distal") |
      startsWith(peaks$annotation, "5'") |
      startsWith(peaks$annotation, "3'")))] <- "ENHANCER"
  peaks$cat[((startsWith(peaks$annotation, "Exon")))] <- "EXON"

  if (!is.null(dfs$ses)) {
    ses <- makeGRangesFromDataFrame(dfs$ses, keep.extra.columns = TRUE)

    peaks$inSE <- "FALSE"
    hits <- findOverlaps(peaks, ses)
    peaks$inSE[queryHits(hits)] <- "TRUE"

    # Add a column for the peak category.
    peaks$cat[((startsWith(peaks$annotation, "Promoter")) & 
      peaks$inSE == "TRUE")] <- "PROMOTER.SE"
    peaks$cat[((startsWith(peaks$annotation, "Intron") | 
      startsWith(peaks$annotation, "Downstream") | 
      startsWith(peaks$annotation, "Distal") |
      startsWith(peaks$annotation, "5'") |
      startsWith(peaks$annotation, "3'")) & 
        peaks$inSE == "TRUE")] <- "ENHANCER.SE"
    peaks$cat[((startsWith(peaks$annotation, "Exon")) & 
      peaks$inSE == "TRUE")] <- "EXON.SE"
  }
  
  return(as.data.frame(peaks))
}
