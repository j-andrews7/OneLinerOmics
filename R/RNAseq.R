#' Run DESeq2 analyses
#'
#' This function performs a high-level differential gene expression analysis 
#' with \code{DESeq2}.
#'
#'
#'
#' @param outpath Path to directory to be used for output. Additional 
#'   directories will be generated within this folder.
#' @param quants.path Path to directory containing a directory for each sample
#'   with a salmon-generated \code{quant.sf} file inside.
#' @param samplesheet Path to samplesheet containing sample metadata.
#' @param tx2gene Path to file with transcript IDs in first column and gene 
#'   identifiers in second. Used for gene-level summarization.
#' @param level String defining variable of interest.
#' @param padj.thresh Number or numeric scalar indicating the adjusted p-value 
#'   cutoff(s) to be used for determining "significant" differential expression.
#'   If multiple are given, multiple tables/plots will be generated using all 
#'   combinations of \code{padj.thresh} and \code{fc.thresh}.
#' @param fc.thresh Number or numeric scalar indicating the log2 fold-change 
#'   cutoff(s) to be used for determining "significant" differential expression.
#'   If multiple are given, multiple tables/plots will be generated using all 
#'   combinations of \code{padj.thresh} and \code{fc.thresh}.
#' @param plot.annos String or character vector defining the column(s) in 
#'   \code{samplesheet} to use to annotate figures.
#' @param plot.box Boolean indicating whether box plots for DEGs should be 
#'   created for each comparison. If so, the \code{top.n} genes will be plotted.
#'   This step is quite time-consuming with many genes.
#' @param top.n Number of differentially expressed genes to create boxplots for, 
#'   ranked by adj. p-value after applying \code{padj.thresh} and 
#'   \code{fc.thresh} thresholds. If multiple thresholds are provided, the 
#'   lowest fold-change and highest adj. p-value thresholds will be used. 
#' @param block String or character vector defining the column(s) in 
#'   \code{samplesheet} to use to block for unwanted variance, e.g. batch or 
#'   technical effects.
#' @param count.filt Number indicating read threshold. Genes with fewer counts 
#'   than this number summed across all samples will be removed from the
#'   analysis.
#' @return Named List containing a \linkS4class{DESeqDataSet} object from 
#'   running \code{\link[DESeq2]{DESeq}}, a 
#'   \linkS4class{RangedSummarizedExperiment} object from running 
#'   \code{\link[DESeq2]{rlog}}, and a \linkS4class{RangedSummarizedExperiment} 
#'   object from running \code{\link[DESeq2]{vst}}.
#'   
#' @import DESeq2
#' @import ggplot2
#' @importFrom pheatmap pheatmap
#' @importFrom tximport tximport
#' @importFrom utils read.table write.table
#'
#' @export
#'
#' @author Jared Andrews
#'
#' @seealso
#' \code{\link[DESeq2]{DESeq}}, for more about differential expression analysis.
#' \code{\link{ProcessDEGs}}, for analyzing and visualizing the results.
#'
RunDESeq2 <- function(outpath, quants.path, samplesheet, tx2gene, level,  
  padj.thresh = 0.05, fc.thresh = 2, plot.annos = NULL, plot.box = TRUE, 
  top.n = 100, block = NULL, count.filt = 10) {
    
  message("### EXPLORATORY DATA ANALYSIS ###\n")
  message("# SET DIRECTORY STRUCTURE AND MODEL DESIGN #\n")
  # Create directory structure and set design formula.
  base <- outpath
  setup <- CreateOutputStructure(block, level, base)
  base <- setup$base
  design <- setup$design
  
  if (is.null(plot.annos)) {
      plot.annos <- level
  }
  
  ### FILE LOADING & PRE-FILTERING ###
  message("# FILE LOADING & PRE-FILTERING #\n")
  tx2gene <- read.table(tx2gene, sep = "\t")    
  samples <- read.table(samplesheet, header = TRUE)
  rownames(samples) <- samples$name
  files <- file.path(quants.path, samples$name, "quant.sf")
  names(files) <- samples$name

  # Read in our actual count files now.
  txi <- tximport(files, type = "salmon", tx2gene = tx2gene)

  message(paste0("\nDesign is: ", design, "\n"))
  # Create the DESeqDataSet object.
  dds <- DESeqDataSetFromTximport(txi, colData = samples, design = design)

  # Pre-filter transcripts with really low read counts.
  message(paste0("\ndds starting with", nrow(dds), "genes."))
  dds <- dds[rowSums(counts(dds)) >= count.filt,]
  message(paste0("\ndds has", nrow(dds), 
    "genes after removing genes with under ", count.filt, " counts total.\n"))

  ### VARIANCE STABILIZATION COMPARISONS ###
  message("\n# VARIANCE STABILIZATION COMPARISONS #\n")
  vst.out <- paste0(base,"/EDAFigures/VarianceTransformations.pdf")
  message(vst.out)
  trans <- PlotVarianceTransformations(dds, vst.out)
  rld <- trans$rld
  vsd <- trans$vsd

  ### SAMPLE DISTANCES ###
  message("\n# PLOTTING SAMPLE DISTANCES #\n")
  dists.out <- paste0(base, "/EDAFigures/SampleDistances.pdf")
  message(dists.out)
  PlotSampleDistances(rld, vsd, dists.out, level, plot.annos)

  ### PCA PLOTS ###
  message("\n# PCA PLOTS #\n")
  pca.out <- paste0(base, "/EDAFigures/pca.pdf")
  message(pca.out)
  PlotEDAPCAs(rld, vsd, pca.out, level, plot.annos)
  
  #======================================#
  ### DIFFERENTIAL EXPRESSION ANALYSIS ###
  message("\n### DIFFERENTIAL EXPRESSION ANALYSIS ###\n")
  dds <- DESeq(dds)

  message("\n# SAVING ROBJECTS #\n")
  message(paste0("DESeq2: ", paste0(base, "/Robjects/dds.rds")))
  message(paste0("Regularized log transformation: ", 
    paste0(base, "/Robjects/rld.rds")))
  message(paste0("Variance stabilized transformation: ", 
    paste0(base, "/Robjects/vsd.rds")))
  saveRDS(dds, file=paste0(base, "/Robjects/dds.rds"))
  saveRDS(rld, file=paste0(base, "/Robjects/rld.rds"))
  saveRDS(vsd, file=paste0(base, "/Robjects/vsd.rds"))

  return(list(dds = dds, rld = rld, vsd = vsd))
}


#' Extract, visualize, and save DEG lists
#'
#' 
#' @param dds A \linkS4class{DESeqDataSet} object as returned by 
#'   \link{RunDESeq2}.
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
#' @param padj.thresh Number or numeric scalar indicating the adjusted p-value 
#'   cutoff(s) to be used for determining "significant" differential expression.
#'   If multiple are given, multiple tables/plots will be generated using all 
#'   combinations of \code{padj.thresh} and \code{fc.thresh}.
#' @param fc.thresh Number or numeric scalar indicating the log2 fold-change 
#'   cutoff(s) to be used for determining "significant" differential expression.
#'   If multiple are given, multiple tables/plots will be generated using all 
#'   combinations of \code{padj.thresh} and \code{fc.thresh}.
#' @param plot.box Boolean indicating whether box plots for DEGs should be 
#'   created for each comparison. If so, the \code{top.n} genes will be plotted.
#' @param top.n Number of differentially expressed genes to create boxplots for, 
#'   ranked by adj. p-value after applying \code{padj.thresh} and 
#'   \code{fc.thresh} thresholds. If multiple thresholds are provided, the 
#'   lowest fold-change and highest adj. p-value thresholds will be used. 
#'
#' @importFrom DESeq2 lfcShrink counts plotCounts
#' @importFrom magrittr %>%
#' @importFrom htmlwidgets saveWidget
#'
#' @author Jared Andrews
#'
ProcessDEGs <- function(dds, rld, vsd, outpath, level, plot.annos, 
  padj.thresh = 0.05, fc.thresh = 2, plot.box = TRUE, top.n = 100) {

  message("\n# COLLECTING RESULTS #\n")
  # Get all possible sample comparisons.
  combs <- combn(colData(rld)[,level], 2)
  combs.seq <- seq(1, length(combs), by = 2)

  res.list <- list()
  for (samp in combs.seq) {
    g1 <- combs[samp]
    g2 <- combs[samp + 1]
    res <- lfcShrink(dds, contrast=c(level, g1, g2), type='ashr')
    res.list[[paste0(g1, "-v-", g2)]] <- res
  }
  
  ##### PLOTTING #####

  ### DEG BOXPLOTS, PCAS, VOLCANOES, HEATMAPS ###
  message("\n# GENERATING PLOTS #\n")
  for (p in padj.thresh) {
    for (fc in fc.thresh) {
      if (plot.box) {
        PlotBoxplots(res.list, dds, rld, outpath, p, fc, top.n, level)
      }
      PlotDEGPCAs(res.list, rld, vsd, outpath, level, plot.annos, p, fc)
      PlotVolcanoes(res.list, dds, outpath, p, fc)
      PlotHeatmaps(res.list, rld, vsd, outpath, p, fc, plot.annos)
      PlotCombinedHeatmaps(res.list, rld, vsd, outpath, p, fc, plot.annos)
      PlotEnrichments(res.list, outpath, p, fc)
    }
  }
  
  ### MA PLOTs ###
  for (r in seq_along(res.list)) {
    res <- res.list[r]
    comp <- names(res.list[r])
    message("Comparison set: ", comp)
    for (p in padj.thresh) {
      message("Generating MA plots with p.adj cutoff of ", p)
      pdf(paste0(outpath,"/MAPlots/", comp, ".padj.", p, ".MA_plot.pdf"))
      plotMA(res, ylim = c(-6, 6), alpha = p)
      dev.off()
    }
  }

  ### SAVING TABLES ###
  message("\n# SAVING RESULTS TABLES #\n")

  write.table(res, file=paste0(outpath, "/", g1, "-v-", g2, ".ALL.csv"), 
    row.names = FALSE, sep = "\t")

  # Make table for just DEGs too.
  write.csv(res, file=paste0(outpath, "/", g1, "-v-", g2, ".DEGs.padj0.1.csv"), 
    row.names = FALSE, sep = "\t")
  
  # This fpm version allows for accurate intergene comparisons of counts.
  write.csv(fpm(dds), file=paste0(outpath, "/", g1, "-v-", g2,".ALL.fpm.csv"), 
    row.names = FALSE, sep = "\t")

  return(list(res.list = res.list, dds = dds, rld = rld, vsd = vsd))
}