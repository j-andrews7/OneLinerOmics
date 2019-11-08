#' Run DiffBind analysis on a samplesheet
#'
#' \code{RunDiffBind} performs a high-level differential binding analysis 
#' with \code{DiffBind}. It, along with \link{ProcessDBRs}, form the crux of the
#' ChIP-seq portion of this package. 
#'
#' The default parameters should be an adequate starting place for most users,
#' but lazy folks can provide multiple thresholds to \code{padj.thresh}
#' and/or \code{fc.thresh} if they aren't sure how stringent or lenient they 
#' need to be with their data.
#'
#' It's generally best to provide an empty directory as the output path, as
#' several directories will be generated. 
#'
#' 
#'
#' @param outpath Path to directory to be used for output. Additional 
#'   directories will be generated within this folder.
#' @param samplesheet Path to samplesheet containing sample metadata.
#' @param level String defining variable of interest from \code{samplesheet}. 
#'   Must be one of: "Treatment", "Condition", "Tissue", or "Factor".
#' @param se Path to file containing consensus SEs, which will be used to
#'   to annotate whether individual peaks fall within an SE or not.
#' @param padj.thresh Number or numeric scalar indicating the adjusted p-value 
#'   cutoff(s) to be used for determining "significant" differential binding.
#'   If multiple are given, multiple tables/plots will be generated using all 
#'   combinations of \code{padj.thresh} and \code{fc.thresh}.
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
#' @param plot.enrich Boolean indicating whether enrichment analyses for DEGs 
#'   should be run and plotted for each comparison. 
#' @param enrich.libs A vector of valid \code{enrichR} libraries to test the 
#'   genes against.
#'
#' Available libraries can be viewed with \code{\link[enrichR]{listEnrichrDbs}} 
#' from the \code{enrichR} package.
#'
#' @importFrom utils read.table
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
RunDiffBind <- function(outpath, samplesheet, 
  level = c("Treatment", "Condition", "Tissue", "Factor"), 
  se = NULL, padj.thresh = 0.05, fc.thresh = 2, 
  block = NULL, n.consensus = 2, plot.enrich = TRUE, 
  enrich.libs = c("GO_Molecular_Function_2018", 
  "GO_Cellular_Component_2018", "GO_Biological_Process_2018", "KEGG_2019_Human",
  "Reactome_2016", "BioCarta_2016", "Panther_2016")) {

  # Check args, as only certain column names are appropriate.
  level <- match.arg(level)

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

  # Set up colors/breaks for heatmaps later.
  breaks = c(seq(-3, -1, length = 250), seq(-1, -0.1, length = 250),
    seq(-0.1, 0.1,length = 1), seq(0.1, 1, length = 250), 
    seq(1, 3, length = 250))
  if (mark=="H3AC") {
    color = colorRampPalette(c("#b35806", "#e08214", "#f5f5f5", "#8073ac",
      "#542788"))(n = 1000)
  } else {
    color = colorRampPalette(c("#8c510a", "#d8b365", "#f5f5f5", "#5ab4ac",
      "#01665e"))(n = 1000)
  }
  
  ### SET DIRECTORY STRUCTURE ###
  message("### SETTING DIRECTORY STRUCTURE ###\n\n")

  # Determine if there is a blocking factor or not and set design.
  if (!is.null(block)) {
    if (!dir.exists(file.path(outpath,paste(block,"Block", sep='')))) {
      dir.create(file.path(outpath,paste(block,"Block", sep='')))
    }
    if ((process_se==TRUE) & (!dir.exists(file.path(se_outpath, 
      paste(block, "Block", sep=''))))) {
        dir.create(file.path(se_outpath,paste(block,"Block", sep='')))
    }
    base = file.path(outpath, paste(block,"Block",sep=''))
    se_base = file.path(se_outpath, paste(block,"Block",sep=''))
  } else {
    if (!dir.exists(file.path(outpath,"NoBlock"))) {
      dir.create(file.path(outpath,"NoBlock"))
    }
    if ((process_se==TRUE) & (!dir.exists(file.path(se_outpath,
      "NoBlock")))) {
        dir.create(file.path(se_outpath,"NoBlock"))
    }
    base = file.path(outpath, "NoBlock")
    se_base = file.path(se_outpath, "NoBlock")
  }
  
  ### PROCESS SEs ###
  if (process_se==TRUE) {
    message("### PROCESSING SEs ###\n\n")
    se <- ProcessSEs(se_samplesheet, rblock, g1, g2, se_base, breaks, 
      color)
  } else {
    se <- read.delim(paste0(se_base, '/', g1, '-v-', g2, 
      '.H3K27AC.Consensus.SEs.txt'), sep = '\t', header = TRUE)
  }  

  ### FINDING DIFFERENTIALLY BOUND REGIONS ###
  message("### FINDING DIFFERENTIALLY BOUND REGIONS ###\n\n")
  samps = dba(sampleSheet = samplesheet)
  count = dba.count(samps)
  if (!is.null(rblock)) {
    cont = dba.contrast(count, count$masks[[g1]], count$masks[[g2]], g1, g2, 
      block = rblock)
  } else {
    cont = dba.contrast(count, count$masks[[g1]], count$masks[[g2]], g1, g2)
  }
  results = dba.analyze(cont)
  
  # CONSENSUS PEAKS #
  message('# DEALING WITH CONSENSUS PEAKS #\n\n')
  txdb=TxDb.Hsapiens.UCSC.hg19.knownGene
  if (!is.null(rblock)){
    report = dba.report(results, th = 1, bCalled = TRUE, 
      bCalledDetail = TRUE, bCounts = TRUE, method = DBA_DESEQ2_BLOCK)
  } else {
    report = dba.report(results, th = 1, bCalled = TRUE, 
      bCalledDetail = TRUE, bCounts = TRUE, method = DBA_DESEQ2)
  }

  # ANNOTATION #
  message('# ANNOTATING PEAKS & MAKING GENERIC PLOTS #\n\n')
  peakAnno <- annotatePeak(report, tssRegion = c(-2000, 2000),
                       TxDb = txdb, annoDb = "org.Hs.eg.db")
  
  df = data.frame(peakAnno)
  
  # Remove garbage chromosomes.
  df <- df[!grepl("chrM", df$seqnames),]
  df <- df[!grepl("rand", df$seqnames),]
  df <- df[!grepl("_", df$seqnames),]

  pdf(paste(base, '/', g1, '-v-', g2, ".", mark, 
    ".ConsensusPeaks.Figures.pdf", sep = ''))
  if (!is.null(rblock)){
    plotAnnoPie(peakAnno)
    plotAnnoBar(peakAnno)
    vennpie(peakAnno)
    upsetplot(peakAnno)
    plotDistToTSS(peakAnno, title = "Distribution of Peaks Relative to TSS")
    dba.plotPCA(results, th = 1, contrast = 1, sub = paste0(g1, '-v-', g2, 
      " Consensus\n", mark, " Peaks"), method = DBA_DESEQ2_BLOCK)

    dba.plotHeatmap(results, method = DBA_DESEQ2_BLOCK, th = 1, 
      contrast = 1, margin = 15, correlations = FALSE, scale = "row", 
      density.info = "none", colScheme = color, breaks = breaks, 
      main = paste0(g1, '-v-', g2, " Consensus\n", mark, " Top 1000"))
    dba.plotHeatmap(results, method = DBA_DESEQ2_BLOCK, th = 1, contrast = 1, 
      margin = 15, correlations = FALSE, scale = "row", density.info = "none", 
      colScheme = color, breaks = breaks, Colv = NULL, 
      main = paste0(g1, '-v-', g2, " Consensus\n", mark, " Top 1000"))
    dba.plotHeatmap(results, method = DBA_DESEQ2_BLOCK, th = 1, contrast = 1, 
      margin = 15, correlations = FALSE, scale = "row", density.info = "none", 
      colScheme = color, breaks = breaks, maxSites = 10000, 
      main = paste0(g1, '-v-', g2, " Consensus\n", mark, " Top 10000"))
    dba.plotHeatmap(results, method = DBA_DESEQ2_BLOCK, th = 1, contrast = 1, 
      margin = 15, correlations = FALSE, scale = "row", density.info = "none", 
      colScheme = color, breaks = breaks, maxSites = 10000, Colv = NULL, 
      main = paste0(g1, '-v-', g2, " Consensus\n", mark, " Top 10000"))
    dba.plotHeatmap(results, method = DBA_DESEQ2_BLOCK, th = 1, contrast = 1, 
      margin = 15, density.info = "none", 
      main = paste0(g1, "-v-", g1, " Consensus\n", mark, 
        " Peaks - Correlation"))
  } else {
    plotAnnoPie(peakAnno)
    plotAnnoBar(peakAnno)
    vennpie(peakAnno)
    upsetplot(peakAnno)
    plotDistToTSS(peakAnno, title = "Distribution of Peaks Relative to TSS")
    dba.plotPCA(results, th = 1, contrast = 1, 
      sub = paste0(g1, '-v-', g2, " Consensus\n", mark, " Peaks"), 
      method = DBA_DESEQ2)

    dba.plotHeatmap(results, method = DBA_DESEQ2, th = 1, contrast = 1, 
      margin = 15, correlations = FALSE, scale = "row", density.info = "none", 
      colScheme = color, breaks = breaks, 
      main = paste0(g1, '-v-', g2, " Consensus\n", mark, " Top 1000"))
    dba.plotHeatmap(results, method = DBA_DESEQ2, th = 1, contrast = 1, 
      margin = 15, correlations = FALSE, scale = "row", density.info = "none", 
      colScheme = color, breaks = breaks, Colv = NULL, 
      main = paste0(g1, '-v-', g2, " Consensus\n", mark, " Top 1000"))
    dba.plotHeatmap(results, method = DBA_DESEQ2, th = 1, contrast = 1, 
      margin = 15, correlations = FALSE, scale = "row", density.info = "none", 
      colScheme = color, breaks = breaks, maxSites = 10000, 
      main = paste0(g1, '-v-', g2, " Consensus\n", mark, " Top 10000"))
    dba.plotHeatmap(results, method = DBA_DESEQ2, th = 1, contrast = 1, 
      margin = 15, correlations = FALSE, scale = "row", density.info = "none", 
      colScheme = color, breaks = breaks, maxSites = 10000, Colv = NULL, 
      main = paste0(g1, '-v-', g2, " Consensus\n", mark, " Top 10000"))
    dba.plotHeatmap(results, method = DBA_DESEQ2, th = 1, contrast = 1, 
      margin = 15, density.info = "none", 
      main = paste0(g1, "-v-", g1, " Consensus\n", mark, 
        " Peaks - Correlation"))
  }
  dev.off()
  
  # DB PEAKS #
  message('# DEALING WITH DB PEAKS #\n\n')
  # Subset so that we can compare localization and pathways between these groups as well.
  if (!is.null(rblock)){
    reportdb = dba.report(results, th = 0.05, bCalled = TRUE, 
      bCalledDetail = TRUE, bCounts = TRUE, method = DBA_DESEQ2_BLOCK)
  } else {
    reportdb = dba.report(results, th = 0.05, bCalled = TRUE, 
      bCalledDetail = TRUE, bCounts = TRUE, method = DBA_DESEQ2)
  }
  
  g1up = reportdb[(reportdb$Fold > 0)]
  g2up = reportdb[(reportdb$Fold < 0)]

  # Create a named list of our subsets.
  files = GRangesList(reportdb, g1up, g2up)
  names(files) = c("All_DB", paste0(g1, ".up"), paste0(g2, ".up"))
  peakAnnoList <- lapply(files, annotatePeak, TxDb = txdb,
    tssRegion=c(-2000, 2000), verbose = FALSE, annoDb = "org.Hs.eg.db")
  ProcessDBRs(peakAnnoList, reportdb, results, rblock, g1, g2, base, mark, 
    color, breaks)
  
  # Categorize and save all peaks.
  dfs <- list("peaks" = df, "ses" = se)
  x <- CategorizePeaks(dfs)
  
  write.table(x, file = paste0(base, '/', g1, '-v-', g2, ".", mark, 
    ".ConsensusPeaks.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
  
  x.db <- x[x$FDR <= 0.05,]
  write.table(x.db, file = paste0(base, '/', g1, '-v-', g2, ".", mark, 
    ".DBRs.FDR0.05.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
  
  return(x)
}

#' Extract, visualize, and save DBRs
#'
#'
#'
#'
#'
#' @export
#'
#' @author Jared Andrews
#'
ProcessDBRs <- function(peakAnnoList, report, results, rblock, g1, g2, base, 
  mark, color, breaks) {
    # peakAnnoList = A GRangesList of all DB peaks, those up in g1, and those up in g2.
    # report = DiffBind object containing only the DB sites.
    # results = The original DiffBind results object.
    # rblock = Used to determine whether or not blocking is used.
    # g1 = Name of group1.
    # g2 = Name of group2.
    # base = Path to use for saving files.
    # mark = Histone mark - used for file naming and heatmap colors.
    # color = Color scheme to be used for the heatmaps.
    # breaks = Breaks to be used in the heatmaps.

  pdf(paste0(base ,"/", g1, '-v-', g2, ".", mark, ".DBRs.Figures.pdf"))

  if (!is.null(rblock)) {
    plotAnnoBar(peakAnnoList)
    plotDistToTSS(peakAnnoList, 
      title = "Distribution of Regions Relative to TSS")
    dba.plotPCA(results, th = 0.05, contrast = 1, 
      sub = paste0(g1, '-v-', g2, " DBRs\n", mark), method = DBA_DESEQ2_BLOCK)
    dba.plotMA(results, th = 0.05, contrast = 1, method = DBA_DESEQ2_BLOCK)
    dba.plotMA(results, th = 0.05, contrast = 1, bXY = TRUE, 
      method = DBA_DESEQ2_BLOCK)
    dba.plotVolcano(results, th = 0.05, contrast = 1, method = DBA_DESEQ2_BLOCK)
    dba.plotBox(results, th = 0.05, contrast = 1, method = DBA_DESEQ2_BLOCK)

    dba.plotHeatmap(results, method = DBA_DESEQ2_BLOCK, th = 0.05, contrast = 1, 
      margin = 15, correlations = FALSE, scale = "row", density.info = "none", 
      colScheme = color, breaks = breaks, main = paste0(g1, '-v-', g2, 
        " DBRs\n", mark, " Top 1000"))
    dba.plotHeatmap(results, method = DBA_DESEQ2_BLOCK, th = 0.05, contrast = 1, 
      margin = 15, correlations = FALSE, scale = "row", density.info = "none", 
      colScheme = color, breaks = breaks, Colv = NULL, 
      main = paste0(g1, '-v-', g2, " DBRs\n", mark, " Top 1000"))
    dba.plotHeatmap(results, method = DBA_DESEQ2_BLOCK, th = 0.05, contrast = 1, 
      margin = 15, correlations = FALSE, scale = "row", density.info = "none", 
      colScheme = color, breaks = breaks, maxSites = 20000, 
      main = paste0(g1, '-v-', g2, " DBRs\n", mark, " Top 20000"))
    dba.plotHeatmap(results, method = DBA_DESEQ2_BLOCK, th = 0.05, contrast = 1, 
      margin = 15, correlations = FALSE, scale = "row", density.info = "none", 
      colScheme = color, breaks = breaks, maxSites = 20000, Colv = NULL, 
      main = paste0(g1, '-v-', g2, " DBRs\n", mark, " Top 20000"))
    dba.plotHeatmap(results, method = DBA_DESEQ2_BLOCK, th = 0.05, contrast = 1, 
      margin = 15, density.info = "none", 
      main = paste0(g1, "-v-", g1, " DBRs\n", mark, " - Correlation"))
  } else {
    plotAnnoBar(peakAnnoList)
    plotDistToTSS(peakAnnoList, 
      title = "Distribution of Regions Relative to TSS")
    dba.plotPCA(results, th = 0.05, contrast = 1, 
      sub = paste0(g1, '-v-', g2, " DBRs\n", mark), method = DBA_DESEQ2)
    dba.plotMA(results, th = 0.05, contrast = 1, method = DBA_DESEQ2)
    dba.plotMA(results, th = 0.05, contrast = 1, bXY = TRUE, 
      method = DBA_DESEQ2)
    dba.plotVolcano(results, th = 0.05, contrast = 1, method = DBA_DESEQ2)
    dba.plotBox(results, th = 0.05, contrast = 1, method = DBA_DESEQ2)

    dba.plotHeatmap(results, method = DBA_DESEQ2, th = 0.05, contrast = 1, 
      margin = 15, correlations = FALSE, scale = "row", density.info = "none", 
      colScheme = color, breaks = breaks, 
      main = paste0(g1, '-v-', g2, " DBRs\n", mark, " Top 1000"))
    dba.plotHeatmap(results, method = DBA_DESEQ2, th = 0.05, contrast = 1, 
      margin = 15, correlations = FALSE, scale = "row", density.info = "none", 
      colScheme = color, breaks = breaks, Colv = NULL, 
      main = paste0(g1, '-v-', g2, " DBRs\n", mark, " Top 1000"))
    dba.plotHeatmap(results, method = DBA_DESEQ2, th = 0.05, contrast = 1, 
      margin = 15, correlations = FALSE, scale = "row", density.info = "none", 
      colScheme = color, breaks = breaks, maxSites = 20000, 
      main = paste0(g1, '-v-', g2, " DBRs\n", mark, " Top 20000"))
    dba.plotHeatmap(results, method = DBA_DESEQ2, th = 0.05, contrast = 1, 
      margin = 15, correlations = FALSE, scale = "row", density.info = "none", 
      colScheme = color, breaks = breaks, maxSites = 20000, Colv = NULL, 
      main = paste0(g1, '-v-', g2, " DBRs\n", mark, " Top 20000"))
    dba.plotHeatmap(results, method = DBA_DESEQ2, th = 0.05, contrast = 1, 
      margin = 15, density.info = "none", 
      main = paste0(g1, "-v-", g1, " DBRs\n", mark, " - Correlation"))
  }
  dev.off()
  
  RunEnrichments(peakAnnoList, g1, g2, base)
}


ProcessSEs <- function(samplesheet, rblock, g1, g2, base, breaks, color) {
  # samplesheet = File path to SE samplesheet.
  # rblock = Blocking factor.
  # g1 = Name of group 1.
  # g2 = Name of group 2.
  # base = File path to output folder.
  # breaks = Breaks to be used for heatmaps.
  # color = Color scheme to be used for heatmaps.
  
  samps = dba(sampleSheet = samplesheet)
  
  count = dba.count(samps)
  if (!is.null(rblock)) {
    cont = dba.contrast(count, count$masks[[g1]], count$masks[[g2]], g1, g2, 
      block=rblock, minMembers = 2)
  } else {
    cont = dba.contrast(count, count$masks[[g1]], count$masks[[g2]], g1, g2, 
      minMembers = 2)
  }
  results = dba.analyze(cont)
  
  # CONSENSUS PEAKS #
  message('\n# DEALING WITH CONSENSUS SEs #\n\n')
  txdb = TxDb.Hsapiens.UCSC.hg19.knownGene
  if (!is.null(rblock)) {
    report = dba.report(results, th = 1, bCalled = TRUE, bCalledDetail = TRUE, 
      bCounts = TRUE, method = DBA_DESEQ2_BLOCK)
  } else {
    report = dba.report(results, th = 1, bCalled = TRUE, bCalledDetail = TRUE, 
      bCounts = TRUE, method = DBA_DESEQ2)
  }

  # ANNOTATION #
  message('\n# ANNOTATING SEs & MAKING GENERIC PLOTS #\n\n')
  peakAnno <- annotatePeak(report, tssRegion = c(-2000, 2000),
                       TxDb = txdb, annoDb = "org.Hs.eg.db")
  
  df = data.frame(peakAnno)
  
  # Clean out the garbage chromosomes.
  df <- df[!grepl("chrM", df$seqnames), ]
  df <- df[!grepl("rand", df$seqnames), ]
  df <- df[!grepl("_", df$seqnames), ]
  
  write.table(df, row.names = FALSE, quote = FALSE, sep = "\t", 
    file = paste0(base, "/", g1, '-v-', g2, ".", mark, ".Consensus.SEs.txt"))

  pdf(paste0(base, "/", g1, '-v-', g2, ".", mark, ".Consensus.SEs.Figures.pdf"))

  if (!is.null(rblock)) {
    plotAnnoPie(peakAnno)
    plotAnnoBar(peakAnno)
    vennpie(peakAnno)
    upsetplot(peakAnno)
    plotDistToTSS(peakAnno, title = "Distribution of SEs Relative to TSS")
    dba.plotPCA(results, report = report, 
      sub = paste0(g1, '-v-', g2, " Consensus\n", mark, " SEs"), 
      method = DBA_DESEQ2_BLOCK)
    
    dba.plotHeatmap(results, method = DBA_DESEQ2_BLOCK, report = report, 
      margin = 15, correlations = FALSE, scale = "row", density.info = "none", 
      colScheme = color, breaks = breaks, 
      main = paste0(g1, '-v-', g2, " Consensus\n", mark, " Top 1000"))
    dba.plotHeatmap(results, method = DBA_DESEQ2_BLOCK, report = report, 
      margin = 15, correlations = FALSE, scale = "row", density.info = "none", 
      colScheme = color, breaks = breaks, Colv = NULL, 
      main = paste0(g1, '-v-', g2, " Consensus\n", mark, " Top 1000"))
    dba.plotHeatmap(results, method = DBA_DESEQ2_BLOCK, report = report, 
      margin = 15, density.info = "none", 
      main = paste0(g1, "-v-", g1, " Consensus\n", mark, " SEs - Correlation"))
  } else {
    plotAnnoPie(peakAnno)
    plotAnnoBar(peakAnno)
    vennpie(peakAnno)
    upsetplot(peakAnno)
    plotDistToTSS(peakAnno, title = "Distribution of SEs Relative to TSS")
    dba.plotPCA(results, report = report, 
      sub = paste0(g1, '-v-', g2, " Consensus\n", mark, " SEs"), 
      method = DBA_DESEQ2)
    
    dba.plotHeatmap(results, method = DBA_DESEQ2, report = report, margin = 15, 
      correlations = FALSE, scale = "row", density.info = "none", 
      colScheme = color, breaks = breaks, 
      main = paste0(g1, '-v-', g2, " Consensus\n", mark, " Top 1000"))
    dba.plotHeatmap(results, method = DBA_DESEQ2, report = report, margin = 15, 
      correlations = FALSE, scale = "row", density.info = "none", 
      colScheme = color, breaks = breaks, Colv = NULL, 
      main = paste0(g1, '-v-', g2, " Consensus\n", mark, " Top 1000"))
    dba.plotHeatmap(results, method = DBA_DESEQ2, report = report, margin = 15, 
      density.info = "none", 
      main = paste0(g1, "-v-", g1, " Consensus\n", mark, " SEs - Correlation"))
  }
  dev.off()
  
  return(df)
}


CategorizePeaks <- function(dfs, groups=FALSE) {
  # dfs = Named list of two dataframes - first for peaks, second for SEs.
  
  # Make them GRanges objects for easy overlap.
  peaks <- makeGRangesFromDataFrame(dfs$peaks, keep.extra.columns = TRUE)
  ses <- makeGRangesFromDataFrame(dfs$ses, keep.extra.columns = TRUE)
  
  peaks$inSE <- "FALSE"
  hits <- findOverlaps(peaks, ses)
  peaks$inSE[queryHits(hits)] <- "TRUE"
  
  # Add a column for the peak category.
  peaks$cat <- "PROMOTER"
  peaks$cat[((startsWith(peaks$annotation, "Promoter")) & 
    peaks$inSE=="TRUE")] <- "PROM.SE"
  peaks$cat[((startsWith(peaks$annotation, "Intron") | 
    startsWith(peaks$annotation, "Downstream") | 
    startsWith(peaks$annotation, "Distal") ) & peaks$inSE=="TRUE") ] <- "CE.SE"
  peaks$cat[((startsWith(peaks$annotation, "Intron") | 
    startsWith(peaks$annotation, "Downstream") | 
    startsWith(peaks$annotation, "Distal") ) & peaks$inSE=="FALSE") ] <- "CE"
  peaks$cat[((startsWith(peaks$annotation, "Exon")) & 
    peaks$inSE=="FALSE") ] <- "EXON"
  peaks$cat[((startsWith(peaks$annotation, "Exon")) & 
    peaks$inSE=="TRUE") ] <- "EXON.SE"
  
  
  # Rearrange some columns for easy plotting later.
  final <- as.data.frame(peaks)
  if (groups == FALSE){
    nsamps <- (NCOL(final) - 27)/2
    final <- subset(final, select=c(1:13,(14+nsamps):NCOL(final), 
      14:(NCOL(final)-(14+nsamps))))
  }
  
  return(final)
}


RunEnrichments <- function(peakAnnoList, g1, g2, outpath) {
  message("# GO ENRICHMENT #\n\n")
  genes = lapply(peakAnnoList, function(i) as.data.frame(i)$geneId)
  genes2 = lapply(peakAnnoList[2:3], function(i) as.data.frame(i)$geneId)
}

