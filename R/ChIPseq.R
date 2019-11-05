#' Run ChIPQC on a sample sheet
#'
#' This function simply runs \code{\link[ChIPQC]{ChIPQC}} on a sample sheet and
#' returns QC metrics. It generates a consensus peakset between all sets before
#' defining the fraction of reads in peaks and other metrics.
#'
#' @param outpath Path to directory to be used for output.
#' @param samplesheet Path to samplesheet containing sample metadata.
#' @param chromosomes String or character vector indicating chromosomes to be 
#'   used for QC. 
#' @return Dataframe containing QC metrics for each sample, including 
#'   cross-correlation scores, fraction of reads in peaks, fragment length, etc.
#'
#' @importFrom ChIPQC ChIPQC QCmetrics
#' @importFrom BiocParallel SerialParam register
#' @importFrom utils write.table
#' 
#' @export
#'
#' @author Jared Andrews
#'
#' @seealso \code{\link[ChIPQC]{ChIPQC}}
#'
RunChIPQC <- function(outpath, samplesheet, chromosomes = "chr10") {

  register(SerialParam())
  exp <- suppressMessages(suppressWarnings(ChIPQC(samplesheet, 
    chromosomes = chromosomes, consensus = TRUE, bCount = FALSE)))
  metrics <- as.data.frame(QCmetrics(exp))
  
  write.table(metrics, file = paste0(outpath, "/QCmetrics.txt"), quote = FALSE,
    sep = "\t")
  
  return(metrics)
}


RunDiffBind <- function(base, samplesheet, se_samplesheet = NULL, g1, g2, mark, 
  block = NULL) {
  # base = Base name of file path.
  # samplesheet = File path to sample sheet.
  # se_samplesheet = File path to the super enhancer sample sheet.
  # g1 & g2 = Strings for groups. Used in file names and paths.
  # mark = Histone mark. Again, file names, colors, and all. "H3AC" or "H3K27AC".
  # block = One of "TREATMENT","CONDITION", or "TISSUE" to account for confounding factors.

  outpath = paste0(base, mark)
  se_outpath = paste0(base, "SEs")
  process_se = FALSE
  se = NULL

  if (mark == "H3K27AC" & !is.null(se_samplesheet)) {
    process_se = TRUE
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

  # Set block properly.
  if (!is.null(block)) {
    if (block=="TREATMENT") {
      rblock=DBA_TREATMENT
    } else if (block=="CONDITION") {
      rblock=DBA_CONDITION
    } else if (block=="TISSUE") {
      rblock=DBA_TISSUE
    } else if (block=="FACTOR") {
      rblock=DBA_FACTOR
    } else {
      stop("Block set to improper value, must be one of 'TREATMENT', 'CONDITION', 'FACTOR' or 'TISSUE'\n\n")
    }
  } else {
    rblock=NULL
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
      block = rblock, minMembers = 2)
  } else {
    cont = dba.contrast(count, count$masks[[g1]], count$masks[[g2]], g1, g2, 
      minMembers = 2)
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
  
  # COMPARING SE # AND WIDTH BETWEEN GROUPS #
  message('\n# COMPARING SE # AND WIDTH BETWEEN GROUPS #\n\n')
  
  pdf(paste0(base, "/SE-Width.Height.Comparisons.pdf"))
  compare.peak_widths(g1, g2, samps, base, se=NULL, process_se=TRUE)
  compare.peak_heights(g1, g2, df, base, process_se=TRUE)
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
  peaks$cat <- "PROM"
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

  ### DB Peaks/SEs ###
  ego.cc.clusters = NULL
  ego.bp.clusters = NULL
  ego.mf.clusters = NULL
  ego.all.clusters = NULL
  
  ego.cc.clusters2 = NULL
  ego.bp.clusters2 = NULL
  ego.mf.clusters2 = NULL
  ego.all.clusters2 = NULL

  # Remove redundant terms
  try({
    ego.bp.clusters <- compareCluster(geneCluster = genes, fun="enrichGO", 
      ont="BP", OrgDb = org.Hs.eg.db, pAdjustMethod = "BH", 
      pvalueCutoff = 0.05, qvalueCutoff = 0.2);
    ego.bp.clusters <- simplify(ego.bp.clusters, cutoff = 0.7, by = "p.adjust", 
      select_fun=min)
  })
  
  try({
    ego.mf.clusters <- compareCluster(geneCluster = genes, fun = "enrichGO", 
      ont = "MF", OrgDb = org.Hs.eg.db, pAdjustMethod = "BH", 
      pvalueCutoff = 0.05, qvalueCutoff = 0.2);
    ego.mf.clusters <- simplify(ego.mf.clusters, cutoff = 0.7, by = "p.adjust", 
      select_fun = min)
  })
  
  try({
    ego.cc.clusters <- compareCluster(geneCluster = genes, fun = "enrichGO", 
      ont = "CC", OrgDb = org.Hs.eg.db,
      pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.2);
    ego.cc.clusters <- simplify(ego.cc.clusters, cutoff = 0.7, by = "p.adjust", 
      select_fun=min)
  })

  try(
    ego.all.clusters <- compareCluster(geneCluster = genes, fun="enrichGO", 
      ont = "ALL", OrgDb = org.Hs.eg.db, pAdjustMethod = "BH", 
      pvalueCutoff = 0.05, qvalueCutoff = 0.2, pool=TRUE)
  )
  
  # Now with only the up/down sets.
  try({
    ego.bp.clusters2 <- compareCluster(geneCluster = genes2, fun = "enrichGO", 
      ont = "BP", OrgDb = org.Hs.eg.db, pAdjustMethod = "BH", 
      pvalueCutoff = 0.05, qvalueCutoff = 0.2);
    ego.bp.clusters2 <- simplify(ego.bp.clusters2, cutoff = 0.7, 
      by = "p.adjust", select_fun = min)
  })
  
  try({
    ego.mf.clusters2 <- compareCluster(geneCluster = genes2, fun = "enrichGO",
      ont = "MF", OrgDb = org.Hs.eg.db, pAdjustMethod = "BH", 
      pvalueCutoff = 0.05, qvalueCutoff = 0.2);
    ego.mf.clusters2 <- simplify(ego.mf.clusters2, cutoff = 0.7, 
      by = "p.adjust", select_fun = min)
  })
  
  try({
    ego.cc.clusters2 <- compareCluster(geneCluster = genes2, fun = "enrichGO",
      ont = "CC", OrgDb = org.Hs.eg.db, pAdjustMethod = "BH", 
      pvalueCutoff = 0.05, qvalueCutoff = 0.2);
    ego.cc.clusters2 <- simplify(ego.cc.clusters2, cutoff = 0.7, 
      by = "p.adjust", select_fun = min)
  })

  try(
    ego.all.clusters2 <- compareCluster(geneCluster = genes2, fun = "enrichGO",
     ont="ALL", OrgDb = org.Hs.eg.db, pAdjustMethod = "BH", 
     pvalueCutoff = 0.05, qvalueCutoff = 0.2, pool=TRUE)
  )

  pdf(paste(outpath,"/GO_Enrichments.DB.pdf", sep=''), height=12, width=12)

  if (!is.null(ego.all.clusters)) {
    write.table(ego.all.clusters, 
      file = paste0(outpath,"/GO_Enrichments.DB.txt"), quote = FALSE, 
      sep = "\t", row.names = FALSE)
  }
  if (!(is.null(ego.bp.clusters))) {
    p = dotplot(ego.bp.clusters, showCategory = 40, 
      title = "GO BP (padj < 0.05) - Top 40\nDBR (FDR <= 0.05)")
    plot(p)
  }
  if (!(is.null(ego.mf.clusters))) {
    p = dotplot(ego.mf.clusters, showCategory = 40, 
      title = "GO MF (padj < 0.05) - Top 40\nDBR (FDR <= 0.05)")
    plot(p)
  }
  if (!(is.null(ego.cc.clusters))) {
    p = dotplot(ego.cc.clusters, showCategory = 40, 
      title = "GO CC (padj < 0.05) - Top 40\nDBR (FDR <= 0.05)")
    plot(p)
  }
  if (!(is.null(ego.all.clusters))) {
    p = dotplot(ego.all.clusters, showCategory=40, 
      title = "GO BP/MF/CC (padj < 0.05) - Top 40\nDBR (FDR <= 0.05)")
    plot(p)
  }
                  
  if (!(is.null(ego.bp.clusters2))) {
    p = dotplot(ego.bp.clusters2, showCategory=40, 
      title = "GO BP (padj < 0.05) - Top 40\nDBR (FDR <= 0.05)")
    plot(p)
  }
  if (!(is.null(ego.mf.clusters2))) {
    p = dotplot(ego.mf.clusters2, showCategory=40, 
      title = "GO MF (padj < 0.05) - Top 40\nDBR (FDR <= 0.05)")
    plot(p)
  }
  if (!(is.null(ego.cc.clusters2))) {
    p = dotplot(ego.cc.clusters2, showCategory=40, 
      title = "GO CC (padj < 0.05) - Top 40\nDBR (FDR <= 0.05)")
    plot(p)
  }
  if (!(is.null(ego.all.clusters2))) {
    p = dotplot(ego.all.clusters2, showCategory=40, 
      title = "GO BP/MF/CC (padj < 0.05) - Top 40\nDBR (FDR <= 0.05)")
    plot(p)
  }
  dev.off()

  ### KEGG Enrichment ###
  kk.clusters = NULL
  kk.clusters2 = NULL
  message("\n\n# PATHWAY ENRICHMENTS #\n\n")
  pdf(paste0(outpath,"/Pathway_Enrichments.DB.pdf"), height = 12, width = 12)

  try(kk.clusters <- compareCluster(geneCluster = genes, fun = "enrichKEGG"))
  if (!(is.null(kk.clusters))) {
    write.table(kk.clusters, file = paste0(outpath,"/KEGG_Enrichments.DB.txt"), 
      quote = FALSE, sep = "\t", row.names = FALSE)
  }
  try(kk.clusters2 <- compareCluster(geneCluster = genes2, fun = "enrichKEGG"))
      
  if (!(is.null(kk.clusters))) {
    p = dotplot(kk.clusters, showCategory = 40, 
      title = "KEGG Pathways (padj < 0.05) - Top 40\nDBR (FDR <= 0.05)")
    plot(p)
  }
  if (!(is.null(kk.clusters2))) {
    p = dotplot(kk.clusters2, showCategory = 40, 
      title = "KEGG Pathways (padj < 0.05) - Top 40\nDBR (FDR <= 0.05)")
    plot(p)
  }

  # KEGG Modules - Sometimes easier to interpret #
  kk.clusters = NULL
  kk.clusters2 = NULL
                  
  try(kk.clusters <- compareCluster(geneCluster = genes, fun = "enrichMKEGG"))
  if (!(is.null(kk.clusters))) {
    write.table(kk.clusters, 
      file = paste0(outpath,"/KEGG_Modules_Enrichments.DB.txt"), quote = FALSE, 
      sep = "\t", row.names = FALSE)
  }
                  
  try(kk.clusters2 <- compareCluster(geneCluster = genes2, fun="enrichMKEGG"))
                  
  if (!(is.null(kk.clusters))) {
    p = dotplot(kk.clusters, showCategory = 40, 
      title = "KEGG Modules (padj < 0.05) - Top 40\nDBR (FDR <= 0.05)")
    plot(p)
  }
  if (!(is.null(kk.clusters2))) {
    p = dotplot(kk.clusters2, showCategory = 40, 
      title = "KEGG Modules (padj < 0.05) - Top 40\nDBR (FDR <= 0.05)")
    plot(p)
  }


  # Reactome Pathways #
  re.clusters = NULL
  re.clusters2 = NULL

  try(re.clusters <- compareCluster(geneCluster = genes, fun="enrichPathway"))
  try(re.clusters2 <- compareCluster(geneCluster = genes2, fun="enrichPathway"))
  if (!(is.null(re.clusters))) {
    write.table(re.clusters, 
      file = paste0(outpath, "/Reactome_Enrichments.DB.txt"), 
      quote = FALSE, sep = "\t", row.names = FALSE)
  }

  if (!(is.null(re.clusters))) {
    p = dotplot(re.clusters, showCategory = 40, 
      title = "Reactome Pathways (padj < 0.05) - Top 40\nDBR (FDR <= 0.05)")
    plot(p)
  }
  
  if (!(is.null(re.clusters2))) {
    p = dotplot(re.clusters2, showCategory = 40, 
      title = "Reactome Pathways (padj < 0.05) - Top 40\nDBR (FDR <= 0.05)")
    plot(p)
  }

  dev.off()

  # Disease Associations #
  message("\n\n# DISEASE ASSOCIATIONS \n\n")
  do.clusters = NULL
  do.clusters2 = NULL
  pdf(paste0(outpath,"/Disease_Enrichments.DB.pdf"), height=12, width=12)

  try(do.clusters <- compareCluster(geneCluster = genes, fun = "enrichDO", 
    ont  = "DO", pvalueCutoff = 0.05))
  try(do.clusters2 <- compareCluster(geneCluster = genes2, fun = "enrichDO", 
    ont  = "DO", pvalueCutoff = 0.05))

  if (!(is.null(do.clusters))) {
    write.table(do.clusters, 
      file = paste0(outpath,"/Disease_Enrichments.DB.txt"), quote = FALSE, 
      sep = "\t", row.names = FALSE)
  } 
  if (!(is.null(do.clusters))) {
    p = dotplot(do.clusters, showCategory = 40, 
      title = "Disease Associations (padj < 0.05) - Top 40\nDBR (FDR <= 0.05)")
    plot(p)
  }
  if (!(is.null(do.clusters2))) {
    p = dotplot(do.clusters2, showCategory = 40, 
      title = "Disease Associations (padj < 0.05) - Top 40\nDBR (FDR <= 0.05)")
    plot(p)
  }
  dev.off()
}

