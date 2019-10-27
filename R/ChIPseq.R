
RunChIPQC <- function(base, samplesheet, facet = NULL) {
  # Base will just be the base path to create the output report in.
  # samplesheet is just a string containing the path to the samplesheet.
  # facet is a string or vector containing the columns to group samples by for plots.
  
  data(blacklist_hg19)
  exp = ChIPQC(samplesheet, consensus=TRUE, bCount=TRUE, bParallel=FALSE)
  print(exp)
  ChIPQCreport(exp, reportFolder=base, facetBy=facet)
  
  return(exp)
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
      message("Block set to improper value, must be one of 'TREATMENT', 'CONDITION', 'FACTOR' or 'TISSUE'\n\n")
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
  
  #message('# DEALING WITH DB SEs #\n\n')
  # Subset so that we can compare localization and pathways between these groups as well.
  #if (!is.null(rblock)) {
    #reportdb = dba.report(results, th = 0.05, bCalled = TRUE, 
      #bCalledDetail = TRUE, bCounts = TRUE, method = DBA_DESEQ2_BLOCK)
  #} else {
    #reportdb = dba.report(results, th = 0.05, bCalled = TRUE, 
      #bCalledDetail = TRUE, bCounts = TRUE, method = DBA_DESEQ2)
  #}
  #g1up = reportdb[(reportdb$Fold > 0)]
  #g2up = reportdb[(reportdb$Fold < 0)]
  
  #if (nrow(as.data.frame(reportdb)) > 1 & !(is.null(reportdb))) {

      # Create a named list of our subsets.
      #files = GRangesList(reportdb,g1up,g2up)
      #names(files) = c("All_DB", paste0(g1,".up"), paste0(g2,".up"))
      #peakAnnoList <- lapply(files, annotatePeak, TxDb = txdb,
        #tssRegion=c(-2000, 2000), verbose = FALSE, annoDb = "org.Hs.eg.db")
      #x = as.data.frame(peakAnnoList[1])
      #write.table(x, 
        #file = paste0(base, "/", g1, '-v-', g2, ".", mark, ".AllDB.SEs.txt"), 
        #sep = "\t", quote = FALSE, row.names = FALSE)
      #process.dbs(peakAnnoList, reportdb, results, rblock, g1, g2, base, mark, 
        #color, breaks)
  #}
  
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
    ego.cc.clusters2 <- simplify(ego.cc.clusters2, cutoff = 0.7, by = "p.adjust", 
      select_fun = min)
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


GetOverlaps <- function(h3, k27, g1, g2, base, block){
  # h3 = dataframe of H3ac peaks.
  # k27 = dataframe of H3K27ac peaks.
  # g1 = Name of group 1 for file naming, etc.
  # g2 = Name of group 2 for file naming, etc.
  # base = Base filepath to save plots.
  # block = Block
  
  # Create output folder if necessary.
  if (!dir.exists(file.path(base,"H3ac-K27ac-Overlap", block))) {
    dir.create(file.path(base,"H3ac-K27ac-Overlap", block), recursive=TRUE)
  }
  base <- file.path(base,"H3ac-K27ac-Overlap",block)
  
  # Convert each to GRanges objects for easy overlap.
  h3 <- makeGRangesFromDataFrame(h3, keep.extra.columns=TRUE)
  k27 <- makeGRangesFromDataFrame(k27, keep.extra.columns=TRUE)
  
  # Additional subsetting.
  # Differentially bound peaks
  h3.db <- h3[h3$FDR <= 0.05]
  k27.db <- k27[k27$FDR <= 0.05]
  
  h3.dbup <- h3.db[h3.db$Fold > 0]
  k27.dbup <- k27.db[k27.db$Fold > 0]
  
  h3.dbdown <- h3.db[h3.db$Fold < 0]
  k27.dbdown <- k27.db[k27.db$Fold < 0]
  
  # All promoters.
  h3.prom <- h3[(h3$cat=="PROM" | h3$cat=="PROM.SE")]
  k27.prom <- k27[(k27$cat=="PROM" | k27$cat=="PROM.SE")]
  
  h3.promup <- h3.prom[h3.prom$Fold > 0]
  k27.promup <- k27.prom[k27.prom$Fold > 0]
  
  h3.promdown <- h3.prom[h3.prom$Fold < 0]
  k27.promdown <- k27.prom[k27.prom$Fold < 0]
  
  # Promoters in SEs.
  h3.prom.se <- h3[(h3$cat=="PROM.SE")]
  k27.prom.se <- k27[(k27$cat=="PROM.SE")]
  
  h3.promup.se <- h3.prom.se[h3.prom.se$Fold > 0]
  k27.promup.se <- k27.prom.se[k27.prom.se$Fold > 0]
  
  h3.promdown.se <- h3.prom.se[h3.prom.se$Fold < 0]
  k27.promdown.se <- k27.prom.se[k27.prom.se$Fold < 0]
  
  # Promoters outside SEs.
  h3.prom.nose <- h3[(h3$cat=="PROM")]
  k27.prom.nose <- k27[(k27$cat=="PROM")]
  
  h3.promup.nose <- h3.prom.nose[h3.prom.nose$Fold > 0]
  k27.promup.nose <- k27.prom.nose[k27.prom.nose$Fold > 0]
  
  h3.promdown.nose <- h3.prom.nose[h3.prom.nose$Fold < 0]
  k27.promdown.nose <- k27.prom.nose[k27.prom.nose$Fold < 0]
  
  # All DB promoters.
  h3.db.prom <- h3.db[(h3.db$cat=="PROM" | h3.db$cat=="PROM.SE")]
  k27.db.prom <- k27.db[(k27.db$cat=="PROM" | k27.db$cat=="PROM.SE")]
  
  h3.db.promup <- h3.db.prom[h3.db.prom$Fold > 0]
  k27.db.promup <- k27.db.prom[k27.db.prom$Fold > 0]
  
  h3.db.promdown <- h3.db.prom[h3.db.prom$Fold < 0]
  k27.db.promdown <- k27.db.prom[k27.db.prom$Fold < 0]
  
  # DB Promoters in SEs.
  h3.db.prom.se <- h3.db[(h3.db$cat=="PROM.SE")]
  k27.db.prom.se <- k27.db[(k27.db$cat=="PROM.SE")]
  
  h3.db.promup.se <- h3.db.prom.se[h3.db.prom.se$Fold > 0]
  k27.db.promup.se <- k27.db.prom.se[k27.db.prom.se$Fold > 0]
  
  h3.db.promdown.se <- h3.db.prom.se[h3.db.prom.se$Fold < 0]
  k27.db.promdown.se <- k27.db.prom.se[k27.db.prom.se$Fold < 0]
  
  # DB Promoters outside SEs.
  h3.db.prom.nose <- h3.db[(h3.db$cat=="PROM")]
  k27.db.prom.nose <- k27.db[(k27.db$cat=="PROM")]
  
  h3.db.promup.nose <- h3.db.prom.nose[h3.db.prom.nose$Fold > 0]
  k27.db.promup.nose <- k27.db.prom.nose[k27.db.prom.nose$Fold > 0]
  
  h3.db.promdown.nose <- h3.db.prom.nose[h3.db.prom.nose$Fold < 0]
  k27.db.promdown.nose <- k27.db.prom.nose[k27.db.prom.nose$Fold < 0]
  
  # All enhancers.
  h3.ce <- h3[(h3$cat=="CE" | h3$cat=="CE.SE")]
  k27.ce <- k27[(k27$cat=="CE" | k27$cat=="CE.SE")]
  
  h3.ceup <- h3.ce[h3.ce$Fold > 0]
  k27.ceup <- k27.ce[k27.ce$Fold > 0]
  
  h3.cedown <- h3.ce[h3.ce$Fold < 0]
  k27.cedown <- k27.ce[k27.ce$Fold < 0]
  
  # Enhancers in SEs.
  h3.ce.se <- h3[(h3$cat=="PROM.SE")]
  k27.ce.se <- k27[(k27$cat=="PROM.SE")]
  
  h3.ceup.se <- h3.ce.se[h3.ce.se$Fold > 0]
  k27.ceup.se <- k27.ce.se[k27.ce.se$Fold > 0]
  
  h3.cedown.se <- h3.ce.se[h3.ce.se$Fold < 0]
  k27.cedown.se <- k27.ce.se[k27.ce.se$Fold < 0]
  
  # enhancers outside SEs.
  h3.ce.nose <- h3[(h3$cat=="PROM")]
  k27.ce.nose <- k27[(k27$cat=="PROM")]
  
  h3.ceup.nose <- h3.ce.nose[h3.ce.nose$Fold > 0]
  k27.ceup.nose <- k27.ce.nose[k27.ce.nose$Fold > 0]
  
  h3.cedown.nose <- h3.ce.nose[h3.ce.nose$Fold < 0]
  k27.cedown.nose <- k27.ce.nose[k27.ce.nose$Fold < 0]
  
  # All DB enhancers.
  h3.db.ce <- h3.db[(h3.db$cat=="PROM" | h3.db$cat=="PROM.SE")]
  k27.db.ce <- k27.db[(k27.db$cat=="PROM" | k27.db$cat=="PROM.SE")]
  
  h3.db.ceup <- h3.db.ce[h3.db.ce$Fold > 0]
  k27.db.ceup <- k27.db.ce[k27.db.ce$Fold > 0]
  
  h3.db.cedown <- h3.db.ce[h3.db.ce$Fold < 0]
  k27.db.cedown <- k27.db.ce[k27.db.ce$Fold < 0]
  
  # DB enhancers in SEs.
  h3.db.ce.se <- h3.db[(h3.db$cat=="PROM.SE")]
  k27.db.ce.se <- k27.db[(k27.db$cat=="PROM.SE")]
  
  h3.db.ceup.se <- h3.db.ce.se[h3.db.ce.se$Fold > 0]
  k27.db.ceup.se <- k27.db.ce.se[k27.db.ce.se$Fold > 0]
  
  h3.db.cedown.se <- h3.db.ce.se[h3.db.ce.se$Fold < 0]
  k27.db.cedown.se <- k27.db.ce.se[k27.db.ce.se$Fold < 0]
  
  # DB enhancers outside SEs.
  h3.db.ce.nose <- h3.db[(h3.db$cat=="PROM")]
  k27.db.ce.nose <- k27.db[(k27.db$cat=="PROM")]
  
  h3.db.ceup.nose <- h3.db.ce.nose[h3.db.ce.nose$Fold > 0]
  k27.db.ceup.nose <- k27.db.ce.nose[k27.db.ce.nose$Fold > 0]
  
  h3.db.cedown.nose <- h3.db.ce.nose[h3.db.ce.nose$Fold < 0]
  k27.db.cedown.nose <- k27.db.ce.nose[k27.db.ce.nose$Fold < 0]
  
  # Overlaps.
  hits <- findOverlaps(h3, k27)
  hits.prom <- findOverlaps(h3.prom, k27.prom)
  hits.promup <- findOverlaps(h3.promup, k27.promup)
  hits.promdown <- findOverlaps(h3.promdown, k27.promdown)
  hits.prom.se <- findOverlaps(h3.prom.se, k27.prom.se)
  hits.promup.se <- findOverlaps(h3.promup.se, k27.promup.se)
  hits.promdown.se <- findOverlaps(h3.promdown.se, k27.promdown.se)
  hits.prom.nose <- findOverlaps(h3.prom.nose, k27.prom.nose)
  hits.promup.nose <- findOverlaps(h3.promup.nose, k27.promup.nose)
  hits.promdown.nose <- findOverlaps(h3.promdown.nose, k27.promdown.nose)
  hits.ce <- findOverlaps(h3.ce, k27.ce)
  hits.ceup <- findOverlaps(h3.ceup, k27.ceup)
  hits.cedown <- findOverlaps(h3.cedown, k27.cedown)
  hits.ce.se <- findOverlaps(h3.ce.se, k27.ce.se)
  hits.ceup.se <- findOverlaps(h3.ceup.se, k27.ceup.se)
  hits.cedown.se <- findOverlaps(h3.cedown.se, k27.cedown.se)
  hits.ce.nose <- findOverlaps(h3.ce.nose, k27.ce.nose)
  hits.ceup.nose <- findOverlaps(h3.ceup.nose, k27.ceup.nose)
  hits.cedown.nose <- findOverlaps(h3.cedown.nose, k27.cedown.nose)
  
  hits.db <- findOverlaps(h3.db, k27.db)
  hits.db.prom <- findOverlaps(h3.db.prom, k27.db.prom)
  hits.db.promup <- findOverlaps(h3.db.promup, k27.db.promup)
  hits.db.promdown <- findOverlaps(h3.db.promdown, k27.db.promdown)
  hits.db.prom.se <- findOverlaps(h3.db.prom.se, k27.db.prom.se)
  hits.db.promup.se <- findOverlaps(h3.db.promup.se, k27.db.promup.se)
  hits.db.promdown.se <- findOverlaps(h3.db.promdown.se, k27.db.promdown.se)
  hits.db.prom.nose <- findOverlaps(h3.db.prom.nose, k27.db.prom.nose)
  hits.db.promup.nose <- findOverlaps(h3.db.promup.nose, k27.db.promup.nose)
  hits.db.promdown.nose <- findOverlaps(h3.db.promdown.nose, 
    k27.db.promdown.nose)
  hits.db.ce <- findOverlaps(h3.db.ce, k27.db.ce)
  hits.db.ceup <- findOverlaps(h3.db.ceup, k27.db.ceup)
  hits.db.cedown <- findOverlaps(h3.db.cedown, k27.db.cedown)
  hits.db.ce.se <- findOverlaps(h3.db.ce.se, k27.db.ce.se)
  hits.db.ceup.se <- findOverlaps(h3.db.ceup.se, k27.db.ceup.se)
  hits.db.cedown.se <- findOverlaps(h3.db.cedown.se, k27.db.cedown.se)
  hits.db.ce.nose <- findOverlaps(h3.db.ce.nose, k27.db.ce.nose)
  hits.db.ceup.nose <- findOverlaps(h3.db.ceup.nose, k27.db.ceup.nose)
  hits.db.cedown.nose <- findOverlaps(h3.db.cedown.nose, k27.db.cedown.nose)
  
  all.overlaps <- list(
    'All.Peaks' = hits, 'Promoters' = hits.prom, 'Promoters.Up' = hits.promup, 
    'Promoters.Down' = hits.promdown, 'Promoters.in.SEs' = hits.prom.se, 
    'Promoters.Up.in.SEs' = hits.promup.se, 
    'Promoters.Down.in.SEs' = hits.promdown.se, 
    'Promoters.Outside.SEs' = hits.prom.nose, 
    'Promoters.Up.Outside.SEs' = hits.promup.nose, 
    'Promoters.Down.outside.SEs' = hits.promdown.nose, 'Enhancers' = hits.ce, 
    'Enhancers.Up' = hits.ceup, 'Enhancers.Down' = hits.cedown, 
    'Enhancers.in.SEs' = hits.ce.se, 'Enhancers.Up.in.SEs' = hits.ceup.se, 
    'Enhancers.Down.in.SEs' = hits.cedown.se, 
    'Enhancers.Outside.SEs' = hits.ce.nose, 
    'Enhancers.Up.Outside.SEs' = hits.ceup.nose, 
    'Enhancers.Down.Outside.SEs' = hits.cedown.nose)
  
  db.overlaps <- list(
    'All.DB.Peaks' = hits.db, 'DB.Promoters' = hits.db.prom, 
    'DB.Promoters.Up' = hits.db.promup, 'DB.Promoters.Down' = hits.db.promdown, 
    'DB.Promoters.in.SEs' = hits.db.prom.se, 
    'DB.Promoters.Up.in.SEs' = hits.db.promup.se, 
    'DB.Promoters.Down.in.SEs' = hits.db.promdown.se, 
    'DB.Promoters.Outside.SEs' = hits.db.prom.nose, 
    'DB.Promoters.Up.Outside.SEs' = hits.db.promup.nose, 
    'DB.Promoters.Down.Outside.SEs' = hits.db.promdown.nose, 
    'DB.Enhancers' = hits.db.ce, 'DB.Enhancers.Up' = hits.db.ceup, 
    'DB.Enhancers.Down' = hits.db.cedown, 'DB.Enhancers.in.SEs' = hits.db.ce.se, 
    'DB.Enhancers.Up.in.SEs' = hits.db.ceup.se, 
    'DB.Enhancers.Down.in.SEs' = hits.db.cedown.se, 
    'DB.Enhancers.Outside.SEs' = hits.db.ce.nose, 
    'DB.Enhancers.Up.Outside.SEs' = hits.db.ceup.nose, 
    'DB.Enhancers.Down.Outside.SEs' = hits.db.cedown.nose)
  
  # Plotting.
  pdf(paste0(base,'/AllPeaks.Overlaps.pdf'))
  plot.hits(all.overlaps)
  dev.off()
  
  pdf(paste0(base,'/DBPeaks.Overlaps.pdf'))
  plot.hits(db.overlaps)
  dev.off()
  
  return(list(all.overlaps, db.overlaps))
}


PlotHits <- function(hits){
  # hits = List of Hits object from findOverlaps().

  for (i in 1:length(hits)) {
    name <- names(hits)[i]
    blips <- hits[[i]]

    # If there aren't any overlaps, just skip it.
    if (length(queryHits(blips)) > 0 & length(subjectHits(blips) > 0)) {
      fit <- euler(c(
        H3ac = (queryLength(blips) - length(unique(queryHits(blips)))), 
        H3K27ac = (subjectLength(blips) - length(unique(subjectHits(blips)))), 
                     "H3ac&H3K27ac" = length(blips)))
      p1 <- plot(fit, quantities = TRUE, fill = c("#009999", "#009900"), 
        alpha = 0.5)

      h3.total <- queryLength(blips)
      h3.uniq <- (queryLength(blips) - length(unique(queryHits(blips))))
      h3.uniq.p <- (h3.uniq / h3.total) * 100
      h3.ovlp <- length(unique(queryHits(blips)))
      h3.ovlp.p <- (h3.ovlp / h3.total) * 100

      k27.total <- subjectLength(blips)
      k27.uniq <- (subjectLength(blips) - 
        length(unique(subjectHits(blips))))
      k27.uniq.p <- (k27.uniq / k27.total) * 100
      k27.ovlp <- length(unique(subjectHits(blips)))
      k27.ovlp.p <- (k27.ovlp / k27.total) * 100

      count <- c(h3.total, h3.uniq, h3.ovlp, k27.total, k27.uniq, k27.ovlp)
      pct <- c(h3.uniq.p, h3.ovlp.p, k27.uniq.p, k27.ovlp.p)
      type <- c("Total", "Unique", "Overlap", "Total", "Unique", "Overlap")
      type.p <- c("Unique", "Overlap", "Unique", "Overlap")
      mark <- c("H3ac", "H3ac", "H3ac", "H3K27ac", "H3K27ac", "H3K27ac")
      mark.p <- c("H3ac", "H3ac", "H3K27ac", "H3K27ac")

      df2 <- data.frame(count = count, type = type, mark = mark)
      df2$type.ordered <- factor(df2$type, levels = c("Total", "Unique", 
        "Overlap"))
      p2 <- ggplot(data = df2, aes(x = mark, y = count, fill = type.ordered)) +
        geom_bar(stat = "identity", position = position_dodge(), 
          colour = "black") + scale_fill_manual(values=c("#00B7EC", "#FAA200", 
            "#E27DAC")) + ylab("Counts") + xlab("") + 
        theme(axis.text = element_text(size = 10),
          axis.text.x = element_text(angle = 315), 
          legend.title = element_blank())

      df3 <- data.frame(pct = pct, type = type.p, mark = mark.p)
      df3$type.ordered <- factor(df3$type, levels=c("Unique", "Overlap"))
      p3 <- ggplot(data = df3, aes(x = mark, y = pct, fill = type.ordered)) +
          geom_bar(stat = "identity", position = position_dodge(), 
            colour = "black") 
          + ylim(0, 100) + scale_fill_manual(values = c("#FAA200", "#E27DAC")) +
          ylab("% Peaks") + xlab("") + 
          theme(axis.text = element_text(size = 10), 
            axis.text.x = element_text(angle = 315), 
            legend.title = element_blank())

      g1 <- arrangeGrob(grobs = list(p2, p3), ncol = 2)
      grid.arrange(grobs = list(p1, g1), ncol = 1, 
        heights = unit(c(3, 3), c("in", "in")), top = name)
    }
  }
}
