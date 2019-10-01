#' Run DESeq2 analyses
#'
#' This function performs a high-level differential gene expression analysis 
#' 
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
#' @param block String or character vector defining the column(s) in 
#'   \code{samplesheet} to use to block for unwanted variance, e.g. batch or 
#'   technical effects.
#' @param read.filt Number indicating read threshold. Genes with fewer counts 
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
  padj.thresh = 0.05, fc.thresh = 2, plot.annos = NULL, block = NULL, 
  read.filt = 10) {
    
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
  dds <- dds[rowSums(counts(dds)) >= read.filt,]
  message(paste0("\ndds has", nrow(dds), 
    "genes after filtering rows with under ", read.filt, " reads total.\n"))

  ### VARIANCE STABILIZATION COMPARISONS ###
  message("\n# VARIANCE STABILIZATION COMPARISONS #\n")
  vst.out <- paste0(base,"/GenericFigures/VarianceTransformations.pdf")
  message(vst.out)
  trans <- RunVarianceTransformations(dds, vst.out)
  rld <- trans$rld
  vsd <- trans$vsd

  ### SAMPLE DISTANCES ###
  message("\n# PLOTTING SAMPLE DISTANCES #\n")
  dists.out <- paste0(base, "/GenericFigures/SampleDistances.pdf")
  message(dists.out)
  RunSampleDistances(rld, vsd, dists.out, level, plot.annos)

  ### PCA PLOTS ###
  message("\n# PCA PLOTS #\n")
  pca.out <- paste0(base, "/GenericFigures/pca.pdf")
  message(pca.out)

  RunPCA(rld, vsd, pca.out, level, plot.annos)
  
  #======================================#
  ### DIFFERENTIAL EXPRESSION ANALYSIS ###
  message("\n### DIFFERENTIAL EXPRESSION ANALYSIS ###\n")
  dds <- DESeq(dds)

  message("\n# SAVING ROBJECTS #\n")
  message(paste0("DESeq2: ", paste0(base, "/Robjects/dds.rds")))
  message(paste0("Regularized log transformation: ", 
    paste0(base, "/Robjects/rld.rds")))
  message(paste0("Variance stabilized transformation: ", 
    paste0(base, "/Robjects/vld.rds")))
  saveRDS(dds, file=paste0(base, "/Robjects/dds.rds"))
  saveRDS(rld, file=paste0(base, "/Robjects/rld.rds"))
  saveRDS(vld, file=paste0(base, "/Robjects/vld.rds"))

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
#' @param padj.thresh Number or numeric scalar indicating the adjusted p-value 
#'   cutoff(s) to be used for determining "significant" differential expression.
#'   If multiple are given, multiple tables/plots will be generated using all 
#'   combinations of \code{padj.thresh} and \code{fc.thresh}.
#' @param fc.thresh Number or numeric scalar indicating the log2 fold-change 
#'   cutoff(s) to be used for determining "significant" differential expression.
#'   If multiple are given, multiple tables/plots will be generated using all 
#'   combinations of \code{padj.thresh} and \code{fc.thresh}.
#' @param top.n Number of differentially expressed genes to create boxplots for, 
#'   ranked by adj. p-value after applying \code{padj.thresh} and 
#'   \code{fc.thresh} thresholds. If multiple thresholds are provided, the 
#'   lowest fold-change and highest adj. p-value thresholds will be used. 
#'
#'   Plots will be created for each comparison. 
#'
#' @importFrom DESeq2 lfcShrink counts plotCounts
#' @importFrom magrittr %>%
#'
#' @export
#'
#' @author Jared Andrews
#'
ProcessDEGs <- function(dds, rld, vsd, outpath, level, padj.thresh = 0.05, 
  fc.thresh = 2, top.n = 100) {

  # Get all possible sample comparisons.
  combs <- combn(colData(rld)[,level], 2)
  combs.seq <- seq(1, length(combs), by = 2)

  res.list <- list()
  for (samp in combs.seq) {
    g1 <- combs[samp]
    g2 <- combs[samp + 1]
    res <- lfcShrink(dds, contrast=c(level, g1, g2), type='ashr')
    res.list[[paste0(g1, "v", g2)]] <- res
  }
  
  ### PLOTTING THE RESULTS ###
  # This gets all the normalized counts for all genes/samples.
  resdata <- merge(as.data.frame(res), 
    as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
  names(resdata)[1] <- "Gene"

  # Make table for just DEGs.
  resSig <- subset(resdata, padj <= 0.1)

  ### DEG BOXPLOTS ###
  message("\n# DEG BOXPLOTS #\n")
  message("Creating boxplots for all genes with padj <= 0.1")
  # If you have >7 levels for your contrast, you need to add colors here.
  fill = c("#A6CEE3", "#B2DF8A", "#FDBF6F", "#CAB2D6", 
    "#f4cae4", "#f1e2cc", "#b3e2cd")
  line = c("#1F78B4", "#33A02C", "#FF7F00", "#6A3D9A", 
    "#e7298a", "#a6761d", "#1b9e77")
  if (nrow(resSig) > 5) {
    for (i in 1:nrow(resSig)) {
      if (!file.exists(paste0(base, "/GeneBoxPlots/", 
        gsub('/','-',resSig$Gene[i]),".BoxPlot.pdf"))) {
        pdf(paste0(base, "/GeneBoxPlots/", gsub('/','-',resSig$Gene[i]), 
          ".BoxPlot.pdf"))
        d <- plotCounts(dds, gene = resSig$Gene[i], intgroup = level, 
          returnData = T)
        p <- ggplot(d, aes(x = d[,level], y = count)) + 
          geom_boxplot(fill = fill[1:length(levels(colData(rld)[,level]))], 
            colour = line[1:length(levels(colData(rld)[,level]))]) + 
          ggtitle(resSig$Gene[i]) + coord_trans(y = "log10")
        print(p)

        d <- plotCounts(dds, gene = resSig$Gene[i], intgroup = level, 
          returnData = T)
        e <- subset(d, (get(level) == g1 | get(level) == g2))
        p <- ggplot(e, aes(x = e[,level], y = count)) + 
          geom_boxplot(fill=fill[1:2],colour=line[1:2]) + 
          ggtitle(resSig$Gene[i]) + coord_trans(y = "log10")
        print(p)
        dev.off()
      }
    }
  }
  
  ### DEG PCA PLOTS ###
  # TODO - Make this a generic function.
  message(paste0("\n# DEG PCA PLOTS #\n", base, "/DEG_pca.pdf\n"))
  pdf(paste0(base, "/DEG_pca.pdf"))
  resdata <- merge(as.data.frame(res), as.data.frame(assay(rld)), 
    by = "row.names", sort = FALSE)
  names(resdata)[1] <- "Gene"
  resSig <- subset(resdata, padj <= 0.1)
  # Show only DEG genes and only the wanted groups.
  if (nrow(resSig) > 5) {
      rld.sub <- rld[resSig$Gene, colData(rld)[,level] %in% c(g1, g2)]
      p <- DESeq2::plotPCA(rld.sub, intgroup = plot.annos) + 
        ggtitle("DEGs - padj <= 0.1")
      print(p)
      p <- DESeq2::plotPCA(rld.sub, intgroup = level) + 
        ggtitle("DEGs - padj <= 0.1")
      print(p)
  }
  
  # Now with lfc threshold on genes as well.
  resSig <- subset(resSig, log2FoldChange >= 2 | log2FoldChange <= -2)
  # Show only DEG genes and only the wanted groups.
  if (nrow(resSig) > 5) {
      rld.sub <- rld[resSig$Gene, colData(rld)[,level] %in% c(g1, g2)]
      p <- DESeq2::plotPCA(rld.sub, intgroup = plot.annos) + 
        ggtitle("DEGs - padj <= 0.1 & LFC >|< 2")
      print(p)
      p <- DESeq2::plotPCA(rld.sub, intgroup = level) + 
        ggtitle("DEGs - padj <= 0.1 & LFC >|< 2")
      print(p)
  }
  
  resSig <- subset(resdata, padj <= 0.05)
  if (nrow(resSig) > 5) {
      # Show only DEG genes and only the wanted groups.
      rld.sub <- rld[resSig$Gene, colData(rld)[,level] %in% c(g1, g2)]
      p <- DESeq2::plotPCA(rld.sub, intgroup = plot.annos) + 
        ggtitle("DEGs - padj <= 0.05")
      print(p)
      p <- DESeq2::plotPCA(rld.sub, intgroup = level) + 
        ggtitle("DEGs - padj <= 0.05")
      print(p)
  }
  
  # Now with lfc threshold on genes as well.
  resSig <- subset(resdata, padj <= 0.05)
  resSig <- subset(resSig, log2FoldChange >= 2 | log2FoldChange <= -2)
  # Show only DEG genes and only the wanted groups.
  if (nrow(resSig) > 5) {
      rld.sub <- rld[resSig$Gene, colData(rld)[,level] %in% c(g1, g2)]
      p <- DESeq2::plotPCA(rld.sub, intgroup = plot.annos) + 
        ggtitle("DEGs - padj <= 0.05 & LFC >|< 2")
      print(p)
      p <- DESeq2::plotPCA(rld.sub, intgroup = level) + 
        ggtitle("DEGs - padj <= 0.05 & LFC >|< 2")
      print(p)
  }
  dev.off()
  
  ### MA PLOT ###
  message(paste0("\n# MA PLOT #\n", base, "/MA_plot.pdf\n"))
  pdf(paste0(base,"/MA_plot.pdf"))
  plotMA(res, ylim = c(-5, 5))
  dev.off()
  
  ### VOLCANO PLOT ###
  message(paste0("# VOLCANO PLOT #\n", getwd(), "/", substr(base, 2, 500), "/", 
    g1, "-v-", g2, ".DEGs.meanCounts50.Volcano.padj0.1.html\n"))
  
  # Set up results table.
  res <- res[order(res$padj),]
  resdata <- merge(as.data.frame(res), 
    as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
  resvol <- subset(resdata, baseMean >= 50)
  names(resvol)[1] <- "Gene"
  
  diff_df <- resvol[c("Gene","log2FoldChange","padj")]

  # add a grouping column; default value is "not significant"
  diff_df["group"] <- "NotSignificant"

  # for our plot, we want to highlight 
  # padj < 0.1 (significance level)
  # log2FoldChange > 2

  # Change groupings as necessary.
  diff_df[which(diff_df['padj'] <= 0.1 & abs(diff_df['log2FoldChange']) <= 2 ), 
    "group"] <- "Significant"
  diff_df[which(diff_df['padj'] >= 0.1 & abs(diff_df['log2FoldChange']) >= 2 ), 
    "group"] <- "FoldChange"
  diff_df[which(diff_df['padj'] <= 0.1 & abs(diff_df['log2FoldChange']) >= 2 ), 
    "group"] <- "Significant&FoldChange"


  p <- plot_ly(
    diff_df, 
    x = ~log2FoldChange, 
    y = -log10(diff_df$padj), 
    text = ~Gene, 
    color = ~group,
    mode = "markers", 
    type = "scatter") %>%
    layout(title = paste(g1,"v",g2, sep='-'), 
      xaxis = list(title = "log2 Fold Change"), 
      yaxis = list(title = "-log10(adjusted p-value)")) 

  # Save plot - has to use absolute path for some reason.
  htmlwidgets::saveWidget(as_widget(p), paste0(getwd(), "/", 
    substr(base, 2, 500), "/", g1, "-v-", g2, 
    ".DEGs.meanCounts50.Volcano.padj0.1.html"))
  
  # PADJ <= 0.05 #
  # Set up results table.
  res <- res[order(res$padj),]
  resdata <- merge(as.data.frame(res), 
    as.data.frame(counts(dds, normalized = TRUE)), by = "row.names", 
    sort = FALSE)
  names(resdata)[1] <- "Gene"
  resvol <- subset(resdata, baseMean >= 50)
  
  diff_df <- resvol[c("Gene","log2FoldChange","padj")]

  # add a grouping column; default value is "not significant"
  diff_df["group"] <- "NotSignificant"

  # for our plot, we want to highlight 
  # padj < 0.05 (significance level)
  # log2FoldChange > 2

  # Change groupings as necessary.
  diff_df[which(diff_df['padj'] <= 0.05 & abs(diff_df['log2FoldChange']) <= 2 ),
    "group"] <- "Significant"
  diff_df[which(diff_df['padj'] >= 0.05 & abs(diff_df['log2FoldChange']) >= 2 ),
    "group"] <- "FoldChange"
  diff_df[which(diff_df['padj'] <= 0.05 & abs(diff_df['log2FoldChange']) >= 2 ),
    "group"] <- "Significant&FoldChange"


  p <- plot_ly(
    diff_df, 
    x = ~log2FoldChange, 
    y = -log10(diff_df$padj), 
    text = ~Gene, color = ~group, 
    mode = "markers", type = "scatter") %>%
      layout(title = paste(g1, "v", g2, sep='-'), 
        xaxis = list(title = "log2 Fold Change"), 
        yaxis = list(title = "-log10(adjusted p-value)")) 

  # Save plot - has to use absolute path for some reason.
  htmlwidgets::saveWidget(as_widget(p), paste0(getwd(), "/", 
    substr(base, 2, 500), "/", g1, "-v-", g2, 
    ".DEGs.meanCounts50.Volcano.padj0.05.html"))
  
  ### HEATMAPS ###
  message(paste0("\n# DEG HEATMAPS #\n"))

  pdf(paste0(base, "/DEG_Heatmaps.pdf"))
  MakeHeatmaps(dds, res, rld, level, g1, g2, plot.annos)
  dev.off()

  ### SAVING TABLES ###
  message("\n# SAVING RESULTS TABLES #\n")
  saveRDS(dds, file=paste0(base, "/dds.rds"))
  saveRDS(res, file=paste0(base, "/res.rds"))
  saveRDS(rld, file=paste0(base, "/rld.rds"))

  write.csv(resdata, file=paste0(base, "/", g1, "-v-", g2, ".ALL.csv"), 
    row.names = F)

  # Make table for just DEGs too.
  resSig <- subset(resdata, padj <= 0.1)
  write.csv(resSig, file=paste0(base, "/", g1, "-v-", g2, ".DEGs.padj0.1.csv"), 
    row.names = F)

  resSig <- subset(resdata, padj <= 0.05)
  write.csv(resSig, file=paste0(base, "/", g1, "-v-", g2, ".DEGs.padj0.05.csv"), 
    row.names = F)
  
  # This fpm version allows for accurate intergene comparisons of counts.
  resdata <- merge(as.data.frame(res), as.data.frame(fpm(dds)), 
    by = "row.names", sort = FALSE)
  names(resdata)[1] <- "Gene"
  write.csv(fpm(dds), file=paste0(base, "/", g1, "-v-", g2,".ALL.fpm.csv"), 
    row.names = F)

  return(list(res=res, dds=dds, rld=rld))
}