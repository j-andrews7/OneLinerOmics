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
#'
#'
#'
#'
#'
#'
RunDESeq2 <- function(outpath, quantspath, samplesheet, tx2gene, level, g1, g2, 
  plot_annos=NULL, block=NULL) {
# plot_annos should be a vector of at least two columns from the sample sheet to use for labeling plots. 
# g1 and g2 should be strings for the groups to be compared for the level (eg. "Progression", "Response")
# block can either be a single string ("cells") or a vector of multiple (c("cells", "disease"))
    
  message("### EXPLORATORY DATA ANALYSIS ###\n")
  message("# SET DIRECTORY STRUCTURE AND MODEL DESIGN #\n")
  ### SET DIRECTORY STRUCTURE AND MODEL DESIGN ###
  # Base folder to create output folders in. Create output folders.
  base <- outpath
  
  # Determine if there is a blocking factor or not and set design.
  if (!is.null(block)) {
    design <- formula(paste("~", paste(c(block, level), sep = "", 
      collapse = " + ")))
    if (!dir.exists(file.path(base, paste0(block, "Block")))) {
      dir.create(file.path(base, paste0(block, "Block")))
    }
    if (!dir.exists(file.path(base, paste0(block, "Block/GenericFigures")))) {
      dir.create(file.path(base, paste0(block, "Block/GenericFigures")))
    }
    if (!dir.exists(file.path(base, paste0(block, "Block/GeneBoxPlots")))) {
      dir.create(file.path(base, paste0(block, "Block/GeneBoxPlots")))
    }
    base <- file.path(base, paste0(block,"Block"))
  } else {
    design <- formula(paste("~", level))
    if (!dir.exists(file.path(base, "NoBlock"))) {
      dir.create(file.path(base, "NoBlock"))
    }
    if (!dir.exists(file.path(base, "NoBlock/GenericFigures"))) {
      dir.create(file.path(base, "NoBlock/GenericFigures"))
    }
    if (!dir.exists(file.path(base, "NoBlock/GeneBoxPlots"))) {
      dir.create(file.path(base, "NoBlock/GeneBoxPlots"))
    }
    base <- file.path(base, "NoBlock")
  }
  
  if (is.null(plot_annos)) {
      plot_annos <- level
  }
  
  message("# FILE LOADING & PRE-FILTERING #\n")
  ### FILE LOADING & PRE-FILTERING ###
  tx2gene <- read.table(tx2gene, sep = "\t")    
  samples <- read.table(samplesheet, header = TRUE)
  rownames(samples) <- samples$name
  files <- file.path(quantspath, samples$name, "quant.sf")
  names(files) <- samples$name

  # Read in our actual count files now.
  txi <- tximport(files, type = "salmon", tx2gene = tx2gene)

  message(paste0("\n\nDesign is: ", design, "\n\n"))
  # Create the DESeqDataSet object.
  dds <- DESeqDataSetFromTximport(txi, colData = samples, design = design)

  # Pre-filter transcripts with really low read counts (<100 across all samples)
  message(paste("\ndds starting with", nrow(dds), "genes."))
  dds <- dds[rowSums(counts(dds)) >= 100,]
  message(paste("\ndds has", nrow(dds), 
    "genes after filtering rows with under 100 reads total.\n"))

  message("\n# VARIANCE STABILIZATION COMPARISONS #\n")
  message(paste0(base, "/GenericFigures/transformation.pdf\n"))
  ### VARIANCE STABILIZATION COMPARISONS ###
  rld <- rlog(dds, blind = FALSE)
  vsd <- vst(dds, blind = FALSE)

  pdf(paste(base,"/GenericFigures/transformation.pdf", sep = ""))

  dds <- DESeq2::estimateSizeFactors(dds)

  df <- bind_rows(
    as_data_frame(log2(counts(dds, normalized=TRUE)[, 1:2]+1)) %>%
           mutate(transformation = "log2(x + 1)"),
    as_data_frame(assay(rld)[, 1:2]) %>% mutate(transformation = "rlog"),
    as_data_frame(assay(vsd)[, 1:2]) %>% mutate(transformation = "vst"))
  
  colnames(df)[1:2] <- c("x", "y")  

  p <- ggplot(df, aes(x = x, y = y)) + geom_hex(bins = 80) +
    coord_fixed() + facet_grid( . ~ transformation) 
  print(p)

  dev.off()

  ### SAMPLE DISTANCES ###
  message(paste0("\n# PLOTTING SAMPLE DISTANCES #\n", base, 
    "/GenericFigures/samp_dist.pdf\n"))
  sampleDists <- dist(t(assay(rld)))

  pdf(paste0(base, "/GenericFigures/samp_dist.pdf"))
  ## ----distheatmap, fig.width = 6.1, fig.height = 4.5----------------------
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
  dev.off()

  ### PCA PLOTS ###
  message(paste0("# PCA PLOTS #\n", base, "/GenericFigures/pca.pdf\n"))
  pcaData <- DESeq2::plotPCA(rld, intgroup = plot_annos, returnData = TRUE)
  percentVar <- round(100 * attr(pcaData, "percentVar"))
  pdf(paste0(base, "/GenericFigures/pca.pdf"))
  p <- DESeq2::plotPCA(rld, intgroup = plot_annos)
  print(p)
  ## ----ggplotpca, fig.width=6, fig.height=4.5------------------------------
  p <- ggplot(pcaData, aes(x = PC1, y = PC2, 
    color = colData(rld)[,plot_annos[1]], 
    shape = colData(rld)[,plot_annos[2]])) +
    geom_point(size =3) +
    xlab(paste0("PC1: ", percentVar[1], "% variance")) +
    ylab(paste0("PC2: ", percentVar[2], "% variance")) +
    coord_fixed()
  print(p)
  # Now with only the groups we want to compare.
  pcaData <- DESeq2::plotPCA(rld.sub, intgroup = plot_annos, returnData = TRUE)
  percentVar <- round(100 * attr(pcaData, "percentVar"))
  p <- DESeq2::plotPCA(rld.sub, intgroup = plot_annos)
  print(p)
  ## ----ggplotpca, fig.width=6, fig.height=4.5------------------------------
  p <- ggplot(pcaData, aes(x = PC1, y = PC2, 
    color = colData(rld.sub)[,plot_annos[1]], 
    shape = colData(rld.sub)[,plot_annos[2]])) +
    geom_point(size =3) +
    xlab(paste0("PC1: ", percentVar[1], "% variance")) +
    ylab(paste0("PC2: ", percentVar[2], "% variance")) +
    coord_fixed()
  print(p)
  dev.off()
  
  #======================================#
  ### DIFFERENTIAL EXPRESSION ANALYSIS ###
  message("\n### DIFFERENTIAL EXPRESSION ANALYSIS ###\n")
  dds <- DESeq(dds)
  res <- lfcShrink(dds, contrast=c(level, g1, g2), type='ashr')
  
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
    for (i in 1:nrow(resSig)){
      if (!file.exists(paste0(base, "/GeneBoxPlots/", 
        gsub('/','-',resSig$Gene[i]),".BoxPlot.pdf"))) {
        pdf(paste0(base, "/GeneBoxPlots/", gsub('/','-',resSig$Gene[i]), 
          ".BoxPlot.pdf"))
        d <- plotCounts(dds, gene = resSig$Gene[i], intgroup = level)
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
      p <- DESeq2::plotPCA(rld.sub, intgroup = plot_annos) + 
        ggtitle("DEGs - padj <= 0.1")
      print(p)
      p <- DESeq2::plotPCA(rld.sub, intgroup = level) + 
        ggtitle("DEGs - padj <= 0.1")
      print(p)
  }
  
  # Now with lfc threshold on genes as well.
  resSig <- subset(resSig, log2FoldChange >= 2 | log2FoldChange <= -2)
  # Show only DEG genes and only the wanted groups.
  if (nrow(resSig) > 5){
      rld.sub <- rld[resSig$Gene, colData(rld)[,level] %in% c(g1, g2)]
      p <- DESeq2::plotPCA(rld.sub, intgroup = plot_annos) + 
        ggtitle("DEGs - padj <= 0.1 & LFC >|< 2")
      print(p)
      p <- DESeq2::plotPCA(rld.sub, intgroup = level) + 
        ggtitle("DEGs - padj <= 0.1 & LFC >|< 2")
      print(p)
  }
  
  resSig <- subset(resdata, padj <= 0.05)
  if (nrow(resSig) > 5){
      # Show only DEG genes and only the wanted groups.
      rld.sub <- rld[resSig$Gene, colData(rld)[,level] %in% c(g1, g2)]
      p <- DESeq2::plotPCA(rld.sub, intgroup = plot_annos) + 
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
  if (nrow(resSig) > 5){
      rld.sub <- rld[resSig$Gene, colData(rld)[,level] %in% c(g1, g2)]
      p <- DESeq2::plotPCA(rld.sub, intgroup = plot_annos) + 
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
  
  ### GENE CLUSTERING - VARIABLE ###
  # First with all samples included.
  message(paste0("\n# VARIABLE GENE CLUSTERING #\n", base, 
    "/Top100_VariableGenes_HeatmapDistances.pdf\n"))
  topVarGenes <- head(order(rowVars(assay(rld)), decreasing = TRUE), 100)
  mat  <- assay(rld)[ topVarGenes, ]
  mat  <- mat - rowMeans(mat)
  anno <- as.data.frame(colData(rld)[, plot_annos])
  pdf(paste0(base,"/Top100_VariableGenes_HeatmapDistances.pdf"))
  p <- pheatmap(mat, annotation_col = anno, main="Top 100 Variable Genes")
  print(p)

  #Then with only the wanted groups.
  rld.sub <- rld[, colData(rld)[,level] %in% c(g1, g2)]
  topVarGenes <- head(order(rowVars(assay(rld.sub)), decreasing = TRUE), 100)
  mat  <- assay(rld.sub)[ topVarGenes, ]
  mat  <- mat - rowMeans(mat)
  anno <- as.data.frame(colData(rld.sub)[, plot_annos])
  p <- pheatmap(mat, annotation_col = anno, main = "Top 100 Variable Genes")
  print(p)
  
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
  MakeHeatmaps(dds, res, rld, level, g1, g2, plot_annos)
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