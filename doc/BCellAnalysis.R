## ---- echo=FALSE, results="hide", message=FALSE----------------------------
knitr::opts_chunk$set(error=FALSE, message=FALSE, warning=FALSE, eval=FALSE)
library(BiocStyle)

## --------------------------------------------------------------------------
#  # R
#  suppressPackageStartupMessages(library(DESeq2))
#  suppressPackageStartupMessages(library(tximport))
#  suppressPackageStartupMessages(library(dplyr))
#  suppressPackageStartupMessages(library(ggplot2))
#  suppressPackageStartupMessages(library(pheatmap))
#  suppressPackageStartupMessages(library(RColorBrewer))
#  suppressPackageStartupMessages(library(repr))
#  suppressPackageStartupMessages(library(enrichR))
#  suppressPackageStartupMessages(library(cowplot))
#  suppressPackageStartupMessages(library(scales))

## --------------------------------------------------------------------------
#  # R
#  # Basic sample sheet.
#  samples <- read.table("./RNA_Seq/Samples.txt", header=TRUE, sep = "\t")
#  rownames(samples) <- samples$name
#  files <- file.path("./RNA_Seq/quants", samples$name, "quant.sf")
#  names(files) <- samples$name
#  
#  # This will be used to map transcript IDs back to gene symbols and collapse to gene-level counts.
#  tx2gene <- read.table("./RNA_Seq/ref/tx2gene.txt", sep = "\t", as.is = TRUE)
#  tx2gene <- tx2gene[,c(1,3)]
#  
#  # Read in our actual count files now.
#  txi <- tximport(files, type="salmon", tx2gene=tx2gene)
#  
#  # Create dds object using disease and ignore effects due to paired/single-end differences.
#  dds <- DESeqDataSetFromTximport(txi,
#                                     colData = samples,
#                                     design = ~ disease)

## --------------------------------------------------------------------------
#  # R
#  # Compare variance stabilization methods.
#  vsd <- vst(dds, blind = FALSE)
#  rld <- rlog(dds, blind = FALSE)
#  
#  suppressPackageStartupMessages(library("dplyr"))
#  suppressPackageStartupMessages(library("ggplot2"))
#  
#  dds <- estimateSizeFactors(dds)
#  
#  df <- bind_rows(
#    as_data_frame(log2(counts(dds, normalized=TRUE)[, 1:2]+1)) %>%
#           mutate(transformation = "log2(x + 1)"),
#    as_data_frame(assay(vsd)[, 1:2]) %>% mutate(transformation = "vst"),
#    as_data_frame(assay(rld)[, 1:2]) %>% mutate(transformation = "rlog"))
#  
#  colnames(df)[1:2] <- c("x", "y")
#  
#  options(repr.plot.width=9, repr.plot.height=7)
#  df <- bind_rows(
#    as_data_frame(log2(counts(dds, normalized=TRUE)[, 1:2]+1)) %>%
#           mutate(transformation = "log2(x + 1)"),
#    as_data_frame(assay(vsd)[, 1:2]) %>% mutate(transformation = "vst"),
#    as_data_frame(assay(rld)[, 1:2]) %>% mutate(transformation = "rlog"))
#  
#  colnames(df)[1:2] <- c("x", "y")
#  ggplot(df, aes(x = x, y = y)) + geom_hex(bins = 80) +
#    coord_fixed() + facet_grid( . ~ transformation)

## --------------------------------------------------------------------------
#  # R
#  # Sample distances heatmap
#  suppressPackageStartupMessages(library("pheatmap"))
#  suppressPackageStartupMessages(library("RColorBrewer"))
#  suppressPackageStartupMessages(library("repr"))
#  options(repr.plot.width=5, repr.plot.height=6)
#  
#  sampleDists <- dist(t(assay(rld)))
#  sampleDistMatrix <- as.matrix(sampleDists)
#  rownames(sampleDistMatrix) <- paste(rld$name, rld$disease.sub, sep = " - " )
#  colnames(sampleDistMatrix) <- NULL
#  colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
#  pheatmap(sampleDistMatrix,
#           clustering_distance_rows = sampleDists,
#           clustering_distance_cols = sampleDists,
#           col = colors)

## --------------------------------------------------------------------------
#  # R
#  # PCA
#  plotPCA(vsd, intgroup=c("disease"))

## --------------------------------------------------------------------------
#  # R
#  dim(dds)
#  keep <- rowSums(counts(dds)) >= 100
#  dds <- dds[keep,]
#  
#  # See how many were removed.
#  dim(dds)

## --------------------------------------------------------------------------
#  # R
#  dds <- DESeq(dds)
#  res <- results(dds)

## --------------------------------------------------------------------------
#  # R
#  CLLvNORM <- results(dds, contrast=c("disease", "CLL", "NORMAL"))
#  CLLvNORM <- lfcShrink(dds, contrast=c("disease", "CLL", "NORMAL"), res = CLLvNORM, type = "ashr")
#  
#  CLLvDL <- results(dds, contrast=c("disease", "CLL", "DL"))
#  CLLvDL <- lfcShrink(dds, contrast=c("disease", "CLL", "DL"), res = CLLvDL, type = "ashr")
#  
#  CLLvFL <- results(dds, contrast=c("disease", "CLL", "FL"))
#  CLLvFL <- lfcShrink(dds, contrast=c("disease", "CLL", "FL"), res = CLLvFL, type = "ashr")
#  
#  DLvFL <- results(dds, contrast=c("disease", "DL", "FL"))
#  DLvDL <- lfcShrink(dds, contrast=c("disease", "DL", "FL"), res = DLvFL, type = "ashr")
#  
#  DLvNORM <- results(dds, contrast=c("disease", "DL", "NORMAL"))
#  DLvNORM <- lfcShrink(dds, contrast=c("disease", "DL", "NORMAL"), res = DLvNORM, type = "ashr")
#  
#  FLvNORM <- results(dds, contrast=c("disease", "FL", "NORMAL"))
#  FLvNORM <- lfcShrink(dds, contrast=c("disease", "FL", "NORMAL"), res = FLvNORM, type = "ashr")

## --------------------------------------------------------------------------
#  # R
#  log2cutoff <- 2
#  qvaluecutoff <- 0.01
#  
#  sigGenes <- unique(c(
#    rownames(subset(CLLvNORM, padj<=qvaluecutoff & abs(log2FoldChange)>=log2cutoff)),
#    rownames(subset(CLLvDL, padj<=qvaluecutoff & abs(log2FoldChange)>=log2cutoff)),
#    rownames(subset(CLLvFL, padj<=qvaluecutoff & abs(log2FoldChange)>=log2cutoff)),
#    rownames(subset(DLvFL, padj<=qvaluecutoff & abs(log2FoldChange)>=log2cutoff)),
#    rownames(subset(DLvNORM, padj<=qvaluecutoff & abs(log2FoldChange)>=log2cutoff)),
#    rownames(subset(FLvNORM, padj<=qvaluecutoff & abs(log2FoldChange)>=log2cutoff))
#  ))
#  
#  heat <- assay(rld)[sigGenes,]
#  
#  options(repr.plot.width=5, repr.plot.height=7)
#  myCol <- colorRampPalette(c("darkblue", "snow", "darkred"))(1000)
#  colors <- c(seq(-2.5,-.11,length=500),seq(-.1,.1,length=1),seq(.11,2.5,length=500))
#  df <- as.data.frame(colData(dds)[,c("disease", "disease.sub")])
#  df$disease.sub <- NULL
#  p <- pheatmap(heat, cluster_rows=TRUE, show_rownames=FALSE,
#           cluster_cols=FALSE, annotation_col=df, fontsize_col = 6, fontsize_row = 6,
#          fontsize = 6, scale = "row", color=myCol, breaks = colors)
#  
#  pdf("./RNA_Seq/Final_Analyses/Figures/padj.05.lfc.2.AllComparisons.DEGs.pdf", height = 7, width = 5)
#  print(p)
#  dev.off()
#  
#  print(p)

## --------------------------------------------------------------------------
#  # R
#  
#  # Quick look at number of DEGs in each comparison.
#  length(rownames(subset(CLLvNORM, padj<=qvaluecutoff & abs(log2FoldChange)>=log2cutoff)))
#  length(rownames(subset(CLLvDL, padj<=qvaluecutoff & abs(log2FoldChange)>=log2cutoff)))
#  length(rownames(subset(CLLvFL, padj<=qvaluecutoff & abs(log2FoldChange)>=log2cutoff)))
#  length(rownames(subset(DLvFL, padj<=qvaluecutoff & abs(log2FoldChange)>=log2cutoff)))
#  length(rownames(subset(DLvNORM, padj<=qvaluecutoff & abs(log2FoldChange)>=log2cutoff)))
#  length(rownames(subset(FLvNORM, padj<=qvaluecutoff & abs(log2FoldChange)>=log2cutoff)))

## ---- results="hide"-------------------------------------------------------
#  # R
#  # Top 100 for each category.
#  DLvFL.100 <- subset(DLvFL, padj <= qvaluecutoff & abs(log2FoldChange) >= log2cutoff)
#  DLvFL.100 <- DLvFL.100[order(-(abs(DLvFL.100$log2FoldChange))),]
#  DLvFL.100 <- DLvFL.100[1:100,]

## --------------------------------------------------------------------------
#  # R
#  RunEnrichr <- function(genes, libraries = c("GO_Molecular_Function_2018",
#  	"GO_Cellular_Component_2018", "GO_Biological_Process_2018", "KEGG_2019_Human",
#  	"Reactome_2016", "BioCarta_2016", "Panther_2016"), outdir = NULL) {
#  
#  	# Required or package will throw an error saying the website is not responding
#  	# due to poor checking.
#  	options(base_address = "http://amp.pharm.mssm.edu/Enrichr/")
#  	options(enrichRLive = TRUE)
#  
#  	# Run enrichments.
#  	message("Submitting sets to the Enrichr server. Please be patient, this can",
#  		" sometimes take a few minutes if the server is under heavy load or down.")
#  	results <- enrichr(genes, libraries)
#  	res.names <- names(results)
#  
#  	if (!is.null(outdir)) {
#  		y = 1
#  		for (i in results) {
#  			df <- as.data.frame(i, sep = "\t")
#  			df$Old.P.value <- NULL
#  			df$Old.Adjusted.P.value <- NULL
#  			write.table(df, file = sprintf("%s/%s.Results.txt", outdir, res.names[y]),
#  				sep = "\t", quote = FALSE, row.names = FALSE)
#  			y <- y + 1
#  		}
#  	}
#  
#  	return(results)
#  }

## --------------------------------------------------------------------------
#  # R
#  VizEnrichments <- function(enrichments, outdir = NULL,
#  	n.terms = 10, remove.insig = TRUE, adj.p.thresh = 0.05,
#  	colors = c("grey", "darkred")) {
#  
#  	if (length(colors) > 2) {
#  		stop("Only two colors can be provided, get that fancy stuff outta here.")
#  	}
#  
#  	ind = 1
#  	for (i in enrichments) {
#  		lib <- names(enrichments[ind])
#  		ind <- ind + 1
#  		# Significance filter.
#  		if (isTRUE(remove.insig)) {
#  			i <- i[which(i$Adjusted.P.value <= adj.p.thresh), ]
#  			if (nrow(i) == 0) {
#  				message(paste0("Skipping ", lib, " due to no terms meeting the ",
#  					"significance threshold."))
#  				next
#  			}
#  		}
#  
#  		message(paste0("Plotting ", lib))
#  
#  		i$log.Adj.p <- -log10(as.numeric(i$Adjusted.P.value))
#  		i.p <- i[order(-i$log.Adj.p),]
#  		i.s <- i[order(-i$Combined.Score),]
#  
#  		# Limit number of terms.
#  		if (!is.null(n.terms)) {
#  			if (nrow(i) > n.terms) {
#  				i.p <- i.p[1:n.terms, ]
#  				i.s <- i.s[1:n.terms, ]
#  			}
#  		}
#  		
#  		p1 <- ggplot(i.p, aes(x = reorder(Term, log.Adj.p),
#  			log.Adj.p, fill = log.Adj.p)) + geom_col() + coord_flip() +
#  			cowplot::theme_cowplot(12) + theme(axis.text.y = element_text(size = 7),
#  				legend.title = element_text(size=10),
#  				legend.text = element_text(size = 10)) +
#  			scale_x_discrete(labels = wrap_format(40)) +
#  			scale_fill_gradient(low = colors[1], high = colors[2]) + xlab("Term") +
#  			ylab("-log10(Adjusted p-value)") + ylim(0, NA)
#  
#  		p2 <- ggplot(i.s, aes(x = reorder(Term, Combined.Score),
#  			Combined.Score, fill = log.Adj.p)) + geom_col() + coord_flip() +
#  			cowplot::theme_cowplot(12) + theme(axis.text.y = element_text(size = 7),
#  				legend.title = element_text(size=10),
#  				legend.text = element_text(size = 10)) +
#  			scale_x_discrete(labels = wrap_format(40)) +
#  			scale_fill_gradient(low = colors[1], high = colors[2]) + xlab("Term") +
#  			ylab("Combined Score (Enrichr)")
#  
#      c <- cowplot::plot_grid(
#        p1, p2, rel_widths = c(1, 1),
#        nrow = 1
#      )
#  
#  		h <- 0.8 + (0.35 * nrow(i.p))
#  
#      if (!is.null(outdir)) {
#  		  pdf(sprintf("%s/%s.Enrichments.pdf", outdir, lib), height = h, width = 12)
#  		  print(c)
#  		  dev.off()
#      } else {
#        print(c)
#      }
#  	}
#  }

## --------------------------------------------------------------------------
#  # R
#  DoEnrichments <- function(samp1, samp2, results, baseout, padj.thresh = 0.01, lfc.thresh = 2, top.n = NULL) {
#      base <- paste0(baseout, samp1, "v", samp2, "/padj.01.lfc2.read100/")
#      dir.create(file.path(base, "AllGenes"), showWarnings = FALSE, recursive = TRUE)
#      dir.create(file.path(base, paste0(samp1, "up")), showWarnings = FALSE, recursive = TRUE)
#      dir.create(file.path(base, paste0(samp2, "up")), showWarnings = FALSE, recursive = TRUE)
#  
#      One.Two <- subset(results, padj <= padj.thresh & abs(log2FoldChange) >= lfc.thresh)
#      One.up <- subset(results, padj <= padj.thresh & log2FoldChange >= lfc.thresh)
#      Two.up <- subset(results, padj <= padj.thresh & log2FoldChange <= -lfc.thresh)
#  
#      One.Two.terms <- RunEnrichr(rownames(One.Two), outdir = paste0(base, "AllGenes"))
#      VizEnrichments(One.Two.terms, outdir = paste0(base, "AllGenes"))
#  
#      One.up.terms <- RunEnrichr(rownames(One.up), outdir = paste0(base, samp1, "up"))
#      VizEnrichments(One.up.terms, outdir = paste0(base, samp1, "up"))
#  
#      Two.up.terms <- RunEnrichr(rownames(Two.up), outdir = paste0(base, samp2, "up"))
#      VizEnrichments(Two.up.terms, outdir = paste0(base, samp2, "up"))
#  
#      # If top hits set, do it for them as well.
#      if (!is.null(top.n)) {
#          dir.create(file.path(base, "AllGenes", paste0("Top", top.n)), showWarnings = FALSE, recursive = TRUE)
#          dir.create(file.path(base, paste0(samp1, "up/", "Top", top.n)), showWarnings = FALSE, recursive = TRUE)
#          dir.create(file.path(base, paste0(samp2, "up/", "Top", top.n)), showWarnings = FALSE, recursive = TRUE)
#  
#          # Top 100 for each category.
#          One.Two.top.n <- One.Two[order(-(abs(One.Two$log2FoldChange))),]
#          One.Two.top.n <- One.Two.top.n[1:top.n,]
#  
#          One.Two.top.n.terms <- RunEnrichr(rownames(One.Two.top.n), outdir = paste0(base, "AllGenes/Top", top.n))
#          VizEnrichments(One.Two.top.n.terms, outdir = paste0(base, "AllGenes/Top", top.n))
#  
#          One.up.top.n <- One.Two[order(-One.Two$log2FoldChange),]
#          One.up.top.n <- One.up.top.n[1:top.n,]
#  
#          One.up.top.n.terms <- RunEnrichr(rownames(One.up.top.n), outdir = paste0(base, samp1, "up/Top", top.n))
#          VizEnrichments(One.up.top.n.terms, outdir = paste0(base, samp1, "up/Top", top.n))
#  
#          Two.up.top.n <- One.Two[order(One.Two$log2FoldChange),]
#          Two.up.top.n <- Two.up.top.n[1:top.n,]
#  
#          Two.up.top.n.terms <- RunEnrichr(rownames(Two.up.top.n), outdir = paste0(base, samp2, "up/Top", top.n))
#          VizEnrichments(Two.up.top.n.terms, outdir = paste0(base, samp2, "up/Top", top.n))
#      }
#  }

## --------------------------------------------------------------------------
#  # R
#  # Run above function for each comparison.
#  baseout <- "./RNA_Seq/Final_Analyses/EnrichmentTablesFigures/"
#  DoEnrichments("CLL", "NORM", CLLvNORM, baseout, top.n = 500)
#  DoEnrichments("CLL", "DL", CLLvDL, baseout, top.n = 500)
#  DoEnrichments("CLL", "FL", CLLvFL, baseout, top.n = 500)
#  DoEnrichments("DL", "FL", DLvFL, baseout, top.n = 500)
#  DoEnrichments("DL", "NORM", DLvNORM, baseout, top.n = 500)
#  DoEnrichments("FL", "NORM", FLvNORM, baseout, top.n = 500)

## --------------------------------------------------------------------------
#  # R
#  write.table(assay(vsd), file = "./DESeq2_out/vsd.counts.txt", sep = "\t", quote = FALSE)
#  write.table(assay(rld), file = "./DESeq2_out/rld.counts.txt", sep = "\t", quote = FALSE)
#  
#  # Retrieve the count/size factor normalized counts.
#  counts <- counts(dds, normalized = TRUE)
#  
#  write.table(counts, file = "./DESeq2_out/normalized.counts.txt", sep = "\t", quote = FALSE)
#  # And log2 transformed.
#  write.table(assay(normTransform(dds)), file = "./DESeq2_out/log2.normalized.counts.txt", sep = "\t", quote = FALSE)

## --------------------------------------------------------------------------
#  # R
#  suppressPackageStartupMessages(library(ChIPQC))
#  suppressPackageStartupMessages(library(DiffBind))
#  suppressPackageStartupMessages(library(BiocParallel))
#  
#  register(SerialParam())

## ---- results="hide"-------------------------------------------------------
#  # R
#  h3ac.exp <- suppressMessages(suppressWarnings(ChIPQC("./ChIP_Seq/ChIPQC/SampleSheet_H3AC_QC_narrow.NoCTRL.csv", chromosomes = "chr14", consensus = TRUE)))
#  k27ac.exp <- suppressMessages(suppressWarnings(ChIPQC("./ChIP_Seq/ChIPQC/SampleSheet_H3K27AC_QC_narrow.NoCTRL.csv", chromosomes = "chr14", consensus = TRUE)))
#  k4me1.exp <- suppressMessages(suppressWarnings(ChIPQC("./ChIP_Seq/ChIPQC/SampleSheet_H3K4ME1_QC_narrow.NoCTRL.csv", chromosomes = "chr14", consensus = TRUE)))
#  faire.exp <- suppressMessages(suppressWarnings(ChIPQC("./ChIP_Seq/ChIPQC/SampleSheet_FAIRE_QC_narrow.NoCTRL.csv", chromosomes = "chr14", consensus = TRUE)))

## --------------------------------------------------------------------------
#  # R
#  # Save metrics.
#  write.table(as.data.frame(QCmetrics(h3ac.exp)), file = "./ChIP_Seq/ChIPQC/H3AC.QCmetrics.txt", quote = FALSE, sep = "\t")
#  write.table(as.data.frame(QCmetrics(k27ac.exp)), file = "./ChIP_Seq/ChIPQC/H3K27AC.QCmetrics.txt", quote = FALSE, sep = "\t")
#  write.table(as.data.frame(QCmetrics(k4me1.exp)), file = "./ChIP_Seq/ChIPQC/H3K4ME1.QCmetrics.txt", quote = FALSE, sep = "\t")
#  write.table(as.data.frame(QCmetrics(faire.exp)), file = "./ChIP_Seq/ChIPQC/FAIRE.QCmetrics.txt", quote = FALSE, sep = "\t")

## --------------------------------------------------------------------------
#  sessionInfo()

