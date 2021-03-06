% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RNAseqViz.R
\name{PlotRNADEGPCAs}
\alias{PlotRNADEGPCAs}
\title{Plot PCAs from variance stabilized counts for differentially expressed genes}
\usage{
PlotRNADEGPCAs(
  res.list,
  rld,
  vsd,
  outpath,
  level,
  plot.annos,
  padj.thresh,
  fc.thresh
)
}
\arguments{
\item{res.list}{Named List containing \linkS4class{DESeqResults} objects for 
all comparisons generated by \code{\link{ProcessDEGs}}.}

\item{rld}{A \linkS4class{RangedSummarizedExperiment} object of 
\code{\link[DESeq2]{rlog}} transformed counts as returned by
\link{RunDESeq2}.}

\item{vsd}{A \linkS4class{RangedSummarizedExperiment} object of 
\code{\link[DESeq2]{vst}} transformed counts as returned by
\link{RunDESeq2}.}

\item{outpath}{Path to directory to be used for output.}

\item{level}{String defining variable of interest.}

\item{plot.annos}{String or character vector defining the column(s) in 
\code{samplesheet} to use to annotate figures.}

\item{padj.thresh}{Number indicating the adjusted p-value 
cutoff to be used for determining "significant" differential expression.}

\item{fc.thresh}{Number indicating the log2 fold-change 
cutoff to be used for determining "significant" differential expression.}
}
\description{
Plot PCAs from variance stabilized counts for differentially expressed genes
}
\author{
Jared Andrews
}
