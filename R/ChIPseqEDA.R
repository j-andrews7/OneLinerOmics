#' Run ChIPQC on a sample sheet
#'
#' This function simply runs \code{\link[ChIPQC]{ChIPQC}} on a sample sheet and
#' returns QC metrics. It generates a consensus peakset between all sets by 
#' retaining and merging all peaks that overlap in at least two samples before
#' defining the fraction of reads in peaks and other metrics.
#'
#' A basic plot showing overlap between peaksets will also be generated, which 
#' can be useful for defining how many samples a peak should be found in to be
#' included in the consensus peakset as defined in \code{\link{RunDiffBind}}.
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
#' @importFrom DiffBind dba.overlap 
#' @importFrom grDevices pdf dev.off
#' @importFrom graphics plot
#' 
#' @export
#'
#' @author Jared Andrews
#'
#' @seealso \code{\link[ChIPQC]{ChIPQC}}
#'
RunChIPQC <- function(outpath, samplesheet, chromosomes = "chr18") {

  register(SerialParam())
  exp <- suppressMessages(suppressWarnings(ChIPQC(samplesheet, 
    chromosomes = chromosomes, consensus = TRUE, bCount = FALSE)))

  exp.all <- suppressMessages(suppressWarnings(ChIPQC(samplesheet, 
    chromosomes = chromosomes, consensus = FALSE, bCount = FALSE)))
  metrics <- as.data.frame(QCmetrics(exp))

  db <- exp.all@DBA

  olap.rate <- dba.overlap(db, mode = DBA_OLAP_RATE)
  pdf(paste0(outpath, "/Overlaps.pdf"))
  p <- plot(olap.rate, type = "b", ylab='# peaks', 
    xlab='Overlap at least this many peaksets')
  print(p)
  dev.off()
  
  write.table(metrics, file = paste0(outpath, "/QCmetrics.txt"), quote = FALSE,
    sep = "\t")
  
  return(metrics)
}