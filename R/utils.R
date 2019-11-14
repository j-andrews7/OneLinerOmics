#' Create output directory structure
#'
#' \code{CreateRNAOutputStructure} generates the output directories used by
#' \code{\link{RunDESeq2}}, \code{\link{ProcessDEGs}}, 
#' \code{\link{RunDiffBind}}, and \code{\link{ProcessDBRs}}.
#'
#' @param block String or character vector defining the variable(s) to use to
#'   block for unwanted variance, e.g. batch or technical effects.
#' @param level String defining variable of interest.
#' @param base String defining the base output path.
#' @param chip Boolean indicating whether output is being created for ChIP-seq.
#' @return A named List containing the modified base output path and design
#'   formula.
#'
#' @importFrom stats formula
#'
#' @author Jared Andrews
#'
CreateOutputStructure <- function(block, level, base, chip = FALSE) {

  design = NULL
  # Generate folders and craft design.
  if (!is.null(block)) {
    if (!dir.exists(file.path(base, paste0(block, ".Block")))) {
      dir.create(file.path(base, paste0(block, ".Block")))
    }

    if (!chip) {
      design <- formula(paste("~", paste(c(block, level), sep = "",
        collapse = " + ")))

      if (!dir.exists(file.path(base, paste0(block, ".Block/GeneBoxPlots")))) {
        dir.create(file.path(base, paste0(block, ".Block/GeneBoxPlots")))
      }
      if (!dir.exists(file.path(base, paste0(block, ".Block/DEGFigures")))) {
        dir.create(file.path(base, paste0(block, ".Block/DEGFigures")))
      }
      if (!dir.exists(file.path(base, paste0(block, ".Block/MAPlots")))) {
        dir.create(file.path(base, paste0(block, ".Block/MAPlots")))
      }
      if (!dir.exists(file.path(base, paste0(block, ".Block/EDAFigures")))) {
        dir.create(file.path(base, paste0(block, ".Block/EDAFigures")))
      }
      if (!dir.exists(file.path(base, paste0(block, ".Block/Heatmaps")))) {
        dir.create(file.path(base, paste0(block, ".Block/Heatmaps")))
      }

    } else {
      for (i in c("DBRFigures", "ConsensusFigures")) {
        if (!dir.exists(file.path(base, paste0(block, ".Block/", i)))) {
          dir.create(file.path(base, paste0(block, ".Block/", i)))
        }
        if (!dir.exists(file.path(base, paste0(block, 
          ".Block/", i, "/MAPlots")))) {

          dir.create(file.path(base, paste0(block, 
            ".Block/", i, "/MAPlots")))
        }
        if (!dir.exists(file.path(base, paste0(block, 
          ".Block/", i, "/PeakAnnotations")))) {

          dir.create(file.path(base, paste0(block, 
            ".Block/", i, "/PeakAnnotations")))
        }
        if (!dir.exists(file.path(base, paste0(block, 
          ".Block/", i, "/SignalBoxPlots")))) {

          dir.create(file.path(base, paste0(block, 
            ".Block/", i, "/SignalBoxPlots")))
        }
        if (!dir.exists(file.path(base, paste0(block, 
          ".Block/", i, "/PCAPlots")))) {

          dir.create(file.path(base, paste0(block, 
            ".Block/", i, "/PCAPlots")))
        }
        if (!dir.exists(file.path(base, paste0(block, ".Block/", i, 
          "/Heatmaps")))) {

          dir.create(file.path(base, paste0(block, ".Block/", i, "/Heatmaps")))
        }
        if (!dir.exists(file.path(base, paste0(block, ".Block/", i, 
          "/VolcanoPlots")))) {

          dir.create(file.path(base, paste0(block, ".Block/", i, 
            "/VolcanoPlots")))
        }
      }
    }

    if (!dir.exists(file.path(base, paste0(block, ".Block/Robjects")))) {
      dir.create(file.path(base, paste0(block, ".Block/Robjects")))
    }
    if (!dir.exists(file.path(base, paste0(block, ".Block/ResultsTables")))) {
      dir.create(file.path(base, paste0(block, ".Block/ResultsTables")))
    }
    if (!dir.exists(file.path(base, paste0(block, ".Block/Enrichments")))) {
      dir.create(file.path(base, paste0(block, ".Block/Enrichments")))
    }
    base <- file.path(base, paste0(block,".Block"))
  } else {
    if (!dir.exists(file.path(base, "NoBlock"))) {
      dir.create(file.path(base, "NoBlock"))
    }

    if (!chip) {
      design <- formula(paste("~", level))

      if (!dir.exists(file.path(base, "NoBlock/GeneBoxPlots"))) {
        dir.create(file.path(base, "NoBlock/GeneBoxPlots"))
      }
      if (!dir.exists(file.path(base, "NoBlock/DEGFigures"))) {
        dir.create(file.path(base, "NoBlock/DEGFigures"))
      }
      if (!dir.exists(file.path(base, "NoBlock/MAPlots"))) {
        dir.create(file.path(base, "NoBlock/MAPlots"))
      }
      if (!dir.exists(file.path(base, "NoBlock/EDAFigures"))) {
        dir.create(file.path(base, "NoBlock/EDAFigures"))
      }
      if (!dir.exists(file.path(base, "NoBlock/Heatmaps"))) {
        dir.create(file.path(base, "NoBlock/Heatmaps"))
      }

    } else {
      for (i in c("DBRFigures", "ConsensusFigures")) {
        if (!dir.exists(file.path(base, "NoBlock", i))) {
          dir.create(file.path(base, "NoBlock", i))
        }
        if (!dir.exists(file.path(base, "NoBlock", i, "SignalBoxPlots"))) {
          dir.create(file.path(base, "NoBlock", i, "SignalBoxPlots"))
        }
        if (!dir.exists(file.path(base, "NoBlock", i, "MAPlots"))) {
          dir.create(file.path(base, "NoBlock", i, "MAPlots"))
        }
        if (!dir.exists(file.path(base, "NoBlock", i, "PeakAnnotations"))) {
          dir.create(file.path(base, "NoBlock", i, "PeakAnnotations"))
        }
        if (!dir.exists(file.path(base, "NoBlock", i, "PCAPlots"))) {
          dir.create(file.path(base, "NoBlock", i, "PCAPlots"))
        }
        if (!dir.exists(file.path(base, "NoBlock", i, "Heatmaps"))) {
          dir.create(file.path(base, "NoBlock", i, "Heatmaps"))
        }
        if (!dir.exists(file.path(base, "NoBlock/", i, "VolcanoPlots"))) {
          dir.create(file.path(base, "NoBlock", i, "VolcanoPlots"))
        }
      }
    }

    if (!dir.exists(file.path(base, "NoBlock/Robjects"))) {
      dir.create(file.path(base, "NoBlock/Robjects"))
    }
    if (!dir.exists(file.path(base, "NoBlock/ResultsTables"))) {
      dir.create(file.path(base, "NoBlock/ResultsTables"))
    }
    if (!dir.exists(file.path(base, "NoBlock/Enrichments"))) {
      dir.create(file.path(base, "NoBlock/Enrichments"))
    }
    base <- file.path(base, "NoBlock")
  }

  return(list(base = base, design = design))
}


#' Save results
#'
#' \code{SaveResults} saves the results for each comparison in \code{results}. 
#' If \code{chip = FALSE}, it will also save normalized gene counts.
#'
#' @param results Either:
#'   \itemize{
#'     \item a named list containing \linkS4class{DESeqResults} objects for 
#'       all comparisons generated by \code{\link{ProcessDEGs}}, or
#'     \item a \code{DBA} object from \code{\link[DiffBind]{dba.analyze}}.
#'   }
#' @param outpath Path to directory to be used for output.
#' @param dds A \linkS4class{DESeqDataSet} object as returned by 
#'   \code{\link[DESeq2]{DESeq}}.
#' @param chip Boolean indicating whether \code{results} are from ChIP-seq 
#'   analysis.
#' @param method Method used for differential binding e.g. \code{DBA_DESEQ2}.
#'
#'   Do not quote this parameter, it is a global variable from \code{DiffBind}.
#'   If a block was used, be sure to include that 
#'   (e.g. \code{DBA_DESEQ2_BLOCK}).
#' @param promoters Scalar vector containing how many basepairs up and 
#'   downstream of the TSS should be used to define gene promoters.
#' @param se Dataframe containing location informatin for super enhancers.
#' @param txdb \code{TxDb} object to use for annotation.
#'
#' @importFrom utils write.table
#' 
#' @export
#'
#' @author Jared Andrews
#'
SaveResults <- function(results, outpath, dds = NULL, chip = FALSE, 
  method = NULL, promoters = NULL, se = NULL, txdb = NULL) {

  if (!chip) {
    out <- paste0(outpath, "/ResultsTables/")

    write.table(counts(dds, normalized = TRUE), file = paste0(out, 
      "NormalizedGeneCounts.txt"), sep = "\t", quote = FALSE)

    write.table(assay(normTransform(dds)), file = paste0(out, 
      "NormalizedGeneCounts.log2.txt"), sep = "\t", quote = FALSE)

    write.table(fpm(dds), file = paste0(out, "GeneCounts.fpm.txt"), 
      sep = "\t", quote = FALSE)

    for (r in seq_along(results)) {
      res <- results[r]
      comp <- names(results[r])
      resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, 
        normalized = TRUE)), by = "row.names", sort = FALSE)
      names(resdata)[1] <- "Gene"


      write.table(resdata, file = paste0(out, comp, 
        ".Results.NormalizedGeneCounts.txt"), row.names = FALSE, 
        sep = "\t", quote = FALSE)
    }
  } else {
    for (i in seq_along(results$contrasts)) {
      out <- paste0(outpath, "/ResultsTables/")
      report = dba.report(results, th = 1, bCalled = TRUE, 
          bCounts = TRUE, method = method, contrast = i)

      peak.anno <- annotatePeak(report, tssRegion = promoters, TxDb = txdb, 
        annoDb = "org.Hs.eg.db")
  
      df <- data.frame(peak.anno)
      dfs <- list("peaks" = df, "ses" = se)
      final <- .categorize_peaks(dfs)

      write.table(final, file = paste0(out, results$contrasts[[i]]$name1, "-v-",
        results$contrasts[[i]]$name2, ".AllPeaks.txt"), sep = "\t", 
        quote = FALSE)
    }
  }
}
