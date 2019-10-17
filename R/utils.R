#' Create output directory structure
#'
#' \code{CreateOutputStructure} generates the output directories used by
#' \code{\link{RunDESeq2}} and \code{\link{ProcessDEGs}}.
#'
#'
#' @param block String or character vector defining the variable(s) to use to
#'   block for unwanted variance, e.g. batch or technical effects.
#' @param level String defining variable of interest.
#' @param base String defining the base output path.
#' @return A named List containing the modified base output path and design
#"   formula.
#'
#' @importFrom stats formula
#'
#' @author Jared Andrews
#'
CreateOutputStructure <- function(block, level, base) {

  # Generate folders and craft design.
  if (!is.null(block)) {
    design <- formula(paste("~", paste(c(block, level), sep = "",
      collapse = " + ")))
    if (!dir.exists(file.path(base, paste0(block, ".Block")))) {
      dir.create(file.path(base, paste0(block, ".Block")))
    }
    if (!dir.exists(file.path(base, paste0(block, ".Block/EDAFigures")))) {
      dir.create(file.path(base, paste0(block, ".Block/EDAFigures")))
    }
    if (!dir.exists(file.path(base, paste0(block, ".Block/Robjects")))) {
      dir.create(file.path(base, paste0(block, ".Block/Robjects")))
    }
    if (!dir.exists(file.path(base, paste0(block, ".Block/GeneBoxPlots")))) {
      dir.create(file.path(base, paste0(block, ".Block/GeneBoxPlots")))
    }
    if (!dir.exists(file.path(base, paste0(block, ".Block/MAPlots")))) {
      dir.create(file.path(base, paste0(block, ".Block/MAPlots")))
    }
    if (!dir.exists(file.path(base, paste0(block, ".Block/Heatmaps")))) {
      dir.create(file.path(base, paste0(block, ".Block/Heatmaps")))
    }
    if (!dir.exists(file.path(base, paste0(block, ".Block/DEGFigures")))) {
      dir.create(file.path(base, paste0(block, ".Block/DEGFigures")))
    }
    if (!dir.exists(file.path(base, paste0(block, ".Block/ResultsTables")))) {
      dir.create(file.path(base, paste0(block, ".Block/ResultsTables")))
    }
    base <- file.path(base, paste0(block,".Block"))
  } else {
    design <- formula(paste("~", level))
    if (!dir.exists(file.path(base, "NoBlock"))) {
      dir.create(file.path(base, "NoBlock"))
    }
    if (!dir.exists(file.path(base, "NoBlock/EDAFigures"))) {
      dir.create(file.path(base, "NoBlock/EDAFigures"))
    }
    if (!dir.exists(file.path(base, "NoBlock/Robjects"))) {
      dir.create(file.path(base, "NoBlock/Robjects"))
    }
    if (!dir.exists(file.path(base, "NoBlock/GeneBoxPlots"))) {
      dir.create(file.path(base, "NoBlock/GeneBoxPlots"))
    }
    if (!dir.exists(file.path(base, paste0(block, "NoBlock/MAPlots")))) {
      dir.create(file.path(base, paste0(block, "NoBlock/MAPlots")))
    }
    if (!dir.exists(file.path(base, paste0(block, "NoBlock/Heatmaps")))) {
      dir.create(file.path(base, paste0(block, "NoBlock/Heatmaps")))
    }
    if (!dir.exists(file.path(base, "NoBlock/DEGFigures"))) {
      dir.create(file.path(base, "NoBlock/DEGFigures"))
    }
    if (!dir.exists(file.path(base, "NoBlock/ResultsTables"))) {
      dir.create(file.path(base, "NoBlock/ResultsTables"))
    }
    base <- file.path(base, "NoBlock")
  }

  return(list(base = base, design = design))
}
