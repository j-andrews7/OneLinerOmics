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
CreateOutputStructure <- function(block, level, base) {

  # Generate folders and craft design.
  if (!is.null(block)) {
    design <- formula(paste("~", paste(c(block, level), sep = "", 
      collapse = " + ")))
    if (!dir.exists(file.path(base, paste0(block, ".Block")))) {
      dir.create(file.path(base, paste0(block, ".Block")))
    }
    if (!dir.exists(file.path(base, paste0(block, ".Block/GenericFigures")))) {
      dir.create(file.path(base, paste0(block, ".Block/GenericFigures")))
    }
    if (!dir.exists(file.path(base, paste0(block, ".Block/Robjects")))) {
      dir.create(file.path(base, paste0(block, ".Block/Robjects")))
    }
    if (!dir.exists(file.path(base, paste0(block, ".Block/GeneBoxPlots")))) {
      dir.create(file.path(base, paste0(block, ".Block/GeneBoxPlots")))
    }
    base <- file.path(base, paste0(block,".Block"))
  } else {
    design <- formula(paste("~", level))
    if (!dir.exists(file.path(base, "NoBlock"))) {
      dir.create(file.path(base, "NoBlock"))
    }
    if (!dir.exists(file.path(base, "NoBlock/GenericFigures"))) {
      dir.create(file.path(base, "NoBlock/GenericFigures"))
    }
    if (!dir.exists(file.path(base, "NoBlock/Robjects"))) {
      dir.create(file.path(base, "NoBlock/Robjects"))
    }
    if (!dir.exists(file.path(base, "NoBlock/GeneBoxPlots"))) {
      dir.create(file.path(base, "NoBlock/GeneBoxPlots"))
    }
    base <- file.path(base, "NoBlock")
  }

  return(list(base = base, design = design))
}