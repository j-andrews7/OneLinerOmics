#' @importFrom grDevices colorRampPalette
.heatmap_colors <- function(breaks, preset = NULL, custom.colors = NULL) {
  preset <- match.arg(preset, c("BuRd", "OrPu", "BrTe", "PuGr", 4))
  preset.colors <- list(
    c("#053061", "#2166ac", "#f5f5f5", "#b2182b", "#67001f"),
    c("#b35806", "#e08214", "#f5f5f5", "#8073ac", "#542788"), 
    c("#8c510a", "#d8b365", "#f5f5f5", "#5ab4ac", "#01665e"),
    c("#542788", "#8073ac", "#f5f5f5", "#5aae61", "#00441b")
    )

  if (is.null(custom.colors)) {
      custom.colors <- preset.colors[[preset]]
    }

  colors <- colorRampPalette(custom.colors)(n = length(breaks) - 1)

  return(colors)
}