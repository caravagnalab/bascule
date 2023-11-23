#' Report of the fit
#'
#' @param x
#'
#' @import ggplot2
#'
#' @return assembled plots

plot_fit = function(x) {
  omega = plot_beta_weights(x)
  if (is.null(omega)) ncols = 1 else ncols = 2
  plots = list(
    expos = plot_exposures(x),
    sigs = plot_signatures(x),
    muts = (plot_data(x, reconstructed=TRUE)+labs(title="Reconstructed")) %>%
      patchwork::wrap_plots(plot_data(x, reconstructed=FALSE)+labs(title="GT"), ncol=ncols)
  )
  if (!is.null(omega)) {
    plots[["omega"]] = omega
    design = "AABB\nCCBB\nDDBB"
  } else { design = "AABB\nCCBB\nCCBB" }

  return(patchwork::wrap_plots(plots, design=design))
}
