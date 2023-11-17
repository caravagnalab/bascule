#' Report of the fit
#'
#' @param x
#'
#' @import ggplot2
#'
#' @return assembled plots

plot_fit = function(x) {
  plots = list(
    expos = plot_exposures(x),
    sigs = plot_signatures(x),
    muts = (plot_data(x, reconstructed=TRUE)+labs(title="Reconstructed")) %>%
      patchwork::wrap_plots(plot_data(x, reconstructed=FALSE)+labs(title="GT"))
  )

  omega = plot_beta_weights(x)
  if (!is.null(omega)) plots[["omega"]] = omega

  design = "AABB\nCCBB\nDDBB"
  return(patchwork::wrap_plots(plots, design=design))
}
