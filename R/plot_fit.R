#' Report of the fit
#'
#' @param x
#'
#' @import ggplot2
#'
#' @return assembled plots

plot_fit = function(x) {
  omega = plot_beta_weights(x)
  centroids = plot_centroids(x)
  mixing_prop = plot_mixture_weights(x)
  if (is.null(omega)) ncols = 1 else ncols = 2
  plots = list(
    expos = plot_exposures(x),
    sigs = plot_signatures(x),
    muts = plot_data(x, reconstructed=FALSE)
  )
  if (!is.null(omega)) {
    plots[["omega"]] = omega
    design = "AABB\nCCBB\nDDBB"
  } else if (have_groups(x)) {
    plots[["centroids"]] = centroids
    plots[["mixing_prop"]] = mixing_prop
    design = "AADE\nBBBB\nCCCC\nCCCC" } else {
      design = "AABB\nBBBB\nCCCC" }

  return(patchwork::wrap_plots(plots, design=design, guides="collect"))
}
