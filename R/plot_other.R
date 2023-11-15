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


plot_beta_star = function(x, types=get_types(x)) {
  lapply(types, function(tid) {
    get_params(x, what="nmf", type=tid)[[1]]$beta_star %>%
      as.data.frame() %>%
      wide_to_long(what="beta") %>%
      reformat_contexts(what=tid) %>%
      dplyr::mutate(type=tid)
  }) %>% do.call(rbind, .) %>%
    plot_signatures_aux()
}


plot_alpha_star = function(x, types=get_types(x)) {
  lapply(types, function(tid) {
    get_params(x, what="nmf", type=tid)[[1]]$alpha_star %>%
      as.data.frame() %>%
      wide_to_long(what="exposures") %>%
      dplyr::mutate(type=tid)
  }) %>% do.call(rbind, .) %>%
    plot_exposures_aux()
}


plot_beta_weights = function(x, types=get_types(x)) {
  beta_w = get_beta_weights(x, types=types)
  if (is.null(beta_w)) return(NULL)

  beta_w %>%
    ggplot() +
    geom_point(aes(x=sigs, y=sigid, size=value)) +
    facet_grid(~type, scales="free") +
    theme_bw()
}


