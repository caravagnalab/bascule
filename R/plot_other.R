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


