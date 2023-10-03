my_ggplot_theme = function() {
  ggplot2::theme_light(base_size = 10) +
    ggplot2::theme(legend.position = "bottom",
                   legend.key.size = ggplot2::unit(0.3, "cm"),
                   panel.background = ggplot2::element_rect(fill = "white")
    )
}


get_signature_colors = function(x) {
  dn = x$fit$denovo_signatures %>% rownames
  sc = x$fit$catalogue_signatures %>% rownames

  if (!is.null(dn)) {
    dn_c = ggsci::pal_nejm()(x$n_denovo)
    names(dn_c) = dn
  } else { dn_c = NULL }

  n_cat = x$fit$catalogue_signatures %>% nrow

  sc_c = ggsci::pal_simpsons()(n_cat)
  names(sc_c) = sc

  c(dn_c, sc_c)
}



merge_colors_palette = function(fit1, fit2, ref_catalogue) {
  ref_names = c(get_signames(fit1), get_signames(fit2)) %>% unique()
  ref_colors = COSMIC_color_palette(catalogue=ref_catalogue)[ref_names] %>% purrr::discard(is.na)
  dn_names = setdiff(c(get_dn_signames(fit1), get_dn_signames(fit2)) %>% unique(), names(ref_colors))
  dn_colors = gen_palette(n=length(dn_names)) %>% setNames(dn_names)
  cls = c(ref_colors, dn_colors)
  return(cls)
}



