



my_ggplot_theme = function ()
{
  ggplot2::theme_light(base_size = 10) +
    ggplot2::theme(legend.position = "bottom",
                   legend.key.size = ggplot2::unit(0.3, "cm"),
                   panel.background = ggplot2::element_rect(fill = "white")
                   )
}


get_signature_colors = function(x)
{
  dn = x$fit$denovo_signatures %>% rownames
  sc = x$fit$catalogue_signatures %>% rownames

  dn_c = ggsci::pal_nejm()(x$n_denovo)
  names(dn_c) = dn

  sc_c = ggsci::pal_simpsons()(x$n_catalogue)
  names(sc_c) = sc

  c(dn_c, sc_c)
}


# ------------------------------------------------------------------------------



# ------------------------------------------------------------------------------

plot_similarity <- function(denovo, reference, limit) {

  cosine_matrix <- cosine.matrix(reference, denovo)

  cosine_matrix <- cosine_matrix %>% dplyr::filter_all(dplyr::any_vars(. > limit))
  cosine_matrix <- tibble::rownames_to_column(cosine_matrix, var="reference")
  cosine_long <- cosine_matrix %>% tidyr::pivot_longer(cols = -c(reference), names_to = 'denovo', values_to = 'cosine')

  p <- ggplot(cosine_long, aes(x=reference, y=cosine)) +
    geom_point() +
    facet_grid(denovo ~ .) +
    theme(axis.text.x=element_text(angle=-90, vjust=0.5, hjust=0)) +
    ggtitle("Denovo vs Reference (Similarity)") +
    ylab("Cosine Similarity") +
    geom_hline(yintercept=c(limit), linetype='dashed', color='red') +
    theme(axis.title.x=element_blank())

  return(p)
}



