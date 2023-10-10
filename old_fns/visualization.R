
#




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



