
#' plot exposure matrix
#'
#' @description creates bar plot of relative exposure matrix, where x-axis are samples and y-axis are their relative contribution.
#' @param x basilica object
#'
#' @return plot
#' @export plot_exposure
#'
#' @examples
plot_exposure <- function(x) {

  alpha <- get_exposure(x, long = FALSE)

  plt <- basilica:::.plot_exposure(x = alpha)

  return(plt)
}

# ------------------------------------------------------------------------------

#' plot signatures
#'
#' @description creates bar plot of inferred signature profiles, where x-axis are 96 substitution bases and y-axis are their relative contribution.
#'
#' @param x basilica object
#' @param denovoSignature if TRUE, plots inferred de-novo signatures, otherwise plots inferred catalogue signatures
#'
#' @return plot
#' @export plot_signatures
#' @examples
plot_signatures <- function(x, denovo = TRUE ) {

  if (denovo==TRUE) {
    beta <- get_denovo_signatures(x)
  } else {
    beta <- get_catalogue_signatures(x)
  }

  plt <- basilica:::.plot_signatures(beta)

  return(plt)
}

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



