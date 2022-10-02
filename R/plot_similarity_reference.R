#' Title
#'
#' @param x
#' @param similarity_cutoff
#'
#' @return
#' @export
#'
#' @examples
plot_similarity_reference <- function(x, similarity_cutoff = 0.4) {

  # Similarity to the reference
  cosine_matrix <- cosine.matrix(
    x$input$reference_catalogue,
    x %>% get_signatures()
  )

  # cosine_matrix_ldf = reshape2::melt(cosine_matrix %>% as.matrix()) %>%
  #   dplyr::rename(
  #     Target = Var1,
  #     Reference = Var2,
  #     Cosine = value
  #   ) %>%
  #   dplyr::as_tibble()

  # Nice colors
  color_gradient = (RColorBrewer::brewer.pal(10, 'Spectral')) %>% rev
  color_gradient[1:5] = color_gradient[1:5] %>% ggplot2::alpha(0.7)
  color_breaks = seq(0, 1, 0.1)

  # Numbers where worth
  numbers = cosine_matrix %>% round(2)
  numbers[numbers < similarity_cutoff] = ''

  # The world is a better place now that I can pheatmap -> ggplot
  ggp = pheatmap::pheatmap(
    mat = cosine_matrix,
    color = color_gradient,
    breaks = color_breaks,
    border_color = 'white',
    cellwidth = 25,
    cellheight = 15,
    display_numbers = numbers
  ) %>% ggplotify::as.ggplot()

  # De novo comparisons
  if(x$n_denovo > 0)
  {
    cosine_matrix_dn = cosine_matrix[, x %>% get_denovo_signatures() %>% rownames()] %>%
      apply(2, which.max)

    cosine_matrix_dnm = cosine_matrix[, x %>% get_denovo_signatures() %>% rownames()] %>%
      apply(2, max)

    extra_plots = NULL

    for(i in (cosine_matrix_dn %>% seq()))
    {
      target = x %>% get_reference_signatures(long = TRUE) %>%
        dplyr::filter(Signature == rownames(cosine_matrix)[cosine_matrix_dn[i]])

      reference = x %>%
        get_denovo_signatures(long = TRUE) %>%
        dplyr::filter(Signature == names(cosine_matrix_dn)[i]) %>%
        dplyr::mutate(Value = -1 * Value)

      sigs = bind_rows(target, reference)

      # Remove parenthesis
      sigs$substitution = stringr::str_extract_all(sigs$Feature, "\\[[^\\]\\[]*]") %>% unlist()
      sigs$substitution = gsub(pattern = '\\[', replacement = '', x = sigs$substitution)
      sigs$substitution = gsub(pattern = '\\]', replacement = '', x = sigs$substitution)

      sigs = sigs %>%
        dplyr::rowwise() %>%
        dplyr::mutate(
          context = paste0(substr(Feature, 1, 1), '_', substr(Feature, 7, 7)),
        )

      max_range = sigs$Value %>% abs %>% max
      brange = seq(- max_range, max_range, max_range/5) %>% round(3)

      col = get_signature_colors(x)
      col = c(col[names(cosine_matrix_dn)[i]], 'black')
      names(col)[2] = rownames(cosine_matrix)[cosine_matrix_dn[i]]

      plt = ggplot2::ggplot(sigs) +
        ggplot2::geom_bar(
          ggplot2::aes(x = context, y = Value, fill = Signature),
          stat = "identity",
          position = "identity") +
        facet_wrap(~substitution, nrow = 1) +
        my_ggplot_theme() +
        theme(
          axis.text.x = ggplot2::element_text(angle = 90, hjust = 0),
          axis.text.y = ggplot2::element_blank()
        ) +
        ylim(-max_range, max_range) +
        scale_fill_manual(values = col) +
        labs(
          title = bquote(
            .(names(col)[1])~'vs'~
              .(names(col)[2])~"cosine similarity"~theta~ '='~.(cosine_matrix_dnm[i]))
        )

      extra_plots = append(extra_plots, list(plt))
    }

    extra_plots = ggpubr::ggarrange(
      plotlist = extra_plots,
      ncol = 1
    )

    ggp = ggpubr::ggarrange(
      ggp,
      extra_plots,
      ncol = 2
    )
  }

  return(ggp)
}
