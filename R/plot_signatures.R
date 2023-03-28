#' plot signatures
#'
#' @description creates bar plot of inferred signature profiles, where x-axis are 96 substitution bases and y-axis are their relative contribution.
#'
#' @param Type
#' @param highlight
#' @param x basilica object
#'
#' @return plot
#' @export plot_signatures
#' @examples
plot_signatures <- function(x, Type = c("De novo", "Catalogue"), highlight = 0.05)
{
  sigs = x %>% get_signatures(long = TRUE) %>%
    dplyr::filter(Type %in% !!Type)

  return(plot_signatures_sbs(x, sigs, highlight))
}

plot_signatures_sbs = function(x, sigs, highlight = 0.05)
{
  # Remove parenthesis
  sigs$substitution = stringr::str_extract_all(sigs$Feature, "\\[[^\\]\\[]*]") %>% unlist()
  sigs$substitution = gsub(pattern = '\\[', replacement = '', x = sigs$substitution)
  sigs$substitution = gsub(pattern = '\\]', replacement = '', x = sigs$substitution)

  # Detect highlights
  sigs = sigs %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      context = paste0(substr(Feature, 1, 1), '_', substr(Feature, 7, 7)),
    ) %>%
    dplyr::group_by(Signature) %>%
    dplyr::mutate(
      cutoff = !!highlight * sum(Value),
      highlight = ifelse(
        (Value > !!highlight * sum(Value)),
        paste0("probability mass above ", !!highlight * 100, '%'),
        NA
      )
    )

  # Extract a crazy map to colour the reference nucleotide
  sigs$ref = substr(sigs$substitution, 1, 1)

  sigs = sigs %>%
    rowwise() %>%
    mutate(crazy_map =
             gsub(
               '_',
               paste0("<span style='color:indianred3'>", ref, "</span>"),
               context
             ))

  # Nice plot
  plt = ggplot2::ggplot(sigs) +
    ggplot2::geom_hline(
      data = sigs %>% dplyr::distinct(Signature, cutoff),
      ggplot2::aes(yintercept = cutoff),
      size = .2,
      linetype = 'dashed',
      color = 'black'
    )  +
    ggplot2::geom_bar(
      ggplot2::aes(x = crazy_map, y = Value, fill = Signature, color = highlight),
      stat="identity",
      position="identity") +
    ggplot2::facet_grid(as.character(Signature)~substitution, scales="free", drop = FALSE, space="free_x") +
    my_ggplot_theme() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 90, hjust = 1)
    ) +
    ggplot2::labs(
      y = "Frequency",
      title = paste0(x$cohort, ' (n = ', x$n_samples, ')')
    ) +
    # ggplot2::scale_fill_manual(values = get_signature_colors(x) %>% ggplot2::alpha(0.7)) +
    ggplot2::scale_color_manual(values = "black", na.value = NA, na.translate = F) +
    ggplot2::guides(fill = ggplot2::guide_legend(
      paste0("p >", highlight * 100, '%')
    )) +
    # ggplot2::guides(
    #   fill = ggplot2::guide_legend(
    #     nrow = ifelse(x$n_catalogue + x$n_denovo > 8, 2, 1)),
    #   color = ggplot2::guide_legend("", override.aes = ggplot2::aes(fill = NA))
    # ) +
    ggplot2::theme(axis.text.x = ggtext::element_markdown(angle = 90, hjust = 0))

  return(plt)
}
