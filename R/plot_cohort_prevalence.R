#' Title
#'
#' @param x
#'
#' @return
#' @export
#'
#' @examples
plot_cohort_prevalence <- function(x) {

  alpha <- get_exposure(x, long = TRUE) %>%
    dplyr::group_by(Signature) %>%
    dplyr::summarise(Exposure = sum(Exposure)) %>%
    dplyr::mutate(Exposure = Exposure/sum(Exposure)) %>%
    dplyr::arrange(dplyr::desc(Exposure))

  alpha$Signature = factor(alpha$Signature,
                           levels = alpha$Signature)


  ggplot2::ggplot(
    data = alpha,
    ggplot2::aes(x="Overall", y=Exposure, fill=Signature)
  ) +
    ggplot2::geom_bar(stat = "identity") +
    my_ggplot_theme() +
    ggplot2::scale_y_continuous(labels=scales::percent) +
    ggplot2::scale_fill_manual(values = get_signature_colors(x)) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 90, hjust = 1)
    ) +
    ggplot2::labs(
      title = paste0(x$cohort, ' (n = ', x$n_samples, ')'),
      y = "Overall exposure",
      x = ''
    ) +
    ggplot2::guides(
      fill = ggplot2::guide_legend(
        nrow = ifelse(x$n_catalogue + x$n_denovo > 8, 2, 1))
    )

  return(plt)
}
