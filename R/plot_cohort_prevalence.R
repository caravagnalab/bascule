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


  plt = ggplot2::ggplot(
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

  alpha <- get_exposure(x, long = TRUE)

  plt2 = ggplot(alpha,
         aes(x = Signature, y = Exposure, color = Signature))+
    geom_jitter(size = .5, alpha = .6) +
    geom_boxplot(color = 'black', width = .3) +
    # scale_fill_manual(values = get_signature_colors(x)) +
    scale_color_manual(values = get_signature_colors(x)) +
    my_ggplot_theme() +
    geom_hline(
      yintercept = x$params$phi,
      color = 'indianred3',
      linetype = 'dashed'
    ) +
    labs(title = "Per sample exposure") +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 90, hjust = 1)
    )

  alpha <- get_exposure(x, long = TRUE)

  seq_cuts = seq(0, 1, by = 0.01)

  regression_data = lapply(seq_cuts,
         function(s)
         {
           alpha %>%
             filter(Exposure >= s) %>%
             group_by(Signature) %>%
             summarise(n = n(), cut = s)
         }
  ) %>%
    Reduce(f = bind_rows)

  pl3 = regression_data %>%
    ggplot() +
    geom_line(
      aes(x = cut, y = n, color = Signature)
    ) +
    geom_point(
      data = regression_data %>% group_by(Signature) %>% filter(n == min(n)) %>% filter(row_number() == n()),
      aes(x = cut, y = n, color = Signature)
    ) +
    scale_color_manual(values = get_signature_colors(x)) +
    my_ggplot_theme()

  pl4 = pl3 +
    scale_y_log10()

  plt = ggpubr::ggarrange(
    plt,
    plt2,
    ncol = 2,
    widths = c(1, 2),
    common.legend = TRUE,
    legend = 'bottom'
  )

  plt34 = ggpubr::ggarrange(pl3, pl4, ncol = 2)

  plt = ggpubr::ggarrange(plt, plt34, ncol = 1,
                    common.legend = TRUE,
                    legend = 'bottom')


  return(plt)
}
