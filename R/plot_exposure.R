#' plot exposure matrix
#'
#' @description creates bar plot of relative exposure matrix, where x-axis are samples and y-axis are their relative contribution.
#' @param x basilica object
#'
#' @return plot
#' @export plot_exposure
#'
#' @examples
plot_exposure <- function(x, sort_by = NULL, thr=0.1) {

  alpha <- get_exposure(x, long = TRUE)

  # plt <- .plot_exposure(x = alpha)
  if ("groups" %in% colnames(alpha))
    samples_order = alpha %>% dplyr::arrange(groups) %>% dplyr::pull(Sample) %>% unique() else
    samples_order = alpha$Sample %>% unique %>% gtools::mixedsort()

  if(!is.null(sort_by))
  {
    samples_order = alpha %>%
      dplyr::filter(Signature == sort_by) %>%
      dplyr::arrange(dplyr::desc(Exposure)) %>%
      dplyr::pull(Sample)
  }

  alpha$Sample = factor(
    alpha$Sample,
    levels = samples_order
  )

  alpha = alpha %>%
    dplyr::mutate(Signature=ifelse(Exposure < thr, "Other", Signature))

  caption = paste0("Sorted by ", sort_by)
  if(is.null(sort_by)) caption = "Sorted by sample"

  keep = alpha$Signature %>% unique()

  other_col = list("Other"="gainsboro") %>% unlist()

  plt = ggplot2::ggplot(
    data = alpha,
    ggplot2::aes(x=Sample, y=Exposure, fill=Signature)
  ) +
    ggplot2::geom_bar(stat = "identity") +
    my_ggplot_theme() +
    ggplot2::scale_y_continuous(labels=scales::percent) +
    ggplot2::scale_fill_manual(values = c(get_signature_colors(x),other_col),
                               breaks = keep) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 90, hjust = 1)
    ) +
    ggplot2::labs(
      title = paste0(x$cohort, ' (n = ', x$n_samples, ')'),
      caption = caption
    )
  # ggplot2::guides(
  #   fill = ggplot2::guide_legend(
  #     nrow = ifelse(x$n_catalogue + x$n_denovo > 8, 2, 1))
  #   )

  if ("groups" %in% names(x$fit))
    plt = plt + ggplot2::facet_grid(groups~., scales="free_x")


  return(plt)
}
