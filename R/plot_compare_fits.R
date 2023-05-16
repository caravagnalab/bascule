#' Title
#'
#' @param x add
#' @param y add
#' @param similarity_cutoff add
#'
#' @return add
#' @export plot_compare_fits

plot_compare_fits <- function(x, y, similarity_cutoff = 0.4) {
  # Similarity to the reference
  cosine_matrix <- cosine.matrix(x %>% get_signatures(),
                                 y %>% get_signatures())

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

  # Closest matches: x to y and viceverse
  x_y = cosine_matrix %>% as.matrix() %>% reshape2::melt()

  x_y_best = x_y %>% group_by(Var1) %>% filter(value == max(value))
  y_x_best = x_y %>% group_by(Var2) %>% filter(value == max(value))

  extra_plots = NULL

  x_y_best = easypar::run(
    FUN = function(i){
      comparative_plot_sbs(
        x %>% get_signatures(long = TRUE) %>% dplyr::filter(Signature == x_y_best$Var1[i]),
        y %>% get_signatures(long = TRUE) %>% dplyr::filter(Signature == x_y_best$Var2[i])
      ) +
        labs(title = paste0("x (", x_y_best$Var1[i], ')',
                            "vs y (", x_y_best$Var2[i], ') = ',
                            x_y_best$value[i] %>% round(3)))
    },
    PARAMS = lapply(1:nrow(x_y_best), list),
    parallel = FALSE
  )

  y_x_best = easypar::run(
    FUN = function(i){
      comparative_plot_sbs(
        x %>% get_signatures(long = TRUE) %>% dplyr::filter(Signature == y_x_best$Var1[i]),
        y %>% get_signatures(long = TRUE) %>% dplyr::filter(Signature == y_x_best$Var2[i])
      ) +
        labs(title = paste0("x (", y_x_best$Var1[i], ')',
                            "vs y (", y_x_best$Var2[i], ') = ',
                            y_x_best$value[i] %>% round(3)))
    },
    PARAMS = lapply(1:nrow(y_x_best), list),
    parallel = FALSE
  )

  # np = extra_plots %>% length()
  # nc = np %>% sqrt() %>% ceiling()
  # nr = nc
  # if ((nr - 1) * nc >= np)
  #   nr = nr - 1
  #
  # extra_plots = ggpubr::ggarrange(
  #   plotlist = extra_plots,
  #   ncol = nc,
  #   nrow = nr,
  #   common.legend = TRUE,
  #   legend = 'bottom'
  # )

  # ggp = ggpubr::ggarrange(ggp, extra_plots, ncol = 2, widths = c(1, nc))

  return(list(cosine = ggp, x_y = x_y_best, y_x = y_x_best))
}


comparative_plot_sbs = function(x, y)
{
  x = x %>% mutate(Object = 'x')
  y = y %>% dplyr::mutate(Value = -1 * Value) %>% mutate(Object = 'y')

  sigs = bind_rows(x, y)

  # Remove parenthesis
  sigs$substitution = stringr::str_extract_all(sigs$Feature, "\\[[^\\]\\[]*]") %>% unlist()
  sigs$substitution = gsub(
    pattern = '\\[',
    replacement = '',
    x = sigs$substitution
  )
  sigs$substitution = gsub(
    pattern = '\\]',
    replacement = '',
    x = sigs$substitution
  )

  sigs = sigs %>%
    dplyr::rowwise() %>%
    dplyr::mutate(context = paste0(substr(Feature, 1, 1), '_', substr(Feature, 7, 7)), )

  max_range = sigs$Value %>% abs %>% max
  brange = seq(-max_range, max_range, max_range / 5) %>% round(3)

  ggplot2::ggplot(sigs) +
    ggplot2::geom_bar(
      ggplot2::aes(x = context, y = Value, fill = Object),
      stat = "identity",
      position = "identity"
    ) +
    facet_wrap(~ substitution, nrow = 1) +
    my_ggplot_theme() +
    theme(
      axis.text.x = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank(),
      axis.title.x = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_blank()
    ) +
    ylim(-max_range, max_range)
}
