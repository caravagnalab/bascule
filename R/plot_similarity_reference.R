#' Function to visualize the similarity of denovo and reference signatures.
#'
#' @param x Basilica object.
#' @param similarity_cutoff add
#' @param reference External reference catalogue to compare the denovo with.
#' @param by_subs Logical. If set to \code{TRUE}, the similarity is computed
#'                separately for each substitution.
#' @param context Logical. If set to \code{TRUE}, the context labels are
#'                reported on the x axis.
#' @param add_pheatmap Logical. If set to \code{TRUE}, the heatmap with the
#'                     similarity values among signatures will be reported.
#'
#' @return add
#' @export plot_similarity_reference

plot_similarity_reference <- function(x, reference = NULL, similarity_cutoff = 0.4,
                                      by_subs = FALSE, context = T, add_pheatmap = T,
                                      similarity = "cosine") {

  if (!is.null(reference))
    x$input$reference_catalogue = reference

  if (similarity == "cosine") {
    if (!by_subs)
      # Similarity to the reference
      cosine_matrix <- cosine.matrix(
        x$input$reference_catalogue,
        x %>% get_signatures()
      ) else
        cosine_matrix = cosine.matrix(x$input$reference_catalogue,
                                      x %>% get_signatures(),
                                      substitutions = get_contexts(x)$subs %>% unique())
  } else if (similarity == "KL") {
    cosine_matrix = sapply(rownames(x$input$reference_catalogue), function(s1) {
      sapply(rownames(get_signatures(x)), function(s2) {
        sigs = as.matrix(rbind(x$input$reference_catalogue[s1,],
                               get_signatures(x)[s2,]))
        suppressMessages(philentropy::KL(sigs, unit="log"))
      }) %>% setNames(rownames(get_signatures(x)))
    })
  } else { cli::cli_alert_danger("Parameter `similarity` must be either `cosine` or `KL`"); stop() }

  cosine_matrix = cosine_matrix[sort(rownames(cosine_matrix)), sort(colnames(cosine_matrix))]

  # Nice colors
  color_gradient = (RColorBrewer::brewer.pal(10, 'Spectral')) %>% rev
  color_gradient[1:5] = color_gradient[1:5] %>% ggplot2::alpha(0.7)
  color_breaks = seq(0, 1, 0.1)

  if (similarity == "KL") color_breaks = seq(0, max(cosine_matrix), length.out=10)

  # Numbers where worth
  numbers = cosine_matrix %>% round(2)
  numbers[numbers < similarity_cutoff] = ''

  # Blacklist
  if(!is.null(x$iterations$blacklist) && (x$iterations$blacklist %>% sapply(length) %>% sum() > 0))
  {
    BList = x$iterations$blacklist %>%
      seq_along() %>%
      lapply(function(i){
        x$iterations$blacklist[[i]] %>%
          as.data.frame() %>%
          mutate(iteration = i)
      }) %>%
      Reduce(f = bind_rows)

    colnames(BList)[1] = "Signature"

    BList = BList %>%
      group_by(Signature) %>%
      arrange(iteration) %>%
      filter(row_number() == 1) %>%
      as.data.frame()

    rownames(BList) = BList$Signature
    BList = BList[, -2, drop = FALSE]
    colnames(BList) = 'blacklist'

    color_BList = colorRampPalette(c("gray", "black"))
    color_BList = color_BList(BList$blacklist %>% unique %>% length())
    names(color_BList) = BList$blacklist %>% unique %>% sort
  } else {
    BList = color_BList = NULL
  }

  if (!is.null(x$iterations)) {
    # Entry list -- exists
    EList = x$iterations$ICS %>%
      seq_along() %>%
      lapply(function(i){
        x$iterations$ICS[[i]] %>%
          rownames() %>%
          as.data.frame() %>%
          mutate(iteration = i)
      }) %>%
      Reduce(f = bind_rows)

    colnames(EList)[1] = "Signature"

    EList = EList %>%
      group_by(Signature) %>%
      arrange(iteration) %>%
      filter(row_number() == 1) %>%
      as.data.frame()

    rownames(EList) = EList$Signature
    EList = EList[, -2, drop = FALSE]
    colnames(EList) = 'detection'
  } else EList = NULL

  if(
    !is.null(EList) &&
    RColorBrewer::brewer.pal.info["Greens", ]$maxcolors >=
    EList$detection %>% unique %>% length()
    )
  {
    color_EList = RColorBrewer::brewer.pal(
      EList$detection %>% unique %>% length(),
      "Greens"
    )
    names(color_EList) = EList$detection %>% unique %>% sort
  } else if (!is.null(EList)) {
    tmp_g = RColorBrewer::brewer.pal(3, "Greens")

    color_EList = colorRampPalette(c(tmp_g[1], tmp_g[3]))
    color_EList = color_EList(EList$detection %>% unique %>% length())
    names(color_EList) = EList$detection %>% unique %>% sort
  }

  cluster_rows = cluster_cols = TRUE
  if (dim(cosine_matrix)[1]==1)
    cluster_rows = FALSE
  if (dim(cosine_matrix)[2]==1)
    cluster_cols = FALSE

  # The world is a better place now that I can pheatmap -> ggplot
  ggp = pheatmap::pheatmap(
    mat = cosine_matrix,
    cluster_rows = cluster_rows,
    cluster_cols = cluster_cols,
    color = color_gradient,
    breaks = color_breaks,
    border_color = 'white',
    cellwidth = 25,
    cellheight = 15,
    annotation_row = BList,
    annotation_col = EList,
    annotation_colors = list(blacklist = color_BList, detection = color_EList),
    display_numbers = numbers
  ) %>% ggplotify::as.ggplot()

  # De novo comparisons
  if(x$n_denovo > 0)
  {
    cosine_matrix_dn = tibble::as_tibble(cosine_matrix)[, x %>% get_denovo_signatures() %>% rownames()] %>%
      apply(2, which.max)

    cosine_matrix_dnm = tibble::as_tibble(cosine_matrix[, x %>% get_denovo_signatures() %>% rownames()]) %>%
      apply(2, max)

    extra_plots = NULL

    for(i in 1:length(cosine_matrix_dn))
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
        ggplot2::facet_wrap(~substitution, nrow = 1) +
        my_ggplot_theme() +
        ggplot2::theme(
          axis.text.x = ggplot2::element_text(angle = 90, hjust = 0),
          axis.text.y = ggplot2::element_blank()
        ) +
        ggplot2::ylim(-max_range, max_range) +
        ggplot2::scale_fill_manual(values = col) +
        ggplot2::labs(
          title = bquote(
            .(names(col)[1])~'vs'~
              .(names(col)[2])~"cosine similarity"~theta~ '='~.(cosine_matrix_dnm[i]))
        )

      if(context == F){ plt = plt + theme(axis.ticks.x = element_blank(),axis.text.x = element_blank()) + labs(x = "")}

      extra_plots = append(extra_plots, list(plt))
    }

    extra_plots = ggpubr::ggarrange(
      plotlist = extra_plots,
      ncol = 1
    )

    if (add_pheatmap) {
    plot = ggpubr::ggarrange(
        plotlist = list(ggp, extra_plots),
        nrow = 1,
        ncol = 2
      )
    } else{
     plot = extra_plots
    }

  } else
    plot = ggp

  return(plot)
}

