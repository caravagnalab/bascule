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

plot_similarity_reference = function(x, reference=NULL, type="SBS", similarity_cutoff=0.8,
                                     context=T, add_pheatmap=T) {

  if (is.null(reference)) reference = get_fixed_signatures(x, types=type, matrix=T)[[type]]
  reference = reference[!rownames(reference) %in% get_signames(x, types=type)[[type]], ]
  if (nrow(reference) == 0) return(NULL)
  signatures = get_signatures(x, types=type, matrix=T)[[type]]
  denovo_sigs = get_denovo_signatures(x, types=type, matrix=F)[[type]]
  reference_sigs = get_fixed_signatures(x, types=type, matrix=F)[[type]] %>%
    dplyr::add_row(reference %>% wide_to_long(what="beta"))

  # Similarity to the reference
  cosine_matrix = lsa::cosine(t(rbind(reference, signatures)))
  cosine_matrix = cosine_matrix[rownames(reference), rownames(signatures)]

  # Nice colors
  color_gradient = (RColorBrewer::brewer.pal(10, "Spectral")) %>% rev
  color_gradient[1:5] = color_gradient[1:5] %>% ggplot2::alpha(0.7)
  color_breaks = seq(0, 1, 0.1)

  # Numbers where worth
  numbers = cosine_matrix %>% round(2)
  numbers[numbers < 0.5] = ''

  cluster_rows = cluster_cols = TRUE
  if (dim(cosine_matrix)[1]==1) cluster_rows = FALSE
  if (dim(cosine_matrix)[2]==1) cluster_cols = FALSE

  # The world is a better place now that I can pheatmap -> ggplot
  ggp = pheatmap::pheatmap(
    mat = cosine_matrix,
    cluster_rows = cluster_rows,
    cluster_cols = cluster_cols,
    color = color_gradient,
    breaks = color_breaks,
    border_color = "white",
    cellwidth = 25,
    cellheight = 15,
    display_numbers = numbers
  ) %>% ggplotify::as.ggplot()

  # De novo comparisons
  if(nrow(denovo_sigs) > 0) {
    denovo_signames = denovo_sigs$sigs %>% unique()
    cosine_matrix_dn = tibble::as_tibble(cosine_matrix)[, denovo_signames] %>%
      apply(2, which.max)

    cosine_matrix_dnm = tibble::as_tibble(cosine_matrix[, denovo_signames]) %>%
      apply(2, max)

    extra_plots = NULL

    colpalette = gen_palette_aux(signames=unique(c(rownames(cosine_matrix), colnames(cosine_matrix))))
    for(i in 1:length(cosine_matrix_dn)) {
      ref = reference_sigs %>%
        dplyr::filter(sigs == rownames(cosine_matrix)[cosine_matrix_dn[i]])

      dn = denovo_sigs %>%
        dplyr::filter(sigs == names(cosine_matrix_dn)[i]) %>%
        dplyr::mutate(value = -1 * value)

      sigs = dplyr::bind_rows(ref, dn) %>% reformat_contexts(what="SBS")

      max_range = sigs$value %>% abs %>% max
      brange = seq(- max_range, max_range, max_range/5) %>% round(3)

      plt = ggplot2::ggplot(sigs) +
        ggplot2::geom_bar(ggplot2::aes(x=context, y=value, fill=sigs),
                          stat="identity", position="identity") +
        ggplot2::facet_wrap(~ variant, nrow=1) +
        ggplot2::theme_bw() +
        ggplot2::theme(axis.text.x=ggplot2::element_text(angle=90, hjust=0),
                       axis.text.y=ggplot2::element_blank()) +
        ggplot2::ylim(-max_range, max_range) +
        ggplot2::scale_fill_manual(values=colpalette) +
        ggplot2::labs(
          title=bquote(
            .(unique(ref$sigs)) ~ "vs" ~
              .(unique(dn$sigs))~"cosine similarity" ~ theta ~ "=" ~ .(cosine_matrix_dnm[i]))
        )

      if(context == F) { plt = plt + theme(axis.ticks.x=element_blank(), axis.text.x=element_blank()) + labs(x="")}

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





# if by_context is TRUE, it computes the cosine similarity by substitution type
cosine.matrix <- function(a, b, substitutions=NULL) {
  # a and b are data.frame

  df <- data.frame(matrix(0, nrow(a), nrow(b)))
  rownames(df) <- rownames(a)
  colnames(df) <- rownames(b)

  cmp = nrow(a) * nrow(b)
  pb <- progress::progress_bar$new(
    format = paste0("  Cosine similarity (n = ", cmp, ") [:bar] :percent eta: :eta"),
    total = cmp,
    clear = FALSE,
    width= 90
  )


  for (i in 1:nrow(a)) {
    denovo <- a[i, ]
    for (j in 1:nrow(b)) {
      ref <- b[j, ]
      pb$tick()

      score <- cosine.vector(denovo, ref, substitutions)
      df[i,j] <- score
    }
  }

  return(df)
}




cosine.vector <- function(a, b, substitutions = NULL) {

  if (is.matrix(a) && nrow(a)>1) a = t(a)
  if (is.matrix(b) && nrow(b)>1) b = t(b)

  if (is.null(substitutions)) {
    if (!identical(colnames(a), colnames(b))) {
      a = a[names(b)]
    }

    numerator <- sum(a * b)
    denominator <- sqrt(sum(a^2)) * sqrt(sum(b^2))
    return(numerator / denominator)
  }

  cosine.tot = 0
  keep_subs = length(substitutions)
  for (ss in substitutions) {
    keep_cols.tmp = grep(ss, colnames(b), value = T)

    if (all(c(a[,keep_cols.tmp], b[,keep_cols.tmp])==0)) {
      keep_subs = keep_subs - 1
      next
    }

    num.tmp = sum(a[,keep_cols.tmp] * b[,keep_cols.tmp])
    denomin.tmp = sqrt(sum(a[,keep_cols.tmp]^2)) * sqrt(sum(b[,keep_cols.tmp]^2))

    if (num.tmp != 0 && denomin.tmp != 0)
      cosine.tot = cosine.tot + (num.tmp / denomin.tmp)
  }

  return(cosine.tot / keep_subs)
}
