#' Plot exposure matrix
#'
#' @description Creates bar plot of relative exposure matrix, where x-axis are
#'              samples and y-axis are their relative contribution.
#' @param x Basilica object
#' @param sample_name Logical. If set to \code{TRUE}, the sample IDs will be
#'                    reported on the x axis.
#' @param levels add
#' @param cls add
#' @param flip_coord add
#'
#' @return plot
#' @export plot_exposures

plot_exposures = function(x, sample_name=F, sigs_levels=NULL, cls=NULL,
                          flip_coord=F, muts=FALSE, sampleIDs=NULL, sort_by=NULL) {

  if (is.null(sampleIDs)) sampleIDs = rownames(x$fit$exposure)

  b = x$fit$exposure

  if (muts) b = b * rowSums(x$fit$x)

  b = b[sampleIDs,]

  if(is.null(cls) && !have_color_palette(x)) cls = gen_palette(ncol(b)) %>% setNames(colnames(b))
  else if (is.null(cls) && have_color_palette(x)) cls = get_color_palette(x)

  if(is.null(sigs_levels)) sigs_levels = sort(colnames(b))

  idcols = c("sample")
  if (have_groups(x)) {
    idcols = c("sample","groups")
    b$groups = x$groups[rownames(x$fit$exposure) == sampleIDs]
  }

  if (!is.null(sort_by))
    sample_levels = b %>% as.data.frame() %>%
      tibble::rownames_to_column(var="sample") %>%
      reshape2::melt(cols=idcols, variable.name="Signature", value.name="alpha") %>%
      dplyr::filter(Signature==sort_by) %>%
      dplyr::arrange(desc(alpha)) %>% dplyr::pull(sample) else sample_levels = rownames(b)

  p = b %>% as.data.frame() %>%
    tibble::rownames_to_column(var="sample") %>%
    reshape2::melt(id=idcols, variable.name="Signature", value.name="alpha") %>%
    dplyr::mutate(sample=factor(sample, levels=sample_levels)) %>%

    ggplot(aes(x=sample, y=alpha, fill=factor(Signature, levels=sigs_levels))) +
    geom_bar(stat = "identity") +
    ggplot2::scale_fill_manual(values = cls) +
    labs(title = "Exposure", x = "") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90)) +
    guides(fill=guide_legend(title="Signatures")) + ylab("")

  if (!sample_name)
    p = p + theme(axis.ticks.x = element_blank(), axis.text.x = element_blank()) + labs(x = "")

  if (flip_coord)
    p = p + coord_flip()

  if (have_groups(x))
    p = p + facet_grid(~groups, scales="free_x")

  return(p)
}



plot_sigs_prevalence = function(x) {
  expos = get_exposure(x, long = T)

  expos %>%
    dplyr::group_by(Signature) %>%
    dplyr::arrange(Signature, Exposure) %>%
    dplyr::mutate(nnn=1, cs=cumsum(nnn)/x$n_samples) %>%
    ggplot(aes(x=Exposure, y=cs, color=Signature)) +
    geom_point(size=.5) +
    geom_line() +
    scale_color_manual(values=x$color_palette) +
    theme_bw() + ylab("Fraction of samples") +
    xlim(0,1) + ylim(0,1)

}



# plot_exposure <- function(x, labels = NULL,sort_by = NULL, thr=0.1){
#
#
#   alpha <- get_exposure(x, long = TRUE)
#
#   # plt <- .plot_exposure(x = alpha)
#   if ("groups" %in% colnames(alpha))
#     samples_order = alpha %>% dplyr::arrange(groups) %>% dplyr::pull(Sample) %>% unique() else
#     samples_order = alpha$Sample %>% unique %>% gtools::mixedsort()
#
#   if(!is.null(sort_by))
#   {
#     samples_order = alpha %>%
#       dplyr::filter(Signature == sort_by) %>%
#       dplyr::arrange(dplyr::desc(Exposure)) %>%
#       dplyr::pull(Sample)
#   }
#
#
#  alpha = alpha %>%
#     dplyr::mutate(Signature=ifelse(Exposure < thr, "Other", Signature))
#
#   caption = paste0("Sorted by ", sort_by)
#   if(is.null(sort_by)) caption = "Sorted by sample"
#
#   if(!is.null(labels)){
#
#     alpha = alpha %>% dplyr::select(-groups) %>% full_join(labels, by = "Sample")
#
# }
#
#
#   keep = alpha$Signature %>% unique()
#
#   other_col = list("Other"="gainsboro") %>% unlist()
#
#   plt = ggplot2::ggplot(
#     data = alpha,
#     ggplot2::aes(x=factor(Sample,levels = samples_order), y=Exposure, fill=Signature)
#   ) +
#     ggplot2::geom_bar(stat = "identity") +
#     my_ggplot_theme() +
#     ggplot2::scale_y_continuous(labels=scales::percent) +
#     ggplot2::scale_fill_manual(values = c(get_signature_colors(x),other_col),
#                                breaks = keep) +
#     ggplot2::theme(
#       axis.text.x = ggplot2::element_text(angle = 90, hjust = 1)
#     ) +
#     ggplot2::labs(
#       title = paste0(x$cohort, ' (n = ', x$n_samples, ')'),
#       caption = caption,
#       x = "Sample"
#     )
#
#   # ggplot2::guides(
#   #   fill = ggplot2::guide_legend(
#   #     nrow = ifelse(x$n_catalogue + x$n_denovo > 8, 2, 1))
#   #   )
#
#  if ("groups" %in% names(x$fit))
#     plt = plt + ggplot2::facet_grid(~groups)
#
#
#   return(plt)
# }
#
#
#
#
#
# plot_exposure = function(x,sample_name = T,levels= NULL, flip_coord = F){
#
#   b = x$fit$exposure
#
#   if(is.null(cls)){ cls = ggsci::pal_simpsons()(ncol(b))
#   names(cls) = colnames(b)
#   }
#
#   if(is.null(levels)){ levels =   colnames(b) }
#
#   p = ggplot(data = b %>% as.data.frame() %>% mutate(sample = rownames(b)) %>%
#                reshape2::melt() %>% dplyr::rename(Signature = variable),
#              aes(x = sample, y  = value,
#                  fill = factor(Signature,levels = levels))) +
#     geom_bar(stat = "identity")  + ggplot2::scale_fill_manual(values = cls) + labs(title = "Expsosure", x = "") +
#     theme(axis.text.x = element_text(angle = 90)) +
#     guides(fill=guide_legend(title="Signatures"))
#
#   if (!sample_name) {
#     p =  p +  theme(axis.ticks.x = element_blank(), axis.text.x = element_blank()) + labs(x = "")
#
#   }
#
#   if(flip_coord){
#
#     p =  p + coord_flip()
#   }
#
#   p
# }
#




