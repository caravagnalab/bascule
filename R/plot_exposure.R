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
                          add_centroid=FALSE, flip_coord=FALSE, facet=TRUE,
                          muts=FALSE, sampleIDs=NULL, sort_by=NULL) {

  if (!have_groups(x)) add_centroid = FALSE
  if (is.null(sampleIDs)) sampleIDs = rownames(get_exposure(x))

  b = get_exposure(x, add_groups=TRUE)

  if (muts) b = b * rowSums(x$fit$x)

  b = b[sampleIDs,]

  if (is.null(cls) && !have_color_palette(x))
    cls = gen_palette(ncol(b)) %>% setNames(colnames(b))
  else if (is.null(cls) && have_color_palette(x))
    cls = get_color_palette(x)

  if (is.null(sigs_levels))
    sigs_levels = sort(colnames(b %>% dplyr::select(-dplyr::contains("groups"))))

  idcols = c("sample")
  if (have_groups(x)) idcols = c("sample","groups")

  b = b %>% as.data.frame() %>%
    tibble::rownames_to_column(var="sample")

  if (add_centroid)
    p_centr = plot_centroids(x, cls=cls, sigs_levels=sigs_levels,
                             flip_coord=flip_coord, sort_by=sort_by)

  b = b %>%
    reshape2::melt(id=idcols, variable.name="Signature", value.name="alpha")

  p = plot_exposures_aux(expos=b, cls=cls, titlee="Exposure", sigs_level=sigs_levels,
                         sample_name=sample_name, flip_coor=flip_coord, facet=facet,
                         sort_by=sort_by)

  if (add_centroid)
    p = patchwork::wrap_plots(p, p_centr, ncol=2, widths=c(9,1), guides="collect")

  return(p)
}


plot_centroids = function(x, cls=NULL, sigs_levels=NULL, flip_coord=FALSE, sort_by=NULL) {
  a_pr = get_centroids(x, normalize=T)
  grps = rownames(a_pr)
  a_pr = a_pr[as.character(grps),]
  rownames(a_pr) = paste0("G", grps %>% stringr::str_replace_all("G",""))

  if (is.null(cls) && !have_color_palette(x))
    cls = gen_palette(ncol(a_pr)) %>% setNames(colnames(a_pr))
  else if (is.null(cls) && have_color_palette(x))
    cls = get_color_palette(x)

  if (is.null(sigs_levels))
    sigs_levels = sort(colnames(a_pr %>% dplyr::select(-dplyr::contains("groups"))))

  a_pr = a_pr %>% as.data.frame() %>% tibble::rownames_to_column(var="sample") %>%
    reshape2::melt(id="sample", variable.name="Signature", value.name="alpha")

  return(
    plot_exposures_aux(expos=a_pr, cls=cls, titlee="Centroids", sigs_level=sigs_levels,
                       sample_name=TRUE, flip_coor=flip_coord, facet=FALSE,
                       sample_levels=paste0("G", sort(grps) %>% stringr::str_replace_all("G","")))
  )
}


plot_exposures_aux = function(expos, cls=NULL, titlee="", sigs_levels=NULL,
                              sample_name=FALSE, flip_coord=FALSE,
                              facet=TRUE, sort_by=NULL, sample_levels=NULL) {

  if (!is.null(sort_by))
    sample_levels = expos %>%
      dplyr::filter(Signature==sort_by) %>%
      dplyr::arrange(desc(alpha)) %>%
      dplyr::pull(sample) else if (is.null(sample_levels))
        sample_levels = expos$sample %>% unique()

  p = expos %>%
    ggplot(aes(x=factor(sample, levels=sample_levels), y=alpha, fill=Signature)) +
    geom_bar(stat="identity", position="stack") +
    labs(title=titlee) +
    theme_bw() + theme(axis.text.x=element_text(angle=90)) +
    guides(fill=guide_legend(title="Signatures")) + ylab("") + xlab("")

  if (!is.null(cls))
    p = p + ggplot2::scale_fill_manual(values=cls) + ggplot2::scale_color_manual(values=cls)

  if ("alpha_centr" %in% colnames(expos))
    p = p + geom_hline(aes(yintercept=alpha_centr, color=Signature),
                       linetype="dashed") + guides(color="none")

  if (!sample_name)
    p = p + theme(axis.ticks.x=element_blank(), axis.text.x=element_blank()) + labs(x="")

  if (flip_coord)
    p = p + coord_flip()

  if ("groups" %in% colnames(expos) && facet)
    p = p + facet_grid(~ groups, scales="free_x", space="free_x")

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



plot_exposures_real = function(x, groups_true, cls=NULL, sort_by=NULL) {
  if (is.null(cls)) cls = get_color_palette(x)

  p = get_exposure(x, long=F, add_groups=T) %>%
    dplyr::mutate(groups_true=groups_true) %>%
    tibble::rownames_to_column(var="sample") %>%
    reshape2::melt(id=c("sample","groups","groups_true"),
                   variable.name="Signature", value.name="alpha") %>%
    plot_exposures_aux(facet=FALSE, cls=cls, sort_by=sort_by)

  p + ggh4x::facet_nested(~groups_true+groups, scale="free_x", space="free_x")
}



