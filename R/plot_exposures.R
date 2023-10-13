plot_exposures = function(x, types=get_types(x), samples=get_samples(x),
                          clusters=get_cluster_labels(x),
                          sample_name=F, sigs_levels=NULL, cls=NULL,
                          add_centroid=FALSE, muts=FALSE, sort_by=NULL) {

  exposures = lapply(types, function(t)
    get_exposure(x, types=types, samples=samples,
                 clusters=clusters, add_groups=TRUE)[[t]] %>%
      dplyr::mutate(type=t)) %>%
    do.call(rbind, .)

  # if (!have_groups(x)) add_centroid = FALSE

  if (is.null(cls)) cls = gen_palette(x, types=types)
  if (is.null(sigs_levels)) sigs_levels = sort(get_signames(x, types=types) %>% unlist(use.names=F))

  # if (add_centroid)
  #   p_centr = plot_centroids(x, cls=cls, sigs_levels=sigs_levels,
  #                            flip_coord=flip_coord, sort_by=sort_by)

  p = plot_exposures_aux(exposures=exposures, cls=cls, titlee="Exposure",
                         sigs_level=sigs_levels, sample_name=sample_name,
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


plot_exposures_aux = function(exposures, cls=NULL, titlee="", sigs_levels=NULL,
                              sample_name=FALSE, sort_by=NULL, sample_levels=NULL) {

  if (!is.null(sort_by))
    sample_levels = exposures %>%
      dplyr::filter(sigs==sort_by) %>%
      dplyr::arrange(desc(value)) %>%
      dplyr::pull(samples) else if (is.null(sample_levels))
        sample_levels = exposures$samples %>% unique()

    p = exposures %>%
      ggplot(aes(x=factor(samples, levels=sample_levels), y=value, fill=sigs)) +
      geom_bar(stat="identity", position="stack") +
      theme_bw() + theme(axis.text.x=element_text(angle=90)) +
      guides(fill=guide_legend(title="Signatures")) + ylab("") + xlab("")

    if (!is.null(cls))
      p = p + ggplot2::scale_fill_manual(values=cls) + ggplot2::scale_color_manual(values=cls)

    if (!sample_name)
      p = p + theme(axis.ticks.x=element_blank(), axis.text.x=element_blank()) + labs(x="")

    if ("clusters" %in% colnames(exposures))
      p = p + facet_grid(type ~ clusters, scales="free_x", space="free_x") else
        p = p + facet_grid(type ~ ., scales="free_x", space="free_x")

    return(p)
}
