#' Function to visualize the estimated exposures
#'
#' @param x basilica object.
#' @param types List of variant types to visualize.
#' @param samples List of samples to visualize.
#' @param clusters List of clusters to visualize.
#' @param sample_name Logical. If `TRUE`, sample names will be reported on the x axis.
#' @param sigs_levels Deprecated.
#' @param cls Custom color palette for signatures.
#' @param add_centroid Logical. If `TRUE`, also clustering's centroids will be plotted.
#' @param muts Deprecated.
#' @param sort_by Signature to sort patients' exposures by.
#'
#' @return ggplot object.
#' @export
#'
#' @examples
#' plot_exposures(example_fit)
plot_exposures = function(x, types=get_types(x), samples=get_samples(x),
                          clusters=get_cluster_labels(x),
                          sample_name=F, sigs_levels=NULL, cls=NULL,
                          add_centroid=FALSE, muts=FALSE, sort_by=NULL) {

  exposures = lapply(types, function(t)
    get_exposure(x, types=types, samples=samples,
                 clusters=clusters, add_groups=TRUE)[[t]] %>%
      dplyr::mutate(type=t)) %>%
    do.call(rbind, .)

  if (is.null(cls)) cls = gen_palette(x, types=types)
  if (is.null(sigs_levels)) sigs_levels = sort(get_signames(x, types=types) %>% unlist(use.names=F))

  p = plot_exposures_aux(exposures=exposures, cls=cls, titlee="Exposure",
                         sigs_level=sigs_levels, sample_name=sample_name,
                         sort_by=sort_by)

  p_centr = plot_centroids(x, types=types, cls=cls, sigs_levels=sigs_levels, sort_by=sort_by)
  if (add_centroid)
    p = patchwork::wrap_plots(p, p_centr, ncol=2, widths=c(9,1), guides="collect")

  return(p)
}

match_type = function(types, sigs) {
  sapply(sigs, function(sid) {
    matches = sapply(types, function(tid) grepl(tid, x=sid))
    which_m = names(matches[matches==TRUE])
    if (length(which_m)==1) return(which_m)
    return(NA)
  }) %>% setNames(NULL)
}


plot_centroids = function(x, types=get_types(x), cls=NULL, sigs_levels=NULL, flip_coord=FALSE, sort_by=NULL) {
  centr = get_centroids(x)
  if (!have_groups(x) || is.null(centr)) return(NULL)
  a_pr = centr %>%
    dplyr::mutate(type=match_type(types, sigs)) %>%
    dplyr::filter(!is.na(type)) %>%
    # dplyr::mutate(clusters=paste0("G",stringr::str_replace_all(clusters,"G",""))) %>%
    # tidyr::separate(col="sigs", into=c("else","sigs"), sep="_") %>% dplyr::mutate("else"=NULL) %>%
    # dplyr::mutate(type=match_type(types, sigs)) %>%
    # dplyr::filter(!is.na(type)) %>%
    dplyr::rename(samples=clusters)

  if (is.null(cls)) cls = gen_palette(x, types=types)

  return(
    plot_exposures_aux(exposures=a_pr, cls=cls, titlee="Centroids",
                       sample_name=TRUE, sample_levels=NULL)
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
