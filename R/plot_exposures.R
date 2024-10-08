#' Function to visualize the estimated exposures
#'
#' @param x bascule object.
#' @param types List of variant types to visualize.
#' @param samples List of samples to visualize.
#' @param clusters List of clusters to visualize.
#' @param sample_name Logical. If `TRUE`, sample names will be reported on the x axis.
#' @param color_palette Custom color palette for signatures.
#' @param add_centroid Logical. If `TRUE`, also clustering's centroids will be plotted.
#' @param sort_by Signature to sort patients' exposures by.
#' @param exposure_thr Only signatures with exposures greater than `exposure_thr` in all samples will be highlighted.
#' @param quantile_thr add
#' @param signatures_list add
#'
#' @return ggplot object.
#' @export plot_exposures
plot_exposures = function(x,
                          types=get_types(x),
                          samples=get_samples(x),
                          clusters=get_cluster_labels(x),
                          sample_name=FALSE,
                          color_palette=NULL,
                          add_centroid=FALSE,
                          sort_by=NULL,
                          exposure_thr=0,
                          quantile_thr=0,
                          signatures_list=get_signames(x)) {

  exposures = lapply(types, function(t)
    get_exposure(x, types=types, samples=samples,
                 clusters=clusters, add_groups=TRUE)[[t]] %>%
      dplyr::mutate(type=t)) %>%
    do.call(rbind, .)

  if (is.null(color_palette)) {
    cls = gen_palette(x, types=sort(types))
  } else {
    cls = color_palette
  }

  # merging signatures where their exposure in all the samples are below the threshold
  to_keep = exposures$sigs %>% unique(); sigs_levels = NULL
  if (exposure_thr > 0 | quantile_thr > 0) {
    scores = get_clusters_score(x, types=types, exposure_thr=exposure_thr,
                                quantile_thr=quantile_thr) %>%
      dplyr::rename("sigs"="signature", "clusters"="cluster") %>% tibble::as_tibble()

    exposures = exposures %>% dplyr::left_join(scores) %>%
      dplyr::mutate(sigs=ifelse(significance, sigs, "Other")) %>%
      dplyr::select(samples, clusters, sigs, value, type)
  }

  to_keep = intersect(unlist(signatures_list), exposures$sigs) %>% as.character()
  sigs_levels = c(sort(to_keep) %>% purrr::discard(.p=function(x) x=="Other"), "Other")
  cls["Other"] = "gainsboro"

  # keep only signatures passed as input
  exposures = exposures %>%
    dplyr::mutate(sigs=ifelse(sigs%in%to_keep, sigs, "Other"))

  p = plot_exposures_aux(exposures=exposures, cls=cls,
                         titlee="Exposure",
                         sample_name=sample_name,
                         sort_by=sort_by,
                         sigs_levels=sigs_levels) +
    scale_fill_manual(values=cls)

  p_centr = plot_centroids(x,
                           types=types,
                           cls=cls,
                           sort_by=sort_by,
                           exposure_thr=exposure_thr,
                           quantile_thr=quantile_thr,
                           sigs_levels=sigs_levels) +
    scale_fill_manual(values=cls)

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


#' Plot clustering centroids
#'
#' @param x bascule object.
#' @param types List of variant types to visualize.
#' @param clusters List of clusters to visualize.
#' @param cls Custom color palette for signatures.
#' @param sort_by Signature to sort patients' exposures by.
#' @param exposure_thr Only signatures with exposures greater than `exposure_thr` in all samples will be highlighted.
#' @param quantile_thr add
#' @param signatures_list add
#' @param ... Additional arguments
#'
#' @import gtools
#'
#' @return ggplot object.
#' @export plot_centroids
plot_centroids = function(x,
                          types = get_types(x),
                          clusters=get_cluster_labels(x),
                          cls = NULL,
                          sort_by = NULL,
                          exposure_thr = 0,
                          quantile_thr = 0,
                          signatures_list = get_signames(x),
                          ...) {

  centr = get_centroids(x)
  if ("sigs_levels" %in% names(list(...))) { sigs_levels = list(...)$sigs_levels } else { sigs_levels = NULL }

  if (!have_groups(x) || is.null(centr)) return(NULL)

  a_pr = centr %>%
    dplyr::mutate(type=match_type(types, sigs)) %>%
    dplyr::filter(!is.na(type)) %>%
    dplyr::filter(clusters %in% !!clusters) %>%
    dplyr::rename(samples=clusters)

  if (is.null(cls)) cls = gen_palette(x, types=sort(types))

  # Just plot significant signatures in each cluster [concise=TRUE]
  if (exposure_thr > 0 | quantile_thr > 0) {
    scores = get_clusters_score(x, types=types, exposure_thr=exposure_thr,
                                quantile_thr=quantile_thr) %>%
      dplyr::rename("sigs"="signature", "samples"="cluster") %>% tibble::as_tibble()

    a_pr = a_pr %>% dplyr::left_join(scores) %>%
      dplyr::mutate(sigs=ifelse(significance, sigs, "Other")) %>%
      dplyr::select(samples, sigs, value, type)
  }

  to_keep = intersect(unlist(signatures_list), a_pr$sigs) %>% as.character()
  sigs_levels = c(gtools::mixedsort(to_keep) %>% purrr::discard(.p=function(x) x=="Other"), "Other")
  cls["Other"] = "gainsboro"

  # keep only signatures passed as input
  a_pr = a_pr %>%
    dplyr::mutate(sigs=ifelse(sigs%in%to_keep, sigs, "Other"))

  sample_levels = gtools::mixedsort(unique(a_pr$samples))

  return(
    plot_exposures_aux(exposures=a_pr,
                       cls=cls,
                       titlee="Centroids",
                       sample_name=TRUE,
                       sample_levels=sample_levels,
                       sigs_levels=sigs_levels) +
      scale_fill_manual(values=cls) + theme(axis.text.x=element_text(angle=0))
  )
}


plot_exposures_aux = function(exposures,
                              cls=NULL,
                              titlee="",
                              sample_name=FALSE,
                              sigs_levels=NULL,
                              sort_by=NULL,
                              sample_levels=NULL) {

  if (!is.null(sigs_levels))
    exposures = exposures %>% dplyr::mutate(sigs=factor(sigs, levels=sigs_levels))

  if (!is.null(sort_by))
    sample_levels = exposures %>%
      dplyr::filter(sigs==sort_by) %>%
      dplyr::group_by(samples) %>%
      dplyr::summarise(value=sum(value)) %>%
      dplyr::arrange(desc(value)) %>%
      dplyr::pull(samples) else if (is.null(sample_levels))
        sample_levels = exposures$samples %>% unique()

  p = exposures %>%
    ggplot(aes(x=factor(samples, levels=sample_levels), y=value, fill=sigs)) +
    geom_bar(stat="identity", position="stack", show.legend=T) +
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
