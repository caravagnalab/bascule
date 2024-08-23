#' Function to visualize the input and reconstructed mutation counts.
#'
#' @param x bascule object.
#' @param types List of variant types to visualize.
#' @param samples List of samples to visualize.
#' @param clusters List of clusters to visualize.
#' @param reconstructed Logical. If `TRUE`, the reconstructed counts, i.e., computed as alpha x beta, will be visualized.
#' @param color_by_sigs Logical. If `TRUE`, each context will report the number of mutations split by signature.
#'
#' @return ggplot object.
#' @export
plot_data = function(x, types=get_types(x), samples=get_samples(x),
                     clusters=get_cluster_labels(x), reconstructed=TRUE,
                     color_by_sigs=FALSE) {

  counts_list = get_input(x, types=types, samples=samples, clusters=clusters,
                          reconstructed=reconstructed, by_sigs=color_by_sigs)
  counts_list = lapply(types, function(tid)
    counts_list[[tid]] %>% dplyr::bind_rows()) %>%
    setNames(types)

  facet_vars = "type ~ variant"
  mapping = aes(x=context, y=value)
  if (have_groups(x)) facet_vars = "type ~ clusters + variant"
  if (color_by_sigs) mapping = aes(x=context, y=value, fill=sigs)

  patch = lapply(types, function(tid) {
    counts_list[[tid]] %>% reformat_contexts(what=tid) %>%
      dplyr::mutate(type=!!tid) %>%
      dplyr::filter(samples %in% !!samples) %>%
      ggplot() +
      geom_histogram(mapping=mapping, stat="identity", position="stack") +
      ggh4x::facet_nested(as.formula(facet_vars), scales="free", space="free_x") +
      scale_fill_manual(values=gen_palette(x)) +
      theme_bw()
    }) %>%
    patchwork::wrap_plots(nrow=length(types)) & ylab(NULL) & xlab(NULL) & theme(plot.margin=margin(5.5, 5.5, 5.5, 0))

  patchwork::wrap_elements(patch) +
    labs(tag="Context") +
    theme(plot.tag=element_text(size=rel(1), angle=0),
          plot.tag.position="bottom")

}


plot_data_differences = function(x, types=get_types(x), samples=get_samples(x),
                                 clusters=get_cluster_labels(x)) {

  counts_reconstructed = get_input(x, types=types, samples=samples, clusters=clusters, reconstructed=TRUE)
  counts_true = get_input(x, types=types, samples=samples, clusters=clusters, reconstructed=FALSE)

  counts_diff = lapply(types, function(tid) {
    counts_reconstructed[[tid]] %>% dplyr::inner_join(counts_true[[tid]] %>% dplyr::rename(value_true=value)) %>%
      dplyr::mutate(difference=abs(value - value_true))
  }) %>% setNames(types)

  facet_vars = "type ~ variant"
  mapping = aes(x=context, y=difference)
  if (have_groups(x)) facet_vars = "type ~ clusters + variant"

  patch = lapply(types, function(tid) {
    counts_diff[[tid]] %>% reformat_contexts(what=tid) %>%
      dplyr::mutate(type=!!tid) %>%
      dplyr::filter(samples %in% !!samples) %>%
      ggplot() +
      geom_histogram(mapping=mapping, stat="identity", position="stack") +
      ggh4x::facet_nested(as.formula(facet_vars), scales="free", space="free_x") +
      scale_fill_manual(values=gen_palette(x)) +
      theme_bw()
  }) %>%
    patchwork::wrap_plots(nrow=length(types)) & ylab(NULL) & xlab(NULL) & theme(plot.margin=margin(5.5, 5.5, 5.5, 0))

  patchwork::wrap_elements(patch) +
    labs(tag="Context") +
    theme(plot.tag=element_text(size=rel(1), angle=0),
          plot.tag.position="bottom")

}
