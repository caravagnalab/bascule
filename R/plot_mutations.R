#' Function to visualize the input and reconstructed mutation counts.
#'
#' @param x bascule object.
#' @param types List of variant types to visualize.
#' @param samples List of samples to visualize.
#' @param clusters List of clusters to visualize.
#' @param reconstructed Logical. If `TRUE`, the reconstructed counts, i.e., computed as alpha x beta, will be visualized.
#' @param color_by_sigs Logical. If `TRUE`, each context will report the number of mutations split by signature.
#' @param plot_by_sample Logical. If `TRUE`, the function will return a list of ggplot objects with the counts for each sample.
#'
#' @return one or a list of ggplot objects.
#' @export plot_data
plot_data = function(x, types=get_types(x), samples=get_samples(x),
                     clusters=get_cluster_labels(x), reconstructed=TRUE,
                     color_by_sigs=FALSE, plot_by_sample=FALSE) {

  counts_list = get_input(x, types=types, samples=samples, clusters=clusters,
                          reconstructed=reconstructed, by_sigs=color_by_sigs)
  counts_list = lapply(types, function(tid)
    counts_list[[tid]] %>% dplyr::mutate(type=tid)) %>%
    dplyr::bind_rows()

  facet_vars = "type ~ variant"
  mapping = aes(x=context, y=value)
  if (have_groups(x)) facet_vars = "type ~ clusters + variant"
  if (color_by_sigs) mapping = aes(x=context, y=value, fill=sigs)

  main_plot = function(counts_list, types, samples) {
    lapply(types, function(tid) {
      counts_list %>%
        dplyr::filter(type==tid) %>%
        # dplyr::select(-type) %>%
        dplyr::filter(samples %in% !!samples) %>%
        reformat_contexts(what=tid) %>%
        ggplot() +
        geom_histogram(mapping=mapping, stat="identity", position="stack") +
        ggh4x::facet_nested(as.formula(facet_vars), scales="free", space="free_x") +
        scale_fill_manual(values=gen_palette(x)) +
        theme_bw()
    }) %>%
      patchwork::wrap_plots(nrow=length(types)) & ylab(NULL) & xlab(NULL) & theme(plot.margin=margin(5.5, 5.5, 5.5, 0))
  }

  if (plot_by_sample) {
    plots = lapply(samples, function(sample_id) main_plot(counts_list, types, sample_id) + labs(title=sample_id)) %>% setNames(samples)
  } else {
    plots = main_plot(counts_list, types, samples)
  }

  # patch = main_plot(counts_list, types, samples)

  # patchwork::wrap_elements(patch) +
  #   labs(tag="Context") +
  #   theme(plot.tag=element_text(size=rel(1), angle=0),
  #         plot.tag.position="bottom")

  return(plots)
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
