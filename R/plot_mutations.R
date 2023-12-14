plot_data = function(x, types=get_types(x), samples=get_samples(x),
                     clusters=get_cluster_labels(x), reconstructed=TRUE) {

  counts_list = get_input(x, types=types, samples=samples, clusters=clusters,
                          reconstructed=reconstructed)

  patch = lapply(types, function(tid) {
    counts_list[[tid]] %>% reformat_contexts(what=tid) %>%
      dplyr::mutate(type=!!tid) %>%
      dplyr::filter(samples %in% !!samples) %>%
      ggplot() +
      geom_histogram(aes(x=context, y=value), stat="identity") +
      facet_grid(type ~ variant, scales="free", space="free_x") +
      theme_bw()
    }) %>%
    patchwork::wrap_plots(nrow=length(types)) & ylab(NULL) & xlab(NULL) & theme(plot.margin = margin(5.5, 5.5, 5.5, 0))

  patchwork::wrap_elements(patch) +
    labs(tag="Context") +
    theme(plot.tag=element_text(size=rel(1), angle=0),
          plot.tag.position="bottom")

}
