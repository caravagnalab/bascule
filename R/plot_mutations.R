plot_data = function(x, types=get_types(x), samples=get_samples(x),
                     clusters=get_cluster_labels(x)) {

  counts_list = get_input(x, types=types, samples=samples, clusters=clusters)

  counts = lapply(types, function(tid) {
    counts_list[[tid]] %>% reformat_contexts(what=tid) %>%
      dplyr::mutate(type=!!tid)
    }) %>% do.call(rbind, .) %>%
    dplyr::filter(samples %in% !!samples)

  counts %>%
    ggplot() +
    geom_histogram(aes(x=context, y=value), stat="identity") +
    ggh4x::facet_nested(.~type+variant, scales="free") +
    theme_bw()

}
