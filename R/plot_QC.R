plot_posterior_probs = function(x) {
  if (!have_groups(x)) return(NULL)

  params = get_pyro_stat(x,what="clustering",statname="params")$infered_params
  return(pheatmap::pheatmap(params$post_probs, cluster_rows=T, cluster_cols=F))
}



plot_gradient_norms = function(x, types=get_types(x)) {
  get_gradient_norms(x, types=types) %>%
    ggplot() +
    geom_line(aes(x=as.integer(step), y=value)) +
    facet_wrap(type~parameter, scales="free_y") +
    theme_bw()
}



plot_scores = function(x, types=get_types(x)) {
  scores = get_scores(x, types=types)

  scores %>%
    # dplyr::mutate(score=ifelse(!grepl("llik", score_id), log(score), score)) %>%

    ggplot() +
    geom_point(aes(x=value, y=score, color=factor(seed))) +
    geom_line(aes(x=value, y=score, color=factor(seed))) +
    ggh4x::facet_nested_wrap(type + parname ~score_id, scales="free",
                             nrow=length(unique(scores$type))) +
    theme_bw()

}



