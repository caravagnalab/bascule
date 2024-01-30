plot_posterior_probs = function(x, samples=get_samples(x)) {
  if (!have_groups(x)) return(NULL)

  params = get_pyro_stat(x,what="clustering",statname="params")[[1]]$infered_params
  return(pheatmap::pheatmap(params$post_probs[samples,], cluster_rows=T, cluster_cols=F))
}



plot_gradient_norms = function(x, types=get_types(x)) {
  get_gradient_norms(x, types=types) %>%
    ggplot() +
    geom_line(aes(x=as.integer(step), y=value)) +
    facet_wrap(type~parameter, scales="free_y") +
    theme_bw() + labs(title="Gradient norms") + xlab("Iteration") + ylab("Norm")
}



plot_scores = function(x, types=get_types(x)) {
  scores = get_scores(x, types=types)

  scores %>%
    ggplot() +
    geom_point(aes(x=as.integer(value), y=score, color=factor(seed))) +
    geom_line(aes(x=as.integer(value), y=score, color=factor(seed))) +
    ggh4x::facet_nested_wrap(type + parname ~score_id, scales="free",
                             nrow=length(unique(scores$type))) +
    theme_bw() + labs(title="Scores") + xlab("") + ylab("Score")

}


plot_losses = function(x) {
  losses = get_losses(x)

  losses %>%
    ggplot() +
    geom_line(aes(x=iteration, y=losses)) +
    facet_grid(what ~ type, scales="free") +
    theme_bw() + xlab("Iterations") + ylab("Loss") +
    labs(title="Loss")
}


plot_likelihoods = function(x) {
  likelihoods = get_likelihoods(x, what="nmf")

  likelihoods %>%
    ggplot() +
    geom_line(aes(x=iteration, y=likelihood)) +
    facet_grid(what ~ type, scales="free") +
    theme_bw() + xlab("Iterations") + ylab("Log-likelihood") +
    labs(title="Log-likelihood")
}


plot_penalty = function(x) {
  penalty = get_penalty(x)
  if (is.null(penalty)) return(NULL)

  penalty %>%
    ggplot() +
    geom_line(aes(x=iteration, y=penalty)) +
    facet_grid(what ~ type, scales="free") +
    theme_bw() + xlab("Iterations") + ylab("Penalty") +
    labs(title="Penalty")
}


plot_QC = function(x) {
  loss = plot_losses(x)
  lik = plot_likelihoods(x)
  penalty = plot_penalty(x)
  gradient_norms = plot_gradient_norms(x)
  scores = plot_scores(x)
  return(patchwork::wrap_plots(scores / gradient_norms / (loss + lik + penalty)))
}



