#' Plot posterior probabilities
#'
#' Plot the posterior probabilities over the assignments.
#'
#' @param x bascule object
#' @param samples List of samples to visualise.
#'
#' @return ggplot object
#' @export plot_posterior_probs
plot_posterior_probs = function(x, samples=get_samples(x)) {
  if (!have_groups(x)) return(NULL)

  params = get_pyro_stat(x, what="clustering", statname="params")[[1]]$infered_params

  params$post_probs[samples,] %>%
    dplyr::rename_all(.funs=function(x) paste0("G", x)) %>%
    tibble::rownames_to_column(var="samples") %>%
    reshape2::melt(id="samples", variable.name="clusters") %>%
    tibble::as_tibble() %>%

    ggplot() +
    geom_tile(aes(x=clusters, y=samples, fill=value)) +
    scale_fill_distiller(direction=1) +
    xlab("") + ylab("Samples") +
    theme_bw() + theme(panel.grid=element_blank())
}


#' Plot parameters gradients norms
#'
#' Visualise the norms of the parameters gradients computed during the inference.
#'
#' @param x bascule object.
#' @param types List of types to visualise gradient norms for.
#'
#' @return patchwork object
#' @export plot_gradient_norms
plot_gradient_norms = function(x, types=get_types(x)) {
  norms = get_gradient_norms(x, types=types)

  norms_nmf = norms %>%
    dplyr::filter(type != "Clustering") %>%
    ggplot() +
    geom_line(aes(x=as.integer(step), y=value)) +
    facet_wrap(type~parameter, scales="free_y", nrow=length(types)) +
    theme_bw() + labs(title="Gradient norms NMF") +
    xlab("Iteration") + ylab("Norm") +
    scale_y_continuous(labels=function(x) scales::scientific(x))

  if (!have_groups(x)) return(norms_nmf)

  norms_clst = norms %>%
    dplyr::filter(type == "Clustering") %>%
    ggplot() +
    geom_line(aes(x=as.integer(step), y=value)) +
    facet_wrap(~parameter, scales="free_y", nrow=length(types)) +
    theme_bw() + labs(title="Gradient norms clustering") +
    xlab("Iteration") + ylab("Norm") +
    scale_y_continuous(labels=function(x) scales::scientific(x))

  patchwork::wrap_plots(norms_nmf / norms_clst)
}



#' Plot model selection scores
#'
#' @description
#' Function to plot the model selection scores for NMF on each variant type and clustering.
#' Reported scores are BIC and negative log-likelihood.
#'
#' @param x bascule object.
#' @param types List of variant types to visualize.
#' @param remove_outliers Logical. If `TRUE`, outliers in each score will be removed.
#'
#' @return ggplot2 object.
#' @export plot_scores
plot_scores = function(x, types=get_types(x), remove_outliers=FALSE) {
  best_scores = get_best_scores(x, types=types)

  scores = get_scores(x) %>%
    dplyr::select(-value) %>% unique() %>%

    dplyr::right_join(get_best_scores(x, types=types) %>% dplyr::select(-score)) %>%

    dplyr::group_by(parname, type, score_id) %>%
    dplyr::mutate(label=replace(NA, score==min(score) & score_id=="bic", "Best fit"))

  if (remove_outliers)
    scores = scores %>%
      dplyr::group_by(score_id, parname, type) %>%
      dplyr::mutate(is.outlier=score %in% boxplot.stats(score)$out) %>%
      dplyr::filter(!is.outlier)

  scores_nmf = scores %>%
    dplyr::filter(parname == "K") %>%
    dplyr::filter(type %in% types) %>%
    ggplot(aes(x=as.integer(value_fit), y=score)) +
    # geom_point(aes(color=factor(seed)), shape=16, size=3) +
    # geom_line(aes(color=factor(seed)), linetype="solid") +

    geom_point(shape=16, size=3) +
    geom_line() +
    ggrepel::geom_label_repel(aes(label=label), box.padding=0.05, size=3,
                             na.rm=T, show.legend=F) +
    ggh4x::facet_nested_wrap(type ~ score_id, scales="free", nrow=length(types)) +
    theme_bw() + labs(title="Scores") + xlab("K") + ylab("Scores NMF") +
    # guides(color=guide_legend(title="Seed")) +
    scale_x_continuous(breaks=scores %>% dplyr::filter(parname == "K") %>% dplyr::pull(value_fit) %>% unique())

  if (!have_groups(x)) return(scores_nmf)

  scores_clst = scores %>%
    dplyr::filter(parname == "G") %>%
    ggplot(aes(x=as.integer(value_fit), y=score)) +
    # geom_point(aes(color=factor(seed)), size=3) +
    # geom_line(aes(color=factor(seed))) +
    geom_point(shape=16, size=3) +
    geom_line() +
    ggrepel::geom_label_repel(aes(label=label), box.padding=0.05, size=3, na.rm=T) +
    ggh4x::facet_nested_wrap(~ score_id, scales="free", nrow=1) +
    theme_bw() + labs(title="Scores") + xlab("G") + ylab("Scores clustering") +
    guides(color=guide_legend(title="Seed")) +
    # scale_y_continuous(labels=function(x) scales::scientific(x, digits=1)) +
    scale_x_continuous(breaks=scores %>% dplyr::filter(parname == "G") %>% dplyr::pull(value_fit) %>% unique())

  design = paste0(rep("A",length(types)), collapse="\n") %>% paste0("\nB")
  patchwork::wrap_plots(scores_nmf, scores_clst, design=design)

}


plot_losses = function(x) {
  losses = get_losses(x)

  losses %>%
    ggplot() +
    geom_line(aes(x=iteration, y=value)) +
    ggh4x::facet_nested(what + type ~ ., scales="free") +
    theme_bw() + xlab("Iterations") + ylab("Loss") +
    scale_y_continuous(labels=function(x) scales::scientific(x)) +
    labs(title="Loss")
}


plot_likelihoods = function(x) {
  likelihoods = get_likelihoods(x)

  likelihoods %>%
    ggplot() +
    geom_line(aes(x=iteration, y=value)) +
    facet_grid(what ~ type, scales="free") +
    theme_bw() + xlab("Iterations") + ylab("Log-likelihood") +
    labs(title="Log-likelihood")
}


plot_QC = function(x) {
  loss = plot_losses(x)
  lik = plot_likelihoods(x)
  # penalty = plot_penalty(x)
  gradient_norms = plot_gradient_norms(x)
  scores = plot_scores(x)
  patchwork::wrap_plots(scores / gradient_norms / (loss + lik))
}



