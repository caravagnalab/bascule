#' Plot a fit and compare it to ground truth.
#'
#' @param x fit to plot.
#' @param x.true ground truth fit.
#' @param reconstructed logical. If \code{TRUE}, it will return the recontructed counts.
#' @param cls list of colors to use. It overwrites the fit color palette.
#' @param sort_by signature code to sort by in the exposure plot.
#' @param what type of signature. A value among "SBS", "DBS", "CN".
#' @param epsilon logical. If set to \code{TRUE}, the estimated error will be plotted.
#' @param sample_name logical. If set to \code{TRUE} the sample names on the x axis will be reported.
#'
#' @return \code{patchwork} object
#' @export plot_fit

plot_fit = function(x, x.true=NULL, reconstructed=T, cls=NULL,
                    sort_by=NULL, what="SBS", epsilon=FALSE,
                    sample_name=FALSE, sampleIDs=NULL, name1="fit", name2="simul",
                    single_samples=FALSE, add_centroid=F, cutoff=.8) {
  if (is.null(sampleIDs)) sampleIDs = rownames(x$fit$x)
  if (!is.null(x.true)) catalogue = get_signatures(x.true) else catalogue = NULL

  if (!is.null(x.true)) {
    if (length(intersect(get_dn_signames(x), get_signames(x.true)))==0)
      x = convert_sigs_names(x, x.simul=x.true, cutoff=cutoff)

    if (is.null(cls))
      cls = merge_colors_palette(x, x.true)

  }

  # make_plots_compare()

  mm = plot_mutations(x, reconstructed=reconstructed, what=what,
                      epsilon=epsilon, sampleIDs=sampleIDs) + ylab("") + labs(title=paste0("Mutations ", name1))
  alp = plot_exposures(x, sort_by=sort_by, cls=cls, sample_name=sample_name,
                       sampleIDs=sampleIDs, add_centroid=add_centroid) + labs(title=paste0("Exposure ", name1))
  bet = plot_signatures(x, cls=cls, what=what, catalogue=catalogue)

  if (is.null(x.true)) return(patchwork::wrap_plots(bet + (mm/alp), guides="collect"))

  fn_sigs = setdiff(x.true %>% get_signames(),
                    x %>% get_signames()) %>% sort()
  fp_sigs = setdiff(x %>% get_signames(),
                    x.true %>% get_signames()) %>% sort()

  mm.true = plot_mutations(x.true, reconstructed=F, what=what,
                           epsilon=epsilon, sampleIDs=sampleIDs) + ylab("") + labs(title=paste0("Mutations ", name2))
  alpha.true = plot_exposures(x.true, sort_by=sort_by, cls=cls, sample_name=sample_name,
                              sampleIDs=sampleIDs) + labs(title=paste0("Exposure ", name2))

  subtitle_str = ""
  if (length(fn_sigs) > 0) subtitle_str = paste0("Missing signatures: ",
                                                 paste(fn_sigs, collapse=","),
                                                 ". ")
  if (length(fp_sigs) > 0) subtitle_str = paste0(subtitle_str, "Added signatures: ",
                                                 paste(fp_sigs, collapse=","),
                                                 ".")

  design = "AAABBB
            AAACCC
            AAADDD
            ###EEE"

  return(patchwork::wrap_plots(bet, mm, mm.true, alp, alpha.true, design=design, guides="collect") &
           patchwork::plot_annotation(title=paste0(name1, " vs ", name2), subtitle=subtitle_str))
}


get_false_negative_signatures = function(x, x.true)
  return(setdiff(x.true %>% get_signatures() %>% rownames(),
                 x %>% get_signatures() %>% rownames()) %>% sort())


get_false_positive_signatures = function(x, x.true)
  return(setdiff(x %>% get_signatures() %>% rownames(),
                 x.true %>% get_signatures() %>% rownames()) %>% sort())


get_samples_high_expos = function(x, x.true, min_expos, sigs=NULL) {
  fn_sigs = get_false_negative_signatures(x, x.true)
  fp_sigs = get_false_positive_signatures(x, x.true)

  samples_fp = x$fit$exposure[x$fit$exposure[, fp_sigs] > min_expos, ] %>% rownames()
  samples_fn = x.true$fit$exposure[x.true$fit$exposure[, fn_sigs] > min_expos, ] %>% rownames()
  samples_sigs = x.true$fit$exposure[x.true$fit$exposure[, sigs] > min_expos, ] %>% rownames()

  return(list("false_negatives"=samples_fn,
              "false_positives"=samples_fp,
              "input_sigs"=samples_sigs))
}



plot_scores = function(x, which="") {
  scores = get_K_scores(x)

  if (!"groups" %in% colnames(scores)) scores = scores %>% dplyr::mutate(groups=1)

  scores %>%
    dplyr::mutate(K=stringr::str_replace_all(K, "k_denovo", "") %>% as.integer(),
                  seed=stringr::str_replace_all(seed, "seed_", "") %>% as.integer()) %>%
    dplyr::mutate(score=ifelse(score_id=="bic", log(score), score)) %>%

    ggplot() +
    geom_point(aes(x=K, y=score, color=factor(seed), shape=factor(groups))) +
    geom_line(aes(x=K, y=score, color=factor(seed), linetype=factor(groups))) +
    facet_wrap(~score_id, scales="free_y") +
    theme_bw()

}


plot_gradient_norms = function(x) {
  x$fit$gradient_norms %>%
    as.data.frame() %>%
    tibble::rownames_to_column(var="step") %>%
    reshape2::melt(id="step", variable.name="parameter", value.name="value") %>%
    ggplot() +
    geom_line(aes(x=as.integer(step), y=value)) +
    facet_wrap(~parameter, scales="free_y") +
    theme_bw()
}


plot_posterior_probs = function(x) {
  pheatmap::pheatmap(x$fit$post_probs, cluster_rows=T, cluster_cols=F)
}






