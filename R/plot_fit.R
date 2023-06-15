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
                    sample_name=FALSE, sampleIDs=NULL,
                    single_samples=FALSE) {
  if (is.null(sampleIDs)) sampleIDs = rownames(x$fit$x)

  if (!is.null(x.true)) catalogue = get_signatures(x.true) else catalogue = NULL

  mm = plot_mutations(x, reconstructed=reconstructed, what=what,
                      epsilon=epsilon, sampleIDs=sampleIDs)
  alp = plot_exposures(x, sort_by=sort_by, cls=cls,
                       sample_name=sample_name, sampleIDs=sampleIDs)
  bet = plot_signatures(x, cls=cls, what=what, catalogue=catalogue)

  if (is.null(x.true)) return(patchwork::wrap_plots(bet + (mm/alp), guides="collect"))

  fn_sigs = setdiff(x.true %>% get_signatures() %>% rownames(),
                    x %>% get_signatures() %>% rownames()) %>% sort()
  fp_sigs = setdiff(x %>% get_signatures() %>% rownames(),
                    x.true %>% get_signatures() %>% rownames()) %>% sort()

  mm.true = plot_mutations(x.true, reconstructed=F, what=what,
                           epsilon=epsilon, sampleIDs=sampleIDs)
  alpha.true = plot_exposures(x.true, sort_by=sort_by, cls=cls,
                              sample_name=sample_name, sampleIDs=sampleIDs)

  subtitle_str = ""
  if (length(fn_sigs) > 0) subtitle_str = paste0("Missing signatures: ", fn_sigs, ". ")
  if (length(fp_sigs) > 0) subtitle_str = paste0(subtitle_str, "Added signatures: ", fp_sigs, ".")

  return(patchwork::wrap_plots(bet + (mm/mm.true/alp/alpha.true), guides="collect") &
           patchwork::plot_annotation(subtitle=subtitle_str))
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

  scores %>%
    dplyr::mutate(K=stringr::str_replace_all(K, "K_", "") %>% as.integer(),
                  seed=stringr::str_replace_all(seed, "seed_", "") %>% as.integer()) %>%
    dplyr::mutate(score=ifelse(score_id=="bic", log(score), score)) %>%

    ggplot() +
    geom_point(aes(x=K, y=score, color=factor(seed))) +
    geom_line(aes(x=K, y=score, color=factor(seed))) +
    facet_wrap(~score_id, scales="free") +
    theme_bw()

}
