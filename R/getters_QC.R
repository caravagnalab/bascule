get_gradient_norms = function(x, types=get_types(x)) {
  vname = "gradient_norms"
  qcs_nmf = get_QC(x, what="nmf", types=types)
  qcs_clustering = get_QC(x, what="clustering")
  norms_nmf = lapply(types, function(tid) qcs_nmf[[tid]][[vname]] %>% as.data.frame() %>%
           tibble::rownames_to_column(var="step") %>%
           reshape2::melt(id="step", variable.name="parameter", value.name="value") %>%
           dplyr::mutate(type=tid)) %>% do.call(rbind, .)
  norms_clustering = qcs_clustering[[vname]] %>% as.data.frame() %>%
    tibble::rownames_to_column(var="step") %>%
    reshape2::melt(id="step", variable.name="parameter", value.name="value") %>%
    dplyr::mutate(type="Clustering")

  return(rbind(norms_nmf, norms_clustering))
}


get_scores = function(x, types=get_types(x)) {
  vname = "scores"
  qcs_nmf = get_QC(x, what="nmf", types=types)

  scores_nmf = lapply(types, function(tid) qcs_nmf[[tid]][[vname]] %>% as.data.frame() %>%
                        dplyr::mutate(type=tid)) %>%
    do.call(rbind, .)

  if (have_groups(x)) {
    qcs_clustering = get_QC(x, what="clustering")
    scores_clustering = qcs_clustering[[vname]] %>% dplyr::mutate(type="Clustering")
  } else scores_clustering = NULL

  return(rbind(scores_nmf, scores_clustering) %>%
           dplyr::mutate(dplyr::across(.cols=c("seed","value"), as.integer)) %>%
           tibble::as_tibble())
}


# params = list("K"=NA,"G"=NA,"seed"=NA)
get_alternative_run = function(x, K=NA, G=NA, seed=NA,
                               types=get_types(x)) {

  if (all(is.na(c(K,G,seed)))) return(cli::cli_alert_warning("No input were provided."))

  if (is.na(seed)) seed = list("nmf"=get_seed(x, what="nmf", types=types),
                               "clustering"=get_seed(x, what="clustering"))
  if (is.na(G)) G = get_G(x, input=TRUE)
  if (is.na(K)) K = get_K(x)

  k_vname = paste0("k_denovo:", K)
  g_vname = paste0("cluster:", G)

  x$clustering = get_alternatives(x, what="clustering")$runs_seed[[paste0("seed:", seed$clustering)]]
  alter_nmf = get_alternatives(x, what="nmf", types=types)
  x$nmf = lapply(types, function(tid) alter_nmf[[tid]]$runs_seed[[paste0("seed:", seed$nmf[[tid]])]]) %>%
    setNames(types)

}


## what %in% c("nmf", "clustering")
get_alternatives = function(x, what, types=get_types(x)) {
  return(get_pyro_stat(x=x, what=what, statname="alternatives", types=types))
}


## what %in% c("nmf", "clustering")
get_QC = function(x, what, types=get_types(x)) {
  return(get_pyro_stat(x=x, what=what, statname="QC", types=types))
}


get_pyro_stat = function(x, what, statname, types=get_types(x)) {
  what = tolower(what)
  if (what=="nmf") return(
    lapply(types, function(tid) x[[what]][[tid]][["pyro"]][[statname]]) %>% setNames(types)
  )
  if (what=="clustering") return(x[[what]][["pyro"]][[statname]])
  cli::cli_alert_warning("`what` must be either 'nmf' or 'clustering'")
}


get_params = function(x, what, types=get_types(x)) {
  params = get_pyro_stat(x, what=what, types=types, statname="params")
  if (what == "nmf")
    return(
      lapply(types, function(tid) params[[tid]]$infered_params)
    )
  return(params$infered_params)
}


get_initial_object = function(x, what="clustering") {
  stopifnot(what=="clustering")

  init_params = get_pyro_stat(x, what=what, statname="params")$init_params
  x[[what]]$params = init_params
  x[[what]]$clusters = tibble::tibble(samples=get_samples(x),
                                      clusters=paste0("G",init_params$init_clusters))
  x[[what]]$centroids = init_params$alpha_prior %>%
    wide_to_long(what="exposures") %>%
    dplyr::rename(clusters=samples)
  return(x)
}


get_losses = function(x, what=get_fittypes(x), types=get_types(x)) {
  lapply(what, function(whatid) {
    lapply(types, function(tid) {
      data.frame(loss=get_QC(x, what=whatid, types=tid)[[1]]$losses,
                 type=tid, what=whatid)
    }) %>% do.call(rbind, .)
  }) %>% do.call(rbind, .)

}

