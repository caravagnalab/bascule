# Fit infos #####
get_n_denovo = function(x) {
  lapply(get_types(x), function(tid) {
    get_denovo_signames(x)[[tid]] %>% length()
  }) %>% setNames(get_types(x))
}


get_n_groups = function(x) {
  if (!have_groups(x)) return(1)
  return(get_cluster_labels(x) %>% length())
}


get_seed = function(x) {
  list("clustering"=get_pyro_stat(x, what="clustering", statname="seed")[[1]],
       "nmf"=get_pyro_stat(x, what="nmf", statname="seed"))
}


get_best_seed = function(x, value, type_id, parname) {
  get_scores(x) %>%
    dplyr::filter(type==type_id, score_id=="bic", parname==!!parname, value==!!value) %>%
    dplyr::filter(score==min(score)) %>%
    dplyr::pull(seed)
}


# Scores #####

get_gradient_norms = function(x, types=get_types(x)) {
  vname = "gradient_norms"
  qcs_nmf = get_QC(x, what="nmf", types=types)
  qcs_clustering = get_QC(x, what="clustering")[[1]]
  norms_nmf = lapply(types, function(tid) (qcs_nmf[[tid]] %>%
                       dplyr::filter(stat==!!vname) %>%
                       dplyr::pull(value))[[1]] %>%
                       data.frame() %>%
                       tibble::rownames_to_column(var="step") %>%
                       reshape2::melt(id="step", variable.name="parameter", value.name="value") %>%
                       dplyr::mutate(type=tid)) %>%
    do.call(rbind, .)
  norms_clustering = (qcs_clustering %>%
    dplyr::filter(stat==!!vname) %>%
    dplyr::pull(value))[[1]] %>%
    data.frame() %>%
    tibble::rownames_to_column(var="step") %>%
    reshape2::melt(id="step", variable.name="parameter", value.name="value") %>%
    dplyr::mutate(type="Clustering")

  return(rbind(norms_nmf, norms_clustering))
}


get_scores = function(x, types=get_types(x)) {

  get_alternatives(x) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(bic=ifelse(
      parname == "K",
      pyro_fit$QC %>% dplyr::filter(stat=="bic") %>% dplyr::pull(value) %>% unlist(),
      pyro_fit$pyro$QC %>% dplyr::filter(stat=="bic") %>% dplyr::pull(value) %>% unlist()
      ),

      likelihood=ifelse(
        parname == "K",
        pyro_fit$QC %>% dplyr::filter(stat=="likelihood") %>% dplyr::pull(value) %>% unlist(),
        pyro_fit$pyro$QC %>% dplyr::filter(stat=="likelihood") %>% dplyr::pull(value) %>% unlist()
      )) %>%

    dplyr::ungroup() %>%

    tidyr::pivot_longer(cols=c("bic","likelihood"), values_to="score", names_to="score_id") %>%

    dplyr::select(-pyro_fit)

}


get_losses = function(x, what=get_fittypes(x), types=get_types(x)) {
  losses = get_QC_stat(x, statname="losses")
  lapply(what, function(whatid) {
    losses[[whatid]] %>%
      data.frame() %>%
      reshape2::melt(variable.name="type") %>%
      dplyr::mutate(what=whatid)
  }) %>% do.call(rbind, .) %>%

    dplyr::group_by(type, what) %>%
    dplyr::mutate(iteration=1:dplyr::n()) %>%

    dplyr::ungroup()
}


# get_likelihoods = function(x, types=get_types(x)) {
#   # get_stats(x, what, types, statname="likelihood") %>%
#   get_QC_stat(x, statname=)[["nmf"]] %>%
#
#     dplyr::group_by(type, what) %>%
#     dplyr::mutate(iteration=1:dplyr::n())
# }


get_QC_stat = function(x, statname) {
  values_nmf = lapply(get_types(x), function(tid) {
    (get_QC(x, what="nmf")[[tid]] %>%
       dplyr::filter(stat==statname) %>%
       dplyr::pull(value))[[1]]
  }) %>% setNames(get_types(x))

  values_cls = (get_QC(x, what="clustering")[[1]] %>%
       dplyr::filter(stat==statname) %>%
       dplyr::pull(value))[[1]]

  list("nmf"=values_nmf,
       "clustering"=values_cls)
}


# get_penalty = function(x, what=get_fittypes(x), types=get_types(x)) {
#   penalty = get_stats(x, what, types, statname="penalty")
#   if (nrow(penalty) == 0) return(NULL)
#
#   penalty %>%
#     dplyr::group_by(type, what) %>%
#     dplyr::mutate(iteration=1:dplyr::n())
# }


# Alternative runs #####

# params = list("K"=NA,"G"=NA,"seed"=NA)
get_alternative_run = function(x, K=get_n_denovo(x), G=get_n_groups(x),
                               seed=c()
                               # types=get_types(x)
                               ) {

  alternatives_df = get_alternatives(x)

  if (is.null(alternatives_df)) {
    cli::cli_alert_warning("No alternatives. Returning the original object.")
    return(x)
  }

  best_K = get_n_denovo(x); best_G = get_n_groups(x)#; best_seed = get_seed(x)
  K = c(K, best_K[setdiff(names(best_K), names(K))])
  G = c(G, best_G[setdiff(names(best_G), names(G))])

  grps = G
  sigs = K %>% setNames(names(K))
  if (have_groups(x) & !is.null(G)) {
    if (is.null(seed[["clustering"]]) & have_groups(x)) {
      s_t = get_best_seed(x, value=G, parname="G", type_id="Clustering")
    } else {
      s_t = seed$Clustering
    }
    x$clustering = (alternatives_df %>%
      dplyr::filter(what_fit == "clustering",
                    seed == s_t,
                    value == grps) %>%
      dplyr::pull(basilica_fit))[[1]][[1]]

    # x$clustering = get_alternatives(x, what="clustering")[[1]]$fits[[grps]][[paste0("seed:",seed$clustering)]][[1]]
  }

  # alter_nmf = get_alternatives(x, types=types)
  x$nmf = lapply(types, function(tid) {
    k_t = sigs[[tid]]
    if (is.null(seed[["nmf"]][[tid]])) {
      s_t = get_best_seed(x, value=K[[tid]], parname="K", type_id=tid)
    } else {
      s_t = seed$nmf[[tid]]
    }
    if (is.null(k_t)) k_t = get_n_denovo(x)[[tid]]
    if (is.null(s_t)) s_t = get_seed(x)$nmf[[tid]]

    alter_t = alter_nmf[[tid]]$fits[[k_t]][[s_t]][[1]]

    alter_t = alternatives_df %>%
      dplyr::filter(what_fit == "nmf",
                    type == tid,
                    seed == s_t,
                    value == grps)

    if (nrow(alter_t) == 0) return(x$nmf[[tid]])

    dplyr::pull(alter_t, basilica_fit)[[1]][[1]]

    # list("exposure"=alter_t$exposure,
    #      "beta_fixed"=alter_t$beta_fixed,
    #      "beta_denovo"=alter_t$beta_denovo,
    #      "pyro"=alter_t)
  }) %>% setNames(types)

  return(x)
}


get_initial_object = function(x, what="clustering") {
  if (what!="clustering" || !have_groups(x)) {
    cli::cli_alert_warning("what != 'clustering' or no groups are in the input object, the input object will be returned.")
    return(x)
  }

  init_params = get_pyro_stat(x, what=what, statname="params")[[1]]$init_params
  x[[what]]$pyro$params$infered_params = init_params
  x[[what]]$clusters = tibble::tibble(samples=get_samples(x),
                                      clusters=paste0("G",init_params$init_clusters))
  x[[what]]$centroids = init_params$alpha_prior %>%
    wide_to_long(what="exposures") %>%
    dplyr::rename(clusters=samples)
  return(x)
}


get_nmf_step1 = function(x) {
  for (tid in get_types(x))
    x[["nmf"]][[tid]] = x[["nmf"]][[tid]][["nmf_step1"]]

}



# Parameters #####

get_params = function(x, what, types=get_types(x)) {
  params = get_pyro_stat(x, what=what, types=types, statname="params")
  if (what == "nmf")
    return(
      lapply(types, function(tid) params[[tid]]$infered_params) %>% setNames(types)
    )
  return(params[[1]]$infered_params)
}


get_nmf_initial_parameters = function(x, what, types=get_types(x)) {
  params = get_pyro_stat(x, what=what, types=types, statname="params")
  lapply(types, function(tid) {
    params[[tid]]$init_params
  }) %>% setNames(types)
}


## Fix
# get_train_params = function(x, what, types=get_types(x)) {
#   qc = get_QC(x, what=what, types=types)
#   if (what=="nmf")
#     lapply(types, function(tid) qc[[tid]][["train_params"]]) %>%
#     setNames(types)
#   else
#     qc[[1]][["train_params"]]
# }


# Aux functions #####

## returns a dataframe with the basilica object for each tested configuration
get_alternatives = function(x) {
  alternatives_nmf = get_pyro_stat(x, what="nmf", statname="alternatives") %>%
    purrr::discard(purrr::is_empty)
  alternatives_cls = get_pyro_stat(x, what="clustering", statname="alternatives") %>%
    purrr::discard(purrr::is_empty)

  alt = do.call(rbind, c(alternatives_nmf, alternatives_cls))
  rownames(alt) = NULL
  return(alt)
}


## what %in% c("nmf", "clustering")
get_QC = function(x, what, types=get_types(x)) {
  return(get_pyro_stat(x=x, what=what, statname="QC", types=types))
}


get_pyro_stat = function(x, what, statname, types=get_types(x)) {
  what = tolower(what)
  if (what=="nmf") return(
    lapply(types, function(tid) {
      val = x[[what]][[tid]][["pyro"]][[statname]]
      if (is.null(val)) return(NULL)
      return(val)
    }) %>% setNames(types)
  )
  if (what=="clustering") {
    return(list(x[[what]][["pyro"]][[statname]]))
  }
  cli::cli_alert_warning("`what` must be either 'nmf' or 'clustering'")
}




