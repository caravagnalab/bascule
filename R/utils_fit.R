pyfit = function(counts,
                 k_list,
                 lr = 0.005,
                 optim_gamma = 0.1,
                 n_steps = 2000,
                 stage = "",
                 py = NULL,
                 clusters = NULL,
                 nonparametric = TRUE,
                 reference_cat = NULL,
                 hyperparameters = NULL,
                 CUDA = FALSE,
                 compile = FALSE,
                 seed_list = c(10),
                 store_fits = TRUE,
                 store_parameters = FALSE,
                 filter_dn = FALSE,
                 type = ""
) {

  TIME = as.POSIXct(Sys.time(), format = "%H:%M:%S")

  if (is.null(py)) py = reticulate::import("pybasilica")

  if (length(k_list) > 1) k_list = reticulate::r_to_py(as.integer(k_list)) else
    k_list = reticulate::r_to_py(list(as.integer(k_list)))

  if (length(seed_list) > 1) seed_list = reticulate::r_to_py(as.integer(seed_list)) else
    seed_list = reticulate::r_to_py(list(as.integer(seed_list)))

  if (!is.null(clusters)) clusters = as.integer(clusters)

  obj = py$fit(x = counts, k_list = k_list, lr = lr, optim_gamma = optim_gamma, n_steps = n_steps,
               cluster = clusters, beta_fixed = reference_cat,
               hyperparameters = hyperparameters, nonparametric = nonparametric,
               store_parameters = store_parameters, stage = stage,
               seed_list = seed_list, compile_model = compile,
               CUDA = CUDA, store_fits = store_fits)

  TIME = difftime(as.POSIXct(Sys.time(), format = "%H:%M:%S"), TIME, units = "mins")

  if (is.list(obj)) {
    bestRun = obj[[1]]
    secondBest = obj[[2]]
  } else {
    bestRun = obj
    secondBest = NULL
  }

  # save python object data in a list
  py_obj = get_list_from_py(bestRun, filter_dn=filter_dn, type=type)
  py_obj$alternatives$secondBest = get_list_from_py(secondBest, filter_dn=filter_dn, type=type)
  py_obj$time = TIME

  return(py_obj)
}


filter_sigs_low_expos = function(x, min_exp=0.15, keep_sigs=NULL) {
  if (is.null(x$exposure)) return(x)

  sbs_keep = x$exposure %>%
    dplyr::mutate(value=ifelse(value < min_exp, 0, value)) %>%
    dplyr::filter(value > 0) %>%
    dplyr::pull(sigs) %>% unique() %>% as.character()

  if (!is.null(keep_sigs)) sbs_keep = c(sbs_keep, keep_sigs) %>% unique()
  if (length(sbs_keep) == 0) return(x)

  if (!is.null(x$beta_fixed)) x$beta_fixed = x$beta_fixed %>% dplyr::filter(sigs %in% sbs_keep)
  if (!is.null(x$beta_denovo)) x$beta_denovo = x$beta_denovo %>% dplyr::filter(sigs %in% sbs_keep)
  x$exposure = x$exposure %>% dplyr::filter(sigs %in% sbs_keep)

  return(x)
}


rename_object = function(x, map_names, types=get_types(x)) {
  ## MISSING CONVERSION OF STORED OBJECTS
  mapp = names(map_names) %>% setNames(map_names)
  for (tid in types) {
    alpha_long = get_exposure(x)[[tid]] %>%
      dplyr::mutate(sigs=mapp[sigs])
    dn_long = get_denovo_signatures(x)[[tid]] %>%
      dplyr::mutate(sigs=mapp[sigs])

    x$nmf[[tid]]$exposure = x$nmf[[tid]]$pyro$exposure = alpha_long
    x$nmf[[tid]]$beta_denovo = x$nmf[[tid]]$pyro$beta_denovo = dn_long

    x$nmf[[tid]]$pyro$params$infered_params$alpha = alpha_long %>% long_to_wide(what="exposures")
    x$nmf[[tid]]$pyro$params$infered_params$beta_d = dn_long %>% long_to_wide(what="beta")

    x$nmf[[tid]]$pyro$params$init_params$alpha = x$nmf[[tid]]$pyro$params$init_params$alpha %>%
      wide_to_long(what="exposures") %>% dplyr::mutate(sigs=mapp[sigs]) %>% long_to_wide(what="exposures")
    x$nmf[[tid]]$pyro$params$init_params$beta_dn_param = x$nmf[[tid]]$pyro$params$init_params$beta_dn_param %>%
      wide_to_long(what="beta") %>% dplyr::mutate(sigs=mapp[sigs]) %>% long_to_wide(what="beta")
  }

  if (!have_groups(x)) {
    new_colnames = colnames(x$clustering$pyro$params$infered_params$alpha_prior)
    for (new_name in names(map_names)) {
      old_name = map_names[[new_name]]
      new_colnames = new_colnames %>%
        stringr::str_replace_all(pattern=old_name, replacement=new_name)
    }

    colnames(x$clustering$pyro$params$infered_params$alpha_prior) =
      colnames(x$clustering$pyro$params$init_params$alpha_prior) =
      colnames(x$clustering$pyro$params$init_params$variances) =
      new_colnames

    x$clustering$centroids = x$clustering$pyro$params$infered_params$alpha_prior %>%
      wide_to_long(what="exposures") %>%
      dplyr::rename(clusters=samples) %>%
      dplyr::mutate(clusters=paste0("G",as.integer(clusters)-1))
  }

  return(x)
}


convert_dn_names = function(x, x.simul) {
  assigned_missing = get_assigned_missing(x, x.simul)

  map_names = lapply(names(assigned_missing), function(tid) {
    am_t = assigned_missing[[tid]]
    # names -> reference names; values -> fit names
    c(am_t$assigned_tp, am_t$added_fp %>% setNames(am_t$added_fp))
  }) %>% unlist()

  rename_object(x, map_names)
}


compute_penalized_bic = function(x) {
  N = get_input(x, matrix=T)[[1]] %>% nrow()
  Ks = get_tested_Ks(x)
  seeds = get_tested_seeds(x)
  lapply(get_types(x), function(tid) {
    lapply(Ks[[tid]], function(k_t) {
      lapply(seeds[[tid]], function(sid) {
        x_k = get_alternative_run(x, K=list(k_t) %>% setNames(tid),
                                  seed=list("nmf"=list(sid) %>% setNames(tid)))
        dn_s = get_denovo_signatures(x_k, matrix=T, types=tid)[[tid]]
        data.frame(K=k_t, seed=sid, score=compute_penalty(dn_s), type=tid, score_id="new_penalty")
      })
    })
  })
}


get_tested_Ks = function(x, types=get_types(x)) {
  lapply(types, function(tid) {
    get_scores(x) %>% dplyr::filter(type==tid, parname=="K") %>%
      dplyr::pull(value) %>% unique()
  }) %>% setNames(types)
}


get_tested_seeds = function(x, types=get_types(x)) {
  lapply(types, function(tid) {
    get_scores(x) %>%
      dplyr::filter(type==tid) %>%
      dplyr::pull(seed) %>% unique()
  }) %>% setNames(types)
}


compute_penalty = function(denovo_signatures) {
  if (is.null(denovo_signatures) || nrow(denovo_signatures)<=1) return(0)
  cosine_s = lsa::cosine(t(denovo_signatures %>% as.matrix()))
  penalty = sum(cosine_s[lower.tri(cosine_s, diag=F)]) * N
  return(penalty)
}


compute_cross_entropy_aux = function(p, q) {
  ce = 0
  for (i in 1:length(p))
    ce = ce + (p[i] * log(q[i]))
  return(-ce)
}
