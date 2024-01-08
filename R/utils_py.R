get_list_from_py = function(py_obj, filter_dn=FALSE, type="") {
  if (is.null(py_obj)) return(NULL)

  py_obj = rename_denovo_py(py_obj, type=type)

  x = get_list_from_py_aux(py_obj, fn=get_list_from_py, type=type)
  x$exposure = py_obj$params$alpha %>% wide_to_long(what="exposures")
  x$beta_denovo = py_obj$params$beta_d %>% wide_to_long(what="beta") %>%
    renormalize_denovo_thr(filter_dn=filter_dn)
  x$beta_fixed = py_obj$params$beta_f %>% wide_to_long(what="beta")

  return(x)
}


rename_denovo_py = function(py_obj, type) {
  if (is.null(py_obj$params$beta_d)) return(py_obj)

  new_dn_names = paste0(type, rownames(py_obj$params$beta_d))

  rownames(py_obj$params$beta_d) = new_dn_names
  colnames(py_obj$params$alpha)[grepl("^D[0-9]*$", colnames(py_obj$params$alpha))] =
    new_dn_names

  rownames(py_obj$init_params$beta_dn_param) = new_dn_names
  colnames(py_obj$init_params$alpha)[grepl("^D[0-9]*$", colnames(py_obj$init_params$alpha))] =
    new_dn_names

  return(py_obj)
}


rename_dn_expos = function(x) {
  map_names = c()
  for (tid in get_types(x)) {
    signames = get_signames(x, types=tid)[[tid]]
    expos = get_exposure(x, types=tid, matrix=T)[[tid]]

    if (all(signames == colnames(expos))) {
      map_names = c(map_names, signames) %>% setNames(c(names(map_names), signames))
      next
    }

    map_names = c(map_names, signames) %>% setNames(c(names(map_names), colnames(expos)))
    colnames(expos) = signames

    x$nmf[[tid]]$exposure = wide_to_long(expos, what="exposures")
    x$nmf[[tid]]$pyro$params$infered_params$alpha = expos
    colnames(x$nmf[[tid]]$pyro$params$init_params$alpha) = signames
  }

  alpha_prior_names = colnames(x$clustering$pyro$params$infered_params$alpha_prior)
  new_names = data.frame(alpha_prior_names) %>%
    tidyr::separate(alpha_prior_names, into=c("var_id", "old_sigs")) %>%
    dplyr::mutate(new_sig=ifelse(old_sigs %in% names(map_names), map_names[old_sigs], old_sigs)) %>%
    dplyr::mutate(new_signame=paste(var_id, new_sig, sep="_"), old_signame=paste(var_id, old_sigs, sep="_"))
  map_names2 = new_names$new_signame %>% setNames(new_names$old_signame)
  map_names2 = new_names$new_sig %>% setNames(new_names$old_sigs)

  colnames(x$clustering$pyro$params$infered_params$alpha_prior) = map_names2
  colnames(x$clustering$pyro$params$init_params$alpha_prior) = map_names2
  x$clustering$centroids = x$clustering$centroids %>% dplyr::mutate(sigs=map_names2[sigs])

  return(x)
}


get_list_from_py_clustering = function(py_obj, type="") {
  if (is.null(py_obj)) return(NULL)
  x = list()
  x$pyro = get_list_from_py_aux(py_obj, fn=get_list_from_py_clustering)
  x$clusters = tibble::tibble(samples = rownames(py_obj$alpha),
                              clusters = paste0("G",py_obj$groups))
  x$centroids = py_obj$params$alpha_prior %>%
    wide_to_long(what="exposures") %>%
    dplyr::rename(clusters=samples) %>%
    dplyr::mutate(clusters=paste0("G", as.integer(clusters)-1))
  return(x)
}


get_list_from_py_aux = function(py_obj, fn, type="") {
  x = list()
  x$params = list(infered_params = py_obj$params,
                  init_params = py_obj$init_params,
                  hyperparameters = py_obj$hyperparameters)

  x$QC = get_QC_from_py(py_obj)

  x$alternatives = get_alternatives_from_py(py_obj, fn=fn, type=type)

  try(expr = { x$seed = py_obj$seed })

  return(x)
}


get_train_params_py = function(obj) {
  if (!obj$store_parameters || is.null(obj$params)) return(NULL)
  train_params = obj[["train_params"]]

  params = data.frame()
  lapply(1:length(train_params), function(i) {
    pars_i = train_params[[i]]

    expos = centroid = pi = sigs = data.frame()
    if (!is.null(pars_i[["alpha"]])) {
      expos = wide_to_long(reticulate::py_to_r(pars_i[["alpha"]]), what="exposures") %>%
        dplyr::mutate(parname="alpha") %>%
        dplyr::rename(rowname=samples, columnname=sigs)
    }

    if (!is.null(pars_i[["alpha_prior"]])) {
      centroid = wide_to_long(reticulate::py_to_r(pars_i[["alpha_prior"]]), what="exposures") %>%
        dplyr::mutate(parname="centroid") %>%
        dplyr::rename(rowname=samples, columnname=sigs)
    }

    if (!is.null(pars_i[["pi"]]))
      pi = data.frame("rowname"=sort(unique(centroid$rowname)),
                      "columnname"="",
                      "value"= pars_i[["pi"]] %>% as.numeric() %>% setNames((sort(unique(centroid$rowname)))),
                      "parname"="pi")

    if (!is.null(pars_i[["beta_d"]])) {
      sigs = wide_to_long(reticulate::py_to_r(pars_i[["beta_d"]]), what="beta") %>%
        dplyr::mutate(parname="beta_denovo") %>%
        dplyr::rename(rowname=sigs, columnname=features)
    }


    dplyr::bind_rows(expos, centroid, pi, sigs) %>%
      dplyr::mutate(iteration=i)
  }) %>% do.call(rbind, .) %>% tibble::as_tibble()

}


get_fits_from_py = function(fits_list, fn) {
  return(
    lapply(names(fits_list), function(i) {
      py_obj = fits_list[[i]]
      if ("x" %in% names(py_obj)) {inp = py_obj$x} else {inp = py_obj$alpha}
      py_obj$convert_to_dataframe(inp)
      fn(py_obj)
    }) %>%
      setNames(names(fits_list))
  )
}


get_scores_from_py = function(scores) {
  if (is.null(scores)) return(NULL)

  res = replace_null(scores) %>%
    as.data.frame(optional=TRUE, check_names=FALSE) %>%
    reshape2::melt(value.name="score")
  parname = ifelse(grepl("k_denovo", res$variable[1]), "K", "G")

  res = res %>%
    tidyr::separate("variable", into=c(parname, "seed", "score_id"), sep="[.]") %>%
    tibble::as_tibble() %>%
    dplyr::select_if(dplyr::where(function(i) any(!is.na(i)))) %>%
    dplyr::mutate(dplyr::across(is.character, function(i)
      stringr::str_replace_all(i, "k_denovo:|cluster:|seed:", ""))) %>%
    tidyr::pivot_longer(cols=parname, names_to="parname")

  return(res)
}


replace_null = function(i, value=NA) {
  j = purrr::map(i, ~ replace(.x, is.null(.x), value))
  purrr::map(j, ~ (if(is.list(.x)) replace_null(.x) else .x))
}


get_QC_from_py = function(py_obj) {
  QC = list(lr = py_obj$lr,
            n_steps = py_obj$n_steps,
            bic = py_obj$bic,
            losses = py_obj$losses,
            penalty = py_obj$regs,
            likelihood = py_obj$likelihoods,
            gradient_norms = py_obj$gradient_norms,
            train_params = get_train_params_py(py_obj))

  if ("scores" %in% names(py_obj))
    QC$scores = get_scores_from_py(py_obj$scores)

  return(QC)
}


get_alternatives_from_py = function(py_obj, fn, type="") {
  alt = list()
  alt$runs_seed = alt$runs_scores = alt$all_fits = NULL

  if ("x" %in% names(py_obj)) {inp = py_obj$x} else {inp = py_obj$alpha}
  if ("fits" %in% names(py_obj)) {
    # print(py_obj$fits)
    # fits_tibble = py_obj$fits %>% tibble::as_tibble()
    alt$fits = lapply(names(py_obj$fits), function(i) {
      fits_i = py_obj$fits[[i]]
      lapply(names(fits_i), function(j) {
        fits_i[[j]] # $convert_to_dataframe(inp)
        tmp = tibble::tibble(V1 = list( fn(fits_i[[j]], type=type) ))
        colnames(tmp) = j
        return(tmp)
      }) %>% dplyr::bind_cols()
    }) %>% setNames(names(py_obj$fits))
  }

  return(alt)
}


