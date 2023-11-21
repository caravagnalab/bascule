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
  for (tid in get_types(x)) {
    signames = get_signames(x, types=tid)[[tid]]
    expos = get_exposure(x, types=tid, matrix=T)[[tid]]
    colnames(expos) = signames

    x$nmf[[tid]]$exposure = wide_to_long(expos, what="exposures")
    x$nmf[[tid]]$pyro$params$infered_params$alpha = expos
    colnames(x$nmf[[tid]]$pyro$params$init_params$alpha) = signames
  }

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
    dplyr::rename(clusters=samples)
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
  if (!obj$store_parameters) return(NULL)
  train_params = obj$train_params
  samples_names = obj$params[["alpha"]] %>% rownames()
  bfixed_names = obj$beta_fixed %>% rownames()
  bdenovo_names = obj$params[["beta_d"]] %>% rownames()
  contexts = obj$params[["beta_d"]] %>% colnames()

  params = data.frame()

  for (i in 1:length(train_params)) {
    pars_i = train_params[[i]]
    expos = pars_i[["alpha"]] %>% as.data.frame()
    rownames(expos) = samples_names
    colnames(expos) = c(bfixed_names, bdenovo_names)

    if ("alpha_prior" %in% names(pars_i)) {
      centroids = pars_i[["alpha_prior"]] %>% as.data.frame()
      rownames(centroids) = (1:nrow(centroids)) -1
      colnames(centroids) = c(bfixed_names, bdenovo_names)
      centroids = centroids %>% tibble::rownames_to_column(var="rowname") %>%
        reshape2::melt(id="rowname",variable.name="columnname",value.name="value") %>%
        dplyr::mutate(iteration=i, paramname="centroid")
    } else { centroids = data.frame() }

    if (("beta_w" %in% names(pars_i))) {
      beta_w = pars_i[["beta_w"]] %>% as.data.frame()
      if (nrow(beta_w)>0) {
        rownames(beta_w) = bdenovo_names
        colnames(beta_w) = c(bfixed_names, "DN")
        beta_w = beta_w %>% tibble::rownames_to_column(var="rowname") %>%
          reshape2::melt(id="rowname",variable.name="columnname",value.name="value") %>%
          dplyr::mutate(iteration=i, paramname="beta_w")
      } else { beta_w = data.frame() }
    } else { beta_w = data.frame() }

    if ("pi" %in% names(pars_i)) {
      pi = pars_i[["pi"]] %>% as.numeric() %>% setNames((sort(unique(centroids$rowname))))
      pi = data.frame("rowname"=names(pi),"value"=pi,"iteration"=i,"paramname"="pi")
    } else { pi = data.frame() }

    sigs = pars_i[["beta_d"]] %>% as.data.frame()
    rownames(sigs) = bdenovo_names
    colnames(sigs) = contexts

    params = params %>% dplyr::bind_rows(
      expos %>% tibble::rownames_to_column(var="rowname") %>%
        reshape2::melt(id="rowname",variable.name="columnname",value.name="value") %>%
        dplyr::mutate(iteration=i, paramname="alpha")
    ) %>% dplyr::bind_rows(
      sigs %>% tibble::rownames_to_column(var="rowname") %>%
        reshape2::melt(id="rowname",variable.name="columnname",value.name="value") %>%
        dplyr::mutate(iteration=i, paramname="beta_d")
    ) %>% dplyr::bind_rows(centroids) %>%
      dplyr::bind_rows(pi) %>%
      dplyr::bind_rows(beta_w)
  }

  return(params)
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


replace_null = function(i) {
  j = purrr::map(i, ~ replace(.x, is.null(.x), NA))
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
        fits_i[[j]]$convert_to_dataframe(inp)
        tmp = tibble::tibble(V1 = list( fn(fits_i[[j]], type=type) ))
        colnames(tmp) = j
        return(tmp)
      }) %>% dplyr::bind_cols()
    }) %>% setNames(names(py_obj$fits))
  }

  return(alt)
}


