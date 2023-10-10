get_list_from_py = function(py_obj) {
  stopifnot(!is.null(py_obj))
  if (is.null(py_obj)) return(NULL)

  x = list()
  # x$x = py_obj$x %>% wide_to_long_counts(type=type)
  x$exposure = py_obj$params$alpha %>% wide_to_long_alpha()
  x$beta_denovo = py_obj$params$beta_d %>% wide_to_long_beta()
  x$beta_fixed = py_obj$params$beta_f %>% wide_to_long_beta()
  x$eps_var = py_obj$params$lambda_epsilon
  x$pi = py_obj$params$pi
  x$post_probs = py_obj$params$post_probs
  x$groups = py_obj$groups

  x$params$infered_params = py_obj$params
  x$params$init_params = py_obj$init_params
  x$params$hyperparameters = py_obj$hyperparameters

  x$QC$bic = py_obj$bic
  x$QC$losses = py_obj$losses
  x$QC$gradient_norms = py_obj$gradient_norms
  x$QC$train_params = get_train_params(py_obj)
  try(expr = { x$seed = py_obj$seed })

  x$QC$runs_seed = x$QC$runs_scores = x$QC$all_fits = NULL
  if ("runs_seed" %in% names(py_obj))
    x$QC$runs_seed = py_obj$runs_seed

  if ("scores_K" %in% names(py_obj))
    x$QC$runs_K = get_scores_from_py(py_obj$scores_K)

  if ("scores_CL" %in% names(py_obj))
    x$QC$runs_CL = get_scores_from_py(py_obj$scores_CL) %>% dplyr::rename(G=K)

  if ("all_fits" %in% names(py_obj)) {
    if (length(py_obj$all_fits) > 0) x$all_fits = NULL
    x$QC$all_fits = get_fits_from_py(py_obj$all_fits, x$x, x$input_catalogue, lr, n_steps)
  }

  return(x)
}



get_train_params = function(obj) {
  if (!obj$store_parameters)
    return(NULL)
  train_params = obj$train_params
  samples_names = obj$params[["alpha"]] %>% rownames()
  bfixed_names = obj$beta_fixed %>% rownames()
  bdenovo_names = obj$params[["beta_d"]] %>% rownames()
  contexts = obj$params[["beta_d"]] %>% colnames()

  params = data.frame()

  for (i in 1:length(train_params)) {
    expos = train_params[[i]][["alpha"]] %>% as.data.frame()
    rownames(expos) = samples_names
    colnames(expos) = c(bfixed_names, bdenovo_names)

    if ("alpha_prior" %in% names(train_params[[i]])) {
      centroids = train_params[[i]][["alpha_prior"]] %>% as.data.frame()
      rownames(centroids) = (1:nrow(centroids)) -1
      colnames(centroids) = c(bfixed_names, bdenovo_names)
      centroids = centroids %>% tibble::rownames_to_column(var="rowname") %>%
        reshape2::melt(id="rowname",variable.name="columnname",value.name="value") %>%
        dplyr::mutate(iteration=i, paramname="centroid")
    } else { centroids = data.frame() }

    if ("pi" %in% names(train_params[[i]])) {
      pi = train_params[[i]][["pi"]] %>% as.numeric() %>% setNames((sort(unique(centroids$rowname))))
      pi = data.frame("rowname"=names(pi),"value"=pi,"iteration"=i,"paramname"="pi")
    } else { pi = data.frame() }

    sigs = train_params[[i]][["beta_d"]] %>% as.data.frame()
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
      dplyr::bind_rows(pi)

  }

  return(params)
}


get_fits_from_py = function(fits, counts, beta_fixed, lr, n_steps)
  return(
    lapply(names(fits), function(i) {
      fits[[i]]$convert_to_dataframe(counts)
      get_list_from_py(fits[[i]], counts, beta_fixed, lr, n_steps, save_stats=F)
    }) %>%
      setNames(names(fits))
  )


get_scores_from_py = function(scores) {
  if (is.null(scores)) return(NULL)

  res = replace_null(scores) %>%
    # res = purrr::discard(scores, is.null) %>%
    as.data.frame() %>%
    reshape2::melt(value.name="score") %>%
    tidyr::separate("variable", into=c("K", "seed", "score_id"), sep="[.]") %>%
    dplyr::select_if(function(i) any(!is.na(i))) %>%
    tibble::as_tibble()

  return(res)
}


replace_null = function(i) {
  j = purrr::map(i, ~ replace(.x, is.null(.x), NA))
  purrr::map(j, ~ (if(is.list(.x)) replace_null(.x) else .x))
}


