fit = function(counts, k, cluster) {

}


fit_single_type = function(..., type, k_list, clusters, reference_catalogue, stage,
                           cohort, filtered_catalogue, min_exposure, keep_sigs) {
  TIME = as.POSIXct(Sys.time(), format = "%H:%M:%S")
  call_info = match.call()

  if (!is.null(reference_catalogue)) {
    x_ref = pyfit(...,
                  k_list = 0, clusters = NULL,
                  beta_fixed = reference_catalogue,
                  stage = "random_noise") %>%

      create_basilica_obj(cohort=cohort,
                          filtered_catalogue=filtered_catalogue)

    x_ref_filt = x_ref %>% filter_sigs_low_expos(min_exp=min_exposure, keep_sigs=keep_sigs)
    catalogue2 = get_signatures(x_ref_filt)
  } else {
    x_ref = x_ref_filt = catalogue2 = NULL
    residues = FALSE
  }

  x_dn = pyfit(...,
               k_list = k, clusters = clusters,
               beta_fixed = catalogue2, stage = stage) %>%

    create_basilica_obj(cohort=cohort,
                        filtered_catalogue=filtered_catalogue)

  TIME = difftime(as.POSIXct(Sys.time(), format = "%H:%M:%S"), TIME, units = "mins")

  x_dn$time = TIME
  x_dn$k_list = k
  x_dn$call = call_info

  return(x_dn)
}


pyfit = function(counts,
                 k_list,
                 lr = 0.005,
                 optim_gamma = 0.1,
                 n_steps = 2000,
                 stage = "",
                 py = NULL,
                 clusters = NULL,
                 nonparametric = TRUE,
                 dirichlet_prior = TRUE,
                 beta_fixed = NULL,
                 hyperparameters = NULL,
                 CUDA = FALSE,
                 compile = FALSE,
                 enforce_sparsity = TRUE,
                 store_parameters = FALSE,
                 regularizer = "cosine",
                 regul_compare = NULL,
                 reg_weight = 1,
                 seed_list = c(10),
                 regul_denovo = TRUE,
                 regul_fixed = TRUE,
                 save_all_fits = FALSE
                 ) {

  TIME = as.POSIXct(Sys.time(), format = "%H:%M:%S")

  if (is.null(py))
    py = reticulate::import("pybasilica")

  if (length(k_list) > 1) k_list = reticulate::r_to_py(as.integer(k_list)) else
    k_list = reticulate::r_to_py(list(as.integer(k_list)))

  if (length(seed_list) > 1) seed_list = reticulate::r_to_py(as.integer(seed_list)) else
    seed_list = reticulate::r_to_py(list(as.integer(seed_list)))

  if (!is.null(clusters)) clusters = as.integer(clusters)

  obj = py$fit(x = counts, k_list = k_list, lr = lr, optim_gamma = optim_gamma, n_steps = n_steps,
               cluster = clusters, beta_fixed = beta_fixed,
               hyperparameters = hyperparameters, nonparametric=nonparametric,
               dirichlet_prior = dirichlet_prior, enforce_sparsity = enforce_sparsity,
               store_parameters = store_parameters, regularizer = regularizer,
               reg_weight = reg_weight, regul_compare = regul_compare,
               regul_denovo = regul_denovo, regul_fixed = regul_fixed,
               stage = stage, seed = seed_list, compile_model = compile,
               CUDA = CUDA, save_all_fits = save_all_fits)

  TIME = difftime(as.POSIXct(Sys.time(), format = "%H:%M:%S"), TIME, units = "mins")

  if (is.list(obj)) {
    bestRun = obj[[1]]
    secondBest = obj[[2]]
  } else {
    bestRun = obj
    secondBest = NULL
  }

  # save python object data in a list
  py_obj = get_list_from_py(bestRun, type=type)
  py_obj$runs_seed = lapply(py_obj$runs_seed, function(i) {
    i[["runs_scores"]] = i[["runs_seed"]] = NULL
    i$convert_to_dataframe(counts)
    get_list_from_py(i, type=type)
  })

  py_obj$secondBest = get_list_from_py(secondBest, type=type)
  py_obj$time = TIME

  return(py_obj)
}


create_basilica_obj = function(py_obj, type="T1", cohort="MyCohort",
                               filtered_catalogue=TRUE) {
  obj = list()
  class(obj) = "basilica_obj"
  obj$cohort = cohort

  obj$input[[type]] = list("counts"=py_obj$x,
                           "fixed_signatures"=py_obj$fixed_signatures)

  if (filtered_catalogue)
    py_obj$denovo_signatures = renormalize_denovo_thr(py_obj$denovo_signatures)

  obj$fit[[type]] = py_obj
  obj$groups = py_obj$groups
  obj$color_palette = gen_palette(obj)
  obj$time = py_obj$time

  return(obj)
}


renormalize_denovo_thr = function(denovo, thr=0.02) {
  if (is.null(denovo)) return(NULL)

  return(
    denovo %>%
      dplyr::mutate(value=replace(value, value<thr, 0)) %>%
      dplyr::group_by(sigs) %>%
      dplyr::mutate(value=value/sum(value)) %>%
      dplyr::ungroup()
  )
}


gen_palette = function(x, type) {
  ref = COSMIC_color_palette(catalogue=get_fixed_signatures(x, type=type, wide=T)) %>%
    setNames(get_fixed_signames(x, type=type))
  dn = ggsci::pal_simpsons()(length(get_denovo_signames(x))) %>%
    setNames(get_denovo_signames(x))
  return(c(ref, dn))
}


COSMIC_color_palette = function(catalogue=COSMIC_filt, seed=14) {
  N = nrow(catalogue)
  set.seed(seed)
  colss = Polychrome::createPalette(N, c("#856de3","#9e461c"), target="normal", range=c(15,80), M=1000)[1:N]
  names(colss) = rownames(catalogue)
  return(colss)
}




get_list_from_py = function(py_obj, type=NULL) {
  stopifnot(!is.null(py_obj))
  if (is.null(py_obj)) return(NULL)

  x = list()
  x$x = py_obj$x %>% wide_to_long_counts(type=type)
  x$exposure = py_obj$params$alpha %>% wide_to_long_alpha(type=type)
  x$denovo_signatures = py_obj$params$beta_d %>% wide_to_long_beta(type=type)
  x$fixed_signatures = py_obj$params$beta_f %>% wide_to_long_beta(type=type)
  x$eps_var = py_obj$params$lambda_epsilon
  x$pi = py_obj$params$pi
  x$post_probs = py_obj$params$post_probs
  x$groups = py_obj$groups

  x$params = py_obj$params[!names(py_obj$params) %in% c("alpha","beta_d","beta_f")]
  x$init_params = py_obj$init_params

  x$bic = py_obj$bic
  x$losses = py_obj$losses
  x$gradient_norms = py_obj$gradient_norms
  x$train_params = get_train_params(py_obj)
  x$hyperparameters = py_obj$hyperparameters
  try(expr = { x$seed = py_obj$seed })

  x$runs_seed = x$runs_scores = x$all_fits = NULL
  if ("runs_seed" %in% names(py_obj))
    x$runs_seed = py_obj$runs_seed

  if ("scores_K" %in% names(py_obj))
    x$runs_K = get_scores_from_py(py_obj$scores_K)

  if ("scores_CL" %in% names(py_obj))
    x$runs_CL = get_scores_from_py(py_obj$scores_CL) %>% dplyr::rename(G=K)

  if ("all_fits" %in% names(py_obj)) {
    if (length(py_obj$all_fits) > 0) x$all_fits = NULL
    x$all_fits = get_fits_from_py(py_obj$all_fits, x$x, x$input_catalogue, lr, n_steps)
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


