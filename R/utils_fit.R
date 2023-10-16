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
                 reference_cat = NULL,
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
                 store_fits = TRUE,
                 filter_dn = TRUE
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
               hyperparameters = hyperparameters, nonparametric=nonparametric,
               dirichlet_prior = dirichlet_prior, enforce_sparsity = enforce_sparsity,
               store_parameters = store_parameters, regularizer = regularizer,
               reg_weight = reg_weight, regul_compare = regul_compare,
               regul_denovo = regul_denovo, regul_fixed = regul_fixed,
               stage = stage, seed_list = seed_list, compile_model = compile,
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
  py_obj = get_list_from_py(bestRun, filter_dn=filter_dn)
  py_obj$alternatives$secondBest = get_list_from_py(secondBest, filter_dn=filter_dn)
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




