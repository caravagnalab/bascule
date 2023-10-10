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
  py_obj = get_list_from_py(bestRun)
  py_obj$runs_seed = lapply(py_obj$runs_seed, function(i) {
    i[["runs_scores"]] = i[["runs_seed"]] = NULL
    i$convert_to_dataframe(counts)
    get_list_from_py(i, type=type)
  })

  py_obj$secondBest = get_list_from_py(secondBest, type=type)
  py_obj$time = TIME

  return(py_obj)
}



create_basilica_obj = function(py_obj, cohort="MyCohort",
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
