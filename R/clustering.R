pyro_clustering = function(exposures, cluster, lr=0.005, n_steps=3000,
                           optim_gamma=0.1, enumer="parallel",
                           hyperparameters=NULL, nonparametric=TRUE,
                           store_parameters=FALSE, store_fits=FALSE,
                           seed_list=c(10), CUDA=TRUE, py=NULL) {

  if (is.null(cluster)) return(NULL)

  TIME = as.POSIXct(Sys.time(), format = "%H:%M:%S")

  if (is.null(py)) py = reticulate::import("pybasilica")
  cluster = reticulate::r_to_py(list(as.integer(cluster)))
  seed_list = reticulate::r_to_py(as.integer(seed_list))
  input_expos = reticulate::r_to_py(exposures %>% setNames(NULL))

  obj = py$fit(alpha=input_expos, cluster=cluster, n_steps=n_steps,
               lr=lr, optim_gamma=optim_gamma, hyperparameters=hyperparameters,
               enumer=enumer, nonparametric=nonparametric, seed_list=seed_list,
               CUDA=CUDA, store_parameters=store_parameters, store_fits=store_fits)

  TIME = difftime(as.POSIXct(Sys.time(), format = "%H:%M:%S"), TIME, units = "mins")

  if (is.list(obj)) {
    bestRun = obj[[1]]
    secondBest = obj[[2]]
  } else {
    bestRun = obj
    secondBest = NULL
  }

  # save python object data in a list
  pyro_fit = get_list_from_py_clustering(bestRun)
  pyro_fit$pyro$alternatives$secondBest = get_list_from_py(secondBest, filter_dn=filter_dn)
  pyro_fit$pyro$time = TIME

  # clustering = list(pyro=pyro_fit,
  #                   clusters=pyro_fit$clusters,
  #                   centroids=pyro_fit$centroids)

  return(pyro_fit)

}

