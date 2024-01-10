pyro_clustering = function(exposures, cluster, lr=0.005, n_steps=3000,
                           optim_gamma=0.1, enumer="parallel", autoguide=FALSE,
                           hyperparameters=NULL, nonparametric=TRUE,
                           store_parameters=FALSE, store_fits=FALSE,
                           seed_list=c(10), CUDA=TRUE, py=NULL) {

  if (is.null(cluster)) return(NULL)
  if (length(cluster)==1) cluster = c(cluster)

  TIME = as.POSIXct(Sys.time(), format = "%H:%M:%S")

  if (is.null(py)) py = reticulate::import("pybasilica")
  cluster = reticulate::r_to_py(as.integer(cluster))
  seed_list = reticulate::r_to_py(as.integer(seed_list))
  input_expos = reticulate::r_to_py(exposures %>% setNames(NULL))

  obj = py$fit(alpha=input_expos, cluster=cluster, n_steps=n_steps, lr=lr,
               optim_gamma=optim_gamma, hyperparameters=hyperparameters,
               enumer=enumer, autoguide=autoguide, nonparametric=nonparametric,
               seed_list=seed_list, CUDA=CUDA, store_parameters=store_parameters,
               store_fits=store_fits)

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
  pyro_fit$pyro$alternatives$secondBest = get_list_from_py_clustering(secondBest)
  pyro_fit$pyro$time = TIME

  # clustering = list(pyro=pyro_fit,
  #                   clusters=pyro_fit$clusters,
  #                   centroids=pyro_fit$centroids)

  return(pyro_fit)

}



empirical_centroids = function(x) {
  lapply(get_types(x), function(tid) {
    expos = get_exposure(x, add_groups=TRUE)[[tid]]
    expos %>% dplyr::group_by(clusters, sigs) %>%
      dplyr::reframe(centroid=mean(value)) %>%
      dplyr::ungroup() %>%
      tidyr::pivot_wider(values_from="centroid", names_from="sigs") %>%
      tibble::column_to_rownames(var="clusters")
  }) %>% do.call(cbind, .)
}


merge_clusters = function(x, cutoff=0.8) {

  repeat {
    alpha_prior = empirical_centroids(x)
    if (nrow(alpha_prior) == 1) return(x)
    cosine_simil = lsa::cosine(t(alpha_prior)) %>% as.data.frame()
    cosine_simil[lower.tri(cosine_simil, diag=T)] = 0

    rownames(cosine_simil) = colnames(cosine_simil) = rownames(alpha_prior)

    merging = cosine_simil %>% tibble::rownames_to_column(var="gid1") %>%
      reshape2::melt(id="gid1", variable.name="gid2", value.name="cosine") %>%
      dplyr::mutate(gid1=as.character(gid1), gid2=as.character(gid2)) %>%
      dplyr::filter(cosine > cutoff) %>% dplyr::select(-cosine) %>%
      dplyr::group_by(gid1) %>%
      dplyr::summarise(cl_old=list(c(gid1,gid2) %>% unique())) %>% dplyr::ungroup() %>%
      dplyr::rename(cl_name=gid1)

    if (nrow(merging) == 0) return(x)

    grps = get_cluster_assignments(x)
    for (i in 1:nrow(merging)) {
      old_cl = dplyr::pull(merging[i,], cl_old)[[1]]
      new_cl = dplyr::pull(merging[i,], cl_name)
      grps = grps %>% dplyr::mutate(clusters=ifelse(clusters %in% old_cl,
                                                    new_cl, clusters))
    }
    x$clustering$centroids = alpha_prior %>%
      tibble::rownames_to_column(var="clusters") %>%
      reshape2::melt(id="clusters", variable.name="sigs", value.name="value") %>%
      tibble::as_tibble()
    x$clustering$clusters = grps
  }

  return(x)
}





