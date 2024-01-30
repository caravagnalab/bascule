# counts is a list with input matrices, names are typenames
fit = function(counts, k_list,
               cluster=NULL,
               reference_cat=list("SBS"=COSMIC_filt, "DBS"=COSMIC_dbs),
               keep_sigs = c("SBS1","SBS5"),

               hyperparameters = NULL,

               lr = 0.005,
               optim_gamma = 0.1,
               n_steps = 3000,
               py = NULL,

               enumer = "parallel",
               nonparametric = TRUE,
               autoguide = FALSE,

               filter_dn = FALSE,
               min_exposure = 0.2,
               CUDA = TRUE,
               compile = FALSE,

               store_parameters = FALSE,
               store_fits = FALSE,

               seed_list = c(10)) {
  if (!is.list(counts)) counts = list("T1"=counts)
  if (!is.list(reference_cat)) reference_cat = list("T1"=reference_cat)

  bas = list(); class(bas) = "basilica_obj"

  types = names(counts)

  # input contains a list of counts matrices
  bas$input = lapply(types, function(t) {
    list(counts=wide_to_long(counts[[t]], what="counts"),
         reference=wide_to_long(reference_cat[[t]], what="beta"))
  }) %>% setNames(types)

  # nmf contains the pyro fits
  bas$nmf = lapply(types, function(tid) {
    pyro_nmf(counts = counts[[tid]],
             k_list = k_list,
             reference_cat = reference_cat[[tid]],
             keep_sigs = keep_sigs,

             hyperparameters = hyperparameters,

             lr = lr,
             optim_gamma = optim_gamma,
             n_steps = n_steps,

             filter_dn = filter_dn,
             min_exposure = min_exposure,
             CUDA = CUDA,
             compile = compile,

             store_parameters = store_parameters,
             store_fits = store_fits,

             seed_list = seed_list,
             py = py,
             type = tid)
  }
  ) %>% setNames(types)

  # clustering contains the clustering
  bas = fit_clustering(bas,
                       cluster = cluster,

                       enumer = enumer,
                       nonparametric = nonparametric,
                       autoguide = autoguide,

                       hyperparameters = hyperparameters,

                       lr = lr,
                       optim_gamma = optim_gamma,
                       n_steps = n_steps,

                       CUDA = CUDA,

                       store_parameters = store_parameters,
                       store_fits = store_fits,

                       seed_list = seed_list,
                       py = py)

  # exposures = get_exposure(bas, matrix=TRUE)
  # bas$clustering = pyro_clustering(exposures = exposures,
  #                                  cluster = cluster,
  #
  #                                  enumer = enumer,
  #                                  nonparametric = nonparametric,
  #
  #                                  hyperparameters = hyperparameters,
  #
  #                                  lr = lr,
  #                                  optim_gamma = optim_gamma,
  #                                  n_steps = n_steps,
  #
  #                                  CUDA = CUDA,
  #
  #                                  store_parameters = store_parameters,
  #                                  store_fits = store_fits,
  #
  #                                  seed_list = seed_list,
  #                                  py = py)
  return(bas)
}



fit_clustering = function(x,
                          cluster,
                          hyperparameters = NULL,

                          lr = 0.005,
                          optim_gamma = 0.1,
                          n_steps = 3000,
                          py = NULL,

                          enumer = "parallel",
                          nonparametric = TRUE,
                          autoguide = FALSE,

                          CUDA = TRUE,
                          compile = FALSE,

                          store_parameters = FALSE,
                          store_fits = FALSE,

                          seed_list = c(10)) {
  exposures = get_exposure(x, matrix=TRUE)
  x$clustering = pyro_clustering(exposures = exposures,
                                 cluster = cluster,

                                 enumer = enumer,
                                 nonparametric = nonparametric,
                                 autoguide = autoguide,

                                 hyperparameters = hyperparameters,

                                 lr = lr,
                                 optim_gamma = optim_gamma,
                                 n_steps = n_steps,

                                 CUDA = CUDA,

                                 store_parameters = store_parameters,
                                 store_fits = store_fits,

                                 seed_list = seed_list,
                                 py = py)
  return(x)
}



