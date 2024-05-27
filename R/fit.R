#' Fit a basilica object
#'
#' @param counts List of mutation counts matrices from multiple variant types.
#' @param k_list List of number of denovo signatures to test.
#' @param cluster Maximum number of clusters. If `NULL`, no clustering will be performed.
#' @param reference_cat List of reference catalogues to use for NMF. Names must be the same as input counts.
#' @param keep_sigs List of reference signatures to keep even if found with low exposures.
#' @param hyperparameters List of hyperparameters passed to the NMF and clustering models.
#' @param lr Learning rate used for SVI.
#' @param optim_gamma Deprecated
#' @param n_steps Number of iterations for inference.
#' @param py User-installed version of \code{pybasilica} package
#' @param enumer Enumeration used for clustering (either `parallel` or `sequential`).
#' @param nonparametric Deprecated. The model only works in nonparametric way.
#' @param autoguide Logical. If `TRUE`, the clustering model will use the Pyro autoguide.
#' @param filter_dn Logical. If `TRUE`, all contexts below 0.01 in denovo signatures will be set to 0, provided the filtered signatures remain consistent with the inferred ones.
#' @param min_exposure Reference signatures with an exposures lower than `min_exposure` will be dropped.
#' @param CUDA Logical. If `TRUE` and a GPU is available, the models will run on GPU.
#' @param compile Deprecated.
#' @param store_parameters Logical. If `TRUE`, parameters at every step of inference will be stored in the object.
#' @param store_fits Logical. If `TRUE`, all tested fits, i.e., for every value of `K`, will be stored in the object.
#' @param seed_list List of seeds used for every input configuration.
#'
#' @return Basilica object.
#' @export fit
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
               store_fits = TRUE,

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
  return(bas)
}



#' Fit clustering
#'
#' @param x Basilica object with signatures deconvolution performed.
#' @param cluster Maximum number of clusters.
#' @param hyperparameters List of hyperparameters passed to the NMF and clustering models.
#' @param lr Learning rate for SVI optimizer.
#' @param optim_gamma Deprecated.
#' @param n_steps Number of steps for the inference.
#' @param py User-installed version of \code{pybasilica} package
#' @param enumer Enumeration used for clustering (either `parallel` or `sequential`).
#' @param nonparametric Deprecated. The model only works in nonparametric way.
#' @param autoguide Logical. If `TRUE`, the clustering model will use the Pyro autoguide.
#' @param CUDA Logical. If `TRUE` and a GPU is available, the models will run on GPU.
#' @param compile Deprecated.
#' @param store_parameters Logical. If `TRUE`, parameters at every step of inference will be stored in the object.
#' @param store_fits Logical. If `TRUE`, all tested fits, i.e., for every value of `K`, will be stored in the object.
#' @param seed_list List of seeds used for every input configuration.
#'
#' @return Basilica object.
#' @export fit_clustering
fit_clustering = function(x,
                          cluster,
                          hyperparameters = NULL,

                          lr = 0.005,
                          optim_gamma = 0.1,
                          n_steps = 3000,
                          py = NULL,

                          enumer = "parallel",
                          nonparametric = TRUE,
                          autoguide = TRUE,

                          CUDA = TRUE,
                          compile = FALSE,

                          store_parameters = FALSE,
                          store_fits = TRUE,

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



