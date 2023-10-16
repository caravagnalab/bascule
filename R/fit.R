# counts is a list with input matrices, names are typenames
fit = function(counts, k_list, cluster, reference_cat,
               keep_sigs = c("SBS1","SBS5"),

               hyperparameters = NULL,

               lr = 0.005,
               optim_gamma = 0.1,
               n_steps = 3000,
               py = NULL,

               stage = "",  # remove

               enumer = "parallel",
               nonparametric = TRUE,

               dirichlet_prior = TRUE,  # remove
               enforce_sparsity = TRUE,  # remove

               filter_dn = TRUE,
               min_exposure = 0.2,
               CUDA = TRUE,
               compile = FALSE,

               store_parameters = FALSE,
               store_fits = FALSE,

               reg_weight = 0.,
               regularizer = "cosine",
               regul_compare = NULL,
               regul_denovo = TRUE,
               regul_fixed = TRUE,

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
  bas$nmf = lapply(types, function(t)
    pyro_nmf(counts = counts[[t]],
             k_list = k_list,
             reference_cat = reference_cat[[t]],
             keep_sigs = keep_sigs,

             hyperparameters = hyperparameters,

             lr = lr,
             optim_gamma = optim_gamma,
             n_steps = n_steps,

             stage = stage,  # remove

             dirichlet_prior = dirichlet_prior,  # remove
             enforce_sparsity = enforce_sparsity,  # remove

             filter_dn = filter_dn,
             min_exposure = min_exposure,
             CUDA = CUDA,
             compile = compile,

             store_parameters = store_parameters,
             store_fits = store_fits,

             reg_weight = reg_weight,
             regularizer = regularizer,
             regul_compare = regul_compare,
             regul_denovo = regul_denovo,
             regul_fixed = regul_fixed,

             seed_list = seed_list,
             py = py)
  ) %>% setNames(types)

  # clustering contains the clustering
  exposures = get_exposure(bas, matrix=TRUE)
  bas$clustering = pyro_clustering(exposures = exposures,
                                   cluster = cluster,

                                   enumer = enumer,
                                   nonparametric = nonparametric,

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


gen_palette = function(x=NULL, types=get_types(x), n=NULL) {
  if (!is.null(n)) return(ggsci::pal_simpsons()(n))

  denovo_names = unique(unlist(get_denovo_signames(x, types=types)))

  ref = COSMIC_color_palette(signames=get_fixed_signames(x, types=types))
  dn = ggsci::pal_simpsons()(length(denovo_names)) %>%
    setNames(denovo_names)
  return(c(ref, dn))
}


COSMIC_color_palette = function(signames=get_signames(x), seed=14) {
  N = length(unique(unlist(signames)))
  set.seed(seed)
  colss = Polychrome::createPalette(N, c("#856de3","#9e461c"), target="normal", range=c(15,80), M=1000)[1:N]
  names(colss) = unique(unlist(signames))
  return(colss)
}


