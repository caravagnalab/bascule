# counts is a list with input matrices, names are typenames
fit = function(counts) {
  if (!is.list(counts)) counts = list("T1"=counts)

  bas = list(); class(bas) = "basilica_obj"

  # input contains a list of counts matrices
  bas$input = lapply(counts, wide_to_long, what="counts")

  # nmf contains the pyro fits
  types = names(counts)
  bas$nmf = lapply(types, function(t)
    pyro_nmf(counts=counts[[t]],
             k_list=k_list,
             lr=lr,
             optim_gamma = optim_gamma,
             n_steps = n_steps,
             stage = stage,
             py = py,
             dirichlet_prior = dirichlet_prior,
             beta_fixed = beta_fixed,
             keep_sigs = keep_sigs,
             min_exposure = min_exposure,
             hyperparameters = hyperparameters,
             CUDA = CUDA,
             compile = compile,
             enforce_sparsity = enforce_sparsity,
             store_parameters = store_parameters,
             regularizer = regularizer,
             regul_compare = regul_compare,
             reg_weight = reg_weight,
             seed_list = seed_list,
             regul_denovo = regul_denovo,
             regul_fixed = regul_fixed,
             save_all_fits = save_all_fits,
             filter_dn = filter_dn)
  ) %>% setNames(types)

  # clustering contains the clustering
  expos = get_exposure(bas, matrix=TRUE)
  bas$clustering = pyro_clustering(exposures=expos, type=t, ...)

  return(bas)
}


# add args for pyro fit; returns nmf for type t
pyro_nmf = function(..., k_list, beta_fixed, stage, cohort,
                    filter_dn, min_exposure, keep_sigs) {

  pyro_fit = nmf_single_type(..., k_list=k_list, beta_fixed=beta_fixed,
                             stage=stage, cohort=cohort, filter_dn=filter_dn,
                             min_exposure=min_exposure, keep_sigs=keep_sigs)

  nmf_t = list(pyro=pyro_fit, exposure=pyro_fit$exposure,
               beta_fixed=pyro_fit$beta_fixed, beta_denovo=pyro_fit$beta_denovo)
  # x$nmf[[type]] = nmf_t
  return(nmf_t)
}


nmf_single_type = function(..., k_list, beta_fixed, stage, cohort,
                           filter_dn, min_exposure, keep_sigs) {
  TIME = as.POSIXct(Sys.time(), format = "%H:%M:%S")
  call_info = match.call()

  if (!is.null(beta_fixed)) {
    x_ref = pyfit(..., k_list=0, beta_fixed=beta_fixed,
                  stage="random_noise", filter_dn=FALSE)

    x_ref_filt = x_ref %>% filter_sigs_low_expos(min_exp=min_exposure, keep_sigs=keep_sigs)
    catalogue2 = long_to_wide(dataframe=x_ref_filt$beta_fixed, what="beta")
  } else {
    x_ref = x_ref_filt = catalogue2 = NULL
    residues = FALSE
  }

  x_dn = pyfit(..., k_list=k_list, beta_fixed=catalogue2, stage=stage, filter_dn=filter_dn)

  TIME = difftime(as.POSIXct(Sys.time(), format = "%H:%M:%S"), TIME, units = "mins")

  x_dn$time = TIME
  x_dn$k_list = k_list
  x_dn$call = call_info

  return(x_dn)
}


renormalize_denovo_thr = function(py_obj, thr=0.02, filter_dn=TRUE) {
  if (is.null(py_obj$beta_denovo)) return(NULL)

  py_obj$beta_denovo = py_obj$beta_denovo %>%
      dplyr::mutate(value=replace(value, value<thr, 0)) %>%
      dplyr::group_by(sigs) %>%
      dplyr::mutate(value=value/sum(value)) %>%
      dplyr::ungroup()
  return(py_obj)
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


