# add args for pyro fit
pyro_nmf = function(counts, k, cluster, type) {
  x = list(); class(x) = "basilica_obj"
  pyro_fit = nmf_single_type(...)
  nmf_i = list(pyro=pyro_fit, exposure=pyro_fit$exposure, beta=pyro_fit$beta)
  x$nmf[[type]] = nmf_i
  return(x)
}


nmf_single_type = function(..., k_list, clusters, reference_catalogue, stage,
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


