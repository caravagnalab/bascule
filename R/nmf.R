# add args for pyro fit; returns nmf for type t
pyro_nmf = function(..., k_list, reference_cat, stage, cohort,
                    filter_dn, min_exposure, keep_sigs, type="") {

  pyro_fit = nmf_single_type(..., k_list=k_list, reference_cat=reference_cat,
                             stage=stage, cohort=cohort, filter_dn=filter_dn,
                             min_exposure=min_exposure, keep_sigs=keep_sigs,
                             type=type)

  nmf_t = list(pyro=pyro_fit, exposure=pyro_fit$exposure,
               beta_fixed=pyro_fit$beta_fixed, beta_denovo=pyro_fit$beta_denovo)
  return(nmf_t)
}


nmf_single_type = function(..., k_list, reference_cat, stage, cohort,
                           filter_dn, min_exposure, keep_sigs, type="") {
  TIME = as.POSIXct(Sys.time(), format = "%H:%M:%S")
  call_info = match.call()

  if (!is.null(reference_cat)) {
    x_ref = pyfit(..., k_list=0, clusters=NULL, reference_cat=reference_cat,
                  stage="random_noise", filter_dn=FALSE, type=type)

    x_ref_filt = x_ref %>% filter_sigs_low_expos(min_exp=min_exposure, keep_sigs=keep_sigs)
    catalogue2 = long_to_wide(dataframe=x_ref_filt$beta_fixed, what="beta")
  } else {
    x_ref = x_ref_filt = catalogue2 = NULL
    residues = FALSE
  }

  x_dn = pyfit(..., k_list=k_list, clusters=NULL, stage=stage,
               filter_dn=filter_dn, reference_cat=catalogue2, type=type)

  TIME = difftime(as.POSIXct(Sys.time(), format = "%H:%M:%S"), TIME, units = "mins")

  x_dn$time = TIME
  x_dn$k_list = k_list
  x_dn$call = call_info

  return(x_dn)
}


set_attribute = function(x, what, type, name, value) {
  if (what == "nmf") x[[what]][[type]][[name]] = value
  if (what == "clustering") x[[what]][[name]] = value
  return(x)
}


filter_denovo = function(x, types=get_types(x), thr=0.02) {
  for (tid in types) {
    denovo = get_denovo_signatures(x, types=tid, matrix=F)[[tid]] %>%
      renormalize_denovo_thr(thr=thr, filter_dn=T)
    x = set_attribute(x, what="nmf", type=tid, name="beta_denovo", value=denovo)
  }
  return(x)
}


renormalize_denovo_thr = function(beta_denovo, thr=0.02, filter_dn=FALSE) {
  if (is.null(beta_denovo)) return(NULL)
  if (!filter_dn) return(beta_denovo)

  return(
    beta_denovo %>%
      dplyr::mutate(value=replace(value, value<thr, 0)) %>%
      dplyr::group_by(sigs) %>%
      dplyr::mutate(value=value/(sum(value)+1e-10)) %>%
      dplyr::ungroup()
  )
}




