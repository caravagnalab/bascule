# Set parameters and dataframes #####

set_denovo_signatures = function(x, sigs, type) {
  x$nmf[[type]]$beta_denovo = sigs
  x$nmf[[type]]$pyro$beta_denovo = sigs
  x$nmf[[type]]$pyro$params$infered_params$beta_d = sigs
  return(x)
}

set_fixed_signatures = function(x, sigs, type) {
  x$nmf[[type]]$beta_fixed = sigs
  x$nmf[[type]]$pyro$beta_fixed = sigs
  x$nmf[[type]]$pyro$params$infered_params$beta_f = long_to_wide(sigs, what="beta")
  return(x)
}

set_exposures = function(x, expos, type) {
  x$nmf[[type]]$exposure = expos
  x$nmf[[type]]$pyro$exposure = expos
  x$nmf[[type]]$pyro$params$infered_params$alpha = long_to_wide(expos, what="expos")
  return(x)
}

set_nmf_init_params = function(x, type, denovo=NULL, expos=NULL) {
  if (!is.null(denovo)) x$nmf[[type]]$pyro$params$init_params$beta_dn_param = denovo
  if (!is.null(expos)) x$nmf[[type]]$pyro$params$init_params$alpha = expos
  return(x)
}




