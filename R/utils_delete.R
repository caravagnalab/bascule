# input
#   denovo ----> denovo signatures (wide)
#   exposure --> exposure (wide)
#   coefs -----> linear combination coefficients
#   sigName ---> denovo dignature name to delete
# output
#   list (fixed, denovo, exposure)

delete.signature_aux = function(denovo, exposure, coefs, sigName) {
  if (!(sigName %in% rownames(denovo))) {
    cli::cli_alert_warning("Wrong signature selected!")
    return(NULL)
  }

  if (all(coefs==0)) {
    cli::cli_alert_warning("Can not delete! This signature is not explained by other signatures")
    return(x)
  }

  exp = exposure %>% dplyr::select(-dplyr::all_of(sigName))
  for (refname in names(coefs)) {
    exp[,refname] = exposure[,refname] + exposure[,sigName]*coefs[refname]
  }

  # exp = exposure + sweep(exposure, MARGIN=2, coefs, `*`)
  # exp = exp[, ! names(exp) %in% c(sigName)]

  # normalizing exposure
  # norm_exp = sweep(exp, 1, rowSums(exp), "/")

  denovo = denovo[!rownames(denovo) %in% c(sigName), ]

  return(
    list(denovo=denovo,
         exposure=exp)
    )
}




