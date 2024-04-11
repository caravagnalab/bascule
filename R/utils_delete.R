

# input
#   fixed -----> fixed signatures (wide)
#   denovo ----> denovo signatures (wide)
#   exposure --> exposure (wide)
#   coefs -----> linear combination coefficients
#   sigName ---> denovo dignature name to delete
# output
#   list (fixed, denovo, exposure)
delete.signature_aux <- function(fixed, denovo, exposure, coefs, sigName) {
  
  if (!(sigName %in% rownames(denovo))) {
    warning("wrong signature selected!")
    return(NULL)
  }
  
  if (all(coefs==0)) {
    warning("can not delete! This signature is not explained by other signatures")
    return(x)
  }
  #exp <- basilica:::get_exposure(x, matrix = TRUE)[[1]]
  exp <- exposure + sweep(exposure, MARGIN = 2, coefs, `*`)
  exp <- exp[, ! names(exp) %in% c(sigName)]
  
  # normalizing exposure
  norm_exp <- sweep(exp, 1, rowSums(exp), "/")
  
  #denovo <- basilica:::get_denovo_signatures(x, matrix = TRUE)[[1]]
  denovo <- denovo[!rownames(denovo) %in% c(sigName), ]
  
  return(
    list(
      fixed = fixed, 
      denovo = denovo, 
      exposure = norm_exp
      )
    )
}




