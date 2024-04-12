refinement_aux = function(x, type) {
  while (TRUE) {
    fixed = basilica:::get_fixed_signatures(x, types=type, matrix=TRUE)[[type]]
    denovo = basilica:::get_denovo_signatures(x, types=type, matrix=TRUE)[[type]]

    if (is.null(denovo)) return(x)
    if (nrow(denovo) == 0) return(x)

    exposure = basilica:::get_exposure(x, types=type, matrix=TRUE)[[type]]

    df = qc.linearCombination(fixed=fixed, denovo=denovo, matrix=FALSE)
    a = df[!duplicated(df$denovos), ] %>% dplyr::select(c(denovos, scores))

    if ((length(a$scores) > 0) & (max(a$scores) > 0)) {
      candidate = a[which.max(a$scores), ]$denovos
      print(paste0("Signature ", candidate, " discarded!"))

      coefs = subset(df[df$denovos == candidate, ])$coef
      b = delete.signature_aux(
        fixed=fixed,
        denovo=denovo,
        exposure=exposure,
        coefs=coefs,
        sigName=candidate
      )
      x$nmf[[type]]$beta_fixed = basilica:::wide_to_long(b$fixed, what="beta")
      x$nmf[[type]]$beta_denovo = basilica:::wide_to_long(b$denovo, what="beta")
      x$nmf[[type]]$exposure = basilica:::wide_to_long(b$exposure, what="exposures")
    } else {
      return(x)
    }
  }
}


refinement = function(x, types=get_types(x)) {
  lapply(types, function(tid) x <<- refinement_aux(x, type=tid))
  return(x)
}



