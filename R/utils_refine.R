

refinement_aux <- function( x,  type = c("SBS") ) {
  
  if ( !(type %in% c("SBS", "DBS")) ) {
    warning("type must be 'SBS' or 'DBS'!")
  }
  
  #discarded.signatures <- c()
  while (TRUE) {
    fixed <- basilica:::get_fixed_signatures(x, types = type, matrix = TRUE)[[1]]
    denovo <- basilica:::get_denovo_signatures(x, types = type, matrix = TRUE)[[1]]

    if (is.null(denovo)) return(x)
    if (nrow(denovo) == 0) return(x)

    exposure <- basilica:::get_exposure(x, types = type, matrix = TRUE)[[1]]
    
    df <- qc.linearCombination(fixed = fixed, denovo = denovo, matrix = FALSE)
    a <- df[!duplicated(df$denovos), ] %>% dplyr::select(c(denovos, scores))

    if ( (length(a$scores) > 0) & (max(a$scores) > 0) ) {
      candidate <- a[which.max(a$scores), ]$denovos
      # discarded.signatures[length(discarded.signatures)+1] <- candidate
      print(paste0("Signature ", candidate, " discarded!"))
      coefs <- subset(df[df$denovos == candidate, ])$coef
      b <- delete.signature_aux(
        fixed=fixed,
        denovo=denovo,
        exposure=exposure,
        coefs=coefs,
        sigName=candidate
      )
      x$nmf[[type]]$beta_fixed <- basilica:::wide_to_long(b$fixed, what = "beta")
      x$nmf[[type]]$beta_denovo <- basilica:::wide_to_long(b$denovo, what = "beta")
      x$nmf[[type]]$exposure <- basilica:::wide_to_long(b$exposure, what = "exposures")
    } else {
      return(x)
    }
  }
}

#-------------------------------------------------------------------------------

refinement <- function( x,  type = c("SBS") ) {
  
  # input quality control
  if ( !is.null(type) & is.character(type) ) {
    if ( length(type) %in% c(1,2) ) {
      if ( length(type) == 1 & (sum(type %in% c("SBS", "DBS"))==1) ) {
        #print("perfect")
      } else if ( length(type) == 2 & (sum(type %in% c("SBS", "DBS")) == 2) & (type[1] != type[2]) ) {
        #print("perfect")
      } else { warning("type should be either c('SBS'), c('DBS') or c('SBS', 'DBS')!") }
    } else { warning("must insert one or two values in type arigument!") }
  } else { warning("type is null or its type is not character!") }
  
  if ( "SBS" %in% type ) {
    print("SBS refinement started:")
    x1 <- refinement_aux( x,  type = c("SBS") )
  }
  if ( "DBS" %in% type ) {
    print("DBS refinement started:")
    x1 <- refinement_aux( x,  type = c("DBS") )
  }
  
  return(x1)
}



