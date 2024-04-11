
linux_path <- "/home/azad/Documents/nextcloud/"
mac_path <- "/Users/azadsadr/Nextcloud/"
path <- mac_path

source(paste(path, "basilica/scripts/nmf-analysis/utils_nmf/utils_qc.R", sep = ""))
source(paste(path, "basilica/scripts/nmf-analysis/utils_nmf/utils_delete.R", sep = ""))


refinement <- function(x) {
  #discarded.signatures <- c()
  while (TRUE) {
    fixed <- basilica:::get_fixed_signatures(x, types = "SBS", matrix = TRUE)[[1]]
    denovo <- basilica:::get_denovo_signatures(x, types = "SBS", matrix = TRUE)[[1]]
    exposure <- basilica:::get_exposure(x, types = "SBS", matrix = TRUE)[[1]]
    df <- qc.linearCombination(fixed = fixed, denovo = denovo, matrix = FALSE)
    a <- df[!duplicated(df$denovos), ] %>% select(c(denovos, scores))
    
    if ( (length(a$scores) > 0) & (max(a$scores) > 0) ) {
      candidate <- a[which.max(a$scores), ]$denovos
      #discarded.signatures[length(discarded.signatures)+1] <- candidate
      print(paste0("signature ", candidate, " discarded!"))
      coefs <- subset(df[df$denovos == candidate, ])$coef
      b <- delete.signature_aux(
        fixed=fixed, 
        denovo=denovo, 
        exposure=exposure, 
        coefs=coefs, 
        sigName=candidate
      )
      x$nmf$SBS$beta_fixed <- basilica:::wide_to_long(b$fixed, what = "beta")
      x$nmf$SBS$beta_denovo <- basilica:::wide_to_long(b$denovo, what = "beta")
      x$nmf$SBS$exposure <- basilica:::wide_to_long(b$exposure, what = "exposures")
    } else {
      return(x)
    }
  }
}

