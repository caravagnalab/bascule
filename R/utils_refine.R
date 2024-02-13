
# input: basilica object
# output: input basilica object + refined denovo, fixed and exposure
# comparing similarity between all possible pairs of denovo-fixed and denovo-denovo signatures
# if there is similarity higher than threshold then remove the denovo signature
# and adding its exposure contribution to the corresponding pair signature
refine_signatures_aux <- function(x, type, threshold) {
  
  #types <- basilica:::get_types(x)
  exposure <- basilica:::get_exposure(x, matrix = TRUE)[[type]]
  fixed <- basilica:::get_fixed_signatures(x, matrix = TRUE)[[type]]
  denovo <- basilica:::get_denovo_signatures(x, matrix = TRUE)[[type]]
  
  a <- basilica:::cosine.matrix(fixed, denovo)
  while (max(a) >= threshold) {
    row <- which(a == max(a), arr.ind = T)[1]
    col <- which(a == max(a), arr.ind = T)[2]
    sigF <- rownames(a[row, ])
    sigD <- colnames(a[col])
    exposure[, sigF] <- exposure[, sigF] + exposure[, sigD]
    exposure[, sigD] <- NULL # remove the denovo signature from exposure matrix
    denovo <- denovo[!(rownames(denovo) %in% sigD), ] # remove the denovo signature from denovo matrix
    a[, col] <- 0 # make zero the denovo signature column in cosine matrix
  }
  
  b <- basilica:::cosine.matrix(denovo, denovo)
  for (i in 1:nrow(b)) {
    for (j in 1:i) {
      b[i,j] <- 0
    }
  }
  while (max(b) >= threshold) {
    row <- which(b == max(b), arr.ind = T)[1]
    col <- which(b == max(b), arr.ind = T)[2]
    sig1 <- rownames(b[row, ])
    sig2 <- colnames(b[col])
    exposure[, sig1] <- exposure[, sig1] + exposure[, sig2]
    exposure[, sig2] <- NULL # remove the denovo signature from exposure matrix
    denovo <- denovo[!(rownames(denovo) %in% sig2), ] # remove the denovo signature from denovo matrix
    b[, sig2] <- 0 # make zero column includes the removed denovo signature in cosine matrix
    b[sig2, ] <- 0 # make zero row includes the removed denovo signature in cosine matrix
  }
  
  return(
    list(
      exposure=basilica:::wide_to_long(exposure, what = "exposures"), 
      fixed=basilica:::wide_to_long(fixed, what = "beta"), 
      denovo=basilica:::wide_to_long(denovo, what = "beta")
    )
  )
}


refine_signatures <- function(x, threshold=0.8) {
  xxx <- lapply(basilica:::get_types(x), refine_signatures_aux, x=x, threshold=threshold)
  names(xxx) <- paste(basilica:::get_types(x), "refined", sep = "_")
  x$nmf$refined <- xx
  return(x)
}


