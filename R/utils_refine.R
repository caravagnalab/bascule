
# input: basilica object
# output: input basilica object + refined denovo, fixed and exposure
# comparing similarity between all possible pairs of denovo-fixed and denovo-denovo signatures
# if there is similarity higher than threshold then remove the denovo signature
# and adding its exposure contribution to the corresponding pair signature
refine_signatures <- function(x, threshold=0.8) {
  
  exposure <- basilica:::long_to_wide(x$nmf$SBS$exposure, what = "exposures")
  fixed <- basilica:::long_to_wide(x$nmf$SBS$beta_fixed, what = "beta")
  denovo <- basilica:::long_to_wide(x$nmf$SBS$beta_denovo, what = "beta")
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
    b[, col] <- 0 # make zero the denovo signature column in cosine matrix
  }
  
  x$nmf$SBS_refined$exposure <- basilica:::wide_to_long(exposure, what = "exposures")
  x$nmf$SBS_refined$beta_fixed <- x$nmf$SBS$beta_fixed
  x$nmf$SBS_refined$beta_denovo <- basilica:::wide_to_long(denovo, what = "beta")
  
  return(x)
}

