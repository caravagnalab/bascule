
#library(dplyr)
#library(tidyverse)
# utils_linear.R

# input :
#   fixed ---> wide
#   denovo --> wide
# output:
#   if matrix==TRUE ---> dataframe - [k_denovo X (k_denovo + k_fixed) + 1] (wide format)
#   if matrix==FALSE --> dataframe (long format)
qc.linearCombination <- function(fixed, denovo, matrix=TRUE) {
  df <- data.frame(
    matrix(
      0,
      nrow = nrow(denovo),
      ncol = nrow(denovo) + nrow(fixed)
    )
  )
  rownames(df) <- rownames(denovo)
  colnames(df) <- c(rownames(fixed), rownames(denovo))

  for (i in 1:nrow(denovo)) {
    a <- solve.quadratic.optimization(
      a = denovo[c(i), ],
      b = rbind(fixed, denovo[-c(i), ]),
      filt_pi = 0.05,
      delta = 0.9,
      thr_exposure = 0.05,
      exposures = NULL,
      return_weights = TRUE
    )
    df[rownames(denovo[c(i), ]), names(a[[1]])] <- a[[1]]
  }

  # adding reconstruction score column to dataframe
  ss <- unlist(lapply(
    rownames(df),
    function(x) computeScore_aux(fixed=fixed, denovo=denovo, coefs=as.numeric(df[x, ]), sigName=x)
  ))
  ss[is.nan(ss)] <- 0
  df$scores <- ss

  if (matrix==TRUE) {
    return(df)
  } else {
    df <- tibble::rownames_to_column(df, "denovos")
    df <- df %>% tidyr::gather(key = signature, value = "coef", -c(1, ncol(df)))
    df <- df %>%
      dplyr::mutate(df %>%
               apply(
                 1,
                 function(x) basilica:::cosine.vector(
                   denovo[x[1], ],
                   rbind(fixed, denovo)[x[3], ]
                 )
               )
      )
    colnames(df) <- c("denovos", "scores", "signature", "coef", "cosine")
    df <- df[, c("denovos", "signature", "coef", "cosine", "scores")]
    return(df)
  }
}

#-------------------------------------------------------------------------------

# input :
#   fixed ------> fixed signatures (wide datafram)
#   denovo -----> denovo signatures (wide datafram)
#   coefs ------> linear combination coefficients (numeric vector)
#   <sigName> --> signature of interest (character)
# output:
#   numeric ----> cosine similarity (true signature vs. reconstructed signature)
computeScore_aux <- function(fixed, denovo, coefs, sigName) {
  #coefs <- (df %>% subset(denovos==sigName) %>% select(coef))$coef
  reconstructed_vector <- t(rbind(fixed, denovo)) %*% coefs
  score <- basilica:::cosine.vector(reconstructed_vector, denovo[sigName,] )
  return(round(score, digits = 3))
}




