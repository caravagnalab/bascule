
#-------------------------------------------------------------------------------

cosine.vector <- function(a, b) {

  numerator <- sum(a * b)
  denominator <- sqrt(sum(a^2)) * sqrt(sum(b^2))
  return(numerator / denominator)
}

#-------------------------------------------------------------------------------

cosine.matrix <- function(a, b) {
  # a and b are data.frame

  df <- data.frame(matrix(0, nrow(a), nrow(b)))
  rownames(df) <- rownames(a)
  colnames(df) <- rownames(b)

  cmp = nrow(a) * nrow(b)
  pb <- progress::progress_bar$new(
    format = paste0("  Cosine similarity (n = ", cmp, ") [:bar] :percent eta: :eta"),
    total = cmp,
    clear = FALSE,
    width= 90
    )

  for (i in 1:nrow(a)) {
    denovo <- a[i, ]
    for (j in 1:nrow(b)) {
      ref <- b[j, ]
      pb$tick()

      score <- cosine.vector(denovo, ref)
      df[i,j] <- score
    }
  }

  return(df)
}

#-------------------------------------------------------------------------------


