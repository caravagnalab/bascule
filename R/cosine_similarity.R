
#-------------------------------------------------------------------------------

cosine.vector <- function(a, b, substitutions = NULL) {

  if (is.matrix(a) & nrow(a)>1) a = t(a)
  if (is.matrix(b) & nrow(b)>1) b = t(b)

  if (is.null(substitutions)) {
    if (!identical(colnames(a), colnames(b))) {
      a = a[names(b)]
    }

    numerator <- sum(a * b)
    denominator <- sqrt(sum(a^2)) * sqrt(sum(b^2))
    return(numerator / denominator)
  }

  cosine.tot = 0
  for (ss in substitutions) {
    keep_cols.tmp = grep(ss, colnames(b), value = T)
    num.tmp = sum(a[,keep_cols.tmp] * b[,keep_cols.tmp])
    denomin.tmp = sqrt(sum(a[,keep_cols.tmp]^2)) * sqrt(sum(b[,keep_cols.tmp]^2))
    cosine.tot = cosine.tot + (num.tmp / denomin.tmp)
  }

  return(cosine.tot / length(substitutions))
}

#-------------------------------------------------------------------------------
# if by_context is TRUE, it computes the cosine similarity by substitution type
cosine.matrix <- function(a, b, substitutions=NULL) {
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

      score <- cosine.vector(denovo, ref, substitutions)
      df[i,j] <- score
    }
  }

  print(df)


  return(df)
}

#-------------------------------------------------------------------------------




get_contexts = function(x) {
  tryCatch(expr = {
    context_names = x %>% get_signatures(long=T) %>% dplyr::select(Feature) %>% unique()
  }, error = function(e)
    context_names = data.frame(Feature = x$denovo_signatures %>% colnames()) )

  return(context_names %>%
           dplyr::mutate(Feature=gsub(pattern = '\\[', replacement = '_', x=Feature)) %>%
           dplyr::mutate(Feature=gsub(pattern = '\\]', replacement = '_', x=Feature)) %>%
           tidyr::separate(Feature, into=c("left","subs","right"), sep="_"))
}
