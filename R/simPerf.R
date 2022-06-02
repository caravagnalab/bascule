

#-------------------------------------------------------------------------------

reconstruction_matrix <- function(m, alpha, fixed, denovo) {
  # all args are data.frame
  theta <- diag(rowSums(m)) # matrix
  alpha <- theta %*% as.matrix(alpha) # matrix
  beta <- as.matrix(rbind(fixed, denovo)) # matrix

  mr_matrix <- alpha %*% as.matrix(beta)
  mr <- round(as.data.frame(mr_matrix))
  return(mr)
}

#-------------------------------------------------------------------------------

fixed.quality <- function(input_catalogue, expected_fixed, inferred_fixed) {
  # all args are data.frame
  inp <- rownames(input_catalogue)
  exp <- rownames(expected_fixed)
  inf <- rownames(inferred_fixed)

  TP <- length(intersect(inf, exp))
  FP <- length(setdiff(inf, exp))
  TN <- length(setdiff(setdiff(inp, exp), inf))
  FN <- length(setdiff(exp, inf))
  Accuracy <- (TP + TN) / (TP + TN + FP + FN)

  a <- list(TP = TP, FP = FP, TN = TN, FN = FN, Accuracy = Accuracy)
  return(a)
}

#-------------------------------------------------------------------------------

fitness.quality <- function(m, mr, limit) {
  counter <- 0
  for (i in 1:nrow(m)) {
    cos <- cosine_sim(m[i,], mr[i, ])
    if (cos > limit) {
      counter <- counter + 1
    }
  }
  return(counter / nrow(m))
}

#-------------------------------------------------------------------------------

compute.mae <- function(m , mr) {
  mae <- sum(abs(m - mr)) / (dim(m)[1] * dim(m)[2])
  return(mae)
}

#-------------------------------------------------------------------------------

compute.mse <- function(m , mr) {
  mse <- sum((m - mr)^2) / (dim(m)[1] * dim(m)[2])
  return(mse)
}

#-------------------------------------------------------------------------------

denovo.quality <- function(exp, inf) {

  if (length(exp)==0 | length(inf)==0) {
    return(list(matrix=0, ratio=0, cosine=0))
  } else {
    df <- data.frame(matrix(nrow = nrow(inf), ncol = nrow(exp)))
    colnames(df) <- rownames(exp)
    rownames(df) <- rownames(inf)

    for (i in 1:nrow(inf)) {
      inferred <- inf[i,]
      inferred_name <- rownames(inferred)
      maxScore <- 0
      bestMatch <- NULL
      for (j in 1:nrow(exp)) {
        target <- exp[j, ]
        target_name <- rownames(target)
        score <- cosine_sim(inferred, target)
        df[inferred_name, target_name] <- score
        '
        if (score > maxScore) {
          maxScore <- score
          bestMatch <- target_name
        }
        '
      }
    }
    '
    match_list <- colnames(df)[max.col(df, ties.method = "first")]
    s <- 0
    for (i in 1:length(match_list)) {
      s <- s + df[rownames(df)[i], match_list[i]]
    }
    cosine_score <- (s / length(match_list))
    '
  }

  a <- nrow(inf) / nrow(exp)
  b <- mean(apply(df, 1, max))
  return(list(matrix=df, ratio=a, cosine=b))
}

#-------------------------------------------------------------------------------







