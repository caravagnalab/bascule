

#-------------------------------------------------------------------------------

reconstruction_matrix <- function(m, alpha, fixed, denovo) {
  # all args are data.frame
  theta <- diag(rowSums(m)) # matrix
  alpha <- theta %*% as.matrix(alpha) # matrix
  beta <- as.matrix(rbind(fixed, denovo)) # matrix

  mr_matrix <- alpha %*% as.matrix(beta)
  mr <- round(as.data.frame(mr_matrix))
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
  return(Accuracy)
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

