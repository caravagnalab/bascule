
#------------------------------------------------------------------ [QC: PASSED]

pyfit <- function(
    x,
    k_list,
    lr,
    n_steps,
    groups=NULL,
    input_catalogue=NULL
) {

  py <- reticulate::import("pybasilica")
  obj <- py$fit(x=x, k_list=k_list, lr=lr, n_steps=n_steps, groups=groups, beta_fixed=input_catalogue)


  # save python object data in a list
  data <- list()

  #data$x <- x
  #data$groups <- groups
  #data$input_catalogue <- input_catalogue
  #data$lr <- lr
  #data$steps <- n_steps
  data$exposure <- obj$alpha
  data$denovo_signatures=obj$beta_denovo
  data$bic=obj$bic
  data$losses=obj$losses

  # output
  # ---------------------------------:
  # exposure          --> data.frame
  # denovo_signatures --> data.frame
  # bic               --> numeric
  # losses            --> numeric

  return(data)
}

#-------------------------------------------------------------------------------

filter_fixed <- function(M, alpha, beta_fixed=NULL, phi=0.05) {

  if (!is.data.frame(M)) {
    warning("invalid count matrix (M) !")
  }
  if (!is.data.frame(alpha)) {
    warning("invalid exposure matrix (alpha) !")
  }
  if (!is.numeric(phi)) {
    warning("invalid parameter phi !")
  }

  if (is.null(beta_fixed)) {
    col_names <- colnames(M)
    df = data.frame(matrix(nrow=0, ncol = length(col_names)))
    colnames(df) = col_names
  }
  else if (is.data.frame(beta_fixed)) {
    theta <- matrix(rowSums(M), nrow = 1)
    #print(theta)
    alpha0 <- theta %*% as.matrix(alpha)
    #print(alpha0)
    contribution <- colSums(alpha0) / sum(alpha0)
    #print(contribution)
    dropped <- which(contribution < phi)
    #print(dropped)
    if (sum(dropped)==0) {
      df <- beta_fixed
      #print('nothing to drop')
    } else {
      df <- beta_fixed[-c(dropped), ]
      #print(paste('dropped', length(dropped), 'signatures'))
    }
  }
  else {
    warning("invalid fixed signatures (beta_fixed) !")
  }

  return(df)
}

#-------------------------------------------------------------------------------


cosine_matrix <- function(a, b) {
  # a and b are data.frame

  df <- data.frame(matrix(0, nrow(a), nrow(b)))
  rownames(df) <- rownames(a)
  colnames(df) <- rownames(b)

  for (i in 1:nrow(a)) {
    denovo <- a[i, ]
    for (j in 1:nrow(b)) {
      ref <- b[j, ]

      score <- cosine_sim(denovo, ref)
      df[i,j] <- score
    }
  }

  return(df)

}

#-------------------------------------------------------------------------------


filter_denovo <- function(beta_denovo=NULL, reference_catalogue, delta=0.9) {

  if (!is.data.frame(reference_catalogue)) {
    warning("Invalid reference catalogue!")
  }

  if (!is.numeric(delta)) {
    warning("Invalid delta argument!")
  }

  if (is.null(beta_denovo)) {
    col_names <- colnames(reference_catalogue)
    df = data.frame(matrix(nrow=0, ncol = length(col_names)))
    colnames(df) = col_names
    return(df)
  }
  else if (is.data.frame(beta_denovo)) {
    cos_matrix <- cosine_matrix(beta_denovo, reference_catalogue)
  }

  match_list <- c()
  while (TRUE) {
    max = which(cos_matrix == max(cos_matrix), arr.ind = TRUE)
    if (cos_matrix[max] < delta) {
      break
    }
    row_index <- as.numeric(max)[1]
    col_index <- as.numeric(max)[2]
    match_list[length(match_list) + 1] <- colnames(cos_matrix[col_index])

    cos_matrix <- cos_matrix[-c(row_index), -c(col_index)]
    if (nrow(cos_matrix)==0) {
      break
    }
  }

  if (length(match_list) == 0) {
    col_names <- colnames(reference_catalogue)
    df = data.frame(matrix(nrow=0, ncol = length(col_names)))
    colnames(df) = col_names
    return(df)
  } else {
    return(reference_catalogue[match_list, ])
  }
}

