

fit <- function(
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

  data$m <- x$x           # data.frame
  data$n_samples <- x$n_samples   # integer
  data$beta_fixed <- x$beta_fixed  # data.frame
  data$k_fixed <- x$k_fixed     # integer
  data$groups <- x$groups      # for now NULL

  data$alpha <- x$alpha       # data.frame
  data$beta_denovo <- x$beta_denovo # data.frame
  data$k_denovo <- x$k_denovo    # integer

  data$n_steps <- x$n_steps     # integer
  data$lr <- x$lr          # numeric
  data$losses <- x$losses      # numeric
  data$bic <- x$bic         # numeric

  #x$model
  #x$guide

  return(data)
}

#-------------------------------------------------------------------------------

tryCatch(
  {
    1 + 1
    print("Everything was fine.")
  },

  errore = function(e) {
    print("There was an error message.")
  },

  warning = function(w) {
    print("There was a warning message.")
  },

  finally = {
    print("finally Executed")
  }
)


filter_fixed <- function(M, alpha, beta_fixed, phi) {

  if (is.null(beta_fixed)) {
    return(NULL)
  }

  if (is.data.frame(M)) {
    warning("invalid M!")
  }
  if (is.data.frame(alpha)) {
    warning("invalid alpha!")
  }
  if (is.numeric(phi)) {
    warning("invalid phi!")
  }

  theta <- matrix(rowSums(M), nrow = 1)
  total_mut <- sum(theta)

  alpha <- theta %*% as.matrix(alpha[rownames(beta_fixed)])
  sig_cont <- colSums(alpha)
  selected <- which(sig_cont > phi * total_mut)
  beta_fixed <- beta_fixed[selected, ]

  return(beta_fixed)
}

#-------------------------------------------------------------------------------

filter_denovo <- function(beta_denovo, reference_catalogue, delta) {
  if (is.data.frame(reference_catalogue) & is.numeric(delta)) {

    if (is.null(beta_denovo)) {
      return(NULL)
    }

    df <- data.frame(matrix(0, nrow(beta_denovo), nrow(reference_catalogue)))
    rownames(df) <- rownames(beta_denovo)
    colnames(df) <- rownames(reference_catalogue)

    for (i in 1:nrow(beta_denovo)) {
      denovo <- beta_denovo[i, ]
      for (j in 1:nrow(reference_catalogue)) {
        ref <- reference_catalogue[j, ]
        score <- cosine_sim(denovo, ref)
        df[i,j] <- score
      }
    }
    match_list <- c()
    while (TRUE) {
      max = which(df == max(df), arr.ind = TRUE)
      if (df[max] < delta) {
        break
      }
      row_ind <- as.numeric(max)[1]
      col_ind <- as.numeric(max)[2]

      match_list[length(match_list) + 1] <- colnames(df[col_ind])

      df <- df[-c(row_ind), -c(col_ind)]
      if (!is.data.frame(beta_denovo)) {
        break
      }
    }

    if (length(match_list) == 0) {
      return(NULL)
    } else {
      return(reference_catalogue[match_list, ])
    }
  }
  else {
    warning('invalid beta_fixed, reference_catalogue or delta argument')
  }
}
