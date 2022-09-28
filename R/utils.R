
#------------------------------------------------------------------ [QC: PASSED]

pyfit <- function(
    x,
    k_list,
    lr,
    n_steps,
    groups=NULL,
    input_catalogue=NULL,
    lambda_rate=NULL,
    sigma=FALSE
) {

  py <- reticulate::import("pybasilica")
  k_list <- reticulate::r_to_py(list(as.integer(k_list)))
  obj <- py$fit(x=x, k_list=k_list, lr=lr, n_steps=n_steps, groups=groups, beta_fixed=input_catalogue)
                # lambda_rate=lambda_rate,
                # sigma=sigma)


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
# M ------------> data.frame
# alpha --------> data.frame
# beta_fixed ---> data.frame / NULL
# phi ----------> numeric
filter.fixed <- function(M, alpha, beta_fixed=NULL, phi=0.05) {

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
    #col_names <- colnames(M)
    #df = data.frame(matrix(nrow=0, ncol = length(col_names)))
    #colnames(df) = col_names
    remained_fixed <- NULL
    dropped_fixed <- NULL
  } else if (is.data.frame(beta_fixed)) {
    theta <- matrix(rowSums(M), nrow = 1)
    #print(theta)
    alpha0 <- theta %*% as.matrix(alpha)
    #print(alpha0)
    contribution <- colSums(alpha0) / sum(alpha0)
    #print(contribution)
    dropped <- which(contribution < phi)
    # TEST ------------------
    #print('======')
    #print('dropped')
    #print(dropped)
    #print('======')
    # TEST ------------------
    #print(dropped)
    if (sum(dropped)==0) {
      remained_fixed <- beta_fixed
      dropped_fixed <- NULL
      #print('nothing to drop')
    } else {
      remained_fixed <- beta_fixed[-c(dropped), ]
      if (nrow(remained_fixed)==0) {
        remained_fixed <- NULL
      }
      dropped_fixed <- beta_fixed[c(dropped), ]
    }
  } else {
    warning("invalid fixed signatures (beta_fixed) !")
  }
  return(list(remained_fixed=remained_fixed, dropped_fixed=dropped_fixed))
  # remained_fixed ----> data.frame / NULL
  # dropped_fixed -----> data.frame / NULL
}

#-------------------------------------------------------------------------------

#' @import dplyr
filter.denovo <- function(reference_catalogue, beta_fixed, beta_denovo=NULL, black_list=NULL, delta=0.9) {

  if (!is.data.frame(reference_catalogue)) {
    warning("Invalid reference catalogue!")
  }
  if (!is.numeric(delta)) {
    warning("Invalid delta argument!")
  }

  # (Reference - Beta Fixed) ----------------------
  if (is.data.frame(beta_fixed)) {
    reference <- dplyr::setdiff(reference_catalogue, beta_fixed)
  } else if (is.null(beta_fixed)) {
    reference <- reference_catalogue
  } else {
    warning('invalid fixed signatures (beta_fixed) !')
  }

  # BETA DENOVO ---------------------------------
  if (is.null(beta_denovo)) {
    return(list(new_fixed=NULL, reduced_denovo=NULL))
  } else if (is.data.frame(beta_denovo)) {
    match_list <- c()
    cos_matrix <- cosine.matrix(beta_denovo, reference)
    while (TRUE) {
      max = which(cos_matrix == max(cos_matrix), arr.ind = TRUE)
      if (cos_matrix[max] < delta) {
        break
      }
      row_index <- as.numeric(max)[1]
      col_index <- as.numeric(max)[2]
      match_list[length(match_list) + 1] <- colnames(cos_matrix[col_index])

      if (dim(cos_matrix)[1]==1 | dim(cos_matrix)[2]==1) {
        cos_matrix <- cos_matrix[-c(row_index), -c(col_index)]
        break
      } else {
        cos_matrix <- cos_matrix[-c(row_index), -c(col_index)]
      }
    }
  } else {
    warning("Invalid beta denovo!")
  }

  match_list <- setdiff(match_list, black_list)
  if (length(match_list) == 0) {
    #col_names <- colnames(reference_catalogue)
    #df = data.frame(matrix(nrow=0, ncol = length(col_names)))
    #colnames(df) = col_names
    return(list(new_fixed=NULL, reduced_denovo=beta_denovo))
  } else {
    if (is.null(dim(cos_matrix))) {
      return(list(new_fixed=reference[match_list, ], reduced_denovo=NULL))
    } else {
      reduced_denovo <- beta_denovo[rownames(cos_matrix), ]
      return(list(new_fixed=reference[match_list, ], reduced_denovo=reduced_denovo))
    }
  }
}

#-------------------------------------------------------------------------------

adjust.denovo.fixed <- function(alpha, fixed_signatures, denovo_signatures, limit=0.9) {

  if (is.null(fixed_signatures) | is.null(denovo_signatures)) {
    return(list(exposure=alpha, denovo_signatures=denovo_signatures))
  }

  cos_matrix <- basilica:::cosine.matrix(denovo_signatures, fixed_signatures)
  while (TRUE) {
    max = which(cos_matrix == max(cos_matrix), arr.ind = TRUE)
    if (cos_matrix[max] < limit) {
      break
    } else {
      row_index <- as.numeric(max)[1]
      col_index <- as.numeric(max)[2]
      denovo_name <- rownames(cos_matrix[col_index])[row_index]
      fixed_name <- colnames(cos_matrix[col_index])

      # remove denovo signature
      denovo_signatures <- denovo_signatures[!(rownames(denovo_signatures) %in% denovo_name), ]

      # adjust fixed signature exposure
      alpha[fixed_name] <- alpha[, fixed_name] + alpha[, denovo_name]
      # remove denovo signature exposure
      alpha <- alpha[ , !names(alpha) %in% denovo_name]

      if (dim(cos_matrix)[1]==1 | dim(cos_matrix)[2]==1) {
        break
      } else {
        cos_matrix <- cos_matrix[-c(row_index), -c(col_index)]
      }
    }
  }
  return(list(exposure=alpha, denovo_signatures=denovo_signatures))
}

#-------------------------------------------------------------------------------

adjust.denovo.denovo <- function(alpha, denovo_signatures, limit=0.9) {

  if (is.null(denovo_signatures)) {
    return(list(exposure=alpha, denovo_signatures=denovo_signatures))
  } else if (nrow(denovo_signatures)==1) {
    return(list(exposure=alpha, denovo_signatures=denovo_signatures))
  }
  cos_matrix <- basilica:::cosine.matrix(denovo_signatures, denovo_signatures)
  for (i in 1:nrow(cos_matrix)) {
    cos_matrix[i, i] <- 0
  }
  while (TRUE) {
    max = which(cos_matrix == max(cos_matrix), arr.ind = TRUE)
    if (cos_matrix[max][1] < limit) {
      break
    } else {
      row_index <- as.numeric(max)[1]
      col_index <- as.numeric(max)[2]
      denovo_one <- rownames(cos_matrix[col_index])[row_index]
      denovo_two <- colnames(cos_matrix[col_index])

      denovo_signatures[paste(denovo_one, denovo_two, sep = ''), ] <- denovo_signatures[denovo_one, ] + denovo_signatures[denovo_two, ]

      denovo_signatures <- denovo_signatures[!(rownames(denovo_signatures) %in% c(denovo_one, denovo_two)), ]

      # adjust signature in exposure
      alpha[paste(denovo_one, denovo_two, sep = '')] <- alpha[, denovo_one] + alpha[, denovo_two]
      # remove signature from exposure
      alpha <- alpha[ , !names(alpha) %in% c(denovo_one, denovo_two)]

      if (dim(cos_matrix)[1]==2 | dim(cos_matrix)[2]==2) {
        break
      } else {
        cos_matrix <- cos_matrix[-c(row_index, col_index), -c(col_index, row_index)]
      }
    }
  }
  return(list(exposure=alpha, denovo_signatures=denovo_signatures))
}



