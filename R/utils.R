

#------------------------------------------------------------------ [QC: PASSED]

pyfit <- function(x,
                  k_list,
                  lr,
                  n_steps,
                  py = NULL,
                  groups = NULL,
                  input_catalogue = NULL,
                  lambda_rate = NULL,
                  sigma = FALSE,
                  CUDA = FALSE,
                  compile = FALSE,
                  enforce_sparsity = FALSE,
                  store_parameters = FALSE,
                  regularizer = "cosine",
                  reg_weight = 1,
                  reg_bic = FALSE) {
  if (is.null(py))
    py <- reticulate::import("pybasilica")

  if (length(k_list) > 1)
    k_list <- reticulate::r_to_py(as.integer(k_list))
  else
    k_list <- reticulate::r_to_py(list(as.integer(k_list)))


  obj <-
    py$fit(
      x = x,
      k_list = k_list,
      lr = lr,
      n_steps = n_steps,
      groups = groups,
      beta_fixed = input_catalogue,
      CUDA = CUDA,
      compile_model = compile,
      enforce_sparsity = enforce_sparsity,
      store_parameters = store_parameters,
      regularizer = regularizer,
      reg_weight = reg_weight,
      reg_bic = reg_bic
    )
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
  data$denovo_signatures = obj$beta_denovo
  data$bic = obj$bic
  data$losses = obj$losses
  data$train_params = get_train_params(obj)
  data$groups = obj$groups

  # output
  # ---------------------------------:
  # exposure          --> data.frame
  # denovo_signatures --> data.frame
  # bic               --> numeric
  # losses            --> numeric

  return(data)
}


get_train_params = function(obj) {
  if (!obj$store_parameters)
    return(NULL)
  train_params = obj$train_params
  samples_names = obj$alpha %>% rownames()
  bfixed_names = obj$beta_fixed %>% rownames()
  bdenovo_names = obj$beta_denovo %>% rownames()
  contexts = obj$beta_denovo %>% colnames()

  alpha_all = data.frame() %>% dplyr::mutate(sample_id=as.character(NA), signature=as.character(NA),
                                             iteration=as.integer(NA), alpha=as.numeric(NA))
  beta_d = data.frame() %>% dplyr::mutate(signature=as.character(NA), context=as.character(NA),
                                          iteration=as.integer(NA), beta=as.numeric(NA))

  for (i in 1:length(train_params)) {
    tmp_a = train_params[[i]][["alpha"]]$numpy() %>% as.data.frame()
    rownames(tmp_a) = samples_names
    colnames(tmp_a) = c(bfixed_names, bdenovo_names)

    tmp_a = tmp_a %>% tibble::rownames_to_column(var="sample_id") %>%
      reshape2::melt(id="sample_id",variable.name="signature",value.name="alpha") %>%
      dplyr::mutate(iteration=i)

    alpha_all = alpha_all %>% dplyr::add_row(tmp_a)

    tmp_b = train_params[[i]][["beta_d"]]$numpy() %>% as.data.frame()
    rownames(tmp_b) = bdenovo_names
    colnames(tmp_b) = contexts

    tmp_b = tmp_b %>% tibble::rownames_to_column(var="signature") %>%
      reshape2::melt(id="signature",variable.name="context",value.name="beta") %>%
      dplyr::mutate(iteration=i)

    beta_d = beta_d %>% dplyr::add_row(tmp_b)
  }

  return(tibble::tibble(alpha=list(alpha_all), beta_d=list(beta_d)))
}

# Filter fixed signatures ------------------------------------------------------
# M ------------> data.frame
# alpha --------> data.frame
# beta_fixed ---> data.frame / NULL
# phi ----------> numeric
filter.fixed <- function(M,
                         alpha,
                         beta_fixed = NULL,
                         phi = 0.05) {
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

    atb = alpha0 %>%
      as_tibble() %>%
      reshape2::melt(id = NULL) %>%
      dplyr::rename(Signature = variable, TMB = value)

    ctb = contribution %>%
      as.list() %>%
      as_tibble %>%
      reshape2::melt(id = NULL) %>%
      dplyr::rename(Signature = variable, proportion = value)

    predicate_collapsed = dplyr::full_join(atb, ctb, by = 'Signature') %>%
      dplyr::arrange(dplyr::desc(proportion)) %>%
      mutate(Signature = ifelse(proportion < phi, crayon::red(Signature), Signature))

    predicate_collapsed$TMB = paste0("TMB = ",
                                     predicate_collapsed$TMB %>% round(0))

    predicate_collapsed$proportion = paste0("\u03c0 = ",
                                            predicate_collapsed$proportion %>% round(3))

    predicate_collapsed$Signature = sprintf("%20s", predicate_collapsed$Signature)
    predicate_collapsed$TMB = sprintf("%20s", predicate_collapsed$TMB)
    predicate_collapsed$proportion = sprintf("%20s", predicate_collapsed$proportion)

    predicate_collapsed = apply(predicate_collapsed, 1, function(x)
      paste(x, collapse = ' '))

    cli::boxx(
      predicate_collapsed,
      header = "TMB filter",
      float = 'center',
      footer = paste0("\u03c0 > ", phi)
    ) %>% cat()
    cat('\n')

    #print(contribution)
    dropped <- which(contribution < phi)
    # TEST ------------------
    #print('======')
    #print('dropped')
    #print(dropped)
    #print('======')
    # TEST ------------------
    #print(dropped)
    if (sum(dropped) == 0) {
      remained_fixed <- beta_fixed
      dropped_fixed <- NULL
      #print('nothing to drop')
    } else {
      remained_fixed <- beta_fixed[-c(dropped),]
      if (nrow(remained_fixed) == 0) {
        remained_fixed <- NULL
      }

      if (any(dropped > nrow(beta_fixed))) {
        cli::boxx("AZAD this is a bug") %>% cat()
        dropped = dropped[dropped < nrow(beta_fixed)]
      }

      dropped_fixed <- beta_fixed[c(dropped),]
    }
  } else {
    warning("invalid fixed signatures (beta_fixed) !")
  }
  return(list(remained_fixed = remained_fixed, dropped_fixed = dropped_fixed))
  # remained_fixed ----> data.frame / NULL
  # dropped_fixed -----> data.frame / NULL
}

# Checks if a signature has at least one patient where it's exposure exceeds phi
filter.fixed_minfreq <- function(alpha, beta_fixed, phi = 0.15)
{
  if (is.null(beta_fixed))
    return(list(remained_fixed = NULL, dropped_fixed = NULL))

  if (!is.null(alpha)) {
    alpha_cat = alpha[, rownames(beta_fixed)]
    predicate = apply(alpha_cat, 2, function(x)
      sum(x > phi))

    stays = colnames(alpha_cat)[predicate > 0]
    goes = colnames(alpha_cat)[predicate == 0]

    if (length(stays) == 0)
      remained_fixed = NULL
    else
      remained_fixed = beta_fixed[stays,]

    if (length(goes) == 0)
      dropped_fixed = NULL
    else
      dropped_fixed = beta_fixed[goes,]

    predicate_collapsed = predicate %>% as_tibble()
    predicate_collapsed$Signature = names(predicate)

    predicate_collapsed = predicate_collapsed %>%
      group_by(value) %>%
      mutate(Signature = paste(Signature, collapse = ', ')) %>%
      distinct() %>%
      arrange(value) %>%
      mutate(Signature = ifelse(value == 0, crayon::red(Signature), Signature))

    predicate_collapsed = paste('n =',
                                predicate_collapsed$value,
                                '[',
                                predicate_collapsed$Signature,
                                ']')

    cli::boxx(
      predicate_collapsed,
      header = "Frequency filter",
      float = 'center',
      footer = paste0("\u03C6 > ", phi, " in n samples")
    )  %>% cat()

    cat('\n')

    return(list(remained_fixed = remained_fixed, dropped_fixed = dropped_fixed))
  }
}

# Fake filter implementation
filter.fixed_nofilter <- function(alpha, beta_fixed)
{
  if (is.null(beta_fixed))
    return(list(remained_fixed = NULL, dropped_fixed = NULL))


  if (!is.null(beta_fixed))
    return(list(remained_fixed = beta_fixed, dropped_fixed = NULL))

  if (!is.null(alpha))
    return(list(remained_fixed = beta_fixed[colnames(alpha)[colnames(alpha) %in% rownames(beta_fixed)],],dropped_fixed = NULL))


}


# Filter denovo signatures -----------------------------------------------------

#' @import dplyr
filter.denovo <-
  function(reference_catalogue,
           beta_fixed,
           beta_denovo = NULL,
           black_list = NULL,
           delta = 0.9) {
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
      return(list(new_fixed = NULL, reduced_denovo = NULL))
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
        match_list[length(match_list) + 1] <-
          colnames(cos_matrix[col_index])

        if (dim(cos_matrix)[1] == 1 | dim(cos_matrix)[2] == 1) {
          cos_matrix <- cos_matrix[-c(row_index),-c(col_index)]
          break
        } else {
          cos_matrix <- cos_matrix[-c(row_index),-c(col_index)]
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
      return(list(new_fixed = NULL, reduced_denovo = beta_denovo))
    } else {
      if (is.null(dim(cos_matrix))) {
        return(list(new_fixed = reference[match_list,], reduced_denovo = NULL))
      } else {
        reduced_denovo <- beta_denovo[rownames(cos_matrix),]
        return(list(new_fixed = reference[match_list,], reduced_denovo = reduced_denovo))
      }
    }
  }

#-------------------------------------------------------------------------------

adjust.denovo.fixed <-
  function(alpha,
           fixed_signatures,
           denovo_signatures,
           limit = 0.9) {
    if (is.null(fixed_signatures) | is.null(denovo_signatures)) {
      return(list(exposure = alpha, denovo_signatures = denovo_signatures))
    }

    cos_matrix <-
      basilica:::cosine.matrix(denovo_signatures, fixed_signatures)
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
        denovo_signatures <-
          denovo_signatures[!(rownames(denovo_signatures) %in% denovo_name),]

        # adjust fixed signature exposure
        alpha[fixed_name] <-
          alpha[, fixed_name] + alpha[, denovo_name]
        # remove denovo signature exposure
        alpha <- alpha[,!names(alpha) %in% denovo_name]

        if (dim(cos_matrix)[1] == 1 | dim(cos_matrix)[2] == 1) {
          break
        } else {
          cos_matrix <- cos_matrix[-c(row_index),-c(col_index)]
        }
      }
    }
    return(list(exposure = alpha, denovo_signatures = denovo_signatures))
  }

#-------------------------------------------------------------------------------

adjust.denovo.denovo <-
  function(alpha, denovo_signatures, limit = 0.9) {
    if (is.null(denovo_signatures)) {
      return(list(exposure = alpha, denovo_signatures = denovo_signatures))
    } else if (nrow(denovo_signatures) == 1) {
      return(list(exposure = alpha, denovo_signatures = denovo_signatures))
    }
    cos_matrix <-
      basilica:::cosine.matrix(denovo_signatures, denovo_signatures)
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

        denovo_signatures[paste(denovo_one, denovo_two, sep = ''),] <-
          denovo_signatures[denovo_one,] + denovo_signatures[denovo_two,]

        denovo_signatures <-
          denovo_signatures[!(rownames(denovo_signatures) %in% c(denovo_one, denovo_two)),]

        # adjust signature in exposure
        alpha[paste(denovo_one, denovo_two, sep = '')] <-
          alpha[, denovo_one] + alpha[, denovo_two]
        # remove signature from exposure
        alpha <-
          alpha[,!names(alpha) %in% c(denovo_one, denovo_two)]

        if (dim(cos_matrix)[1] == 2 | dim(cos_matrix)[2] == 2) {
          break
        } else {
          cos_matrix <-
            cos_matrix[-c(row_index, col_index),-c(col_index, row_index)]
        }
      }
    }
    return(list(exposure = alpha, denovo_signatures = denovo_signatures))
  }


# filter based on linear projection with constraints

#' @import dplyr
filter.denovo.QP <-
  function(reference,
           # beta_fixed,
           beta_denovo = NULL,
           black_list = NULL,
           delta = 0.9,
           filt_pi = 0.05,
           thr_exposure = 0.05,
           exposures = NULL,
           substitutions = NULL) {
    ## if denovo = TRUE -> check also if the denovo have exposure > thr in same samples

    if (!is.data.frame(reference)) warning("Invalid reference catalogue!")
    if (!is.numeric(delta)) warning("Invalid delta argument!")

    # (Reference - Beta Fixed) -------------------------------------------------
    # if (is.data.frame(beta_fixed)) {
    #   reference <- dplyr::setdiff(reference_catalogue, beta_fixed)
    # } else if (is.null(beta_fixed)) {
    #   reference <- reference_catalogue
    # } else {
    #   warning('invalid fixed signatures (beta_fixed) !')
    # }

    if (nrow(reference)==0) return(list(new_fixed = NULL, reduced_denovo = beta_denovo))

    # BETA DENOVO --------------------------------------------------------------
    if (is.null(beta_denovo)) return(list(new_fixed = NULL, reduced_denovo = NULL))

    if (is.data.frame(beta_denovo)) {
      if (is.null(exposures)) {
        a = beta_denovo
        b = reference
      } else {
        a = reference
        b = beta_denovo
      }
      ### names of catalogue signatures to include + names de novo to remove
      res_optimization <-
        solve.quadratic.optimization(a,
                                     b,
                                     delta = delta,
                                     filt_pi = filt_pi,
                                     thr_exposure = thr_exposure,
                                     exposures = exposures,
                                     substitutions = substitutions)
      match_list <- res_optimization$catalogue_to_include
      } else warning("Invalid beta denovo!")

    match_list <- setdiff(match_list, black_list)
    if (length(match_list) == 0) {
      # col_names <- colnames(reference_catalogue)
      # df = data.frame(matrix(nrow=0, ncol = length(col_names)))
      # colnames(df) = col_names
      return(list(new_fixed = NULL, reduced_denovo = beta_denovo))
    } else {
      if (is.null(res_optimization$denovo_to_include)) {
        return(list(new_fixed = reference[match_list,], reduced_denovo = NULL))
      } else {
        reduced_denovo <- beta_denovo[res_optimization$denovo_to_include,]
        return(list(new_fixed = reference[match_list,], reduced_denovo = reduced_denovo))
      }
    }
  }



solve.quadratic.optimization <-
  function(a,
           b,
           filt_pi = 0.05,
           delta = 0.9,
           exposures = NULL,
           thr_exposure = 0.05,
           substitutions = NULL) {
    # a and b are data.frame

    df <- data.frame(matrix(0, nrow(a), nrow(b)))
    rownames(df) <- rownames(a)
    colnames(df) <- rownames(b)

    cmp = nrow(a)
    pb <- progress::progress_bar$new(
      format = paste0("  Quadratic programming solver (n = ", cmp, ") [:bar] :percent eta: :eta"),
      total = cmp,
      clear = FALSE,
      width = 90)

    b_m <- as.matrix(b) %>% t
    Rinv <- solve(chol(t(b_m) %*% b_m))

    res <- lapply(
      1:nrow(a),
      FUN = function(i) {
        optim_res <-
          solve.quadratic.optimization.aux(
            v = as.matrix(a)[i, ] %>% t(),
            Z = b_m,
            Rinv = Rinv,
            filt_pi = filt_pi,
            delta = delta,
            a_name = rownames(a)[i],
            exposures = exposures,
            thr_exposure = thr_exposure,
            substitutions = substitutions)
        pb$tick()
        return(optim_res)
      }
    )

    catalogue_to_include <-
      lapply(res, function(x)
        x[[1]]) %>% do.call(c, .) %>% unique()
    denovo_to_include <-
      lapply(res, function(x)
        x[[2]]) %>% do.call(c, .)

    if (!is.null(exposures))
      denovo_to_include = setdiff(rownames(b), denovo_to_include)

    return(
      list(
        catalogue_to_include = catalogue_to_include,
        denovo_to_include = denovo_to_include
      )
    )
  }


solve.quadratic.optimization.aux <-
  function(v,
           Z,
           Rinv,
           a_name,
           filt_pi = 0.05,
           delta = 0.9,
           exposures = NULL,
           thr_exposure = 0.05,
           substitutions = NULL) {
    d <- v %*% Z
    b <- c(1, rep(0, length(d)))
    C <- cbind(rep(1, length(d)), diag(length(d)))
    pis <-
      quadprog::solve.QP(
        Dmat = Rinv,
        factorized = TRUE,
        dvec = d,
        Amat = C,
        bvec = b,
        meq = 1
      )$solution
    pis[pis < filt_pi] <- 0
    reconstructed_vector <- Z %*% pis

    vec1 = matrix(v)
    rownames(vec1) = rownames(reconstructed_vector)

    cos_sim <- cosine.vector(vec1, reconstructed_vector, substitutions = substitutions)

    # return -> first element is catalogue to include
    #           second element is denovo to include

    if (is.null(exposures)) {
      if (!is.na(cos_sim) && cos_sim > delta)
        return(list(colnames(Z)[pis > 0], NULL))
      else
        return(list(NULL, a_name))
    } else {
      if (!is.na(cos_sim) && cos_sim > delta) {
        if (sum(pis>0) <= 1) return(list(NULL, colnames(Z)[pis > 0]))

        keep.tmp = colnames(Z)[pis > 0]
        exposures.tmp = exposures[,keep.tmp]
        exposures.tmp$n_denovo = apply(exposures.tmp > thr_exposure, 1,
                                       function(x) length(unique(x))==1)

        print(a_name)
        print(colnames(Z)[pis > 0])
        print(sum(exposures.tmp$n_denovo))

        if (sum(exposures.tmp$n_denovo) >= nrow(exposures)*0.9)
          # return the reference that can be explained as a linear comb of denovo
          return(list(a_name, colnames(Z)[pis > 0]))
        else
          return(list(NULL, NULL))
        } else {
          return(list(NULL, NULL))
        }
    }


    if (!is.na(cos_sim) && cos_sim > delta) {
      if (!is.null(exposures) && sum(pis > 0)>1) {

        print(a_name)

        keep.tmp = colnames(Z)[pis > 0]
        exposures.tmp = exposures[,keep.tmp]
        exposures.tmp$n_denovo = apply(exposures.tmp > thr_exposure, 1, function(x) length(unique(x))==1)

        ## if both either are present or not in the same patients (more than 90% of the times)
        if (sum(exposures.tmp$n_denovo) >= nrow(exposures)*0.9)
          return(list(colnames(Z)[pis>0],NULL))
      } else if (!is.null(exposures))
        return(list(NULL, NULL))

      return(list(colnames(Z)[pis > 0], NULL))

    } else {
      if (!is.null(exposures)) return(list(NULL, NULL))
      return(list(NULL, a_name))
    }

}

# Renormalize denovo -----------------------------------------------------------

renormalize_denovo_thr = function(denovo, thr=0.02) {
  denovo.tmp = denovo
  denovo.tmp[denovo.tmp < thr] = 0
  return(denovo.tmp / rowSums(denovo.tmp))
}


filter.denovo.phi = function(exposures, phi, denovo) {
  exposures.denovo = exposures[,denovo] %>% as.data.frame()
  colnames(exposures.denovo) = denovo
  sbs_low_exp = (exposures.denovo < phi) %>% colSums()

  rmv = sbs_low_exp[sbs_low_exp > nrow(exposures.denovo)] %>% names()
  if (length(rmv) == 0) return(exposures)

  return(
    exposures %>% dplyr::select(-dplyr::all_of(rmv)) %>% renormalize_denovo_thr(thr=0)
  )
}



