refinement = function(x, types=get_types(x)) {
  lapply(types, function(tid) x <<- refinement_aux(x, type=tid))
  return(x)
}


refinement_aux = function(x, type) {
  while (TRUE) {
    fixed = get_fixed_signatures(x, types=type, matrix=TRUE)[[type]]
    denovo = get_denovo_signatures(x, types=type, matrix=TRUE)[[type]]

    if (is.null(denovo)) return(x)
    if (nrow(denovo) == 0) return(x)

    exposure = get_exposure(x, types=type, matrix=TRUE)[[type]]

    df = qc.linearCombination(fixed=fixed, denovo=denovo, matrix=FALSE)
    a = df[!duplicated(df$denovos), ] %>% dplyr::select(c(denovos, scores))

    if ((length(a$scores) > 0) & (max(a$scores) > 0)) {
      candidate = a[which.max(a$scores), ]$denovos
      cli::cli_process_start(paste0("Deleting signature ", candidate))

      coefs = subset(df[df$denovos == candidate, ])$coef
      updated_dfs = delete.signature_aux(
        denovo=denovo,
        exposure=exposure,
        coefs=coefs,
        sigName=candidate
      )

      updated_init_dfs = delete.signature_aux(
        denovo=get_nmf_initial_parameters(x, what="nmf")[[type]]$beta_dn_param,
        exposure=get_nmf_initial_parameters(x, what="nmf")[[type]]$alpha,
        coefs=coefs,
        sigName=candidate
      )

      x = set_denovo_signatures(x, type=type,
                                sigs=wide_to_long(updated_dfs$denovo, what="beta"))
      x = set_exposures(x, type=type,
                        expos=wide_to_long(updated_dfs$exposure, what="exposures"))
      x = set_nmf_init_params(x, type=type,
                              denovo=updated_init_dfs$denovo,
                              expos=updated_init_dfs$exposure)

      cli::cli_process_done(paste0("Signature ", candidate, " discarded"))
    } else {
      return(x)
    }
  }
}



# Linear combination #####

solve.quadratic.optimization = function(a,
                                        b,
                                        filt_pi = 0.05,
                                        delta = 0.9,
                                        thr_exposure = 0.05,
                                        exposures = NULL,
                                        return_weights = FALSE) {

  df = data.frame(matrix(0, nrow(a), nrow(b)))
  rownames(df) = rownames(a)
  colnames(df) = rownames(b)

  cmp = nrow(a)
  b_m = t(as.matrix(b))
  Rinv = solve(chol(t(b_m) %*% b_m))

  res = lapply(
    1:nrow(a),
    FUN = function(i) {
      optim_res =
        solve.quadratic.optimization.aux(
          v = as.matrix(a)[i, ] %>% t(),
          Z = b_m,
          Rinv = Rinv,
          filt_pi = filt_pi,
          delta = delta,
          exposures = exposures,
          thr_exposure = thr_exposure,
          return_weights = return_weights)
      return(optim_res)
    }
  ) %>% setNames(rownames(a))

  return(res)
}


solve.quadratic.optimization.aux = function(v,
                                            Z,
                                            Rinv,
                                            filt_pi = 0.05,
                                            delta = 0.9,
                                            exposures = NULL,
                                            thr_exposure = 0.05,
                                            return_weights = FALSE) {
  d = v %*% Z
  b = c(1, rep(0, length(d)))
  C = cbind(rep(1, length(d)), diag(length(d)))
  pis =
    quadprog::solve.QP(
      Dmat = Rinv,
      factorized = TRUE,
      dvec = d,
      Amat = C,
      bvec = b,
      meq = 1
    )$solution
  pis[pis < filt_pi] = 0
  reconstructed_vector = Z %*% pis

  vec1 = matrix(v)
  rownames(vec1) = rownames(reconstructed_vector)

  cos_sim = lsa::cosine(as.numeric(v), as.numeric(reconstructed_vector)) %>%
    as.numeric()

  if (is.null(exposures)) {
    if (!is.na(cos_sim) && cos_sim > delta) {
      if (return_weights) return(pis[pis > 0] %>% setNames(colnames(Z)[pis > 0]))
      return(colnames(Z)[pis > 0])
    }

  }
  return(NULL)
}



# input :
#   fixed ---> wide
#   denovo --> wide
# output:
#   if matrix==TRUE ---> dataframe - [k_denovo X (k_denovo + k_fixed) + 1] (wide format)
#   if matrix==FALSE --> dataframe (long format)

qc.linearCombination = function(fixed, denovo, matrix=TRUE) {
  df = data.frame(matrix(0, nrow=nrow(denovo), ncol=nrow(denovo) + nrow(fixed)))
  rownames(df) = rownames(denovo)
  colnames(df) = c(rownames(fixed), rownames(denovo))

  for (i in 1:nrow(denovo)) {
    a = solve.quadratic.optimization(
      a=denovo[c(i), ],
      b=rbind(fixed, denovo[-c(i), ]),
      filt_pi=0.05,
      delta=0.9,
      thr_exposure=0.05,
      exposures=NULL,
      return_weights=TRUE
    )
    df[rownames(denovo[c(i), ]), names(a[[1]])] = a[[1]]
  }

  # adding reconstruction score column to dataframe
  ss = unlist(lapply(
    rownames(df),
    function(x) computeScore_aux(fixed=fixed, denovo=denovo, coefs=as.numeric(df[x, ]), sigName=x)
  ))
  ss[is.nan(ss)] = 0
  df$scores = ss

  if (matrix == TRUE) {
    return(df)
  } else {
    df = tibble::rownames_to_column(df, "denovos")
    df = df %>% tidyr::gather(key=signature, value="coef", -c(1, ncol(df)))
    df = df %>%
      dplyr::mutate(df %>%
                      apply(
                        1,
                        function(x) basilica:::cosine.vector(
                          denovo[x[1], ],
                          rbind(fixed, denovo)[x[3], ]
                        )
                      )
      )
    colnames(df) = c("denovos", "scores", "signature", "coef", "cosine")
    df = df[, c("denovos", "signature", "coef", "cosine", "scores")]
    return(df)
  }
}


# input :
#   fixed ------> fixed signatures (wide dataframe)
#   denovo -----> denovo signatures (wide dataframe)
#   coefs ------> linear combination coefficients (numeric vector)
#   <sigName> --> signature of interest (character)
# output:
#   numeric ----> cosine similarity (true signature vs. reconstructed signature)

computeScore_aux = function(fixed, denovo, coefs, sigName) {
  reconstructed_vector = t(rbind(fixed, denovo)) %*% coefs
  score = basilica:::cosine.vector(reconstructed_vector, denovo[sigName,])
  return(round(score, digits=3))
}



