refinement = function(x, types=get_types(x)) {
  lapply(types, function(tid) x <<- refinement_aux(x, type=tid))
  return(x)
}


refinement_aux = function(x, type) {
  if (!type %in% get_types(x)) {
    cli::cli_alert_warning("Type {type} not present. Returning the initial object.")
    return(x)
  }
  repeat {
    denovo = get_denovo_signatures(x, types=type, matrix=TRUE)[[type]]

    if (is.null(denovo)) break
    if (nrow(denovo) == 0) break

    fixed = get_fixed_signatures(x, types=type, matrix=TRUE)[[type]]
    exposure = get_exposure(x, types=type, matrix=TRUE)[[type]]

    df = qc.linearCombination(fixed=fixed, denovo=denovo, matrix=FALSE)
    a = df %>% dplyr::select(denovos, scores) %>% unique()

    if (nrow(a) == 0 | all(df$coefs==0)) {
      cli::cli_process_done(msg_done="No remaining signature is explained by others.")
      break
    }

    candidate = a[which.max(a$scores), ]$denovos

    cli::cli_process_start("Deleting signature {candidate}")

    coefs = df %>% dplyr::filter(denovos==candidate) %>% dplyr::pull(coefs) %>%
      setNames(df %>% dplyr::filter(denovos==candidate) %>% dplyr::pull(signature))

    updated_dfs = delete.signature_aux(denovo=denovo, exposure=exposure,
                                       coefs=coefs, sigName=candidate)

    updated_init_dfs = delete.signature_aux(
      denovo=get_nmf_initial_parameters(x, what="nmf")[[type]]$beta_dn_param,
      exposure=get_nmf_initial_parameters(x, what="nmf")[[type]]$alpha,
      coefs=coefs, sigName=candidate)

    old_n_denovo = nrow(denovo)
    scores = get_scores(x, types=type) %>%
      dplyr::filter(parname=="K" & value!=old_n_denovo) %>%
      dplyr::select(-type)
    x = set_scores(x, scores, type=type, what="nmf")
    x = set_denovo_signatures(x, type=type,
                              sigs=wide_to_long(updated_dfs$denovo, what="beta"))
    x = set_exposures(x, type=type,
                      expos=wide_to_long(updated_dfs$exposure, what="exposures"))
    x = set_nmf_init_params(x, type=type,
                            denovo=updated_init_dfs$denovo,
                            expos=updated_init_dfs$exposure)
    x = recompute_scores(x, types=type)

    cli::cli_process_done("Signature {candidate} discarded")
  }
  return(x)
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
                                            return_weights = FALSE) {
  d = v %*% Z
  bvec = c(1, rep(0, length(d)))
  C = cbind(rep(1, length(d)), diag(length(d)))
  pis =
    quadprog::solve.QP(
      Dmat = Rinv,
      factorized = TRUE,
      dvec = d,
      Amat = C,
      bvec = bvec,
      meq = 1
    )$solution
  # pis[pis < filt_pi] = 0
  reconstructed_vector = Z %*% pis

  vec1 = matrix(v)
  rownames(vec1) = rownames(reconstructed_vector)

  cos_sim = lsa::cosine(as.numeric(v), as.numeric(reconstructed_vector)) %>%
    as.numeric()

  if (!is.na(cos_sim) & cos_sim > delta) {
    if (return_weights) return(pis[pis > filt_pi] %>% setNames(colnames(Z)[pis > filt_pi]))
    return(colnames(Z)[pis > filt_pi])
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
    qp = solve.quadratic.optimization(
      a=denovo[c(i),],
      b=rbind(fixed, denovo[-c(i),]),
      filt_pi=0,
      delta=0.9,
      return_weights=TRUE
    )
    df[rownames(denovo[c(i), ]), names(qp[[1]])] = qp[[1]]
  }

  # adding reconstruction score column to dataframe
  ss = unlist(lapply(
    rownames(df),
    function(x) compute_refinement_score_aux(fixed=fixed, denovo=denovo, coefs=as.numeric(df[x, ]), sigName=x)
  ))
  ss[is.nan(ss)] = 0
  df$scores = ss

  if (matrix == TRUE) {
    return(df)
  } else {
    return(
      df %>% tibble::rownames_to_column(var="denovos") %>%
        reshape2::melt(id=c("denovos","scores"), variable.name="signature", value.name="coefs") %>%
        dplyr::filter(denovos!=signature) %>%
        dplyr::rowwise() %>%
        dplyr::mutate(cosine=lsa::cosine(as.numeric(denovo[denovos,]),
                                         as.numeric(rbind(fixed, denovo)[signature,]))[[1]])
    )
  }
}


# input :
#   fixed ------> fixed signatures (wide dataframe)
#   denovo -----> denovo signatures (wide dataframe)
#   coefs ------> linear combination coefficients (numeric vector)
#   <sigName> --> signature of interest (character)
# output:
#   numeric ----> cosine similarity (true signature vs. reconstructed signature)

compute_refinement_score_aux = function(fixed, denovo, coefs, sigName) {
  reconstructed_vector = t(rbind(fixed, denovo)) %*% coefs
  score = basilica:::cosine.vector(reconstructed_vector, denovo[sigName,])
  return(round(score, digits=3))
}


recompute_scores = function(x, types=get_types(x)) {
  for (type in types) {
    scores = lapply(c("bic","aic","llik"), function(score_id) {
        tibble::tibble(seed=get_seed(x)$nmf[[type]],
                       score_id=score_id,
                       score=compute_score(x, type=type, score_id=score_id),
                       parname="K",
                       value=get_n_denovo(x)[[type]])
    }) %>% dplyr::bind_rows() %>%
      dplyr::add_row(get_scores(x, types=type) %>% dplyr::select(-type))

    x = set_scores(x, scores=scores, type=type, what="nmf")
  }
  return(x)
}


compute_score = function(x, type, score_id) {
  if (score_id == "bic") compute_bic(x, type=type)
  else if (score_id == "aic") compute_aic(x, type=type)
  else if (score_id == "llik") compute_likelihood(x, type=type)
  else NULL
}


compute_likelihood = function(x, type) {
  betas = get_signatures(x, types=type, matrix=T)[[type]]
  alphas = get_exposure(x, types=type, matrix=T)[[type]]
  counts = get_input(x, types=type, matrix=T)[[type]]
  theta = rowSums(counts)
  alphas_hat = lapply(rownames(alphas), function(sample_id) alphas[sample_id, ] * theta[sample_id]) %>%
    do.call(rbind, .)

  rate = as.matrix(alphas_hat) %*% as.matrix(betas)

  dpois(as.matrix(counts), lambda=rate, log=T) %>% sum()
}


compute_bic = function(x, type) {
  llik = compute_likelihood(x, type=type)
  n_pars = compute_n_parameters(x, type=type)
  n_pars * log(length(get_samples(x))) - (2 * llik)
}

compute_aic = function(x, type) {
  llik = compute_likelihood(x, type=type)
  n_pars = compute_n_parameters(x, type=type)
  2 * n_pars - (2 * llik)
}


compute_n_parameters = function(x, type) {
  n_denovo = get_n_denovo(x)[[type]]
  expos = get_exposure(x, types=type, matrix=T)[[type]]
  n_pars = n_denovo
  if (!is.null(expos)) n_pars = n_pars + dim(expos)[1] * dim(expos)[2]
  return(n_pars)
}






