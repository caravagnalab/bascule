#' De novo signatures refinement.
#'
#' @description
#' Function to refine the inferred de novo signatures.
#' The function performs a linear combination between de novo and reference
#' signatures. If a de novo signature can be explained as a linear combination
#' of one or more reference signatures, it will be removed and its exposures
#' will be distributed among the similar signatures.
#'
#' @param x basilica object.
#' @param types List of variant types to perform de novo refinement on.
#'
#' @return basilica object.
#' @export refine_denovo_signatures

refine_denovo_signatures = function(x, types=get_types(x)) {

  alternatives = get_alternatives(x) %>% dplyr::filter(type %in% types)

  alternatives_refined = alternatives %>%

    dplyr::rename(tid=type) %>%
    dplyr::filter(parname=="K") %>%

    dplyr::rowwise() %>%
    dplyr::mutate(results_refined=list(
      refinement_single_fit(x, pyro_fit=pyro_fit, type=tid)
    )) %>%

    dplyr::mutate(value_refined=results_refined[["value_refined"]],
                  score_refined=results_refined[["score_refined"]],
                  pyro_fit_refined=results_refined[["pyro_fit_refined"]]) %>%
    dplyr::ungroup() %>%

    dplyr::rename(type=tid)

  best_fit = alternatives_refined %>%
    dplyr::group_by(type) %>%
    dplyr::filter(score_refined==min(score_refined))

  ## substitute the best fits as main objects and add the new alternatives
  for (i in 1:nrow(best_fit)) {
    tid = best_fit[[i,"type"]]
    seed = best_fit[[i,"seed"]]
    pyro_fit_refined = best_fit[[i,"pyro_fit_refined"]][[1]]

    x = set_attribute(x, what="nmf", type=tid, name="pyro", value=pyro_fit_refined)
    x = set_exposures(x, expos=pyro_fit_refined$exposure, type=tid)
    x = set_denovo_signatures(x, sigs=pyro_fit_refined$beta_denovo, type=tid)
  }

  new_alternatives = alternatives_refined %>%
    dplyr::select(parname, type, seed, value, value_refined, pyro_fit_refined, seed) %>%
    dplyr::rename(value_fit=value_refined, pyro_fit=pyro_fit_refined) %>%

    dplyr::bind_rows(alternatives %>% dplyr::filter(parname=="G"))

  x = set_alternatives(x, new_alternatives)

  return(x)
}


refinement_single_fit = function(x, type, pyro_fit=NULL) {
  if (!type %in% get_types(x)) {
    cli::cli_alert_warning("Type {type} not present. Returning the initial object.")
    return(x)
  }

  if (!is.null(pyro_fit)) {
    x = set_attribute(x, what="nmf", type=type, name="pyro", value=pyro_fit)
    x = set_exposures(x, expos=pyro_fit$exposure, type=type)
    x = set_denovo_signatures(x, sigs=pyro_fit$beta_denovo, type=type)
  }

  repeat {
    denovo = get_denovo_signatures(x, types=type, matrix=TRUE)[[type]]

    if (is.null(denovo)) break
    if (nrow(denovo) == 0) break

    fixed = get_fixed_signatures(x, types=type, matrix=TRUE)[[type]]
    exposure = get_exposure(x, types=type, matrix=TRUE)[[type]]

    df = qc_linearCombination(fixed=fixed, denovo=denovo, matrix=FALSE)
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

    x = set_denovo_signatures(x, type=type,
                              sigs=wide_to_long(updated_dfs$denovo, what="beta"))
    x = set_exposures(x, type=type,
                      expos=wide_to_long(updated_dfs$exposure, what="exposures"))
    x = set_nmf_init_params(x, type=type,
                            denovo=updated_init_dfs$denovo,
                            expos=updated_init_dfs$exposure)
    x = recompute_scores(x, type=type)

    cli::cli_process_done("Signature {candidate} discarded")
  }

  return(
    list("value_refined"=get_n_denovo(x)[[type]],
         "score_refined"=get_scores(x) %>%
           dplyr::filter(type==!!type, score_id=="bic") %>%
           dplyr::pull(score) %>% unlist(),
         "pyro_fit_refined"=list(x[["nmf"]][[type]][["pyro"]]))
  )
}



# Linear combination #####
solve_quadratic_optimization = function(a,
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
        solve_quadratic_optimization_aux(
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


solve_quadratic_optimization_aux = function(v,
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

qc_linearCombination = function(fixed, denovo, matrix=TRUE) {
  df = data.frame(matrix(0, nrow=nrow(denovo), ncol=nrow(denovo) + nrow(fixed)))
  rownames(df) = rownames(denovo)
  colnames(df) = c(rownames(fixed), rownames(denovo))

  for (i in 1:nrow(denovo)) {
    qp = solve_quadratic_optimization(
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
  score = cosine.vector(reconstructed_vector, denovo[sigName,])
  return(round(score, digits=3))
}


## Delete signatures ####
# input
#   denovo ----> denovo signatures (wide)
#   exposure --> exposure (wide)
#   coefs -----> linear combination coefficients
#   sigName ---> denovo dignature name to delete
# output
#   list (denovo, exposure)
delete.signature_aux = function(denovo, exposure, coefs, sigName) {

  if (!(sigName %in% rownames(denovo))) {
    cli::cli_alert_warning("Wrong signature selected!")
    return(NULL)
  }

  if (all(coefs==0)) {
    cli::cli_alert_warning("Can not delete! This signature is not explained by other signatures")
    return(x)
  }

  exp = exposure %>% dplyr::select(-dplyr::all_of(sigName))
  for (refname in names(coefs)) {
    exp[,refname] = exposure[,refname] + exposure[,sigName]*coefs[refname]
  }

  denovo = denovo[!rownames(denovo) %in% c(sigName), ]

  return(list(denovo=denovo, exposure=exp))
}


# Recompute fit scores + aux ####

set_alternatives = function(x, alternatives) {
  for (tid in unique(alternatives$type)) {
    if (tid %in% get_types(x)) x$nmf[[tid]]$pyro$alternatives = alternatives %>% dplyr::filter(type==tid)
    else x$clustering$pyro$alternatives = alternatives %>% dplyr::filter(type==tid)
  }

  return(x)
}


recompute_scores = function(x, type) {
  x$nmf[[type]]$pyro$QC = x$nmf[[type]]$pyro$QC %>%
    dplyr::mutate(value=replace(value, stat == "bic", compute_bic(x, type=type)),
                  value=replace(value, stat == "likelihood", compute_likelihood(x, type=type)))

  return(x)
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
