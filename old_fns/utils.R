pyfit = function(counts,
                 k_list,
                 lr = 0.005,
                 optim_gamma = 0.1,
                 n_steps = 2000,
                 stage = "",
                 py = NULL,
                 clusters = NULL,
                 nonparametric = TRUE,
                 dirichlet_prior = TRUE,
                 beta_fixed = NULL,
                 hyperparameters = NULL,
                 CUDA = FALSE,
                 compile = FALSE,
                 enforce_sparsity = TRUE,
                 store_parameters = FALSE,
                 regularizer = "cosine",
                 regul_compare = NULL,
                 reg_weight = 1,
                 seed_list = c(10),
                 regul_denovo = TRUE,
                 regul_fixed = TRUE,
                 store_fits = FALSE
                 ) {

  TIME = as.POSIXct(Sys.time(), format = "%H:%M:%S")

  if (is.null(py))
    py = reticulate::import("pybasilica")

  if (length(k_list) > 1) k_list = reticulate::r_to_py(as.integer(k_list)) else
    k_list = reticulate::r_to_py(list(as.integer(k_list)))

  if (length(seed_list) > 1) seed_list = reticulate::r_to_py(as.integer(seed_list)) else
    seed_list = reticulate::r_to_py(list(as.integer(seed_list)))

  if (!is.null(clusters)) clusters = as.integer(clusters)

  obj = py$fit(x = counts, k_list = k_list, lr = lr, optim_gamma = optim_gamma, n_steps = n_steps,
               cluster = clusters, beta_fixed = beta_fixed,
               hyperparameters = hyperparameters, nonparametric=nonparametric,
               dirichlet_prior = dirichlet_prior, enforce_sparsity = enforce_sparsity,
               store_parameters = store_parameters, regularizer = regularizer,
               reg_weight = reg_weight, regul_compare = regul_compare,
               regul_denovo = regul_denovo, regul_fixed = regul_fixed,
               stage = stage, seed = seed_list, compile_model = compile,
               CUDA = CUDA, store_fits = store_fits)

  TIME = difftime(as.POSIXct(Sys.time(), format = "%H:%M:%S"), TIME, units = "mins")

  if (is.list(obj)) {
    bestRun = obj[[1]]
    secondBest = obj[[2]]
  } else {
    bestRun = obj
    secondBest = NULL
  }

  # save python object data in a list
  data = get_list_from_py(bestRun, counts, beta_fixed, lr, n_steps)
  data$runs_seed = lapply(data$runs_seed, function(i) {
    i[["runs_scores"]] = i[["runs_seed"]] = NULL
    i$convert_to_dataframe(counts)
    get_list_from_py(i, counts, beta_fixed, lr, n_steps)
  })

  data$secondBest = get_list_from_py(secondBest, counts, beta_fixed, lr, n_steps)
  data$time = TIME

  return(data)
}


get_list_from_py = function(py_obj, counts, input_catalogue, lr, n_steps, save_stats=T) {
  if (is.null(py_obj)) return(NULL)

  x = list()
  x$x = py_obj$x

  x$exposure = py_obj$params$alpha
  x$denovo_signatures = py_obj$params$beta_d
  x$eps_var = py_obj$params$lambda_epsilon
  x$pi = py_obj$params$pi
  x$post_probs = py_obj$params$post_probs
  x$groups = py_obj$groups

  x$params = py_obj$params
  x$init_params = py_obj$init_params

  x$bic = py_obj$bic
  x$losses = py_obj$losses
  x$gradient_norms = py_obj$gradient_norms
  x$train_params = get_train_params(py_obj)
  x$hyperparameters = py_obj$hyperparameters
  try(expr = { x$seed = py_obj$seed })

  if (!save_stats) return(x)

  x$runs_seed = x$runs_scores = x$all_fits = NULL
  if ("runs_seed" %in% names(py_obj))
    x$runs_seed = py_obj$runs_seed

  if ("scores_K" %in% names(py_obj))
    x$runs_K = get_scores_from_py(py_obj$scores_K)

  if ("scores_CL" %in% names(py_obj))
    x$runs_CL = get_scores_from_py(py_obj$scores_CL) %>% dplyr::rename(G=K)

  if ("all_fits" %in% names(py_obj)) {
    if (length(py_obj$all_fits) > 0) x$all_fits = NULL
    x$all_fits = get_fits_from_py(py_obj$all_fits, x$x, x$input_catalogue, lr, n_steps)
  }

  return(x)
}


get_fits_from_py = function(fits, counts, beta_fixed, lr, n_steps)
  return(
    lapply(names(fits), function(i) {
      fits[[i]]$convert_to_dataframe(counts)
        get_list_from_py(fits[[i]], counts, beta_fixed, lr, n_steps, save_stats=F)
    }) %>%
      setNames(names(fits))
  )


get_scores_from_py = function(scores) {
  if (is.null(scores)) return(NULL)

  print(scores)
  print(replace_null(scores))

  res = replace_null(scores) %>%
  # res = purrr::discard(scores, is.null) %>%
    as.data.frame() %>%
    reshape2::melt(value.name="score") %>%
    tidyr::separate("variable", into=c("K", "seed", "score_id"), sep="[.]") %>%
    tibble::as_tibble()

  return(res)
}


replace_null = function(i) {
  j = purrr::map(i, ~ replace(.x, is.null(.x), NA))
  purrr::map(j, ~ (if(is.list(.x)) replace_null(.x) else .x))
}


get_train_params = function(obj) {
  if (!obj$store_parameters)
    return(NULL)
  train_params = obj$train_params
  samples_names = obj$params[["alpha"]] %>% rownames()
  bfixed_names = obj$beta_fixed %>% rownames()
  bdenovo_names = obj$params[["beta_d"]] %>% rownames()
  contexts = obj$params[["beta_d"]] %>% colnames()

  params = data.frame()

  for (i in 1:length(train_params)) {
    expos = train_params[[i]][["alpha"]] %>% as.data.frame()
    rownames(expos) = samples_names
    colnames(expos) = c(bfixed_names, bdenovo_names)

    if ("alpha_prior" %in% names(train_params[[i]])) {
      centroids = train_params[[i]][["alpha_prior"]] %>% as.data.frame()
      rownames(centroids) = (1:nrow(centroids)) -1
      colnames(centroids) = c(bfixed_names, bdenovo_names)
      centroids = centroids %>% tibble::rownames_to_column(var="rowname") %>%
        reshape2::melt(id="rowname",variable.name="columnname",value.name="value") %>%
        dplyr::mutate(iteration=i, paramname="centroid")
    } else { centroids = data.frame() }

    if ("pi" %in% names(train_params[[i]])) {
      pi = train_params[[i]][["pi"]] %>% as.numeric() %>% setNames((sort(unique(centroids$rowname))))
      pi = data.frame("rowname"=names(pi),"value"=pi,"iteration"=i,"paramname"="pi")
    } else { pi = data.frame() }

    sigs = train_params[[i]][["beta_d"]] %>% as.data.frame()
    rownames(sigs) = bdenovo_names
    colnames(sigs) = contexts

    params = params %>% dplyr::bind_rows(
      expos %>% tibble::rownames_to_column(var="rowname") %>%
        reshape2::melt(id="rowname",variable.name="columnname",value.name="value") %>%
        dplyr::mutate(iteration=i, paramname="alpha")
    ) %>% dplyr::bind_rows(
      sigs %>% tibble::rownames_to_column(var="rowname") %>%
        reshape2::melt(id="rowname",variable.name="columnname",value.name="value") %>%
        dplyr::mutate(iteration=i, paramname="beta_d")
    ) %>% dplyr::bind_rows(centroids) %>%
      dplyr::bind_rows(pi)

    # alpha_all = alpha_all %>% dplyr::add_row(tmp_a)
    # beta_d = beta_d %>% dplyr::add_row(tmp_b)
  }

  return(params)

  # return(tibble::tibble(alpha=list(alpha_all), beta_d=list(beta_d)))
}





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
    # cos_sim = cosine.vector(vec1, reconstructed_vector)

    if (is.null(exposures)) {
      if (!is.na(cos_sim) && cos_sim > delta) {
        if (return_weights) return(pis[pis > 0] %>% setNames(colnames(Z)[pis > 0]))
        return(colnames(Z)[pis > 0])
      }

    }

    # # exposure not null -> makes sense if we test more than one dn vs cosmic
    # } else {
    #   if (!is.na(cos_sim) && cos_sim > delta) {
    #     # if none or one linear combination
    #     if (sum(pis>0) <= 1) return(colnames(Z)[pis > 0])
    #
    #     keep.tmp = colnames(Z)[pis > 0]
    #     exposures.tmp = exposures[, keep.tmp]
    #     exposures.tmp$n_denovo = apply(exposures.tmp > thr_exposure, 1,
    #                                    function(x) length(unique(x))==1)
    #
    #     if (sum(exposures.tmp$n_denovo) >= nrow(exposures)*0.9)
    #       # return the reference that can be explained as a linear comb of denovo
    #       return(colnames(Z)[pis > 0])
    #   }
    # }
  return(NULL)
}

# Renormalize denovo

renormalize_denovo_thr = function(denovo, thr=0.02) {
  if (is.null(denovo)) return(NULL)
  denovo.tmp = denovo
  denovo.tmp[denovo.tmp < thr] = 0
  return(denovo.tmp / rowSums(denovo.tmp))
}



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

  return(df)
}

