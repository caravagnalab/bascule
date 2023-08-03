two_steps_inference = function(x,
                               k,
                               enforce_sparsity1=TRUE,
                               enforce_sparsity2=FALSE,
                               min_exposure=0.2,
                               reference_catalogue=COSMIC_filt_merged,
                               filtered_catalogue=TRUE,
                               keep_sigs=c("SBS1", "SBS5"),
                               lr=0.05,
                               n_steps=500,
                               groups=NULL,
                               clusters=NULL,
                               nonparametric = FALSE,
                               compile=FALSE,
                               cohort="MyCohort",
                               regularizer="cosine",
                               py=NULL,
                               hyperparameters=NULL,
                               reg_weight=1.,
                               reg_bic=TRUE,
                               CUDA=FALSE,
                               verbose=FALSE,
                               seed_list=c(10,27,92),
                               initializ_seed = FALSE,
                               save_runs_seed = TRUE,
                               initializ_pars_fit = TRUE,
                               new_hier = FALSE,
                               regul_denovo = TRUE,
                               regul_fixed = TRUE,
                               save_all_fits = FALSE,
                               do_initial_fit = FALSE) {

  TIME = as.POSIXct(Sys.time(), format = "%H:%M:%S")

  if (!is.null(reference_catalogue)) {
    x_ref = pyfit(
      counts = x,
      py = py,
      lr = lr,
      n_steps = n_steps,
      k_list = 0,
      hyperparameters = hyperparameters,
      groups = groups,
      clusters = NULL,
      input_catalogue = reference_catalogue,
      regularizer = regularizer,
      reg_weight = reg_weight,
      reg_bic = reg_bic,
      compile = compile,
      CUDA = CUDA,
      verbose = verbose,
      enforce_sparsity = enforce_sparsity1,
      stage = "random_noise",
      seed_list = seed_list,
      initializ_seed = initializ_seed,
      save_runs_seed = save_runs_seed,
      initializ_pars_fit = initializ_pars_fit,
      new_hier = new_hier,
      regul_denovo = regul_denovo,
      regul_fixed = regul_fixed,
      do_initial_fit = FALSE) %>%

      create_basilica_obj(input_catalogue=reference_catalogue[keep_sigs, ],
                          reference_catalogue=reference_catalogue,
                          cohort=cohort,
                          filtered_catalogue=filtered_catalogue)

    x_ref_filt = x_ref %>% filter_sigs_low_expos(min_exp=min_exposure, keep_sigs=keep_sigs)
    catalogue2 = get_signatures(x_ref_filt)
  } else {
    x_ref = x_ref_filt = catalogue2 = NULL
    residues = FALSE
  }


  x_dn = pyfit(
    counts = x,
    py = py,
    lr = lr,
    n_steps = n_steps,
    k_list = k,
    hyperparameters = hyperparameters,
    groups = groups,
    clusters = clusters,
    nonparametric = nonparametric,
    input_catalogue = catalogue2,
    regularizer = regularizer,
    reg_weight = reg_weight,
    reg_bic = reg_bic,
    compile = compile,
    CUDA = CUDA,
    verbose = verbose,
    enforce_sparsity = enforce_sparsity2,
    stage = "",
    regul_compare = NULL,
    seed_list = seed_list,
    initializ_seed = initializ_seed,
    save_runs_seed = save_runs_seed,
    initializ_pars_fit = initializ_pars_fit,
    new_hier = new_hier,
    regul_denovo = regul_denovo,
    regul_fixed = regul_fixed,
    save_all_fits = save_all_fits,
    do_initial_fit = do_initial_fit) %>%

    create_basilica_obj(input_catalogue=reference_catalogue[keep_sigs,],
                        reference_catalogue=reference_catalogue,
                        cohort=cohort,
                        filtered_catalogue=filtered_catalogue)

  TIME = difftime(as.POSIXct(Sys.time(), format = "%H:%M:%S"), TIME, units = "mins")

  # merged = merge_fits(x_ref, x_dn, x_ref_filt, min_exposure, keep_sigs)
  merged = x_dn
  merged$time = TIME
  merged$k_list = k

  # return(list("tot"=merged, "step1"=x_ref, "step1_filt"=x_ref_filt, "step2"=x_dn))
  return(merged)
}


merge_fits = function(x1, x2, x1_filt, min_exposure, keep_sigs, residues) {
  if (!residues)
    return(x2)

  merged = x1 %>% filter_sigs_low_expos(min_exp=min_exposure, keep_sigs=keep_sigs)
  merged$n_denovo = x2$n_denovo
  merged$fit$denovo_signatures = NULL

  if (merged$n_denovo == 0)
    merged$fit$exposure = normalize_exposures(x1_filt) %>% get_exposure()
  else {
    resid_expos = 1 - rowSums(merged$fit$exposure)
    denovo_norm = x2$fit$exposure / rowSums(x2$fit$exposure) * resid_expos

    merged$fit$exposure = cbind(merged$fit$exposure, denovo_norm)
    merged$fit$denovo_signatures = x2$fit$denovo_signatures
  }

  merged$color_palette = gen_palette(get_signatures(merged) %>% nrow()) %>%
    setNames(sort(rownames(get_signatures(merged))))

  return(merged)
}


create_basilica_obj_simul = function(simul, cohort="MySimul") {
  ss = list()
  class(ss) = "basilica_obj"
  ss$cohort = cohort
  ss$n_samples = nrow(simul$counts[[1]])
  ss$n_denovo = nrow(simul$alpha[[1]])
  ss$input = list("counts"=simul$counts[[1]],
                  "reference_catalogue"=simul$beta[[1]],
                  "input_catalogue"=simul$beta[[1]])

  ss$fit = list("input_catalogue"=simul$beta[[1]],
                "catalogue_signatures"=simul$beta[[1]],
                "denovo_signatures"=NULL,
                "exposure"=simul$alpha[[1]],
                "x"=simul$counts[[1]])

  centr = simul$alpha_prior[[1]] %>% as.data.frame()
  rownames(centr) = paste0("G",rownames(centr))
  ss$fit$params = list("alpha_prior"=centr)

  ss$groups = simul$groups[[1]]
  if (is.null(ss$groups)) {
    ss$groups = simul$counts[[1]] %>% tibble::rownames_to_column(var="sample") %>%
      tidyr::separate(sample, into=c("gid","sample"), sep="_") %>%
      dplyr::pull(gid) %>% stringr::str_replace_all("G","")
  }
  ss$groups = paste0("G", ss$groups)

  ss$color_palette = gen_palette(get_signatures(ss) %>% nrow()) %>%
    setNames(sort(rownames(get_signatures(ss))))

  if ("private" %in% colnames(simul))
    ss$sigs[["private"]] = simul$private[[1]]

  if ("private_shared" %in% colnames(simul))
    ss$sigs[["private_shared"]] = simul$private_shared[[1]]

  if ("shared" %in% colnames(simul))
    ss$sigs[["shared"]] = simul$shared[[1]]

  return(ss)
}


# create_basilica_obj_simul_old = function(simul, cohort="MySimul") {
#   ss = list()
#   class(ss) = "basilica_obj"
#   ss$cohort = cohort
#   ss$n_samples = nrow(simul$x[[1]])
#   ss$n_denovo = nrow(simul$exp_denovo[[1]])
#   ss$input = list("counts"=simul$x[[1]],
#                   "reference_catalogue"=simul$ref_cat[[1]],
#                   "input_catalogue"=simul$exp_fixed[[1]])
#
#   ss$fit = list("input_catalogue"=simul$exp_fixed[[1]],
#                 "catalogue_signatures"=simul$ref_cat[[1]],
#                 "denovo_signatures"=simul$exp_denovo[[1]],
#                 "exposure"=simul$exp_exposure[[1]],
#                 "x"=simul$x[[1]])
#
#   if (!is.null(simul$groups))
#     ss$groups = simul$groups[[1]]
#
#   ss$color_palette = gen_palette(get_signatures(ss) %>% nrow()) %>%
#     setNames(sort(rownames(get_signatures(ss))))
#
#   ss$private_sigs = list()
#   if ("private_rare" %in% colnames(simul))
#     ss$private_sigs[["private_rare"]] = simul$private_rare[[1]]
#
#   if ("private_common" %in% colnames(simul))
#     ss$private_sigs[["private_common"]] = simul$private_common[[1]]
#
#   return(ss)
# }


normalize_exposures = function(x) {
  x$fit$exposure = x$fit$exposure / rowSums(x$fit$exposure)
  return(x)
}


filter_sigs_low_expos = function(x, min_exp=0.15, keep_sigs=NULL) {
  sbs_keep = get_exposure(x, long=TRUE) %>%
    # x$fit$exposure %>%
    # tibble::rownames_to_column(var="sample") %>%
    # reshape2::melt(id="sample", value.name="alpha", variable.name="sigs") %>%

    dplyr::mutate(Exposure=ifelse(Exposure < min_exp, 0, Exposure)) %>%
    dplyr::filter(Exposure > 0) %>%
    dplyr::pull(Signature) %>% unique() %>% as.character()

  if (!is.null(keep_sigs)) sbs_keep = c(sbs_keep, keep_sigs) %>% unique()
  if (length(sbs_keep) == 0) return(x)

  x$fit$input_catalogue =
    x$fit$input_catalogue[intersect(sbs_keep, get_fixed_signames(x)),]

  x$fit$catalogue_signatures =
    x$fit$catalogue_signatures[intersect(sbs_keep, get_catalogue_signames(x)),]

  x$fit$denovo_signatures =
    x$fit$denovo_signatures[intersect(sbs_keep, get_dn_signames(x)),]

  x$fit$exposure = x$fit$exposure[,sbs_keep]

  # x$input$input_catalogue =
  #   x$input$input_catalogue[intersect(sbs_keep, get_fixed_signames(x)),]

  return(x)
}


map_groups = function(groups) {
  n_gr = groups %>% unique() %>% length()
  grps.int = 1:n_gr %>% setNames(groups %>% unique())

  grps.int = grps.int-1

  return(grps.int[groups] %>% setNames(NULL))
}


create_basilica_obj = function(fit, input_catalogue, reference_catalogue, cohort="MyCohort", filtered_catalogue=TRUE) {
  # fit is the output of "pyfit"
  obj = list()
  class(obj) = "basilica_obj"

  obj$cohort = cohort
  obj$n_samples = nrow(fit$x)

  if ("denovo_signatures" %in% names(fit))
    obj$n_denovo = nrow(fit$denovo_signatures) else
      obj$n_denovo = 0

  if (filtered_catalogue && obj$n_denovo > 0)
    fit$denovo_signatures = renormalize_denovo_thr(fit$denovo_signatures)

  obj$input = list("counts"=fit$x,
                   "reference_catalogue"=reference_catalogue,
                   "input_catalogue"=fit$input_catalogue)

  fit$catalogue_signatures = fit$input_catalogue

  obj$fit = fit
  obj$groups = fit$groups

  obj$color_palette = gen_palette(get_signatures(obj) %>% nrow()) %>%
    setNames(sort(rownames(get_signatures(obj))))

  obj$time = fit$time

  return(obj)
}


compute_residuals = function(x, min_exp=0.2, keep_sigs=NULL) {
  xf = x %>% filter_sigs_low_expos(min_exp=min_exp, keep_sigs=keep_sigs)
  sigs_order = rownames(xf$fit$input_catalogue)

  orig_counts = xf$input$counts
  t1 = (xf$fit$exposure * rowSums(orig_counts))[,sigs_order] %>% as.matrix()
  t2 = xf$fit$input_catalogue[sigs_order,] %>% as.matrix()
  reconstr_counts = t1 %*% t2

  residual_counts = abs(orig_counts - reconstr_counts)

  return(residual_counts)
}



gen_palette = function(n) {
  # library(RColorBrewer)
  # qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
  # col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
  # return(sample(col_vector, n))
  return(ggsci::pal_simpsons()(n))
}



filter_exposures = function(x, min_expos=0.1) {
  expos = get_exposure(x)
  expos[expos < min_expos] = 0

  x$fit$exposure = x$fit$params$alpha = expos / rowSums(expos)

  if (!is.null(x$fit$params$alpha_prior)) {
    x$fit$params$alpha_prior[x$fit$params$alpha_prior < min_expos] = 0
    x$fit$params$alpha_prior = x$fit$params$alpha_prior / rowSums(x$fit$params$alpha_prior)
  }

  return(x)
}


recompute_centroids = function(x.fit) {
  grps = unique(get_groups(x.fit)) %>% as.character()
  expos = get_exposure(x.fit, add_groups=T)
  new_centroids = lapply(grps, function(gid)
    expos %>% dplyr::filter(groups==gid) %>%
      dplyr::select(-groups) %>%
      colMeans()
    ) %>% do.call(rbind, .)
  rownames(new_centroids) = paste0("G",grps)

  return(x.fit %>% set_new_centroids(new_centroids))
}


set_new_centroids = function(x.fit, new_centroids) {
  x.fit$fit$params$alpha_prior_fit = get_centroids(x.fit, normalize=FALSE)
  x.fit$fit$params$alpha_prior = as.data.frame(new_centroids)
  return(x.fit)
}


merge_clusters = function(x.fit, cutoff=0.8) {
  alpha_prior = get_centroids(x.fit, normalize=TRUE)
  if (nrow(alpha_prior) == 1) return(x.fit)
  cosine_simil = lsa::cosine(t(alpha_prior)) %>% as.data.frame()
  cosine_simil[lower.tri(cosine_simil, diag=T)] = 0

  rownames(cosine_simil) = colnames(cosine_simil) = rownames(alpha_prior)

  merging = cosine_simil %>% tibble::rownames_to_column(var="gid1") %>%
    reshape2::melt(id="gid1", variable.name="gid2", value.name="cosine") %>%
    dplyr::mutate(gid1=as.character(gid1), gid2=as.character(gid2)) %>%
    dplyr::filter(cosine > cutoff) %>% dplyr::select(-cosine) %>%
    dplyr::group_by(gid1) %>%
    dplyr::summarise(cl_old=list(c(gid1,gid2) %>% unique())) %>% dplyr::ungroup() %>%
    dplyr::rename(cl_name=gid1)

  if (nrow(merging) == 0) return(x.fit)

  grps = get_groups(x.fit)
  for (i in 1:nrow(merging)) {
    old_cl = dplyr::pull(merging[i,], cl_old)[[1]]
    new_cl = dplyr::pull(merging[i,], cl_name)
    grps[grps %in% old_cl] = new_cl
  }

  x.fit$groups = grps
  x.fit = recompute_centroids(x.fit)

  return(x.fit)
}


fix_assignments = function(x.fit, cutoff=0.8, max_iters=20) {
  init_fit = x.fit
  curr_z = get_groups(x.fit)
  i = 0
  repeat {
    i = i+1
    x.fit = recompute_centroids(x.fit)
    # x.fit = merge_clusters(x.fit, cutoff=cutoff)

    x.fit = recompute_assignments(x.fit)
    new_z = get_groups(x.fit)

    if (setequal(new_z, curr_z) || i > max_iters) break
    curr_z = new_z
  }

  if (i > max_iters)
    return(init_fit %>% recompute_centroids() %>% merge_clusters())

  return(x.fit %>% merge_clusters(cutoff=cutoff))
}


recompute_assignments = function(x.fit) {
  grps = get_groups(x.fit)
  centroids = get_centroids(x.fit, normalize=TRUE)[as.character(unique(grps)),]
  # pi = get_mixture_weights(x.fit)[as.character(unique(grps))]
  pi = table(get_groups(x.fit)) / x.fit$n_samples

  counts = get_data(x.fit, reconstructed=FALSE)
  n_muts = rowSums(counts)
  beta = get_signatures(x.fit)[colnames(centroids), ]

  z = c(); ll_k = data.frame() # K x N

  for (k in unique(grps)) {
    alpha_k = centroids[as.character(k),][rep(1, times=x.fit$n_samples),]
    rownames(alpha_k) = names(n_muts)
    rate = as.matrix(alpha_k * n_muts) %*% as.matrix(beta)

    ll_k = ll_k %>% dplyr::bind_rows(
      log(pi[as.character(k)]) +
        rowSums(dpois(x=as.matrix(counts), lambda=rate, log=TRUE))
    )
  }
  rownames(ll_k) = unique(grps)

  probs = logsumexp(ll_k) %>%
    dplyr::mutate(probs=exp(ll_k-logsumexpdiff))

  new_probs = probs %>% dplyr::select(sampleid, groupid, probs) %>%
    tidyr::pivot_wider(values_from=probs, names_from=groupid) %>%
    tibble::column_to_rownames(var="sampleid")

  x.fit = set_assignment_probs(x.fit, new_probs)

  new_z = probs %>% dplyr::select(sampleid, groupid, probs) %>% unique() %>%
    dplyr::group_by(sampleid) %>%
    dplyr::filter(probs==max(probs)) %>% dplyr::select(sampleid, groupid) %>%
    tibble::column_to_rownames(var="sampleid") %>%
    dplyr::mutate(groupid=as.character(groupid))

  x.fit$groups = new_z[rownames(counts), "groupid"]

  return(x.fit)
}


set_assignment_probs = function(x.fit, new_probs) {
  x.fit$fit$post_probs = new_probs
  return(x.fit)
}



logsumexp = function(ll_k) {
  summed_lk = t(ll_k) %>% as.data.frame() %>%
    tibble::rownames_to_column(var="sampleid") %>%
    reshape2::melt(id="sampleid", variable.name="groupid", value.name="ll_k") %>%

    dplyr::group_by(sampleid) %>%

    # compute max among groups for each sample
    dplyr::mutate(m=max(ll_k)) %>%
    dplyr::mutate(expdiff=exp(ll_k-m)) %>%
    dplyr::reframe(logsumexpdiff=log(sum(expdiff)) + m,
                   ll_k=ll_k, groupid=groupid)

  return(summed_lk)
}





