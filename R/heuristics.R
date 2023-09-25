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


fix_assignments = function(x.fit, cutoff=0.8, max_iters=20, merge_groups=FALSE) {
  cli::cli_warn("The function `fix_assignments` contains undetected bugs. The original object will be returned.")
  return(x.fit)
  if (!have_groups(x.fit)) return(x.fit)

  init_fit = x.fit
  curr_z = get_groups(x.fit)
  i = 0
  repeat {
    i = i+1
    x.fit = recompute_centroids(x.fit)

    x.fit = recompute_assignments(x.fit)
    new_z = get_groups(x.fit)

    if (setequal(new_z, curr_z) || i > max_iters) break
    curr_z = new_z
  }

  if (i > max_iters)
    return(init_fit %>% recompute_centroids() %>% merge_clusters())

  if (!merge_groups) return(x.fit)
  return(x.fit %>% merge_clusters(cutoff=cutoff))
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


recompute_assignments = function(x.fit) {
  cli::cli_warn("The function `recompute_assignments` contains undetected bugs.\nThe original object will be returned.")
  return(x.fit)

  grps = get_groups(x.fit)
  centroids = get_centroids(x.fit, normalize=TRUE)[as.character(unique(grps)),]
  pi = get_mixture_weights(x.fit)[as.character(unique(grps))]
  # pi = table(get_groups(x.fit)) / x.fit$n_samples

  counts = get_data(x.fit, reconstructed=FALSE)
  # n_muts = rowSums(counts)
  # beta = get_signatures(x.fit)[colnames(centroids), ]
  alpha = get_exposure(x.fit)

  # rate = as.matrix(diag(n_muts) %*% as.matrix(alpha)) %*% as.matrix(beta)
  # lprob_pois = rowSums(dpois(x=as.matrix(counts), lambda=rate, log=TRUE))

  z = c(); ll_k = data.frame() # K x N

  for (k in unique(grps)) {
    alpha_k = centroids[as.character(k),]# [rep(1, times=x.fit$n_samples),]
    # rownames(alpha_k) = names(n_muts)

    lprob_alpha = log(gtools::ddirichlet(x=alpha, alpha=as.numeric(centroids[as.character(k),])))
    names(lprob_alpha) = rownames(alpha)

    # lprob = data.frame(log(pi[as.character(k)]) + lprob_alpha) # + lprob_pois
    # colnames(lprob) = as.character(k)

    ll_k = ll_k %>% dplyr::bind_rows(
      as.data.frame(t(data.frame(log(pi[as.character(k)]) + lprob_alpha)))
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



## Signatures linear combination ####
# check if each signature of catalogue sign1 is a linear combination of catalogue sign2
filter_signatures_QP = function(sign1, sign2,
                                delta = 0.9,
                                filt_pi = 0.05,
                                thr_exposure = 0.05,
                                exposures = NULL,
                                return_weights = FALSE) {

  res_optimization =
    solve.quadratic.optimization(a = sign1,
                                 b = sign2,
                                 delta = delta,
                                 filt_pi = filt_pi,
                                 thr_exposure = thr_exposure,
                                 exposures = exposures,
                                 return_weights = return_weights) %>%
    purrr::discard(function(i) is.null(i))

  return(res_optimization)
  }

