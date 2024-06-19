# Convert signatures names based on external object or catalogue #####

convert_dn_names = function(x, x.simul=NULL, reference_cat=NULL, cutoff=0.8) {
  if (is.null(x.simul) & is.null(reference_cat)) {
    cli::cli_alert_warning("No signatures as input. Returning the original object.")
    return(x)
  }
  assigned_missing = get_assigned_missing(x, x.simul=x.simul, reference_cat=reference_cat, cutoff=cutoff)

  map_names = lapply(names(assigned_missing), function(tid) {
    am_t = assigned_missing[[tid]]
    # names -> reference names; values -> fit names
    c(am_t$assigned_tp, am_t$added_fp %>% setNames(am_t$added_fp))
  }) %>% unlist()

  rename_object(x, map_names)
}


get_assigned_missing = function(x, x.simul=NULL, reference_cat=NULL, cutoff=0.8) {
  types = get_types(x)
  lapply(types, function(tid) {
    sigs.fit = get_signatures(x, matrix=T)[[tid]]
    if (!is.null(x.simul)) sigs.simul = get_signatures(x.simul, matrix=T)[[tid]]
    else if (!is.null(reference_cat)) sigs.simul = reference_cat[[tid]]

    assigned = compare_sigs_inf_gt(sigs.fit, sigs.simul, cutoff=cutoff)
    missing = setdiff(rownames(sigs.simul), names(assigned))
    added = setdiff(rownames(sigs.fit), assigned)

    return(list("assigned_tp"=assigned, "missing_fn"=missing, "added_fp"=added))
  }) %>% setNames(types)
}


compare_sigs_inf_gt = function(sigs.fit, sigs.simul, cutoff=0.8) {
  common = intersect(rownames(sigs.fit), rownames(sigs.simul))
  unique_inf = setdiff(rownames(sigs.fit), common)
  unique_gt = setdiff(rownames(sigs.simul), common)

  if (length(unique_inf) == 0 || length(unique_gt) == 0)
    return(common %>% setNames(common))

  total_sigs = rbind(sigs.fit[!rownames(sigs.fit) %in% common,],
                     sigs.simul[!rownames(sigs.simul) %in% common,])
  cosine_matr = as.data.frame(lsa::cosine(t(total_sigs)))[unique_gt, ]
  cosine_matr = cosine_matr[, colnames(cosine_matr) %in% unique_inf, drop=F]

  if (length(unique_inf) == 1 && length(unique_gt) == 1) {
    cosine_matr = as.data.frame(cosine_matr)
    rownames(cosine_matr) = unique_gt
    colnames(cosine_matr) = unique_inf
  }

  assign_similar = cosine_matr %>% as.data.frame() %>%
    tibble::rownames_to_column(var="gt") %>%
    reshape2::melt(id="gt", variable.name="inf", value.name="cosine") %>%
    dplyr::filter(cosine >= cutoff)

  if (nrow(assign_similar) == 0) return(common %>% setNames(common))

  assign_similar = assign_similar %>%
    dplyr::group_by(gt) %>%
    dplyr::mutate(inf=as.character(inf)) %>%
    dplyr::filter(cosine == max(cosine)) %>% dplyr::arrange(gt)

  # if (nrow(sigs.simul) > nrow(sigs.fit))
  if (any(duplicated(assign_similar$inf)))
    assign_similar = assign_similar %>% dplyr::group_by(inf) %>%
    dplyr::filter(cosine == max(cosine)) %>% dplyr::ungroup()

  assigned = c(common, assign_similar$inf) %>% setNames(c(common, assign_similar$gt))

  return(assigned)
}




rename_object = function(x, map_names, types=get_types(x)) {
  ## MISSING CONVERSION OF STORED OBJECTS
  mapp = names(map_names) %>% setNames(map_names)
  for (tid in types) {
    if (is.null(get_denovo_signames(x)[[tid]])) next

    alpha_long = get_exposure(x)[[tid]] %>%
      dplyr::mutate(sigs=mapp[sigs])
    dn_long = get_denovo_signatures(x)[[tid]] %>%
      dplyr::mutate(sigs=mapp[sigs])

    x = set_exposures(x, expos=alpha_long, type=tid)
    x = set_denovo_signatures(x, sigs=dn_long, type=tid)

    init_params = get_nmf_initial_parameters(x, what="nmf")[[tid]]
    x = set_nmf_init_params(x, type=tid,
                            denovo=init_params$beta_dn_param %>%
                              wide_to_long(what="beta") %>%
                              dplyr::mutate(sigs=mapp[sigs]) %>%
                              long_to_wide(what="beta"),
                            expos=init_params$alpha %>%
                              wide_to_long(what="exposures") %>%
                              dplyr::mutate(sigs=mapp[sigs]) %>%
                              long_to_wide(what="exposures"))
  }

  if (have_groups(x)) {
    new_colnames = colnames(x$clustering$pyro$params$infered_params$alpha_prior)
    for (new_name in names(map_names)) {
      old_name = map_names[[new_name]]
      new_colnames = new_colnames %>%
        stringr::str_replace_all(pattern=old_name, replacement=new_name)
    }

    colnames(x$clustering$pyro$params$infered_params$alpha_prior) =
      colnames(x$clustering$pyro$params$init_params$alpha_prior) =
      colnames(x$clustering$pyro$params$init_params$variances) =
      new_colnames

    x$clustering$centroids = x$clustering$pyro$params$infered_params$alpha_prior %>%
      wide_to_long(what="exposures") %>%
      dplyr::rename(clusters=samples) %>%
      dplyr::mutate(clusters=paste0("G",as.integer(clusters)-1))
  }

  return(x)
}



# Set parameters and dataframes #####

set_denovo_signatures = function(x, sigs, type) {
  x$nmf[[type]]$beta_denovo = sigs
  x$nmf[[type]]$pyro$beta_denovo = sigs
  x$nmf[[type]]$pyro$params$infered_params$beta_d = sigs
  return(x)
}

set_fixed_signatures = function(x, sigs, type) {
  x$nmf[[type]]$beta_fixed = sigs
  x$nmf[[type]]$pyro$beta_fixed = sigs
  x$nmf[[type]]$pyro$params$infered_params$beta_f = long_to_wide(sigs, what="beta")
  return(x)
}

set_exposures = function(x, expos, type) {
  x$nmf[[type]]$exposure = expos
  x$nmf[[type]]$pyro$exposure = expos
  x$nmf[[type]]$pyro$params$infered_params$alpha = long_to_wide(expos, what="exposures")
  return(x)
}

set_nmf_init_params = function(x, type, denovo=NULL, expos=NULL) {
  if (!is.null(denovo)) x$nmf[[type]]$pyro$params$init_params$beta_dn_param = denovo
  if (!is.null(expos)) x$nmf[[type]]$pyro$params$init_params$alpha = expos
  return(x)
}


# set_scores = function(x, scores, type, what="nmf") {
#   x[[what]][[type]]$pyro$QC$scores = scores
#   return(x)
# }


remove_alternatives = function(x, types=get_types(x)) {
  for (tid in types)
    x$nmf[[tid]]$pyro$alternatives = NULL
  x
}




