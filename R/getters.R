# Signatures #####
## names ####
get_signames = function(x, types=get_types(x)) {
  lapply(types, function(t)
    get_signatures(x, types=t, matrix=TRUE)[[t]] %>% rownames()) %>%
    setNames(types)
}

get_fixed_signames = function(x, types=get_types(x)) {
  lapply(types, function(t)
    get_fixed_signatures(x, types=t, matrix=TRUE)[[t]] %>% rownames()) %>%
    setNames(types)
}

get_denovo_signames = function(x, types=get_types(x)) {
  lapply(types, function(t)
    get_denovo_signatures(x, types=t, matrix=TRUE)[[t]] %>% rownames()) %>%
    setNames(types)
}

## signature profiles #####
get_signatures = function(x, types=get_types(x), matrix=FALSE) {
  return(get_signatures_aux(x=x, types=types, matrix=matrix, what="all"))
}

get_fixed_signatures = function(x, types=get_types(x), samples=get_samples(x),
                                clusters=get_cluster_labels(x), matrix=FALSE) {
  return(get_signatures_aux(x=x, types=types, matrix=matrix, what="fixed"))
}


get_denovo_signatures = function(x, types=get_types(x), samples=get_samples(x),
                                 clusters=get_cluster_labels(x), matrix=FALSE) {
  return(get_signatures_aux(x=x, types=types, matrix=matrix, what="denovo"))
}

get_signatures_aux = function(x, what, types=get_types(x), matrix=FALSE) {
  ## what %in% c("denovo","fixed","all")
  out = lapply(types, function(t) {
    if (what=="fixed") x$nmf[[t]]$beta_fixed else if (what=="denovo")
      x$nmf[[t]]$beta_denovo else if (what=="all")
        rbind(x$nmf[[t]]$beta_fixed, x$nmf[[t]]$beta_denovo) else NULL
  }) %>% setNames(types)

  if(matrix)
    out = lapply(out, function(df_t) long_to_wide(df_t, what="beta")) %>%
      setNames(types)
  return(out)
}



# Input #####
get_input = function(x, types=get_types(x), samples=get_samples(x),
                     clusters=get_cluster_labels(x), matrix=FALSE,
                     reconstructed=FALSE, add_groups=FALSE, by_sigs=FALSE) {
  out = lapply(types, function(tid) {
    if (reconstructed | by_sigs) {
      expos = get_exposure(x, types=tid, samples=samples, clusters=clusters, matrix=T)[[tid]]
      betas = get_signatures(x, types=tid, matrix=T)[[tid]]
      theta = rowSums(x$input[[tid]]$counts %>% long_to_wide(what="counts"))

      if (by_sigs)
        w = lapply(colnames(expos), function(signame) {
          expos_s = expos[, signame]
          betas_s = betas[signame, ]
          (as.matrix(expos_s*theta) %*% as.matrix(betas_s)) %>%
            wide_to_long(what="counts") %>%
            dplyr::mutate(sigs=signame)
        }) %>% setNames(colnames(expos))
      else
        w = (as.matrix(expos*theta) %*% as.matrix(betas)) %>%
          wide_to_long(what="counts")

    } else {
      w = x$input[[tid]]$counts %>%
        dplyr::filter(samples %in% !!samples)
    }

    if(!is.null(clusters)) {
      clusters_df = x %>%
        get_cluster_assignments(samples=samples, clusters=clusters)
      if (by_sigs)
        w = lapply(w, function(w_i)
          w_i %>% dplyr::right_join(clusters_df, by="samples")) %>%
          setNames(names(w))
      else w = w %>% dplyr::right_join(clusters_df, by="samples")
    }

    return(w)
    }) %>% setNames(types)

  if (matrix) {
    if (by_sigs)
      out = lapply(out, function(df_t) {
        lapply(df_t, function(df_t_s) {
          long_to_wide(df_t_s %>% dplyr::select(-dplyr::contains("clusters")), what="counts")
        }) %>% setNames(names(df_t))
      }) %>% setNames(types)
    else
      out = lapply(out, function(df_t)
        long_to_wide(df_t %>% dplyr::select(-dplyr::contains("clusters")), what="counts")) %>%
        setNames(types)
    }

  return(out)
}

get_input_signatures = function(x, types=get_types(x), matrix=F) {
  lapply(types, function(tid) {
    sigs = x$input[[tid]]$reference
    if (matrix) sigs = long_to_wide(sigs, what="beta")
    sigs
  }) %>% setNames(types)
}

get_input_signames = function(x, types=get_types(x)) {
  sigs = get_input_signatures(x, types=types, matrix=T)
  lapply(sigs, rownames) %>% setNames(names(sigs))
}


get_scores_CL = function(x) {
  if (is.null(x$fit$runs_CL)) return(NULL)
  return(x$fit$runs_CL %>% dplyr::select_if(dplyr::where(function(i) any(!is.na(i)))))
}


# Exposures ####
# Get exposure, it can subset by types, samples and clusters. It can return
# a list of matrices.
get_exposure = function(x, types=get_types(x), samples=get_samples(x),
                        clusters=get_cluster_labels(x), add_groups=FALSE,
                        matrix=FALSE) {
  out = lapply(types, function(t) {
      w = x$nmf[[t]]$exposure %>%
        dplyr::filter(samples %in% !!samples)

      if(!is.null(clusters)) {
        which_selection = x %>%
          get_cluster_assignments(samples=samples, clusters=clusters)

        w = w %>% dplyr::inner_join(which_selection, by="samples")
        if (!add_groups) w = w %>% dplyr::select(-clusters)
      }
      return(w)
    }) %>% setNames(types)

  if(matrix)
    out = lapply(out, function(df_t) long_to_wide(df_t, what="exposures")) %>%
      setNames(types)
  return(out)
}



# Clustering #####
get_centroids = function(x, matrix=F) {
  if (is.null(x$clustering$centroids)) return(NULL)

  unq_labels = get_cluster_labels(x)

  centr = x$clustering$centroids %>%
    dplyr::mutate(clusters=paste0("G",stringr::str_replace_all(clusters,"G","")),
                  sigs=stringr::str_replace_all(sigs,"^[0-9]+_","")) %>%
    dplyr::filter(clusters %in% unq_labels)

  if (matrix)
    return(centr %>%
             tidyr::pivot_wider(names_from="sigs", values_from="value") %>%
             tibble::column_to_rownames(var="clusters"))

  return(centr)
}

get_mixing_proportions = function(x) {
  if (!have_groups(x)) return(NULL)
  pis = get_params(x, what="clustering")[["pi"]]
  if (is.null(pis)) return(NULL)
  cnames = paste0("G",1:length(pis)-1)
  data.frame(value=pis, clusters=factor(cnames, levels=cnames))
}

# Get cluster labels ("C1", "C2")
get_cluster_labels = function(x) {
  if(is.null(x$clustering)) return(NULL)

  x$clustering$cluster$clusters %>% unique()
}

# Get clustering assignments (tibble)
get_cluster_assignments = function(x, samples=get_samples(x), clusters=get_cluster_labels(x)) {
  if(is.null(x$clustering)) return(NULL)

  x$clustering$clusters %>%
    dplyr::filter(samples %in% !!samples, clusters %in% !!clusters)
}



# Fit infos #####

# Get type of data used for signatures: e.g., "SBS", "DBS"
get_types = function(x) {
  if (is.null(x)) return(NULL)
  return(x$input %>% names())
}

get_fittypes = function(x) {
  if (is.null(x)) return(NULL)
  return(names(x)[names(x)!="input"])
}

# Get samples names
get_samples = function(x) {x$input[[1]]$counts$samples %>% unique()}

get_K = function(x, types=get_types(x)) {
  return(lapply(get_signames(x, types=types), length))
}

get_G = function(x, input=FALSE) {
  if (!input) return(length(get_cluster_labels(x)))
  return(
    get_pyro_stat(x, what="clustering",
                  statname="params")[["infered_params"]]$post_probs %>% ncol()
  )
}

get_seed = function(x, what, types=get_types(x)) {
  return(get_pyro_stat(x, what=what, statname="seed"))
}



# Auxiliary fns #####

# what %in% c("beta","exposures","counts")
wide_to_long = function(dataframe, what) {
  if (is.null(dataframe) || nrow(dataframe)==0) return(NULL)

  cols = dplyr::case_when(
    what == "beta" ~ list(variables="features", ids="sigs"),
    what == "exposures" ~ list(variables="sigs", ids="samples"),
    what == "counts" ~ list(variables="features", ids="samples")
  )
  dataframe %>% as.data.frame() %>%
    tibble::rownames_to_column(var=cols$ids) %>%
    reshape2::melt(id=cols$ids, variable.name=cols$variables) %>%
    dplyr::mutate(dplyr::across(is.factor, as.character)) %>%
    tibble::as_tibble()
}


long_to_wide = function(dataframe, what) {
  if (is.null(dataframe) || nrow(dataframe)==0) return(NULL)
  cols = dplyr::case_when(
    what == "beta" ~ list(variables="features", ids="sigs"),
    what == "exposures" ~ list(variables="sigs", ids="samples"),
    what == "counts" ~ list(variables="features", ids="samples")
  )
  dataframe %>%
    dplyr::select(-dplyr::contains("type")) %>%
    tidyr::pivot_wider(names_from=cols$variables, values_from="value") %>%
    tibble::column_to_rownames(var=cols$ids)
}



# Deprecated ####
get_beta_weights = function(x, types=get_types(x)) {
  if (is.null(get_params(x, what="nmf", types=types[1])[[1]]$beta_w))
    return(NULL)
  lapply(types, function(tid) {
    get_params(x, what="nmf", types=tid)[[1]]$beta_w %>%
      tibble::rownames_to_column(var="sigid") %>%
      reshape2::melt(variable.name="sigs") %>%
      dplyr::mutate(type=tid)
  }) %>% do.call(rbind, .)
}
