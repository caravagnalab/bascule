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

get_signames = function(x, types=get_types(x)) {
  lapply(types, function(t)
    get_signatures(x, types=t, matrix=TRUE)[[t]] %>% rownames()) %>%
    setNames(types)
}


get_fixed_signatures = function(x, types=get_types(x), samples=get_samples(x),
                                clusters=get_cluster_labels(x), matrix=FALSE) {
  return(get_signatures_aux(x=x, types=types, matrix=matrix, what="fixed"))
}


get_denovo_signatures = function(x, types=get_types(x), samples=get_samples(x),
                                 clusters=get_cluster_labels(x), matrix=FALSE) {
  return(get_signatures_aux(x=x, types=types, matrix=matrix, what="denovo"))
}


# Get signatures, it can subset by types and return a list of matrices.
get_signatures = function(x, types=get_types(x), matrix=FALSE) {
  return(get_signatures_aux(x=x, types=types, matrix=matrix, what="all"))
}


## what %in% c("denovo","fixed","all")
get_signatures_aux = function(x, what, types=get_types(x), matrix=FALSE) {
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


get_input = function(x, types=get_types(x), samples=get_samples(x),
                     clusters=get_cluster_labels(x), matrix=FALSE) {
  out = lapply(types, function(t) {
    w = x$input[[t]]$counts %>%
      dplyr::filter(samples %in% !!samples)

    if(!is.null(clusters)) {
      which_selection = x %>%
        get_cluster_assignments(samples=samples, clusters=clusters) %>%
        dplyr::pull(samples)

      w = w %>% dplyr::filter(samples %in% which_selection)
    }

    return(w)
    }) %>% setNames(types)

  if(matrix)
    out = lapply(out, function(df_t) long_to_wide(df_t, what="counts")) %>%
      setNames(types)
  return(out)
}

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



# Get cluster labels ("C1", "C2")
get_cluster_labels = function(x) {
  if(is.null(x$clustering)) return(NULL)

  x$clustering$cluster$clusters %>% unique()
}

# Get clustering assignments (tibble)
get_cluster_assignments = function(x, samples=get_samples(x), clusters=get_cluster_labels(x)) {
  if(is.null(x$clustering)) return(NULL)

  return(
    x$clustering$clusters %>%
      dplyr::filter(samples %in% !!samples,
                    clusters %in% !!clusters)
  )
}


# Get type of data used for signatures: e.g., "SBS", "DBS"
get_types = function(x) {
  if (is.null(x)) return(NULL)
  return(x$input %>% names())
}


# Get samples names
get_samples = function(x) {x$input[[1]]$counts$samples %>% unique()}


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




