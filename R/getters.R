get_fixed_signames = function(x, type="T1") {
  return(get_fixed_signatures(x, type=type)$sigs %>% unique())
}


get_denovo_signames = function(x, type="T1") {
  return(get_denovo_signatures(x, type=type)$sigs %>% unique())
}


get_input = function(x, types=get_types(x), samples=get_samples(x),
                     clusters=get_cluster_labels(x), matrix=FALSE) {
  out = lapply(types, function(t) {
    w = x$input[[t]] %>%
      dplyr::filter(samples %in% !!samples)

    if(!is.null(clusters)) {
      which_selection = x %>%
        get_cluster_assignments(samples=samples, clusters=clusters) %>%
        dplyr::pull(samples)

      w = w %>% filter(samples %in% which_selection)
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
                        clusters=get_cluster_labels(x), matrix=FALSE) {
  out = lapply(types, function(t) {
      w = x$nmf[[t]]$exposure %>%
        dplyr::filter(samples %in% !!samples)

      if(!is.null(clusters)) {
        which_selection = x %>%
          get_cluster_assignments(samples=samples, clusters=clusters) %>%
          dplyr::pull(samples)

        w = w %>% filter(samples %in% which_selection)
      }
      return(w)
    }) %>% setNames(types)

  if(matrix)
    out = lapply(out, function(df_t) long_to_wide(df_t, what="exposures")) %>%
      setNames(types)
  return(out)
}

# Get signatures, it can subset by types and return
# a list of matrices.
get_signatures = function(x, types=get_types(x), matrix=FALSE) {
  out = lapply(types, function(t) {
    w = x$nmf[[t]]$signatures %>%
      dplyr::filter(samples %in% !!samples)

    if(!is.null(clusters)) {
      which_selection = x %>%
        get_cluster_assignments(samples = samples, clusters = clusters) %>%
        dplyr::pull(samples)
      w = w %>% filter(samples %in% which_selection)
    }
    return(w)
  }) %>% setNames(types)

  if(matrix)
    out = lapply(out, function(df_t) long_to_wide(df_t, what="beta")) %>%
      setNames(types)
  return(out)
}


# Get cluster labels ("C1", "C2")
get_cluster_labels = function(x) {
  if(is.null(x$clustering)) return(NULL)

  x$clustering$cluster %>% unique()
}

# Get clustering assignments (tibble)
get_cluster_assignments = function(x, samples=get_samples(x), clusters=get_cluster_labels(x)) {
  if(is.null(x$clustering)) return(NULL)

  x$clustering %>%
    dplyr::filter(samples %in% samples,
                  cluster %in% clusters)
}


# Get type of data used for signatures: e.g., "SBS", "DBS"
get_types = function(x) {x$input %>% names()}


# Get samples names
get_samples = function(x) {x$input[[1]]$samples %>% unique()}


# what %in% c("beta","exposures","counts")
wide_to_long = function(dataframe, what) {
  cols = dplyr::case_when(
    what == "beta" ~ list(variables="features", ids="sigs"),
    what == "exposures" ~ list(variables="sigs", ids="samples"),
    what == "counts" ~ list(variables="features", ids="samples")
  )
  dataframe %>% as.data.frame() %>%
    tibble::rownames_to_column(var=cols$ids) %>%
    reshape2::melt(id=cols$ids, variable.name=cols$variables) %>%
    tibble::as_tibble()
}


long_to_wide = function(dataframe, what) {
  if (is.null(dataframe) | nrow(dataframe)==0) return(NULL)
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




