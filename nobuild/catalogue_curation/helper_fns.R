## Functions ####
compute_percentile = function(cat_long, p=.5)
  return(
    cat_long %>%
      dplyr::group_by(sig) %>%
      dplyr::summarise(perc=quantile(dens, p=p)) %>%
      dplyr::ungroup()
  )


filter_catalogue_by_perc = function(cat_long, p=c(.3)) {
  cat_long.final = data.frame()
  for (p_tmp in p)
    if (nrow(cat_long.final) == 0)
      cat_long.final = filter_single_perc(cat_long, p_tmp) else
        cat_long.final = cat_long.final %>%
          dplyr::add_row(filter_single_perc(cat_long, p_tmp))

      return(cat_long.final)
}


filter_single_perc = function(cat_long, p)
  return(
    cat_long %>%
      dplyr::left_join(compute_percentile(cat_long, p=p), by="sig") %>%
      dplyr::mutate(dens=ifelse(dens<=perc, 0, dens)) %>%

      # renormalize
      normalize_cat() %>%
      dplyr::mutate(type=paste0("filt_perc_",p))
  )


normalize_cat = function(cat_long)
  return(
    cat_long %>%
      dplyr::group_by(sig) %>%
      dplyr::mutate(dens=dens / sum(dens)) %>%
      dplyr::ungroup()
  )


compute_similarity = function(cat_wide)
  return(lsa::cosine(t(cat_wide)))


long_to_wide = function(cat_long) {
  if ("type" %in% colnames(cat_long)) {
    final = list()
    for (tt in cat_long$type %>% unique())
      final[[tt]] = cat_long %>%
        dplyr::filter(type==tt) %>%
        tidyr::pivot_wider(id_cols="sig", names_from="context", values_from="dens") %>%
        tibble::column_to_rownames(var="sig")
  } else
    final = cat_long %>%
      tidyr::pivot_wider(id_cols="sig", names_from="context", values_from="dens") %>%
      tibble::column_to_rownames(var="sig")

  return(final)
}


wide_to_long = function(cat_wide) {
  if (is.list(cat_wide)) {
    final = create_empty_long()
    for (i in names(cat_wide))
      final = final %>%
        dplyr::add_row(
          cat_wide %>%
            tibble::rownames_to_column(var="sig") %>%
            reshape2::melt(id="sig", variable.name="context", value.name="dens") %>%
            dplyr::mutate(type=i))

  } else
    final =  cat_wide %>%
      tibble::rownames_to_column(var="sig") %>%
      reshape2::melt(id="sig", variable.name="context", value.name="dens")
  return(final)
}

create_empty_long = function()
  return(
    data.frame() %>%
      dplyr::mutate(sig=as.character(NA),
                    context=as.character(NA),
                    dens=as.numeric(NA),
                    type=as.character(NA))
  )



get_cosine_long = function(cos_wide) {
  cos_wide[upper.tri(cos_wide)] = -1

  return(
    cos_wide %>% as.data.frame() %>%
      tibble::rownames_to_column(var="sig1") %>%
      reshape2::melt(id="sig1", variable.name="sig2", value.name="simil") %>%
      dplyr::filter(simil>-1)
  )
}


add_subs = function(cat_long)
  return(
    cat_long %>%
      dplyr::mutate(subs=gsub("\\[","_",context)) %>%
      dplyr::mutate(subs=gsub("\\]","_",subs)) %>%
      tidyr::separate("subs", into=c("left","subs","right"), sep="_") %>%
      dplyr::mutate(context=paste0(left,"_",right))
  )
