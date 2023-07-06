## Functions ####
filter_cosmic_catalogue = function(catalogue, sbs5=COSMIC_catalogue["SBS5",]) {
  rm_sigs = list(
    not_validated = paste0("SBS", c(8, 9, "10c", "10d", 12, 16, "17a",
                                    19, 23, 29, 32, 33, 34, 37, 38, 39, 40,
                                    41, 84, 85, 86, 89, 91, 92, 93, 94)),
    unk_aetiology = paste0("SBS", c(8, 12, 16, "17a", "17b", 19, 23, 28, 33, 34,
                                    37, 39, 40, 41, 89, 91, 93, 94)),
    seq_artifacts = paste0("SBS", c(27,43:60,95))
  )

  keep = lapply(rownames(catalogue), function(sname) {
    snames = strsplit(sname, " ")[[1]]
    if (any(!snames %in% c(rm_sigs$not_validated, rm_sigs$seq_artifacts, "SBS5"))) return(sname)
  }) %>% unlist()

  # keep = rownames(COSMIC_merged)[!grepl(paste(rm_sigs$not_validated, collapse="|"), rownames(COSMIC_merged))]
  catalogue.tmp = catalogue[keep,]

  thr = c(0.01, 0.02)
  p = c(0.3, 0.5)

  catalogue_long = wide_to_long(catalogue.tmp)
  catalogue_long.filt = catalogue_long %>%
    filter_catalogue(p=p, thr=thr)


  ## Export filtered catalogue - thr = 0.02 #####
  catalogue_filt = (catalogue_long.filt %>% dplyr::filter(type=="filt_thr_0.02") %>% long_to_wide())$filt_thr_0.02
  # COSMIC_filt_merged["SBS5",] = catalogue["SBS5",]
  catalogue_filt["SBS5",] = sbs5
  catalogue_filt = catalogue_filt[!rownames(catalogue_filt) %in% "SBS40 SBS3 SBS5",]

  return(catalogue_filt)
}


compute_percentile = function(cat_long, p=.5)
  return(
    cat_long %>%
      dplyr::group_by(sig) %>%
      dplyr::summarise(thr=quantile(dens, p=p)) %>%
      dplyr::ungroup()
  )


filter_catalogue = function(cat_long, p=c(0.), thr=c(0.)) {
  cat_long.final = create_empty_df_long()
  # if (any(p>0) && any(thr>0)) stop("Both p and thr are greater than 0.")
  if (all(p==0) && all(thr==0)) stop("Not filtering. Both p and thr are 0.")

  for (p_tmp in p)
    cat_long.final = cat_long.final %>%
      dplyr::add_row(filter_single_perc(cat_long, p_tmp))

  for (t_tmp in thr)
    cat_long.final = cat_long.final %>%
      dplyr::add_row(filter_single_thr(cat_long, t_tmp))

  return(cat_long.final)
}


filter_single_thr = function(cat_long, thr) {
  if (thr == 0) return()
  return(
    cat_long %>%
      dplyr::mutate(dens=ifelse(dens<thr,0,dens)) %>%
      normalize_cat() %>%
      dplyr::mutate(type=paste0("filt_thr_",thr),
                    thr=thr)
  )
}


filter_single_perc = function(cat_long, p) {
  if (p == 0) return()
  return(
    cat_long %>%
      dplyr::left_join(compute_percentile(cat_long, p=p), by="sig") %>%
      dplyr::mutate(dens=ifelse(dens<=thr, 0, dens)) %>%

      # renormalize
      normalize_cat() %>%
      dplyr::mutate(type=paste0("filt_perc_",p))
  )
}


normalize_cat = function(cat_long)
  return(
    cat_long %>%
      dplyr::group_by(sig) %>%
      dplyr::mutate(dens=dens / sum(dens)) %>%
      dplyr::ungroup()
  )


compute_similarity = function(cat_wide) {
  if (is.data.frame(cat_wide)) return(lsa::cosine(t(cat_wide)) %>% as.data.frame())

  return(
    lapply(cat_wide, function(i) lsa::cosine(t(i)) %>% as.data.frame)
    )
}


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
  if (!is.data.frame(cat_wide)) {
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


create_empty_df_long = function()
  return(
    data.frame() %>%
      dplyr::mutate(sig=as.character(NA),
                    context=as.character(NA),
                    dens=as.numeric(NA),
                    type=as.character(NA),
                    thr=as.numeric(NA))
  )



get_cosine_long = function(cos_wide) {
  if (is.data.frame(cos_wide))
    cos_wide = list("nofilt"=cos_wide)

  cos_long.list = lapply(names(cos_wide), function(i) {
    cos_wide[[i]][upper.tri(cos_wide[[i]])] = -1
    return(
      cos_wide[[i]] %>% as.data.frame() %>%
        tibble::rownames_to_column(var="sig1") %>%
        reshape2::melt(id="sig1", variable.name="sig2", value.name="simil") %>%
        dplyr::filter(simil>-1) %>% dplyr::mutate(type=i)
    )
  })

  return(do.call(rbind, cos_long.list))
}


add_subs = function(cat_long)
  return(
    cat_long %>%
      dplyr::mutate(subs=gsub("\\[","_",context)) %>%
      dplyr::mutate(subs=gsub("\\]","_",subs)) %>%
      tidyr::separate("subs", into=c("left","subs","right"), sep="_") %>%
      dplyr::mutate(context=paste0(left,"_",right))
  )
