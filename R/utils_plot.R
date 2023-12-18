reformat_contexts = function(signatures, what) {
  if(what == "SBS") {
    signatures = signatures %>%
      dplyr::mutate(variant=paste0(substr(start=3, stop=3, features), ">", substr(start=5, stop=5, features)),
                    context=paste0(substr(start=1, stop=1, features), '_', substr(start=7, stop=7, features)))
  }

  if(what == "DBS") {
    signatures = signatures %>%
      dplyr::mutate(variant=paste0(substr(start=1, stop=2, features),">NN"),
                    context=substr(start=4, stop=5, features))
  }

  if(what == "ID") {
    signatures = signatures %>%
      dplyr::mutate(features=as.character(features),
                    variant=substr(start=1, stop=nchar(features) - 2, features),
                    context=substr(start=nchar(features), stop=nchar(features), features))
  }

  if(what == "CNV") {
    signatures = signatures %>%
      dplyr::mutate(features=as.character(features),
                    variant=paste0(str_split(features, pattern=":")[[1]][1], ":",
                                   str_split(features, pattern=":")[[1]][2]),
                    context= paste0(str_split(features,pattern=":")[[1]][3]))
  }

  return(signatures)
}


gen_palette = function(x=NULL, types=get_types(x), n=NULL) {
  if (!is.null(n)) return(ggsci::pal_simpsons()(n))
  return(gen_palette_aux(signames=get_signames(x, types=types)))
}


COSMIC_color_palette = function(seed=55) {
  catalogues = list(COSMIC_sbs_filt, COSMIC_filt_merged,
                    COSMIC_dbs, COSMIC_cn, COSMIC_indels, COSMIC_sbs)
  sigs = sapply(catalogues, rownames) %>% unlist() %>% unique()
  set.seed(seed)
  return(Polychrome::createPalette(length(sigs), c("#6B8E23","#4169E1"), M=1000,
                                   target="normal", range=c(15,80)) %>%
           setNames(sigs))
}


gen_palette_aux = function(signames=get_signames(x), seed=14) {
  types = names(signames)

  colss = lapply(types, function(tid) {
    cosmic_cols = COSMIC_color_palette()
    ref_sigs = intersect(signames[[tid]], names(cosmic_cols))
    dn_sigs = setdiff(signames[[tid]], names(cosmic_cols))
    n_dn = length(dn_sigs)
    set.seed(seed)
    dn_cols = Polychrome::createPalette(n_dn, c("#483D8B","#9e461c"), M=1000,
                                        target="normal", range=c(15,80))[1:n_dn] %>%
      setNames(dn_sigs)
    c(cosmic_cols[ref_sigs], dn_cols)
  }) %>% unlist()

  return(colss)
}


