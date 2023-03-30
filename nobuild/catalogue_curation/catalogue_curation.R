library(ggplot2)

## Datasets ####

rm_sigs = list(
  not_validated = paste0("SBS", c("7c", "7d", 8, 9, "10c", "10d", 12, 16, "17a",
                                  19, 23, 25, 29, 32, 33, 34, 37, 38, 39, 40,
                                  41, 84, 85, 86, 89, 91, 92, 93, 94)),
  unk_aetiology = paste0("SBS", c(8, 12, 16, "17a", "17b", 19, 23, 28, 33, 34,
                                  37, 39, 40, 41, 89, 91, 93, 94)),
  seq_artifacts = paste0("SBS", c(27,43:60,95))
)

catalogue = COSMIC_catalogue[!rownames(COSMIC_catalogue)%in%rm_sigs$not_validated,]

thr = c(0.01, 0.02)
p = c(0.3, 0.5)

catalogue_long = wide_to_long(catalogue)
catalogue_long.filt_perc = filter_catalogue_by_perc(catalogue_long, p=p)
catalogue_long.filt_thr = catalogue_long %>%
  dplyr::mutate(dens=ifelse(dens<thr,0,dens)) %>% normalize_cat()

cosine.all = compute_similarity(catalogue)
# cosine.filt_perc = compute_similarity(long_to_wide(catalogue_long.filt_perc))
cosine.filt_thr = compute_similarity(long_to_wide(catalogue_long.filt_thr))

cosine_long.all = get_cosine_long(cosine.all) %>%
  dplyr::mutate(type="nofilt") %>%
  dplyr::add_row(
    get_cosine_long(cosine.filt_perc) %>%
      dplyr::mutate(type="filt_perc")
  ) %>%
  dplyr::add_row(
    get_cosine_long(cosine.filt_thr) %>%
      dplyr::mutate(type="filt_thr")
  )



## Visualizations ####

### beta distributions ####
catalogue_long %>%
  dplyr::left_join(compute_percentile(catalogue_long, p=p), by="sig") %>%
  ggplot() +
  geom_histogram(aes(x=dens), bins=50) +
  geom_vline(aes(xintercept=perc)) +
  facet_wrap(~sig, scales="free") +
  theme_bw()


catalogue_long %>%
  dplyr::left_join(compute_percentile(catalogue_long, p=p), by="sig") %>%

  dplyr::group_by(sig) %>%
  dplyr::summarise(cum_fn=list( ecdf(dens) ), perc=unique(perc)) %>%
  dplyr::ungroup() %>%

  dplyr::mutate(x_eval=list(seq(0,1,length.out=200))) %>%
  tidyr::unnest("x_eval") %>%

  dplyr::rowwise() %>%
  dplyr::mutate(y_eval=cum_fn(x_eval)) %>%
  dplyr::ungroup() %>%

  dplyr::group_by(sig) %>%
  dplyr::mutate(max_y=min(max(y_eval),1)) %>%
  dplyr::ungroup() %>%

  dplyr::filter(y_eval<max_y) %>%

  ggplot() +
  # geom_histogram(aes(x=dens), bins=50) +
  geom_point(aes(x=x_eval, y=y_eval), size=.5) +
  geom_vline(aes(xintercept=perc)) +
  facet_wrap(~sig, scales="free") +
  theme_bw()


catalogue_long %>%
  ggplot() +
  geom_boxplot(aes(y=sig, x=dens), outlier.size=.5) +
  theme_bw()


### cosine similarities ####
pheatmap(cosine.all, cluster_rows=F, cluster_cols=F)
pheatmap(cosine.filt_perc, cluster_rows=F, cluster_cols=F)
pheatmap(cosine.filt_thr, cluster_rows=F, cluster_cols=F)


cosine_long.all %>%
  ggplot() +
  geom_histogram(aes(x=simil), bins=50) +
  ggh4x::facet_nested_wrap(~sig1+type) +
  theme_bw()

cosine_long.all %>%
  ggplot() +
  geom_boxplot(aes(x=type, y=simil)) +
  theme_bw()


## Save the report ####
catalogue_long.complete = catalogue_long %>% dplyr::mutate(type="nofilt", perc=0) %>%
  dplyr::add_row(catalogue_long.filt_perc %>% dplyr::mutate(type="filt_perc")) %>%
  dplyr::add_row(catalogue_long.filt_thr %>% dplyr::mutate(type="filt_thr"))


sigsnames = catalogue_long.complete$sig %>% unique()

plots = lapply(sigsnames, function(ss) {
  p1 = cosine_long.all %>%
    dplyr::filter(sig1==ss | sig2==ss, simil!=1) %>%
    ggplot() +
    geom_histogram(aes(x=simil), bins=50) +
    ggh4x::facet_nested_wrap(~type) +
    theme_bw()
  p2 = cosine_long.all %>%
    dplyr::filter(sig1==ss | sig2==ss, simil!=1) %>%
    ggplot() +
    geom_boxplot(aes(x=type, y=simil)) +
    theme_bw()
  p3 = catalogue_long.complete %>%
    add_subs() %>%
    dplyr::filter(sig==ss) %>%
    ggplot() +
    geom_histogram(aes(x=context, y=dens), stat="identity") +
    ggh4x::facet_nested(sig+type~subs) +
    theme_bw()

  patchwork::wrap_plots((p1+p2)/p3) &
    patchwork::plot_annotation(title=paste0("Signature ", ss),
                               subtitle=paste0("thr=", thr, ", p=", p))
  } ) %>% setNames(sigsnames)

pdf(paste0("./nobuild/catalogue_cur.p",p,".pdf"), height = 6, width = 12)
print(plots)
dev.off()

