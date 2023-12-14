library(ggplot2)
source("~/GitHub/basilica/nobuild/catalogue_curation/helper_fns.R")

## Datasets ####

COSMIC_merged = readRDS("~/GitHub/basilica/nobuild/catalogue_curation/COSMIC_merged.rds")
COSMIC_filt_merged = filter_cosmic_catalogue(COSMIC_merged, sbs5=COSMIC_sbs["SBS5",])
usethis::use_data(COSMIC_filt_merged, overwrite=T)
save(COSMIC_filt_merged, file="~/GitHub/basilica/nobuild/catalogue_curation/COSMIC_filt_merged.rda")

catalogue = COSMIC_sbs
COSMIC_sbs_filt = filter_cosmic_catalogue(catalogue, sbs5=catalogue["SBS5",])
usethis::use_data(COSMIC_sbs_filt, overwrite=T)
save(COSMIC_sbs_filt, file="~/GitHub/basilica/nobuild/catalogue_curation/COSMIC_filt.rda")

catalogue = COSMIC_sbs_filt
cosine.all = compute_similarity(catalogue)
# cosine.filt_perc = compute_similarity(long_to_wide(catalogue_long.filt_perc))
cosine.filt = compute_similarity(long_to_wide(catalogue_long.filt))

cosine_long.all = get_cosine_long(cosine.all) %>%
  dplyr::add_row(get_cosine_long(cosine.filt))



## Visualizations ####

### beta distributions ####
catalogue_long %>%
  dplyr::left_join(compute_percentile(catalogue_long, p=.5), by="sig") %>%
  ggplot() +
  geom_histogram(aes(x=dens), bins=50) +
  geom_vline(aes(xintercept=thr)) +
  facet_wrap(~sig, scales="free") +
  theme_bw()


catalogue_long %>%
  dplyr::left_join(compute_percentile(catalogue_long, p=.5), by="sig") %>%

  dplyr::group_by(sig) %>%
  dplyr::summarise(cum_fn=list( ecdf(dens) ), thr=unique(thr)) %>%
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
  geom_vline(aes(xintercept=thr)) +
  facet_wrap(~sig, scales="free") +
  theme_bw()


catalogue_long %>%
  ggplot() +
  geom_boxplot(aes(y=sig, x=dens), outlier.size=.5) +
  theme_bw()


### cosine similarities ####
pheatmap(cosine.all, cluster_rows=F, cluster_cols=F)
filt_heatmap = lapply(cosine.filt, function(i) pheatmap(i, cluster_rows=F, cluster_cols=F))


cosine_long.all %>%
  dplyr::filter(type %in% c("nofilt", "filt_perc_0.05")) %>%
  ggplot() +
  geom_histogram(aes(x=simil), bins=50) +
  ggh4x::facet_nested_wrap(~sig1+type) +
  theme_bw()

cosine_long.all %>%
  dplyr::filter(type %in% c("nofilt", "filt_thr_0.02")) %>%
  ggplot() +
  geom_histogram(aes(x=simil), bins=50) +
  ggh4x::facet_nested_wrap(~sig1+type) +
  theme_bw()

cosine_long.all %>%
  ggplot() +
  geom_boxplot(aes(x=type, y=simil)) +
  theme_bw()


## Save the report ####
catalogue_long.complete = catalogue_long %>% dplyr::mutate(type="nofilt", thr=0) %>%
  dplyr::add_row(catalogue_long.filt)


sigsnames = catalogue_long.complete$sig %>% unique()

plots = lapply(sigsnames, function(ss) {
  p1 = cosine_long.all %>%
    dplyr::filter(sig1==ss | sig2==ss, simil!=1) %>%
    ggplot() +
    geom_histogram(aes(x=simil), bins=50) +
    facet_grid(~type) +
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
    patchwork::plot_annotation(title=paste0("Signature ", ss))
  } ) %>% setNames(sigsnames)

pdf(paste0("./nobuild/catalogue_curation/catalogue_cur.pdf"), height = 10, width = 12)
print(plots)
dev.off()


