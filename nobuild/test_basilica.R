devtools::load_all("~/GitHub/simbascule/")
load_deps()
devtools::load_all("~/GitHub/bascule/")


## SIMULATED ####
counts = get_input(example_dataset, matrix=T)
reference_cat = list("SBS"=COSMIC_sbs_filt, "DBS"=COSMIC_dbs)
max_K = sapply(get_signames(example_dataset), length) %>% max()
K_list = 0:2
x = fit(counts=counts, k_list=K_list, n_steps=3000,
        reference_cat=reference_cat, keep_sigs=c("SBS1","SBS5"), store_fits=TRUE)
x_refined = refine_denovo_signatures(x)
x_refined_cluster = fit_clustering(x_refined, cluster=3)

example_dataset %>% plot_exposures()
x_refined_cluster %>% merge_clusters() %>% plot_exposures()

example_dataset %>% plot_signatures()
x_refined_cluster %>% merge_clusters() %>% plot_signatures()






## REAL ####

## Old fits #####
x = readRDS("~/Google Drive/My Drive/work/bascule_shared/fit_27052024/fit_wcat.Breast.Rds")
for (tid in get_types(x)) {
  alt_t = x$nmf[[tid]]$pyro$alternatives

  x$nmf[[tid]]$pyro$alternatives = alt_t %>%
    dplyr::rowwise() %>%
    dplyr::mutate(seed=pyro_fit[[1]]$seed,
                  pyro_fit=list(pyro_fit[[1]])) %>%
    dplyr::ungroup()
}

x$clustering$pyro$alternatives = x$clustering$pyro$alternatives %>%
  dplyr::rowwise() %>%
  dplyr::mutate(seed=pyro_fit[[1]]$pyro$seed,
                pyro_fit=list(pyro_fit[[1]])) %>%
  dplyr::ungroup()

x_refined = refine_denovo_signatures(x, types="SBS")

saveRDS(x_refined, "~/Desktop/fit_wcat.Breast_refined.Rds")

x %>% plot_signatures()

cls_list = x %>% plot_cluster_scores()
x %>% plot_cls_score_heatmap()


## Last fit ####
x_orig = readRDS("~/Google Drive/My Drive/work/bascule_shared/fit_12062024/fit_refined_cls2.Breast.Rds") %>% merge_clusters()
x_mapped = x_orig %>% convert_dn_names(reference_cat=list("SBS"=COSMIC_sbs_filt, "DBS"=COSMIC_dbs))

get_assigned_missing(x_orig, reference_cat=list("SBS"=COSMIC_sbs_filt, "DBS"=COSMIC_dbs))

x_orig %>% plot_centroids(quantile_thr=0.5)
x_mapped %>% plot_centroids(quantile_thr=0.5)

x_mapped %>% get_centroids() %>% dplyr::filter(sigs=="SBS40a") %>% dplyr::pull(value)
x %>% get_centroids() %>% dplyr::filter(sigs=="SBSD7") %>% dplyr::pull(value)

x %>% plot_exposures()
x %>% convert_dn_names(reference_cat=list("SBS"=COSMIC_sbs_filt, "DBS"=COSMIC_dbs)) %>% plot_exposures()

x %>% plot_centroids()
x %>% plot_exposures()
x %>% plot_centroids(quantile_thr = 0.5)
x %>% plot_exposures(quantile_thr = 0.5)





