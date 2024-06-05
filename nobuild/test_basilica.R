devtools::load_all("~/GitHub/simbasilica/")
load_deps()
devtools::load_all("~/GitHub/basilica/")


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
x = readRDS("~/Google Drive/My Drive/work/basilica_shared/fit_27052024/fit_wcat.Breast.Rds")
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

