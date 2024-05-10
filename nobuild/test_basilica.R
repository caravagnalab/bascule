devtools::load_all("~/GitHub/simbasilica/")
load_deps()
devtools::load_all("~/GitHub/basilica/")

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

