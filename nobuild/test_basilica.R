devtools::load_all("~/GitHub/simbasilica/")
load_deps()

ex_d = example_dataset
ex_d$clustering = NULL

counts = get_input(example_dataset, matrix=T)
reference_cat = list("SBS"=COSMIC_sbs_filt, "DBS"=COSMIC_dbs)
max_K = sapply(get_signames(example_dataset), length) %>% max()
# K_list = 0:(max_K-2+1)
K_list = 0:2

x = fit(counts=counts, k_list=K_list, reference_cat=reference_cat, keep_sigs=c("SBS1","SBS5"), store_fits=TRUE)
x_init = x

# x %>% plot_scores()
# x %>% recompute_scores() %>% plot_scores()
scores = plot_scores(x) + labs(title="Scores")
scores_refined = plot_scores(refinement(x)) + labs(title="Scores post refinement")
data_diff = plot_data_differences(x) + labs(title="Data differences")
data_diff_refined = plot_data_differences(refinement(x)) + labs(title="Data differences post refinement")

pdf("~/Dropbox/dropbox_shared/2022. Basilica/simulations/plots_refinement.pdf", height=16, width=16)
patchwork::wrap_plots((scores+scores_refined), data_diff, data_diff_refined, ncol=1, heights=c(1,2,2), guides="collect")
# patchwork::wrap_plots(data_diff, data_diff_refined, ncol=1)
dev.off()

ex_d %>% plot_data(types="SBS", color_by_sigs=F) %>%
  patchwork::wrap_plots(x %>% plot_data(types="SBS", color_by_sigs=T))

ex_d %>% plot_data(types="SBS", color_by_sigs=F) %>%
  patchwork::wrap_plots(x %>% refinement() %>% plot_data(types="SBS", color_by_sigs=T))

plot_data_differences(x)

x %>% get_alternative_run(K=c("SBS"=0)) %>% plot_exposures() %>%
  patchwork::wrap_plots(
    x %>% refinement() %>% plot_exposures()
  )

compute_likelihood(x %>% refinement(), type="SBS")
compute_likelihood(x %>% get_alternative_run(K=c("SBS"=0)), type="SBS")

x2 = refinement(x)
x2 %>% plot_exposures() %>%
  patchwork::wrap_plots(x %>% refinement() %>% plot_exposures())

x2 %>% plot_exposures()
x %>% plot_exposures()

x2 %>% plot_signatures()
x %>% plot_signatures()

convert_dn_names(x, reference_cat = list("SBS"=COSMIC_sbs_filt))

get_scores(x)

