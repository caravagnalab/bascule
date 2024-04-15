devtools::load_all("~/GitHub/simbasilica/")
load_deps()

counts = get_input(example_dataset, matrix=T)
reference_cat = list("SBS"=COSMIC_sbs_filt, "DBS"=COSMIC_dbs)
max_K = sapply(get_signames(example_dataset), length) %>% max()
# K_list = 2:(max_K-2+1)
K_list = 2:2

x = fit(counts=counts, k_list=K_list, reference_cat=reference_cat, keep_sigs=c("SBS1","SBS5"))
x2 = refinement(x)

x2 %>% plot_exposures()
x %>% plot_exposures()

x2 %>% plot_signatures()
x %>% plot_signatures()

convert_dn_names(x2, reference_cat = list("SBS"=COSMIC_sbs_filt))

