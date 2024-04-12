devtools::load_all("~/GitHub/simbasilica/")
load_deps()

counts = get_input(example_dataset, matrix=T)
reference_cat = list("SBS"=COSMIC_sbs_filt, "DBS"=COSMIC_dbs)
max_K = sapply(get_signames(example_dataset), length) %>% max()
K_list = 2:(max_K-2+1)

x = fit(counts=counts, k_list=K_list, reference_cat=reference_cat, keep_sigs=c("SBS1","SBS5"))

x2 = refinement(x)


