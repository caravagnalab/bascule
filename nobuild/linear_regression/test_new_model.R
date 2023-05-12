if (Sys.getenv("GITHUB_PATH") == "") path=paste0("~/dati_elenab/signatures/") else path=Sys.getenv("GITHUB_PATH")

# header = read.csv("./nobuild/test_groups/counts_all.tsv", header=F, sep="\t", nrows=1) %>% setNames(NULL)
# counts = read.csv("./nobuild/test_groups/counts_all.tsv", header=T, row.names=1, sep="\t")
# colnames(counts) = header

py = reticulate::import_from_path(module="pybasilica", path="~/GitHub/pybasilica/")
devtools::load_all()
library(ggplot2)

simul = readRDS("/Users/elenab/GitHub/simbasilica/nobuild/simulations/simulations_elena/simul.N350.G2.s1.Rds")
x.true = create_basilica_obj_simul(simul)

# fit with the old model
# simul.fit = readRDS("/Users/elenab/GitHub/simbasilica/simulations/fit.N350.G2.s1.Rds")
sigs_true = x.true %>% get_signatures() %>% rownames()
counts = simul$x[[1]]
x = two_steps_inference(counts, k_list=0:7,
                        input_catalogue=COSMIC_filtered,
                        enforce_sparsity1=TRUE,
                        enforce_sparsity2=FALSE)

x$tot %>% plot_exposures()
x$tot %>% plot_signatures()
x$tot %>% plot_similarity_reference()
x$tot %>% plot_mutations()


x.nosparse = two_steps_inference(counts, k_list=0:7,
                                 input_catalogue=COSMIC_filtered,
                                 enforce_sparsity1=FALSE,
                                 enforce_sparsity2=FALSE)
x.nosparse$tot %>% plot_exposures()
x.nosparse$tot %>% plot_signatures()
x.nosparse$tot %>% plot_similarity_reference()



sigs_fit = x$tot %>% get_signatures() %>% rownames()
wrong = sigs_fit[!grepl(paste(sigs_true, collapse="|"), sigs_fit)]
tot_sigs = unique(c(sigs_fit, sigs_true))
sigs_colors = gen_palette(n=length(tot_sigs)) %>% setNames(tot_sigs)

muts_simul = x.true %>% plot_mutations() + labs(title="Mutations")
alpha_simul = x.true %>% plot_exposures(sample_name = F, cls=sigs_colors, sort_by = "SBS2") +
  theme(legend.position="bottom")
sigs_simul = x.true %>% plot_signatures(cls=sigs_colors)
plots_simul = patchwork::wrap_plots(sigs_simul + (muts_simul/alpha_simul))


muts_fit = x$tot %>% plot_mutations(reconstructed=TRUE) + labs(title="Mutations")
alpha_fit = x$tot %>% plot_exposures(sample_name = F, cls=sigs_colors, sort_by = "SBS2") +
  theme(legend.position="bottom")
sigs_fit = x$tot %>% plot_signatures(cls=sigs_colors)
plots_fit = patchwork::wrap_plots(sigs_fit + (muts_fit/alpha_fit))


# patients with missing signature's exposure higher than 0.1
idd = x$tot$fit$exposure[x$tot$fit$exposure[, wrong] > 0.1, ] %>% rownames()
sigs = colnames(x$tot$fit$exposure)[colSums(x$tot$fit$exposure[idd,] > 0.001) > 0]
x$tot$fit$x[idd, ]

p_tot = lapply(idd, function(pid) {
  p1 = plot_mutations(x$tot, sampleIDs=c(pid))
  p2 = plot_exposures(x$tot, sampleIDs=c(pid), muts = T, cls=sigs_colors)
  p2_true = plot_exposures(x.true, sampleIDs=c(pid), muts = T, cls=sigs_colors, sort_by = "SBS2") +
    labs(title="True exposure")
  p3 = plot_signatures(x$tot, cls = sigs_colors)
  patchwork::wrap_plots(p3 + (p1 / (p2+p2_true))) + patchwork::plot_annotation(title=paste0("Patient ", pid))
  })



idd_low = x$tot$fit$exposure[x$tot$fit$exposure[, wrong] < 0.1, ] %>% rownames() %>% sample(1)
sigs_low = colnames(x$tot$fit$exposure)
x$tot$fit$x[idd_low, ]

p_tot_low = lapply(idd_low, function(pid) {
  p1 = plot_mutations(x$tot, sampleIDs=c(pid))
  p2 = plot_exposures(x$tot, sampleIDs=c(pid), muts = T, cls=sigs_colors)
  p2_true = plot_exposures(x.true, sampleIDs=c(pid), muts = T, cls=sigs_colors, sort_by = "SBS2") +
    labs(title="True exposure")
  p3 = plot_signatures(x$tot, cls=sigs_colors)
  patchwork::wrap_plots(p3 + (p1 / (p2+p2_true))) + patchwork::plot_annotation(title=paste0("Patient ", pid))
})


pdf("~/GitHub/basilica/nobuild/linear_regression/simulated_fits.pdf", height = 10, width = 18)
plots_simul + patchwork::plot_annotation(title="Simulated (true) values")
plots_fit + patchwork::plot_annotation(title="Estimated values")
p_tot
p_tot_low
dev.off()



x$tot$fit$exposure * rowSums(x$tot$fit$x)



x.true %>% plot_exposures()
x.true %>% plot_similarity_reference()




library(ggplot2)
col_pal = gen_palette(20)
exp1 = x1$exposure %>%
  tibble::rownames_to_column(var="sample") %>%
  reshape2::melt(id="sample", value.name="alpha", variable.name="sigs") %>%

  dplyr::mutate(alpha=ifelse(alpha < 0.08, 0, alpha)) %>%
  dplyr::filter(alpha>0) %>%

  ggplot() +
  geom_bar(aes(x=sample, y=alpha, fill=sigs), stat="identity") +
  scale_fill_manual(values=col_pal) +
  theme_bw()
  # geom_histogram(aes(x=alpha), bins=100) +
  # scale_y_log10()

exp2 = x1.sparse$exposure %>%
  tibble::rownames_to_column(var="sample") %>%
  reshape2::melt(id="sample", value.name="alpha", variable.name="sigs") %>%

  dplyr::mutate(alpha=ifelse(alpha < 0.1, 0, alpha)) %>%
  dplyr::filter(alpha>0) %>%

  ggplot() +
  geom_bar(aes(x=sample, y=alpha, fill=sigs), stat="identity") +
  scale_fill_manual(values=col_pal) +
  theme_bw()
  # geom_histogram(aes(x=alpha), bins=100) +
  # scale_y_log10()

patchwork::wrap_plots(exp1 / exp2)



# x1 = fit(x = counts_groups$counts, py=py, k=0:7,
#         cohort="simul1", input_catalogue=COSMIC_filtered,
#         # input_catalogue=reference_sub["SBS1.5",],
#         reference_catalogue=reference_sub, delta=0.85, phi=0.1,
#         reg_weight=1., reg_bic=TRUE, filtered_cat = TRUE)


