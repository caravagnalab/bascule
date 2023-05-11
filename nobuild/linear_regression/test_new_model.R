if (Sys.getenv("GITHUB_PATH") == "") path=paste0("~/dati_elenab/signatures/") else path=Sys.getenv("GITHUB_PATH")


# header = read.csv("./nobuild/test_groups/counts_all.tsv", header=F, sep="\t", nrows=1) %>% setNames(NULL)
# counts = read.csv("./nobuild/test_groups/counts_all.tsv", header=T, row.names=1, sep="\t")
# colnames(counts) = header


source("~/GitHub/basilica/nobuild/linear_regression/functions.R")
py = reticulate::import_from_path(module="pybasilica", path="~/GitHub/pybasilica/")
devtools::load_all()
library(ggplot2)

simul = readRDS("/Users/elenab/GitHub/simbasilica/simulations/simul.N350.G2.s1.Rds")
x.true = create_basilica_obj_simul(simul)

# fit with the old model
# simul.fit = readRDS("/Users/elenab/GitHub/simbasilica/simulations/fit.N350.G2.s1.Rds")

counts = simul$x[[1]]
x = two_step_inference(counts, k_list=1:7,
                       input_catalogue=COSMIC_filtered,
                       enforce_sparsity1=TRUE,
                       enforce_sparsity2=FALSE)

x$tot %>% plot_exposures()
x$tot %>% plot_signatures()
x$tot %>% plot_similarity_reference()
x$tot %>% plot_mutations()


x.nosparse = two_step_inference(counts, k_list=1:7,
                                input_catalogue=COSMIC_filtered[rownames(COSMIC_filtered)!="SBS17b",],
                                enforce_sparsity1=FALSE,
                                enforce_sparsity2=FALSE)
x.nosparse$tot %>% plot_exposures()
x.nosparse$tot %>% plot_signatures()
x.nosparse$tot %>% plot_similarity_reference(reference=COSMIC_filtered["SBS17b",])


muts_simul = x.true %>% plot_mutations() + labs(title="Mutations")
alpha_simul = x.true %>% plot_exposures()
sigs_simul = x.true %>% plot_signatures()
plots_simul = patchwork::wrap_plots(sigs_simul + (muts_simul/alpha_simul))


muts_fit = x$tot %>% plot_mutations() + labs(title="Mutations")
alpha_fit = x$tot %>% plot_exposures()
sigs_fit = x$tot %>% plot_signatures()
plots_fit = patchwork::wrap_plots(sigs_fit + (muts_fit/alpha_fit))




idd = x$tot$fit$exposure[x$tot$fit$exposure[,c("SBS6")] > 0.2, ] %>% rownames()
sigs = colnames(x$tot$fit$exposure)[x$tot$fit$exposure[idd,] > 0.001]
x$tot$fit$x[idd, ]

p1 = plot_mutations(x$tot, sampleIDs=c(idd))
plot_mutations(x$tot, by_sig=T)
p2 = plot_exposures(x$tot, sampleIDs=c(idd), muts = T)
p3 = plot_signatures(x$tot, signatures = sigs)
patchwork::wrap_plots(p3 + (p1 / p2))



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


