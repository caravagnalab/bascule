counts_dbs = readRDS("~/Desktop/DBS_counts.rds")
catalogue_dbs = readRDS("~/Desktop/DBS_cosmic_catalogue.rds")

dbs.fit = two_steps_inference(counts_dbs, k_list=0:5,
                              enforce_sparsity1=T,
                              enforce_sparsity2=F,
                              input_catalogue=catalogue_dbs,
                              py=py,
                              regularizer="KL")
dbs1 = dbs.fit$tot
dbs1 %>% plot_mutations(epsilon = T)
make_plots(dbs1, what="DBS", epsilon=TRUE)





counts_sbs = readRDS("~/Desktop/SBS_counts.rds")
# catalogue = readRDS("~/Desktop/DBS_cosmic_catalogue.rds")

sbs.fit1 = two_steps_inference(as.data.frame(counts_sbs),
                              k_list=0:5,
                              enforce_sparsity1=T,
                              enforce_sparsity2=F,
                              input_catalogue=COSMIC_filt_merged,
                              keep_sigs=c("SBS1","SBS40 SBS3 SBS5"),
                              py=py,
                              run_on_resid=TRUE,
                              regularizer="KL")
sbs1 = sbs.fit1$tot
make_plots(sbs1, what="SBS", epsilon=TRUE)
patchwork::wrap_plots(plot_exposures(sbs1, sampleIDs = "UPN06_REL"),
                      plot_mutations(sbs1, sampleIDs="UPN06_REL"),
                      plot_signatures(sbs1), design="AABBBB
                                                      ##CCCC
                                                      ##CCCC")


sbs.fit2 = two_steps_inference(as.data.frame(counts_sbs),
                              k_list=0:5,
                              enforce_sparsity1=T,
                              enforce_sparsity2=F,
                              input_catalogue=COSMIC_filt_merged[c("SBS1","SBS40 SBS3 SBS5"),],
                              keep_sigs=c("SBS1","SBS40 SBS3 SBS5"),
                              py=py,
                              run_on_resid=TRUE,
                              regularizer="KL")
sbs2 = sbs.fit2$tot
make_plots(sbs2, what="SBS", epsilon=TRUE)
patchwork::wrap_plots(plot_exposures(sbs2, sampleIDs = "UPN06_REL"),
                      plot_mutations(sbs2, sampleIDs="UPN06_REL"),
                      plot_signatures(sbs2), design="AABBBB
                                                      ##CCCC
                                                      ##CCCC")



sbs.fit3 = two_steps_inference(as.data.frame(counts_sbs),
                                k_list=0:5,
                                enforce_sparsity1=T,
                                enforce_sparsity2=F,
                                input_catalogue=COSMIC_filt_merged[c("SBS1","SBS40 SBS3 SBS5"),],
                                keep_sigs=c("SBS1","SBS40 SBS3 SBS5"),
                                py=py,
                                run_on_resid=FALSE,
                                regularizer="KL")
sbs3 = sbs.fit3$tot
make_plots(sbs3, what="SBS", epsilon=TRUE)
patchwork::wrap_plots(plot_exposures(sbs3, sampleIDs = "UPN06_REL"),
                      plot_mutations(sbs3, sampleIDs="UPN06_REL"),
                      plot_signatures(sbs3), design="AABBBB
                                                      ##CCCC
                                                      ##CCCC")




