library(lsa)
library(dplyr)
library(ggplot2)
source("~/GitHub/simbasilica/nobuild/script_test/helper_fns.R")
py = reticulate::import_from_path("pybasilica","~/GitHub/pybasilica/")
devtools::load_all("~/GitHub/basilica/")
devtools::load_all("~/GitHub/simbasilica/")

# DBS ####
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




# SBS ####
counts_sbs = readRDS("~/Desktop/SBS_counts.rds")

## Whole catalogue
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
make_plots(sbs1, what="SBS", epsilon=TRUE, sample_name=TRUE)
patchwork::wrap_plots(plot_exposures(sbs1, sampleIDs = "UPN06_REL"),
                      plot_mutations(sbs1, sampleIDs="UPN06_REL"),
                      plot_signatures(sbs1), design="AABBBB
                                                      ##CCCC
                                                      ##CCCC")


sbs.fit1b = two_steps_inference(as.data.frame(counts_sbs),
                                k_list=0:5,
                                enforce_sparsity1=T,
                                enforce_sparsity2=F,
                                input_catalogue=COSMIC_filt_merged[c("SBS1","SBS40 SBS3 SBS5"),],
                                keep_sigs=c("SBS1","SBS40 SBS3 SBS5"),
                                py=py,
                                run_on_resid=FALSE,
                                regularizer="KL")
sbs1b = sbs.fit1b$tot
make_plots(sbs1b, what="SBS", epsilon=TRUE, sample_name=TRUE)
patchwork::wrap_plots(plot_exposures(sbs1, sampleIDs = "UPN06_REL"),
                      plot_mutations(sbs1, sampleIDs="UPN06_REL"),
                      plot_signatures(sbs1), design="AABBBB
                                                      ##CCCC
                                                      ##CCCC")


sbs.fit1c = two_steps_inference(as.data.frame(counts_sbs),
                                k_list=0:5,
                                enforce_sparsity1=T,
                                enforce_sparsity2=F,
                                input_catalogue=COSMIC_filt_merged[c("SBS1","SBS40 SBS3 SBS5"),],
                                keep_sigs=c("SBS1","SBS40 SBS3 SBS5"),
                                py=py,
                                run_on_resid=FALSE,
                                regularizer="cosine")
sbs1c = sbs.fit1c$tot
make_plots(sbs1c, what="SBS", epsilon=TRUE, sample_name=TRUE, sort_by="D1")
patchwork::wrap_plots(plot_exposures(sbs1c, sampleIDs = "UPN06_REL"),
                      plot_mutations(sbs1c, sampleIDs="UPN06_REL"),
                      plot_signatures(sbs1c), design="AABBBB
                                                      ##CCCC
                                                      ##CCCC")


## Only SBS1,5
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

## Second step on all counts
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

plot_sigs_prevalence(sbs3)
make_plots(sbs3, what="SBS", epsilon=TRUE, sample_name=T)
patchwork::wrap_plots(plot_exposures(sbs3, sampleIDs = "UPN06_REL"),
                      plot_mutations(sbs3, sampleIDs="UPN06_REL"),
                      plot_signatures(sbs3), design="AABBBB
                                                      ##CCCC
                                                      ##CCCC")


sbs.fit3b = two_steps_inference(as.data.frame(counts_sbs),
                                k_list=0:5,
                                enforce_sparsity1=T,
                                enforce_sparsity2=F,
                                input_catalogue=COSMIC_filt_merged[c("SBS1","SBS40 SBS3 SBS5"),],
                                keep_sigs=c("SBS1","SBS40 SBS3 SBS5"),
                                py=py,
                                run_on_resid=FALSE,
                                regularizer="cosine")
sbs3b = sbs.fit3b$tot

get_exposure(sbs3b, long=T) %>%
  dplyr::group_by(Signature) %>%
  dplyr::summarise(min_exp=min(Exposure), max_exp=max(Exposure))

sbs3b %>% plot_sigs_prevalence()
make_plots(sbs3b, what="SBS", epsilon=TRUE, sample_name=T)
patchwork::wrap_plots(plot_exposures(sbs3, sampleIDs = "UPN06_REL"),
                      plot_mutations(sbs3, sampleIDs="UPN06_REL"),
                      plot_signatures(sbs3), design="AABBBB
                                                      ##CCCC
                                                      ##CCCC")



x_recon = (get_data(sbs.fit2$step1_filt, reconstructed=FALSE) -
  get_data(sbs.fit2$step1_filt, reconstructed=TRUE)) %>% abs

sbs2.tmp = sbs2
sbs2.tmp$fit$x = compute_residuals(sbs.fit2$step1, min_exp = 0.2, keep_sigs = c("SBS1","SBS40 SBS3 SBS5"))
pl_recon = lapply(rownames(x_recon), function(sid)
  plot_mutations(sbs2.tmp, sampleIDs=sid) + labs(title=sid)) %>%
  patchwork::wrap_plots()

ggsave("~/Desktop/pl_recon2.pdf", plot=pl_recon, height=10, width=14)

pl_true = lapply(rownames(sbs2$fit$x), function(sid)
  plot_mutations(sbs2, sampleIDs=sid) + labs(title=sid)) %>%
  patchwork::wrap_plots()

ggsave("~/Desktop/pl_true.pdf", plot=pl_true, height=10, width=14)


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




## Comparison among llk #####

sbs1.5 = sbs.fit2$tot
sbs.dn = sbs.fit2b$tot

recon1 = get_data(sbs1.5, reconstructed = TRUE)
recon2 = get_data(sbs.dn, reconstructed = TRUE)

ll1 = sapply(1:ncol(counts_sbs), function(ctx)
  sapply(1:nrow(counts_sbs), function(sid)
    dpois(counts_sbs[sid, ctx],
          lambda=recon1[sid, ctx], log=TRUE)
  ))

ll2 = sapply(1:ncol(counts_sbs), function(ctx)
  sapply(1:nrow(counts_sbs), function(sid)
    dpois(counts_sbs[sid, ctx],
          lambda=recon2[sid, ctx], log=TRUE)
  ))

sum(-ll1)
sum(-ll2)

k1 = 22*96
k2 = 1*96 + 22*1


# for fixed in beta_fixed:
#   for denovo in beta_denovo:
#   loss += torch.log(F.kl_div(torch.log(fixed), torch.log(denovo), log_target = True, reduction="batchmean"))

denovo = sbs.dn$fit$denovo_signatures
reference = sbs.dn$fit$catalogue_signatures

denovo[denovo == 0] = 1e-07
reference[reference == 0] = 1e-07

kl = 0
for (i in 1:nrow(reference))
  kl = kl + philentropy::KL(as.matrix(rbind(log(denovo),
                                            log(reference[i,]))), unit="log")

-2*sum(ll1) + k1*log(22)


# _log_like += self.reg_weight * (reg * self.x.shape[0] * self.x.shape[1])
-2*(sum(ll2)+kl*22*96) + k2*log(22)




