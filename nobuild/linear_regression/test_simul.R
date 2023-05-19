library(lsa)
library(dplyr)
library(ggplot2)
source("~/GitHub/simbasilica/nobuild/script_test/helper_fns.R")
py = reticulate::import_from_path("pybasilica","~/GitHub/pybasilica/")
devtools::load_all("~/GitHub/basilica/")
devtools::load_all("~/GitHub/simbasilica/")

path = "~/GitHub/basilica/nobuild/linear_regression/"

make_plots = function(x, x.true=NULL, reconstructed=T, cls=NULL,
                      sort_by=NULL, what="SBS", epsilon=FALSE,
                      sample_name=FALSE) {

  mm = plot_mutations(x, reconstructed=reconstructed, what=what, epsilon=epsilon)
  alp = plot_exposures(x, sort_by=sort_by, cls=cls, sample_name=sample_name)
  bet = plot_signatures(x, cls=cls, what=what)

  if (is.null(x.true)) return(patchwork::wrap_plots(bet + (mm/alp), guides="collect"))

  mm.true = plot_mutations(x.true, reconstructed=F, what=what, epsilon=epsilon)
  alpha.true = plot_exposures(x.true, sort_by=sort_by, cls=cls, sample_name=sample_name)

  return(patchwork::wrap_plots(bet + (mm/mm.true/alp/alpha.true), guides="collect"))
}

## Generate Data ####
cosine_limit = .8
shared = c("SBS1", "SBS40 SBS3 SBS5")
private_common = c("SBS17b", "SBS2", "SBS20")
private_rare = c("SBS8 SBS4")

denovo_catalogue = COSMIC_filt_merged[c(private_common, private_rare),]
reference_catalogue = COSMIC_filt_merged[shared,]
denovo_cosine = lsa::cosine(denovo_catalogue %>% t())
reference_cosine = lsa::cosine(reference_catalogue %>% t())
x = single_dataset(350, 2, 150:350, reference_catalogue, denovo_catalogue,
                   reference_cosine, denovo_cosine, mut_range=10:5000,
                   private_sigs=list("rare"=private_rare,"common"=private_common),
                   private_fracs=list("rare"=0.05,"common"=0.3),
                   cosine_limit, seed=23, cohort_name="test1",
                   out_path=path)

x.simul = readRDS(paste0(path, "simul.N350.G2.s23.test1.Rds"))
xx = create_basilica_obj_simul(x.simul)


## Inference on simulated data ####
x.fit1 = two_steps_inference(x.simul$x[[1]],
                             k=0:7,
                             input_catalogue=COSMIC_filt_merged,
                             enforce_sparsity1=TRUE,
                             enforce_sparsity2=FALSE,
                             regularizer="KL",
                             py=py)
x1 = x.fit1$tot
saveRDS(x1, paste0(path, "fit1.N350.G2.s23.test1.Rds"))

plot_similarity_reference(x1,
                          reference=COSMIC_filt_merged[rownames(get_signatures(xx)),],
                          similarity="KL",
                          similarity_cutoff = 0)


x.fit2 = two_steps_inference(x.simul$x[[1]],
                             k=0:7,
                             input_catalogue=COSMIC_filt_merged[rownames(COSMIC_filt_merged)!="SBS20",],
                             enforce_sparsity1=TRUE,
                             enforce_sparsity2=FALSE,
                             py=py)
x2 = x.fit2$tot
saveRDS(x2, paste0(path, "fit2.N350.G2.s23.test1.Rds"))


x.fit3 = two_steps_inference(x.simul$x[[1]],
                             k=0:7,
                             input_catalogue=COSMIC_filt_merged[
                               !rownames(COSMIC_filt_merged)%in%c("SBS8 SBS4"),],
                             enforce_sparsity1=TRUE,
                             enforce_sparsity2=FALSE,
                             py=py)
x3 = x.fit3$tot
saveRDS(x3, paste0(path, "fit3.N350.G2.s23.test1.Rds"))


sigs = c(rownames(get_signatures(xx)),
         rownames(get_signatures(x1)),
         rownames(get_signatures(x2)),
         rownames(get_signatures(x3))) %>% unique()
cls = gen_palette(length(sigs)) %>% setNames(sigs)

true.plots = make_plots(xx, reconstructed=F, cls=cls) & theme(legend.position="bottom")
x1.plots = make_plots(x1, xx, reconstructed=T, cls=cls) & theme(legend.position="bottom")
x2.plots = make_plots(x2, xx, reconstructed=T, cls=cls) & theme(legend.position="bottom")
x3.plots = make_plots(x3, xx, reconstructed=T, cls=cls) & theme(legend.position="bottom")


pdf(paste0(path, "report.simul_tests.pdf"), height = 10, width = 14)
true.plots & patchwork::plot_annotation(title="True data")
x1.plots & patchwork::plot_annotation(title="Fit 1 - whole catalogue")
x2.plots & patchwork::plot_annotation(title="Fit 2 - removed SBS20")
x3.plots & patchwork::plot_annotation(title="Fit 3 - removed SBS4 SBS8")

plot_similarity_reference(x1, reference=COSMIC_filt_merged[rownames(get_signatures(xx)),])
plot_similarity_reference(x2, reference=COSMIC_filt_merged[rownames(get_signatures(xx)),])
plot_similarity_reference(x3, reference=COSMIC_filt_merged[rownames(get_signatures(xx)),])

# patchwork::wrap_plots(plot_signatures(x1, signatures="SBS20", cls=cls) /
#                         plot_signatures(x2,
#                                         signatures=setdiff(rownames(get_signatures(x1)),
#                                                            rownames(get_signatures(x2))),
#                                         cls=cls))
dev.off()


# pdf("./heatmap.pdf", height = 20, width = 10)
# print(plot_similarity_reference(x.new2, reference=COSMIC_filt_merged))
# dev.off()


cosine_limit = .8
n_fixed = 2

shared = c("SBS1", "SBS40 SBS3 SBS5")
private_common = c("SBS8 SBS4", "SBS2", "SBS20")
private_rare = c("SBS17b")

reference_catalogue = COSMIC_filt_merged[shared,]
denovo_catalogue = COSMIC_filt_merged[c(private_common, private_rare),]

denovo_cosine = lsa::cosine(denovo_catalogue %>% t())
reference_cosine = lsa::cosine(reference_catalogue %>% t())

x.simul2 = single_dataset(500, 2, 150:350, reference_catalogue, denovo_catalogue,
                          reference_cosine, denovo_cosine,
                          private_sigs=list("rare"=private_rare,"common"=private_common),
                          private_fracs=list("rare"=0.05,"common"=0.3),
                          cosine_limit, seed=23, out_path="./nobuild/script_test/simulations/")
xx2 = create_basilica_obj_simul(x.simul2)
x.fit2 = two_steps_inference(x.simul2$x[[1]],
                             k=0:7,
                             input_catalogue=COSMIC_filt_merged,
                             enforce_sparsity1=TRUE,
                             enforce_sparsity2=FALSE,
                             py=py)
x.new2 = x.fit2$tot
saveRDS(x.new2, "~/GitHub/simbasilica/nobuild/test_new1.Rds")




x.new2 %>% plot_exposures()
x.new2 %>% plot_signatures()
x.new2 %>% plot_similarity_reference()
x.new2 %>% plot_mutations()



