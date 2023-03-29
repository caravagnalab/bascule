
map_groups = function(groups) {
  n_gr = groups %>% unique() %>% length()
  grps.int = 1:n_gr %>% setNames(groups %>% unique())

  grps.int = grps.int-1

  return(grps.int[groups] %>% setNames(NULL))
}

devtools::load_all()

reticulate::use_condaenv("basilica-env")

header = read.csv("./nobuild/test_groups/counts_all.tsv", header=F, sep="\t", nrows=1) %>% setNames(NULL)
counts = read.csv("./nobuild/test_groups/counts_all.tsv", header=T, row.names=1, sep="\t")
colnames(counts) = header
counts.gel = counts %>% dplyr::filter(cohort=="GEL") %>%
  dplyr::filter(organ%in%c("Colorectal","Lung","Breast"))
groups = counts.gel$organ

keep = sample(1:nrow(counts.gel), 2000)
counts.gel.sub = counts.gel[keep,]
groups.sub = counts.gel.sub$organ

table(groups)

py = reticulate::import_from_path(module="pybasilica", path="~/dati_elenab/signatures/pybasilicah/")

obj.nogroups = fit(counts.gel.sub %>% dplyr::select(-organ, -cohort), py=py, k=3:7, cohort="GEL")
obj.groups = fit(counts.gel.sub %>% dplyr::select(-organ, -cohort), py=py, k=3:7, cohort="GEL", groups=map_groups(groups.sub))
# obj.cluster = fit(counts.gel.sub %>% dplyr::select(-organ, -cohort), py=py, k=3:7, cohort="GEL", cluster=T)

saveRDS(obj.nogroups, "./nobuild/test_groups/gel.std.fit.Rds")
saveRDS(obj.groups, "./nobuild/test_groups/gel.hier.fit.Rds")
saveRDS(counts.gel.sub %>% dplyr::mutate(groups=groups.sub), "./nobuild/test_groups/input_data.Rds")


obj.nogroups = readRDS("./nobuild/test_groups/gel.std.fit.Rds")
obj.groups = readRDS("./nobuild/test_groups/gel.hier.fit.Rds")


sbs1 = plot_signatures(obj.groups)
sbs2 = plot_signatures(obj.nogroups, Type="De novo")
# patchwork::wrap_plots(sbs1, sbs2)
#
plot_exposure(obj.groups)

w = plot_similarity_reference(obj.groups)
pdf("./test.pdf", width = 44, height = 28)
print(w)
dev.off()


w = plot_similarity_reference(obj.nogroups)
pdf("./test.nogroups.pdf", width = 44, height = 28)
print(w)
dev.off()



