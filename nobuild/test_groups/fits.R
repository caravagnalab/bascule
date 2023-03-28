reticulate::use_condaenv("basilica-env")

header = read.csv("./nobuild/test_groups/counts_all.tsv", header=F, sep="\t", nrows=1) %>% setNames(NULL)
counts = read.csv("./nobuild/test_groups/counts_all.tsv", header=T, row.names=1, sep="\t")
colnames(counts) = header
counts.gel = counts %>% dplyr::filter(cohort=="GEL") %>%
  dplyr::filter(organ%in%c("Colorectal","Lung","Breast"))
keep = sample(1:nrow(counts.gel), 100)

counts.gel.sub = counts.gel[keep,]

groups = counts.gel.sub$organ
table(groups)

py = reticulate::import_from_path(module="pybasilica", path="~/GitHub/pybasilicah/")

obj.nogroups = fit(counts.gel.sub %>% dplyr::select(-organ, -cohort), py=py, k=c(3, 5, 7), cohort="GEL")
obj.groups = fit(counts.gel.sub %>% dplyr::select(-organ, -cohort), py=py, k=c(3, 5, 7), cohort="GEL", groups=map_groups(groups))


sbs1 = plot_signatures(obj.nogroups, Type="De novo")
sbs2 = plot_signatures(obj.groups, Type="De novo")
patchwork::wrap_plots(sbs1, sbs2)

plot_exposure(obj.groups)



map_groups = function(groups) {
  n_gr = groups %>% unique() %>% length()
  grps.int = 1:n_gr %>% setNames(groups %>% unique())

  grps.int = grps.int-1

  return(grps.int[groups] %>% setNames(NULL))
}
