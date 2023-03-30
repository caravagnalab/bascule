if (Sys.getenv("GITHUB_PATH") == "") path=paste0("~/dati_elenab/signatures/") else path=Sys.getenv("GITHUB_PATH")


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

# keep = sample(1:nrow(counts.gel), 2000)
# counts.gel.sub = counts.gel[keep,]
# groups.sub = counts.gel.sub$organ

counts.gel.sub = counts.gel %>% dplyr::filter(organ == "Colorectal")
keep = sample(1:nrow(counts.gel.sub), 100)
counts.gel.sub = counts.gel.sub[keep,]
groups.sub = counts.gel.sub$organ

table(groups)

py_path = paste0(path, "pybasilica")
py = reticulate::import_from_path(module="pybasilica", path=py_path)

reference_crc = COSMIC_catalogue[c("SBS1","SBS5","SBS6","SBS10a","SBS10b"),]

obj.nogroups = fit(counts.gel.sub %>% dplyr::select(-organ, -cohort), py=py, k=1:5,
                   cohort="GEL_crc", input_catalogue=COSMIC_catalogue[c("SBS1","SBS5"),],
                   reference_catalogue=COSMIC_catalogue, cosine_by_subs=FALSE)

obj.nogroups2 = fit(counts.gel.sub %>% dplyr::select(-organ, -cohort), py=py, k=1:5,
                   cohort="GEL_crc", input_catalogue=COSMIC_catalogue[c("SBS1","SBS5"),],
                   reference_catalogue=reference_crc, cosine_by_subs=TRUE, delta=.85)

# obj.groups = fit(counts.gel.sub %>% dplyr::select(-organ, -cohort), py=py, k=0:7,
#                  cohort="GEL_crc", groups=map_groups(groups.sub), reference_catalogue=reference_crc)
# obj.cluster = fit(counts.gel.sub %>% dplyr::select(-organ, -cohort), py=py, k=3:7, cohort="GEL", cluster=T)

saveRDS(obj.nogroups, "./nobuild/test_groups/gel.std.fit.Rds")
# saveRDS(obj.groups, "./nobuild/test_groups/gel.hier.fit.Rds")
# saveRDS(counts.gel.sub %>% dplyr::mutate(groups=groups.sub), "./nobuild/test_groups/input_data.Rds")
#
#
# obj.nogroups = readRDS("./nobuild/test_groups/gel.std.fit.Rds")
# obj.groups = readRDS("./nobuild/test_groups/gel.hier.fit.Rds")


# sbs1 = plot_signatures(obj.groups)
plot_signatures(x)
ggplot2::ggsave("./nobuild/test_groups/sigs.nogroups.pdf", height=8, width=8)
# patchwork::wrap_plots(sbs1, sbs2)
#
plot_exposure(x, sort_by = "D2")
ggplot2::ggsave("./nobuild/test_groups/exposure.nogroups.pdf", height=8, width=8)


w = plot_similarity_reference(x, by_subs = F)
pdf("./test.nogroups.pdf", width = 44, height = 28)
print(w)
dev.off()

w2 = plot_similarity_reference(obj.nogroups2, by_subs = T)
pdf("./test.nogroups2.pdf", width = 44, height = 28)
print(w2)
dev.off()



