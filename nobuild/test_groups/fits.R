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
# counts.gel = counts %>% dplyr::filter(cohort=="GEL") %>%
#   dplyr::filter(organ%in%c("Colorectal","Lung","Breast"))
# groups = counts.gel$organ

# keep = sample(1:nrow(counts.gel), 2000)
# counts.gel.sub = counts.gel[keep,]
# groups.sub = counts.gel.sub$organ

sample_patients = function(counts, N=100, organ_list=c("Colorectal"), cohort_list=c("GEL","ICGC","Hartwig"), seed=5) {
  # organ = c(organ)
  # cohort = c(cohort)
  counts.sub = counts %>% dplyr::filter(organ %in% organ_list, cohort %in% cohort_list)
  set.seed(seed)
  keep = sample(1:nrow(counts.sub), N)
  counts.sub = counts.sub[keep,]
  groups.sub = map_groups(counts.sub$organ)
  return(list("counts"=counts.sub %>% dplyr::select(-cohort, -organ),
              "groups_idx"=groups.sub,
              "organ"=counts.sub$organ,
              "cohort"=counts.sub$cohort))
}

## First try with only one organ -> CRC
counts_groups = sample_patients(counts, cohort_list="GEL", organ_list=c("Colorectal","Lung","Breast"))
table(counts_groups$organ)

py_path = paste0(path, "pybasilica")
py = reticulate::import_from_path(module="pybasilica", path=py_path)
# py = NULL

reference_sub = COSMIC_filtered[c("SBS1","SBS5","SBS6","SBS10a","SBS10b","SBS4",),]

obj.nogroups = fit(counts.gel.sub %>% dplyr::select(-organ, -cohort), py=py, k=0:15,
                   cohort="GEL_crc", input_catalogue=COSMIC_filtered[c("SBS1","SBS5"),],
                   reference_catalogue=reference_crc, reg_weight=1., reg_bic=FALSE)

plot_signatures(obj.nogroups)
plot_exposure(obj.nogroups)
plot_similarity_reference(obj.nogroups, reference = reference_crc)

obj.groups = fit(counts.gel.sub %>% dplyr::select(-organ, -cohort), py=py,
                 k=0:15, groups=groups.sub, cohort="GEL_crc_lung",
                 input_catalogue=COSMIC_filtered[c("SBS1","SBS5"),],
                 reference_catalogue=reference_crc,
                 reg_weight=1., reg_bic=FALSE)

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
plot_signatures(obj.nogroups)
ggplot2::ggsave("./nobuild/test_groups/sigs.nogroups.pdf", height=8, width=8)
# patchwork::wrap_plots(sbs1, sbs2)
#
plot_exposure(obj.nogroups, sort_by = "D4")
ggplot2::ggsave("./nobuild/test_groups/exposure.nogroups.pdf", height=8, width=8)


w = plot_similarity_reference(obj.nogroups, by_subs = F)
pdf("./test.nogroups.pdf", width = 44, height = 28)
print(w)
dev.off()

w2 = plot_similarity_reference(obj.nogroups2, by_subs = T)
pdf("./test.nogroups2.pdf", width = 44, height = 28)
print(w2)
dev.off()



