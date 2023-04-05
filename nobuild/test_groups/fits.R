if (Sys.getenv("GITHUB_PATH") == "") path=paste0("~/dati_elenab/signatures/") else path=Sys.getenv("GITHUB_PATH")


map_groups = function(groups) {
  n_gr = groups %>% unique() %>% length()
  grps.int = 1:n_gr %>% setNames(groups %>% unique())

  grps.int = grps.int-1

  return(grps.int[groups] %>% setNames(NULL))
}

compute_cosine_denovo = function(denovo_list, by_subs=F) {
  denovo.all = lapply(names(denovo_list), function(dd) {
    rownames(denovo_list[[dd]])=paste0(dd,".",rownames(denovo_list[[dd]]))
    return(denovo_list[[dd]])
  })

  denovo.all = do.call(rbind, denovo.all)

  substitutions = NULL
  if (by_subs)
    substitutions = denovo_list[[1]] %>% colnames() %>% lapply(function(i) strsplit(i, "\\[|\\]")[[1]][[2]]) %>% unlist() %>% unique()

  return(cosine.matrix(denovo.all, denovo.all, substitutions = substitutions))
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
counts_groups = sample_patients(counts, N=200, cohort_list="GEL", organ_list=c("Colorectal","Lung"), seed=4)
table(counts_groups$organ)

py_path = paste0(path, "pybasilica")
py = reticulate::import_from_path(module="pybasilica", path=py_path)
# py = NULL

sbs.names = c("SBS1","SBS5",
              "SBS2","SBS6","SBS10a","SBS10b","SBS21","SBS28","SBS44",  # CRC
              "SBS2","SBS3","SBS4","SBS8","SBS13","SBS15","SBS17","SBS18",  # Lung
              "SBS3","SBS8","SBS10","SBS13","SBS17","SBS18","SBS20","SBS26")  # Breast

reference_sub = COSMIC_filtered[intersect(sbs.names, rownames(COSMIC_filtered)),]

# sbs1.5 = (COSMIC_filtered["SBS1",] + COSMIC_filtered["SBS5",]) / rowSums((COSMIC_filtered["SBS1",] + COSMIC_filtered["SBS5",]))
# rownames(sbs1.5) = "SBS1.5"
# reference_sub2 = rbind(reference_sub[!(rownames(reference_sub)%in%c("SBS1","SBS5")),], sbs1.5)

x.nogrps = fit(x = counts_groups$counts, py=py, k=0:7,
               cohort="GEL.crc_lung",
               input_catalogue=reference_sub[c("SBS1","SBS5"),],
               # input_catalogue=reference_sub["SBS1.5",],
               reference_catalogue=reference_sub, delta=0.85, phi=0.1,
               reg_weight=1., reg_bic=TRUE, filtered_cat = TRUE)

saveRDS(x.nogrps, "./nobuild/test_fit.Rds")

x.grps = fit(x=counts_groups$counts, py=py, k=0:7,
             groups=counts_groups$groups_idx,
             cohort="GEL.crc_lung",
             input_catalogue=reference_sub[c("SBS1","SBS5"),],
             # input_catalogue=reference_sub["SBS1.5",],
             reference_catalogue=reference_sub, delta=0.85, phi=0.1,
             reg_weight=1., reg_bic=TRUE, filtered_cat = TRUE)

saveRDS(x.nogrps, "./nobuild/test_fit.Rds")
saveRDS(x.grps, "./nobuild/test_fit.hier.Rds")

# x.nogrps.noreg = fit(counts_groups$counts, py=py, k=0:7,
#                      cohort="GEL.crc_lung",
#                      input_catalogue=reference_sub[c("SBS1","SBS5"),],
#                      # input_catalogue=reference_sub["SBS1.5",],
#                      reference_catalogue=reference_sub, delta=0.85, phi=0.1,
#                      reg_weight=0., reg_bic=FALSE, filtered_cat = TRUE)
#
# x.nogrps.noregbic = fit(counts_groups$counts, py=py, k=0:7,
#                       cohort="GEL.crc_lung",
#                       input_catalogue=reference_sub[c("SBS1","SBS5"),],
#                       # input_catalogue=reference_sub["SBS1.5",],
#                       reference_catalogue=reference_sub, delta=0.85, phi=0.1,
#                       reg_weight=1., reg_bic=FALSE, filtered_cat=TRUE)

pl1 = plot_similarity_reference(x.nogrps, reference = COSMIC_filtered)
ggplot2::ggsave("./nobuild/prova.png", height = 18, width = 20, bg = "white")

pl1.grps = plot_similarity_reference(x.grps, reference = COSMIC_filtered)
ggplot2::ggsave("./nobuild/prova.hier.png", height = 18, width = 20, bg = "white")

pl2 = plot_similarity_reference(x.nogrps.noreg, reference = COSMIC_filtered)
pl3 = plot_similarity_reference(x.nogrps.noregbic, reference = COSMIC_filtered)

pdf("./nobuild/compare.pdf", height = 18, width = 20)
print(pl1)
print(pl2)
print(pl3)
dev.off()

sid = c("GEL-2498766-11", "GEL-2736989-11")
muts = counts_groups$counts[sid,]
plots_muts = function(muts) {
  return(
    muts %>%
      as.data.frame() %>%
      tibble::rownames_to_column(var="pid") %>%
      reshape2::melt(id="pid") %>%
      dplyr::rename(Var1=pid, Var2=variable) %>%
      dplyr::mutate(substitution=paste0(substr(start=3, stop=3, Var2),">",substr(start=5, stop=5, Var2)),
                    context=paste0(substr(start=1, stop=1, Var2), "_", substr(start=7, stop=7, Var2))) %>%
      ggplot() +
      geom_bar(aes(value, x=context, fill=Var1), stat="identity") +
      facet_grid(Var1 ~ substitution, scales="free") +
      my_ggplot_theme() +
      theme(axis.text.x = element_text(angle=90)) +
      # theme(axis.ticks.x=element_blank(),axis.text.x=element_blank()) +
      # scale_fill_manual(values=cls) +
      guides(fill="none")  +
      labs(x="Context", y="", title="Number of muts")
  )
}



denovo_list = list("reg"=x.nogrps$fit$denovo_signatures,
                   "noreg"=x.nogrps.noreg$fit$denovo_signatures,
                   "regbic"=x.nogrps.noregbic$fit$denovo_signatures)

cosine.denovo = compute_cosine_denovo(denovo_list, by_subs = F)
numbers = cosine.denovo %>% round(2)
numbers[numbers < 0.5] = ''
pheatmap(cosine.denovo, cluster_rows = F, cluster_cols = F, display_numbers = numbers)


# run with groups
x.grps = fit(counts_groups$counts, py=py, groups=counts_groups$groups_idx,
             k=0:10, cohort="GEL.crc_lung.Hier",
             input_catalogue=COSMIC_filtered[c("SBS1","SBS5"),],
             reference_catalogue=reference_sub,
             reg_weight=1., reg_bic=TRUE, filtered_cat = TRUE)



pdf("~/GitHub/basilica/nobuild/test_groups/test1.pdf", height = 8, width = 12)
plot_signatures(x.nogrps)
plot_signatures(x.grps)

plot_exposure(x.nogrps, sort_by="D5")
plot_exposure(x.grps)

plot_similarity_reference(x.nogrps, by_subs = T)
plot_similarity_reference(x.grps, reference = reference_sub)
dev.off()


intersect(sbs.names, rownames(x.nogrps$fit$catalogue_signatures))
intersect(sbs.names, rownames(x.grps$fit$catalogue_signatures))





# plot_signatures(obj.nogroups)
# ggplot2::ggsave("./nobuild/test_groups/sigs.nogroups.pdf", height=8, width=8)
#
# plot_exposure(obj.nogroups, sort_by = "D4")
# ggplot2::ggsave("./nobuild/test_groups/exposure.nogroups.pdf", height=8, width=8)
#
#
# w = plot_similarity_reference(obj.nogroups, by_subs = F)
# pdf("./test.nogroups.pdf", width = 44, height = 28)
# print(w)
# dev.off()
#
# w2 = plot_similarity_reference(obj.nogroups2, by_subs = T)
# pdf("./test.nogroups2.pdf", width = 44, height = 28)
# print(w2)
# dev.off()



