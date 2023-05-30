sbs.names = c("SBS1","SBS5",
              "SBS2","SBS6","SBS10a","SBS10b","SBS21",  # CRC
              "SBS2","SBS3","SBS4","SBS8","SBS13","SBS15","SBS17a", "SBS17b", "SBS18",  # Lung
              "SBS3","SBS8","SBS10","SBS13","SBS17","SBS18","SBS20","SBS26")

cosine = lsa::cosine(t(COSMIC_filt_merged))
cosine %>% as.data.frame() %>%
  tibble::rownames_to_column(var="sbs1") %>%
  reshape2::melt(variable.name="sbs2", value.name="cosine") %>%
  dplyr::mutate(cosine = ifelse(cosine < 0.7, 0, cosine)) %>%
  dplyr::filter(sbs1=="SBS11") %>%
  tidyr::pivot_wider(names_from="sbs1", values_from="cosine", values_fill=1) %>%
  tibble::column_to_rownames(var="sbs2") %>% as.matrix() %>% heatmap()


cosine %>% as.data.frame() %>%
  tibble::rownames_to_column(var="sbs1") %>%
  reshape2::melt(variable.name="sbs2", value.name="cosine") %>%
  dplyr::filter(cosine < 0.7)



# private = paste0("SBS", c(13,"10b",14,20,21,42,44,52))
private = select_fixed_sbs(COSMIC_filt_merged[!rownames(COSMIC_filt_merged) %in%
                                                c("SBS1","SBS5","SBS17b",
                                                  "SBS87", "SBS88", private.final),],
                           cosine_limit = 0.3, n_fixed = 5) %>% rownames()

cosine_priv = lsa::cosine(t(COSMIC_filt_merged[private.final,]))
cosine_priv %>% as.data.frame() %>%
  tibble::rownames_to_column(var="sbs1") %>%
  reshape2::melt(variable.name="sbs2", value.name="cosine") %>%
  tidyr::pivot_wider(names_from="sbs1", values_from="cosine", values_fill=1) %>%
  tibble::column_to_rownames(var="sbs2") %>% as.matrix() %>% pheatmap(legend = T)

COSMIC_filt_merged[private.final,] %>%
  tibble::rownames_to_column("sbs") %>%
  reshape2::melt() %>%
  reformat_contexts(what="SBS") %>%
  ggplot() +
  geom_histogram(aes(x=context, y=value), stat="identity") +
  facet_grid(Var1 ~ substitution, scales="free_y")


private.final = c("SBS10b", "SBS28", "SBS56 SBS10a", "SBS90", "SBS2", "SBS13", "SBS20", "SBS22")



