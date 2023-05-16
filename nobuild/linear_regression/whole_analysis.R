
x1.sparse = pyfit(
  x = counts,
  py = py,
  lr = 0.05,
  n_steps = 500,
  k_list = 0,
  groups = NULL,
  input_catalogue = COSMIC_filtered,
  regularizer = TRUE,
  reg_weight = 1,
  reg_bic = TRUE,
  compile = FALSE,
  enforce_sparsity = TRUE,
  stage = "random_noise"
)

xx1.sparse = create_basilica_obj(x1.sparse)
hhp = plot_similarity_reference(xx1.sparse %>% filter_exposures(min_exp=0.15))
ggsave("~/Desktop/heatmap.pdf", plot=hhp, height = 20, width = 20)
xx1.sparse %>% filter_exposures() %>% plot_signatures()
xx1.sparse %>% filter_exposures() %>% plot_exposures()



x1 = pyfit(
  x = counts,
  py = py,
  lr = 0.05,
  n_steps = 500,
  k_list = 0,
  groups = NULL,
  input_catalogue = COSMIC_filtered,
  regularizer = TRUE,
  reg_weight = 1,
  reg_bic = TRUE,
  compile = FALSE,
  enforce_sparsity = FALSE,
  stage = "random_noise"
)

xx1 = create_basilica_obj(x1)
hhp = plot_similarity_reference(xx1 %>% filter_exposures(), by_subs = F)
xx1 %>% filter_exposures() %>% plot_signatures()
xx1 %>% filter_exposures() %>% plot_exposures()




resid_counts = compute_residuals(xx1)

x2 = pyfit(
  x = round(resid_counts),
  py = py,
  lr = 0.05,
  n_steps = 500,
  k_list = 2:7,
  groups = NULL,
  input_catalogue = NULL,
  regularizer = TRUE,
  reg_weight = 1,
  reg_bic = TRUE,
  compile = FALSE,
  enforce_sparsity = FALSE,
  stage = ""
)


xx2 = create_basilica_obj(x2, NULL, n_denovo=x2$denovo_signatures %>% rownames())
xx2 %>% plot_exposures()
xx2 %>% plot_signatures()
xx2 %>% plot_similarity_reference(reference = COSMIC_filtered)
