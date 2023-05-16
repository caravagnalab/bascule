two_steps_inference = function(counts,
                               k_list,
                               enforce_sparsity1=TRUE,
                               enforce_sparsity2=FALSE,
                               min_exposure=0.2,
                               input_catalogue=COSMIC_filtered,
                               cohort="MyCohort") {
  x1 = pyfit(
    x = counts,
    py = py,
    lr = 0.05,
    n_steps = 500,
    k_list = 0,
    groups = NULL,
    input_catalogue = input_catalogue,
    regularizer = "KL",
    reg_weight = 1,
    reg_bic = TRUE,
    compile = FALSE,
    enforce_sparsity = enforce_sparsity1,
    stage = "random_noise"
  )
  xx1 = create_basilica_obj(x1, input_catalogue, cohort=cohort)

  resid_counts = compute_residuals(xx1, min_exp=min_exposure)
  x2 = pyfit(
    x = round(resid_counts),
    py = py,
    lr = 0.05,
    n_steps = 500,
    k_list = k_list,
    groups = NULL,
    input_catalogue = NULL,
    regularizer = "KL",
    reg_weight = 1,
    reg_bic = TRUE,
    compile = FALSE,
    enforce_sparsity = enforce_sparsity2,
    stage = "",
    regul_compare = xx1 %>% filter_exposures(min_exp=min_exposure) %>% get_signatures()
  )
  xx2 = create_basilica_obj(x2, input_catalogue=NULL, cohort=cohort)

  x_tot = xx1 %>% filter_exposures(min_exp=min_exposure)
  x_tot$n_denovo = xx2$n_denovo
  x_tot$fit$denovo_signatures = xx2$fit$denovo_signatures

  if (x_tot$n_denovo == 0)
    x_tot$fit$exposure = normalize_exposures(xx1 %>% filter_exposures(min_exp=min_exposure)) %>% get_exposure()
  else {
    resid_expos = 1 - rowSums(x_tot$fit$exposure)
    denovo_norm = xx2$fit$exposure / rowSums(xx2$fit$exposure) * resid_expos
    x_tot$fit$exposure = cbind(x_tot$fit$exposure, denovo_norm)
  }

  return(list("tot"=x_tot,
              "step1"=xx1,
              "step1_filt"=xx1 %>% filter_exposures(min_exp=min_exposure),
              "step2"=xx2))
}


create_basilica_obj_simul = function(simul, cohort="MySimul") {
  ss = list()
  class(ss) = "basilica_obj"
  ss$cohort = cohort
  ss$n_samples = nrow(simul$x[[1]])
  ss$n_denovo = nrow(simul$exp_denovo[[1]])
  ss$input = list("counts"=simul$x[[1]],
                  "reference_catalogue"=simul$ref_cat[[1]],
                  "input_catalogue"=simul$exp_fixed[[1]])

  ss$fit = list("input_catalogue"=simul$exp_fixed[[1]],
                "catalogue_signatures"=simul$ref_cat[[1]],
                "denovo_signatures"=simul$exp_denovo[[1]],
                "exposure"=simul$exp_exposure[[1]],
                "x"=simul$x[[1]])

  return(ss)
}


normalize_exposures = function(x) {
  x$fit$exposure = x$fit$exposure / rowSums(x$fit$exposure)
  return(x)
}


filter_exposures = function(x, min_exp=0.15) {
  sbs_keep = x$fit$exposure %>%
    tibble::rownames_to_column(var="sample") %>%
    reshape2::melt(id="sample", value.name="alpha", variable.name="sigs") %>%

    dplyr::mutate(alpha=ifelse(alpha < min_exp, 0, alpha)) %>%
    dplyr::filter(alpha > 0) %>%
    dplyr::pull(sigs) %>% unique() %>% as.character()

  x$fit$input_catalogue = x$fit$input_catalogue[sbs_keep,]
  x$fit$catalogue_signatures = x$fit$catalogue_signatures[sbs_keep,]
  x$fit$exposure = x$fit$exposure[,sbs_keep]

  x$input$input_catalogue = x$input$input_catalogue[sbs_keep,]
  x$input$reference_catalogue = x$input$reference_catalogue[sbs_keep,]

  return(x)
}


map_groups = function(groups) {
  n_gr = groups %>% unique() %>% length()
  grps.int = 1:n_gr %>% setNames(groups %>% unique())

  grps.int = grps.int-1

  return(grps.int[groups] %>% setNames(NULL))
}


create_basilica_obj = function(fit, input_catalogue, cohort="MyCohort") {
  obj = list()
  class(obj) = "basilica_obj"

  obj$cohort = cohort
  obj$n_samples = nrow(fit$x)

  if ("denovo_signatures" %in% names(fit))
    obj$n_denovo = nrow(fit$denovo_signatures) else
      obj$n_denovo = 0

  obj$input = list("counts"=fit$x,
                   "reference_catalogue"=COSMIC_filtered,
                   "input_catalogue"=fit$input_catalogue)

  fit$catalogue_signatures = fit$input_catalogue

  obj$fit = fit

  return(obj)
}


compute_residuals = function(x, min_exp=0.2) {
  xf = x %>% filter_exposures(min_exp=min_exp)
  sigs_order = rownames(xf$fit$input_catalogue)

  orig_counts = xf$input$counts
  t1 = (xf$fit$exposure * rowSums(orig_counts))[,sigs_order] %>% as.matrix()
  t2 = xf$fit$input_catalogue[sigs_order,] %>% as.matrix()
  reconstr_counts = t1 %*% t2

  residual_counts = abs(orig_counts - reconstr_counts)

  return(residual_counts)
}



gen_palette = function(n) {
  # library(RColorBrewer)
  # qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
  # col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
  # return(sample(col_vector, n))
  return(ggsci::pal_simpsons()(n))
}



