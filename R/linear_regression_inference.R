two_steps_inference = function(counts,
                               k_list,
                               enforce_sparsity1=TRUE,
                               enforce_sparsity2=FALSE,
                               min_exposure=0.2,
                               input_catalogue=COSMIC_filt_merged,
                               keep_sigs=c("SBS1", "SBS40 SBS3 SBS5"),
                               lr=0.05,
                               n_steps=500,
                               groups=NULL,
                               compile=FALSE,
                               cohort="MyCohort",
                               regularizer="KL",
                               run_on_resid=TRUE,
                               py=NULL) {

  xx1 = NULL; xx1_filt = NULL
  if (!is.null(input_catalogue)) {
    x1 = pyfit(
      x = counts,
      py = py,
      lr = lr,
      n_steps = n_steps,
      k_list = 0,
      groups = NULL,
      input_catalogue = input_catalogue,
      regularizer = regularizer,
      reg_weight = 1,
      reg_bic = TRUE,
      compile = compile,
      enforce_sparsity = enforce_sparsity1,
      stage = "random_noise"
    )
    xx1 = create_basilica_obj(x1, input_catalogue, cohort=cohort)
    xx1_filt = xx1 %>% filter_exposures(min_exp=min_exposure, keep_sigs=keep_sigs)

    catalogue2 = get_signatures(xx1_filt)
  } else catalogue2 = NULL

  resid_counts = counts; regul_compare = NULL

  if (run_on_resid) {
    resid_counts = compute_residuals(xx1,
                                     min_exp=min_exposure,
                                     keep_sigs=keep_sigs)
    regul_compare = catalogue2
    catalogue2 = NULL
  }

  x2 = pyfit(
    x = round(resid_counts),
    py = py,
    lr = lr,
    n_steps = n_steps,
    k_list = k_list,
    groups = NULL,
    input_catalogue = catalogue2,
    regularizer = regularizer,
    reg_weight = 1,
    reg_bic = TRUE,
    compile = compile,
    enforce_sparsity = enforce_sparsity2,
    stage = "",
    regul_compare = regul_compare
  )
  xx2 = create_basilica_obj(x2, input_catalogue=NULL, cohort=cohort)

  if (!run_on_resid)
    return(list("tot"=xx2,
                "step1"=xx1,
                "step1_filt"=xx1_filt,
                "step2"=xx2))

  x_tot = xx1 %>% filter_exposures(min_exp=min_exposure, keep_sigs=keep_sigs)
  x_tot$n_denovo = xx2$n_denovo
  x_tot$fit$denovo_signatures = NULL

  if (x_tot$n_denovo == 0)
    x_tot$fit$exposure = normalize_exposures(xx1_filt) %>% get_exposure()
  else {
    resid_expos = 1 - rowSums(x_tot$fit$exposure)
    denovo_norm = xx2$fit$exposure / rowSums(xx2$fit$exposure) * resid_expos

    # if (max(denovo_norm) > 0.1) {
    x_tot$fit$exposure = cbind(x_tot$fit$exposure, denovo_norm)
    x_tot$fit$denovo_signatures = xx2$fit$denovo_signatures
    # }
  }

  x_tot$color_palette = gen_palette(get_signatures(x_tot) %>% nrow()) %>%
    setNames(sort(rownames(get_signatures(x_tot))))

  return(list("tot"=x_tot,
              "step1"=xx1,
              "step1_filt"=xx1_filt,
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

  if (!is.null(simul$groups))
    ss$groups = simul$groups[[1]]

  ss$color_palette = gen_palette(get_signatures(ss) %>% nrow()) %>%
    setNames(sort(rownames(get_signatures(ss))))

  return(ss)
}


normalize_exposures = function(x) {
  x$fit$exposure = x$fit$exposure / rowSums(x$fit$exposure)
  return(x)
}


filter_exposures = function(x, min_exp=0.15, keep_sigs=NULL) {
  sbs_keep = x$fit$exposure %>%
    tibble::rownames_to_column(var="sample") %>%
    reshape2::melt(id="sample", value.name="alpha", variable.name="sigs") %>%

    dplyr::mutate(alpha=ifelse(alpha < min_exp, 0, alpha)) %>%
    dplyr::filter(alpha > 0) %>%
    dplyr::pull(sigs) %>% unique() %>% as.character()

  if (!is.null(keep_sigs)) sbs_keep = c(sbs_keep, keep_sigs) %>% unique()

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

  obj$color_palette = gen_palette(get_signatures(obj) %>% nrow()) %>%
    setNames(sort(rownames(get_signatures(obj))))

  return(obj)
}


compute_residuals = function(x, min_exp=0.2, keep_sigs=NULL) {
  xf = x %>% filter_exposures(min_exp=min_exp, keep_sigs=keep_sigs)
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



