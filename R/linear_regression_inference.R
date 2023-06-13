two_steps_inference = function(x,
                               k,
                               enforce_sparsity1=TRUE,
                               enforce_sparsity2=FALSE,
                               min_exposure=0.2,
                               reference_catalogue=COSMIC_filt_merged,
                               filtered_catalogue=TRUE,
                               keep_sigs=c("SBS1", "SBS5"),
                               lr=0.05,
                               n_steps=500,
                               groups=NULL,
                               compile=FALSE,
                               cohort="MyCohort",
                               regularizer="cosine",
                               residues=FALSE,
                               py=NULL,
                               reg_weight=0.,
                               reg_bic=TRUE,
                               CUDA=FALSE,
                               verbose=FALSE) {

  TIME = as.POSIXct(Sys.time(), format = "%H:%M:%S")

  if (!is.null(reference_catalogue)) {
    x_ref = pyfit(
      x = x,
      py = py,
      lr = lr,
      n_steps = n_steps,
      k_list = 0,
      groups = groups,
      input_catalogue = reference_catalogue,
      regularizer = regularizer,
      reg_weight = reg_weight,
      reg_bic = reg_bic,
      compile = compile,
      CUDA = CUDA,
      verbose = verbose,
      enforce_sparsity = enforce_sparsity1,
      stage = "random_noise") %>%

      create_basilica_obj(input_catalogue=reference_catalogue[keep_sigs, ],
                          reference_catalogue=reference_catalogue,
                          cohort=cohort,
                          filtered_catalogue=filtered_catalogue)

    x_ref_filt = x_ref %>% filter_exposures(min_exp=min_exposure, keep_sigs=keep_sigs)
    catalogue2 = get_signatures(x_ref_filt)
  } else {
    x_ref = x_ref_filt = catalogue2 = NULL
    residues = FALSE
  }

  if (residues) {
    resid_counts = compute_residuals(x_ref, min_exp=min_exposure, keep_sigs=keep_sigs)
    regul_compare = catalogue2
    catalogue2 = NULL
  } else {
    resid_counts = x
    regul_compare = NULL
  }

  x_dn = pyfit(
    x = round(resid_counts),
    py = py,
    lr = lr,
    n_steps = n_steps,
    k_list = k,
    groups = groups,
    input_catalogue = catalogue2,
    regularizer = regularizer,
    reg_weight = reg_weight,
    reg_bic = reg_bic,
    compile = compile,
    CUDA = CUDA,
    verbose = verbose,
    enforce_sparsity = enforce_sparsity2,
    stage = "",
    regul_compare = regul_compare
  ) %>% create_basilica_obj(input_catalogue=reference_catalogue[keep_sigs,],
                            reference_catalogue=reference_catalogue,
                            cohort=cohort,
                            filtered_catalogue=filtered_catalogue)

  TIME = difftime(as.POSIXct(Sys.time(), format = "%H:%M:%S"), TIME, units = "mins")

  # if (!residues) merged = x_dn
  merged = merge_fits(x_ref, x_dn, x_ref_filt, min_exposure, keep_sigs, residues)
  merged$time = TIME
  merged$k_list = k

  return(list("tot"=merged, "step1"=x_ref, "step1_filt"=x_ref_filt, "step2"=x_dn))
}


merge_fits = function(x1, x2, x1_filt, min_exposure, keep_sigs, residues) {
  if (!residues)
    return(x2)

  merged = x1 %>% filter_exposures(min_exp=min_exposure, keep_sigs=keep_sigs)
  merged$n_denovo = x2$n_denovo
  merged$fit$denovo_signatures = NULL

  if (merged$n_denovo == 0)
    merged$fit$exposure = normalize_exposures(x1_filt) %>% get_exposure()
  else {
    resid_expos = 1 - rowSums(merged$fit$exposure)
    denovo_norm = x2$fit$exposure / rowSums(x2$fit$exposure) * resid_expos

    merged$fit$exposure = cbind(merged$fit$exposure, denovo_norm)
    merged$fit$denovo_signatures = x2$fit$denovo_signatures
  }

  merged$color_palette = gen_palette(get_signatures(merged) %>% nrow()) %>%
    setNames(sort(rownames(get_signatures(merged))))

  return(merged)
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

  ss$private_sigs = list()
  if ("private_rare" %in% colnames(simul))
    ss$private_sigs[["private_rare"]] = simul$private_rare[[1]]

  if ("private_common" %in% colnames(simul))
    ss$private_sigs[["private_common"]] = simul$private_common[[1]]

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
  if (length(sbs_keep) == 0) return(x)

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


create_basilica_obj = function(fit, input_catalogue, reference_catalogue, cohort="MyCohort", filtered_catalogue=TRUE) {
  # fit is the output of "pyfit"
  obj = list()
  class(obj) = "basilica_obj"

  obj$cohort = cohort
  obj$n_samples = nrow(fit$x)

  if ("denovo_signatures" %in% names(fit))
    obj$n_denovo = nrow(fit$denovo_signatures) else
      obj$n_denovo = 0

  if (filtered_catalogue && obj$n_denovo > 0)
    fit$denovo_signatures = renormalize_denovo_thr(fit$denovo_signatures)

  obj$input = list("counts"=fit$x,
                   "reference_catalogue"=reference_catalogue,
                   "input_catalogue"=fit$input_catalogue)

  fit$catalogue_signatures = fit$input_catalogue

  obj$fit = fit
  obj$groups = fit$groups

  obj$color_palette = gen_palette(get_signatures(obj) %>% nrow()) %>%
    setNames(sort(rownames(get_signatures(obj))))

  obj$time = fit$time

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



