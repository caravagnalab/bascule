fit = function(x,
               k,
               enforce_sparsity=FALSE,
               min_exposure=0.2,
               reference_catalogue=COSMIC_filt_merged,
               filtered_catalogue=TRUE,
               keep_sigs=c("SBS1", "SBS5"),
               lr=0.05,
               n_steps=2000,
               groups=NULL,
               clusters=NULL,
               nonparametric=FALSE,
               compile=FALSE,
               cohort="MyCohort",
               regularizer="cosine",
               py=NULL,
               hyperparameters=NULL,
               reg_weight=1.,
               reg_bic=TRUE,
               CUDA=FALSE,
               verbose=FALSE,
               seed_list=c(10,27,92),
               initializ_seed=FALSE,
               save_runs_seed=TRUE,
               initializ_pars_fit=TRUE,
               new_hier=FALSE,
               regul_denovo=TRUE,
               regul_fixed=TRUE,
               save_all_fits=FALSE,
               do_initial_fit=FALSE) {

  TIME = as.POSIXct(Sys.time(), format = "%H:%M:%S")

  if (!is.null(reference_catalogue)) {
    x_ref = pyfit(
      counts = x,
      py = py,
      lr = lr,
      n_steps = n_steps,
      k_list = 0,
      hyperparameters = hyperparameters,
      groups = groups,
      clusters = NULL,
      input_catalogue = reference_catalogue,
      regularizer = regularizer,
      reg_weight = reg_weight,
      reg_bic = reg_bic,
      compile = compile,
      CUDA = CUDA,
      verbose = verbose,
      enforce_sparsity = TRUE,
      stage = "random_noise",
      seed_list = seed_list,
      initializ_seed = initializ_seed,
      save_runs_seed = save_runs_seed,
      initializ_pars_fit = initializ_pars_fit,
      new_hier = new_hier,
      regul_denovo = regul_denovo,
      regul_fixed = regul_fixed,
      do_initial_fit = FALSE) %>%

      create_basilica_obj(input_catalogue=reference_catalogue[keep_sigs, ],
                          reference_catalogue=reference_catalogue,
                          cohort=cohort,
                          filtered_catalogue=filtered_catalogue)

    x_ref_filt = x_ref %>% filter_sigs_low_expos(min_exp=min_exposure, keep_sigs=keep_sigs)
    catalogue2 = get_signatures(x_ref_filt)
  } else {
    x_ref = x_ref_filt = catalogue2 = NULL
    residues = FALSE
  }

  x_dn = pyfit(
    counts = x,
    py = py,
    lr = lr,
    n_steps = n_steps,
    k_list = k,
    hyperparameters = hyperparameters,
    groups = groups,
    clusters = clusters,
    nonparametric = nonparametric,
    input_catalogue = catalogue2,
    regularizer = regularizer,
    reg_weight = reg_weight,
    reg_bic = reg_bic,
    compile = compile,
    CUDA = CUDA,
    verbose = verbose,
    enforce_sparsity = enforce_sparsity,
    stage = "",
    regul_compare = NULL,
    seed_list = seed_list,
    initializ_seed = initializ_seed,
    save_runs_seed = save_runs_seed,
    initializ_pars_fit = initializ_pars_fit,
    new_hier = new_hier,
    regul_denovo = regul_denovo,
    regul_fixed = regul_fixed,
    save_all_fits = save_all_fits,
    do_initial_fit = do_initial_fit) %>%

    create_basilica_obj(input_catalogue=reference_catalogue[keep_sigs,],
                        reference_catalogue=reference_catalogue,
                        cohort=cohort,
                        filtered_catalogue=filtered_catalogue)

  TIME = difftime(as.POSIXct(Sys.time(), format = "%H:%M:%S"), TIME, units = "mins")

  merged = x_dn
  merged$time = TIME
  merged$k_list = k

  return(merged)
}


create_basilica_obj_simul = function(simul, cohort="MySimul") {
  ss = list()
  class(ss) = "basilica_obj"
  ss$cohort = cohort
  ss$n_samples = nrow(simul$counts[[1]])
  ss$n_denovo = nrow(simul$alpha[[1]])
  ss$input = list("counts"=simul$counts[[1]],
                  "reference_catalogue"=simul$beta[[1]],
                  "input_catalogue"=simul$beta[[1]])

  ss$fit = list("input_catalogue"=simul$beta[[1]],
                "catalogue_signatures"=simul$beta[[1]],
                "denovo_signatures"=NULL,
                "exposure"=simul$alpha[[1]],
                "x"=simul$counts[[1]])

  centr = simul$alpha_prior[[1]] %>% as.data.frame()
  rownames(centr) = paste0("G",rownames(centr))
  ss$fit$params = list("alpha_prior"=centr)

  ss$groups = simul$groups[[1]]
  if (is.null(ss$groups)) {
    ss$groups = simul$counts[[1]] %>% tibble::rownames_to_column(var="sample") %>%
      tidyr::separate(sample, into=c("gid","sample"), sep="_") %>%
      dplyr::pull(gid) %>% stringr::str_replace_all("G","")
  }
  ss$groups = paste0("G", ss$groups)

  ss$color_palette = gen_palette(get_signatures(ss) %>% nrow()) %>%
    setNames(sort(rownames(get_signatures(ss))))

  if ("private" %in% colnames(simul))
    ss$sigs[["private"]] = simul$private[[1]]

  if ("private_shared" %in% colnames(simul))
    ss$sigs[["private_shared"]] = simul$private_shared[[1]]

  if ("shared" %in% colnames(simul))
    ss$sigs[["shared"]] = simul$shared[[1]]

  return(ss)
}


filter_sigs_low_expos = function(x, min_exp=0.15, keep_sigs=NULL) {
  sbs_keep = get_exposure(x, long=TRUE) %>%
    # x$fit$exposure %>%
    # tibble::rownames_to_column(var="sample") %>%
    # reshape2::melt(id="sample", value.name="alpha", variable.name="sigs") %>%

    dplyr::mutate(Exposure=ifelse(Exposure < min_exp, 0, Exposure)) %>%
    dplyr::filter(Exposure > 0) %>%
    dplyr::pull(Signature) %>% unique() %>% as.character()

  if (!is.null(keep_sigs)) sbs_keep = c(sbs_keep, keep_sigs) %>% unique()
  if (length(sbs_keep) == 0) return(x)

  x$fit$input_catalogue =
    x$fit$input_catalogue[intersect(sbs_keep, get_fixed_signames(x)),]

  x$fit$catalogue_signatures =
    x$fit$catalogue_signatures[intersect(sbs_keep, get_catalogue_signames(x)),]

  x$fit$denovo_signatures =
    x$fit$denovo_signatures[intersect(sbs_keep, get_dn_signames(x)),]

  x$fit$exposure = x$fit$exposure[,sbs_keep]

  return(x)
}


# map_groups = function(groups) {
#   n_gr = groups %>% unique() %>% length()
#   grps.int = 1:n_gr %>% setNames(groups %>% unique())
#
#   grps.int = grps.int-1
#
#   return(grps.int[groups] %>% setNames(NULL))
# }


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


gen_palette = function(n) {
  return(ggsci::pal_simpsons()(n))
}



