#' Fit function
#'
#' @param x mutations catalogue
#' @param k number or list of values of denovo signatures to test
#' @param enforce_sparsity if TRUE, Beta prior on exposures centroids
#' @param min_exposure filters exposure after the first fit with the reference catalogue
#' @param reference_catalogue signatures to look for in the data.
#' After one fit of the model only those signatures with a strong signal will be kept (see "min_exposure").
#' @param filtered_catalogue whether "reference_catalogue" has been filtered from low signals
#' @param keep_sigs signatures in "reference_catalogue" to keep despite their signals
#' @param lr learning rare
#' @param n_steps number of VI steps
#' @param clusters number or list of values to test for clustering
#' @param nonparametric whether to run the nonparametric clustering
#' @param dirichlet_prior whether to use Dirichlet prior on exposures' centroids (instead of a Normal)
#' @param compile add
#' @param cohort cohort name
#' @param regularizer type of regularizer for the signatures (only considered if "reg_weight" is > 0)
#' @param py python package to use (if NULL, the package will be installed from PyPI)
#' @param hyperparameters list of hyperparameters
#' @param reg_weight weight of the regularization
#' @param CUDA add
#' @param verbose add
#' @param seed_list list of seeds to use for distinct runs
#' @param regul_denovo add
#' @param regul_fixed add
#' @param store_fits add
#' @param do_initial_fit add
#'
#' @return basilica object
#'
#' @importFrom magrittr  %>%
#' @importFrom dplyr mutate group_by select rename as_tibble bind_rows reframe
#' @importFrom dplyr ungroup filter summarise pull setdiff add_row contains
#' @importFrom dplyr arrange full_join inner_join rowwise any_vars filter_all
#'
#' @importFrom tidyr separate gather as_tibble nest pivot_wider pivot_longer
#'
#' @importFrom tibble rownames_to_column column_to_rownames as_tibble
#'
#' @export

fit = function(x,
               k,
               enforce_sparsity=FALSE,
               min_exposure=0.2,
               reference_catalogue=COSMIC_filt_merged,
               filtered_catalogue=TRUE,
               keep_sigs=c("SBS1", "SBS5"),
               lr=0.005,
               optim_gamma=0.1,
               n_steps=2000,
               # groups=NULL,
               clusters=NULL,
               hyperparameters=NULL,
               nonparametric=TRUE,
               dirichlet_prior=TRUE,
               py=NULL,
               reg_weight=0.,
               regularizer="cosine",
               regul_denovo=TRUE,
               regul_fixed=TRUE,
               stage="",
               seed_list=c(10,27,92),
               # initializ_seed=FALSE,
               # save_runs_seed=TRUE,
               # initializ_pars_fit=TRUE,
               store_parameters=FALSE,
               store_fits=FALSE,
               # do_initial_fit=FALSE,
               compile=FALSE,
               CUDA=FALSE,
               # verbose=FALSE,
               cohort="MyCohort") {

  TIME = as.POSIXct(Sys.time(), format = "%H:%M:%S")
  call_info = match.call()

  if (!is.null(reference_catalogue)) {
    x_ref = pyfit(
      counts = x,
      py = py,
      lr = lr,
      optim_gamma = optim_gamma,
      n_steps = n_steps,
      k_list = 0,
      hyperparameters = hyperparameters,
      # groups = groups,
      clusters = NULL,
      dirichlet_prior = FALSE,
      input_catalogue = reference_catalogue,
      regularizer = regularizer,
      reg_weight = reg_weight,
      compile = compile,
      CUDA = CUDA,
      # verbose = verbose,
      enforce_sparsity = TRUE,
      stage = "random_noise",
      seed_list = seed_list,
      # initializ_seed = initializ_seed,
      # save_runs_seed = save_runs_seed,
      # initializ_pars_fit = initializ_pars_fit,
      store_parameters = FALSE,
      regul_denovo = regul_denovo,
      regul_fixed = regul_fixed,
      # do_initial_fit = FALSE
      ) %>%

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
    optim_gamma = optim_gamma,
    n_steps = n_steps,
    k_list = k,
    hyperparameters = hyperparameters,
    # groups = groups,
    clusters = clusters,
    nonparametric = nonparametric,
    dirichlet_prior = dirichlet_prior,
    input_catalogue = catalogue2,
    regularizer = regularizer,
    reg_weight = reg_weight,
    compile = compile,
    CUDA = CUDA,
    # verbose = verbose,
    enforce_sparsity = enforce_sparsity,
    stage = stage,
    regul_compare = NULL,
    seed_list = seed_list,
    # initializ_seed = initializ_seed,
    # save_runs_seed = save_runs_seed,
    # initializ_pars_fit = initializ_pars_fit,
    store_parameters = store_parameters,
    regul_denovo = regul_denovo,
    regul_fixed = regul_fixed,
    store_fits = store_fits
    # do_initial_fit = do_initial_fit
    ) %>%

    create_basilica_obj(input_catalogue=reference_catalogue[keep_sigs,],
                        reference_catalogue=reference_catalogue,
                        cohort=cohort,
                        filtered_catalogue=filtered_catalogue)

  TIME = difftime(as.POSIXct(Sys.time(), format = "%H:%M:%S"), TIME, units = "mins")

  merged = x_dn
  merged$time = TIME
  merged$k_list = k
  merged$call = call_info

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


create_basilica_obj = function(fit, input_catalogue, reference_catalogue,
                               cohort="MyCohort", filtered_catalogue=TRUE) {
  # fit is the output of "pyfit"
  obj = list()
  class(obj) = "basilica_obj"

  obj$cohort = cohort
  obj$n_samples = nrow(fit$x)

  # if ("denovo_signatures" %in% names(fit))
  #   obj$n_denovo = nrow(fit$denovo_signatures) else
  #     obj$n_denovo = 0

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



