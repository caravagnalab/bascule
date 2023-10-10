#' get exposure matrix
#'
#' @param x basilica object
#' @param long if TRUE return the long format exposure matrix (default=FALSE)
#'
#' @return a data.frame where rows are samples and columns are inferred signature profiles
#' @export get_exposure

get_exposure = function(x, long = FALSE, add_groups = FALSE) {

  stopifnot(inherits(x, "basilica_obj"))

  alpha = x$fit$exposure

  if (is.matrix(alpha)) {
    alpha = alpha %>% as.data.frame()
    colnames(alpha) = c(x$fit$catalogue_signatures %>% rownames(),
                        x$fit$denovo_signatures %>% rownames())
  }

  if (add_groups && have_groups(x))
    alpha$groups = x$groups

  if (long) {
    is_denovo = function(n){
      n %in% (x$fit$denovo_signatures %>% rownames())
    }

    alpha$Sample = rownames(alpha)
    if ("groups" %in% colnames(alpha))
      alpha = tidyr::gather(alpha, key="Signature", value="Exposure", c(-Sample, -groups)) else
        alpha = tidyr::gather(alpha, key="Signature", value="Exposure", c(-Sample))

    alpha = alpha %>%
      dplyr::mutate(Type = ifelse(
        is_denovo(Signature),
        "De novo",
        "Catalogue"
      )) %>%
      tidyr::as_tibble()

  }

  return(alpha)
}

#' get catalogue signatures
#'
#' @param x basilica object
#'
#' @return a data.frame where rows are inferred signatures (included in reference catalogue) and columns are 96 substitution bases.
#' @export get_catalogue_signatures

get_catalogue_signatures = function(x, long = FALSE) {
  stopifnot(inherits(x, "basilica_obj"))

  sigs = x$fit$catalogue_signatures

  if(long)
    sigs = reshape2::melt(sigs %>% as.matrix()) %>%
      dplyr::rename(
        Signature = Var1,
        Feature = Var2,
        Value = value
      ) %>%
    dplyr::as_tibble()%>%
    dplyr::mutate(Type = 'Catalogue')

  return(sigs)
}

#' get de novo signatures
#'
#' @param x basilica object
#'
#' @return a data.frame where rows are inferred signatures (not included in reference catalogue) and columns are 96 substitution bases.
#' @export get_denovo_signatures

get_denovo_signatures = function(x,  long = FALSE) {
  stopifnot(inherits(x, "basilica_obj"))

  sigs = x$fit$denovo_signatures

  if(is.null(sigs)) return(NULL)

  if(long)
    sigs = reshape2::melt(sigs %>% as.matrix()) %>%
    dplyr::rename(
      Signature = Var1,
      Feature = Var2,
      Value = value
    ) %>%
    dplyr::as_tibble() %>%
    dplyr::mutate(Type = 'De novo')

  return(sigs)
}


get_fixed_signatures = function(x,  wide=F) {
  stopifnot(inherits(x, "basilica_obj"))

  sigs = x$fit$input_catalogue

  if(is.null(sigs)) return(NULL)

  if(long)
    sigs = reshape2::melt(sigs %>% as.matrix()) %>%
    dplyr::rename(
      Signature = Var1,
      Feature = Var2,
      Value = value
    ) %>%
    dplyr::as_tibble() %>%
    dplyr::mutate(Type = 'De novo')

  return(sigs)
}


#' get de novo and catalouge signatures
#'
#' @param x basilica object
#'
#' @return a data.frame where rows are inferred signatures (not included in reference catalogue) and columns are 96 substitution bases.
#' @export get_denovo_signatures

get_signatures = function(x,  long = FALSE) {
  stopifnot(inherits(x, "basilica_obj"))

  sigs_dn = x %>% get_denovo_signatures(long = !!long)
  sigs_ct = x %>% get_catalogue_signatures(long = !!long)

  sigs = sigs_ct %>%
    dplyr::bind_rows(sigs_dn)

  return(sigs)
}


get_reference_signatures = function(x, long = FALSE) {
  stopifnot(inherits(x, "basilica_obj"))

  sigs = x$input$reference_catalogue

  if(long)
    sigs = reshape2::melt(sigs %>% as.matrix()) %>%
    dplyr::rename(
      Signature = Var1,
      Feature = Var2,
      Value = value
    ) %>%
    dplyr::as_tibble() %>%
    dplyr::mutate(Type = 'Reference')

  return(sigs)
}


have_groups = function(x) {
  if ("groups" %in% names(x) && !is.null(x[["groups"]])) return(TRUE)
  return(FALSE)
}


have_epsilon = function(x) {
  if ("eps_var" %in% names(x$fit)) return(TRUE)
  return(FALSE)
}


have_color_palette = function(x) {
  if ("color_palette" %in% names(x)) return(TRUE)
  return(FALSE)
}


get_color_palette = function(x) {
  if ("color_palette" %in% names(x)) return(x[["color_palette"]])
  return(NULL)
}


subset_catalogue = function(catalogue, rm_sigs=NULL, keep_sigs=NULL) {
  if (is.null(rm_sigs) && is.null(keep_sigs))
    return(catalogue)

  if (is.null(rm_sigs))
    return(catalogue[keep_sigs,])

  if (is.null(keep_sigs))
    return(catalogue[!rownames(catalogue) %in% rm_sigs, ])

  keep_fin = intersect(which(!rownames(catalogue) %in% rm_sigs), which(rownames(catalogue) %in% keep_sigs)) %>% unique()
  return(catalogue[keep_fin])
}


get_groups_with_sigs = function(x, sigs, thr=0) {
  if (!have_groups(x))
    return()

  return(
    get_exposure(x, add_groups=T, long=T) %>%
      dplyr::filter(Signature %in% sigs) %>%
      dplyr::group_by(groups, Signature) %>%
      dplyr::reframe(is_present=any(Exposure > 0)) %>%
      dplyr::ungroup() %>% dplyr::filter(is_present) %>%
      dplyr::select(groups, Signature) %>%
      tidyr::nest(data=groups)
  )
}


get_samples_with_sigs = function(x, sigs, thr=0, return_idx=FALSE) {
  if (!return_idx)
    return( rownames(get_exposure(x))[which(get_exposure(x)[,sigs] > thr)] )

  return( which(get_exposure(x)[,sigs] > thr) )
}


get_sigs_group = function(x, groupID, thr=0) {
  samples_g = rownames(get_group(x, groupID))
  expos = get_exposure(x)[samples_g,]
  expos[expos < thr] = 0
  return(
    colnames(expos)[colSums(expos) > 0]
  )
}


get_group = function(x, groupIDs, return_idx=FALSE) {
  if (!have_groups(x))
    return()
  if (!return_idx)
    return(
      x %>% get_data() %>% dplyr::mutate(groups=x$groups) %>%
        dplyr::filter(groups %in% groupIDs) %>%
        dplyr::select(-groups)
    )

  return(
    x %>% get_data() %>% dplyr::mutate(groups=x$groups) %>%
      dplyr::filter(groups %in% groupIDs) %>%
      dplyr::select(-groups) %>% rownames()
  )
}


add_groups = function(x, groups) {
  x$groups = groups
  return(x)
}

get_dn_signames = function(x) {
  return(rownames(get_denovo_signatures(x)))
}

get_fixed_signames = function(x) {
  return(rownames(get_fixed_signatures(x)))
}

get_catalogue_signames = function(x) {
  return(rownames(get_catalogue_signatures(x)))
}

get_signames = function(x) {
  return(rownames(get_signatures(x)))
}


get_K_scores = function(x) {
  return(x$fit$runs_K)
}



get_secondBest_run = function(x) {
  return(
    x$fit$secondBest %>%
      create_basilica_obj(input_catalogue=get_fixed_signatures(x),
                          reference_catalogue=get_reference_signatures(x),
                          cohort=x$cohort,
                          filtered_catalogue=TRUE)
  )
}


get_fit_by_id = function(x, idd) {
  new_fit = x$fit$all_fits[[idd]]

  return(
    new_fit %>%
      create_basilica_obj(input_catalogue=get_fixed_signatures(x),
                          reference_catalogue=get_reference_signatures(x),
                          cohort=x$cohort,
                          filtered_catalogue=TRUE)
  )
}


get_groups = function(x) {
  if (have_groups(x)) return(x$groups)
  return(rep(0, time=x$n_samples))
}


get_centroids = function(x, normalize=FALSE) {
  if (!have_groups(x)) return(NULL)

  centr = x$fit$params$alpha_prior

  if (any(grepl("G", get_groups(x)))) to_paste = "G" else to_paste = ""
  if (!any(grepl("G", rownames(centr)))) rownames(centr) = 1:nrow(centr) -1
  rownames(centr) = paste0(to_paste, rownames(centr) %>% stringr::str_replace_all("G",""))

  if (normalize) return(centr / rowSums(centr))
  return(centr)
}


get_mixture_weights = function(x) {
  if (!have_groups(x)) return(NULL)

  cl_names = rownames(get_centroids(x))
  return(x$fit$params$pi %>% setNames(cl_names))
}


COSMIC_color_palette = function(catalogue=COSMIC_filt, seed=14) {
  N = nrow(catalogue)
  set.seed(seed)
  colss = Polychrome::createPalette(N, c("#856de3","#9e461c"), target="normal", range=c(15,80), M=1000)[1:N]
  names(colss) = rownames(catalogue)
  return(colss)
}


get_contexts = function(x) {
  tryCatch(expr = {
    context_names = x %>% get_signatures(long=T) %>% dplyr::select(Feature) %>% unique()
  }, error = function(e)
    context_names = data.frame(Feature = x$denovo_signatures %>% colnames()) )

  return(context_names %>%
           dplyr::mutate(Feature=gsub(pattern = '\\[', replacement = '_', x=Feature)) %>%
           dplyr::mutate(Feature=gsub(pattern = '\\]', replacement = '_', x=Feature)) %>%
           tidyr::separate(Feature, into=c("left","subs","right"), sep="_"))
}


get_obj_initial_params = function(x) {
  params = x$fit$init_params
  # x$fit$exposure = x$fit$params$alpha = params$alpha / rowSums(params$alpha)
  x$fit$groups = x$groups = params$init_clusters
  x$fit$pi = x$fit$params$pi = params$pi_param

  x$fit$params$alpha_prior = params$alpha_prior_param

  return(x)
}


get_adjusted_fit_lc = function(x) {
  if ("lc_check" %in% names(x)) return(x$lc_check)
  return(x)
}

