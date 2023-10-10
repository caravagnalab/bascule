get_fixed_signames = function(x, type="T1") {
  return(get_fixed_signatures(x, type=type)$sigs %>% unique())
}


get_denovo_signames = function(x, type="T1") {
  return(get_denovo_signatures(x, type=type)$sigs %>% unique())
}


get_fixed_signatures = function(x, type=NULL, wide=F) {
  stopifnot(inherits(x, "basilica_obj"))
  if (!wide) return(x$fit$fixed_signatures)

  return(x$fit$fixed_signatures %>% long_to_wide_beta(type=type))
}


get_denovo_signatures = function(x, type=NULL, wide=F) {
  stopifnot(inherits(x, "basilica_obj"))
  if (!wide) return(x$fit$denovo_signatures %>% dplyr::filter(type==type))

  return(x$fit$denovo_signatures %>% long_to_wide_beta(type=type))
}


wide_to_long_beta = function(dataframe, type=NULL) {
  stopifnot(!is.null(type))
  return(
    dataframe %>% tibble::rownames_to_column(var="sigs") %>%
      reshape2::melt(id="sigs", variable.name="features") %>%
      dplyr::mutate(type=type) %>% tibble::as_tibble()
  )
}


long_to_wide_beta = function(dataframe, type) {
  dataframe %>%
    dplyr::filter(type==type) %>% dplyr::select(-type) %>%
    tidyr::pivot_wider(names_from="features", values_from="value") %>%
    tibble::column_to_rownames(var="sigs")
}



wide_to_long_alpha = function(dataframe, type=NULL) {
  stopifnot(!is.null(type))
  return(
    dataframe %>% tibble::rownames_to_column(var="samples") %>%
      reshape2::melt(id="samples", variable.name="sigs") %>%
      dplyr::mutate(type=type) %>% tibble::as_tibble()
  )
}


wide_to_long_counts = function(dataframe, type=NULL) {
  stopifnot(!is.null(type))
  return(
    dataframe %>% tibble::rownames_to_column(var="samples") %>%
      reshape2::melt(id="samples", variable.name="features") %>%
      dplyr::mutate(type=type) %>% tibble::as_tibble()
  )
}


