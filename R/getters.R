
#' get exposure matrix
#'
#' @param x basilica object
#' @param long if TRUE return the long format exposure matrix (default=FALSE)
#'
#' @return a data.frame where rows are samples and columns are inferred signature profiles
#' @export get_exposure
#'
#' @examples
get_exposure <- function(x, long = FALSE) {

  stopifnot(inherits(x, "basilica_obj"))

  alpha <- x$fit$exposure

  if (long) {

    is_denovo = function(n){
      n %in% (x$fit$denovo_signatures %>% rownames())
    }

    alpha$Sample <- rownames(alpha)
    alpha <- tidyr::gather(alpha, key="Signature", value="Exposure", c(-Sample))
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
#'
#' @examples
get_catalogue_signatures <- function(x, long = FALSE) {
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
#'
#' @examples
get_denovo_signatures <- function(x,  long = FALSE) {
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


#' get de novo and catalouge signatures
#'
#' @param x basilica object
#'
#' @return a data.frame where rows are inferred signatures (not included in reference catalogue) and columns are 96 substitution bases.
#' @export get_denovo_signatures
#'
#' @examples
get_signatures <- function(x,  long = FALSE) {
  stopifnot(inherits(x, "basilica_obj"))

  sigs_dn = x %>% get_denovo_signatures(long = !!long)
  sigs_ct = x %>% get_catalogue_signatures(long = !!long)

  sigs = sigs_ct %>%
    dplyr::bind_rows(sigs_dn)

  return(sigs)
}

get_reference_signatures <- function(x,  long = FALSE) {
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

