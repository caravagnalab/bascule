
#' get exposure matrix
#'
#' @param x basilica object
#' @param long if TRUE return the long format exposure matrix (default=FALSE)
#'
#' @return a data.frame where rows are samples and columns are inferred signature profiles
#' @export get_exposure
#'
#' @examples
get_exposure <- function(x, long=FALSE) {
  #stopifnot(inherits(x, "basilica"))

  alpha <- x$fit$exposure

  if (long) {
    alpha$samples <- rownames(alpha)
    alpha <- tidyr::gather(alpha, key="Signature", value="Exposure", c(-samples))
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
get_catalogue_signatures <- function(x) {
  #stopifnot(inherits(x, "basilica"))
  return(x$fit$catalogue_signatures)
}

#' get de novo signatures
#'
#' @param x basilica object
#'
#' @return a data.frame where rows are inferred signatures (not included in reference catalogue) and columns are 96 substitution bases.
#' @export get_denovo_signatures
#'
#' @examples
get_denovo_signatures <- function(x) {
  #stopifnot(inherits(x, "basilica"))
  return(x$fit$denovo_signatures)
}

