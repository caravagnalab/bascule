
#' get exposure matrix
#'
#' @param x basilica object
#' @param long if TRUE return the long format exposure matrix (default=FALSE)
#'
#' @return
#' @export
#'
#' @examples
get_exposure <- function(x, long=FALSE) {
  stopifnot(inherits(x, "basilica"))

  alpha <- x$fit$exposure

  if (long) {
    alpha$Sample <- rownames(alpha)
    alpha <- tidyr::gather(alpha, key="Signature", value="Exposure", c(-Sample))
  }

  return(alpha)
}

#' get catalogue signatures
#'
#' @param x basilica object
#'
#' @return inferred signatures which are from reference catalogue
#' @export
#'
#' @examples
get_catalogue_signatures <- function(x) {
  stopifnot(inherits(x, "basilica"))
  return(x$fit$catalogue_signatures)
}

#' get de novo signatures
#'
#' @param x basilica object
#'
#' @return inferred signatures which are not included in reference catalogue
#' @export
#'
#' @examples
get_denovo_signatures <- function(x) {
  stopifnot(inherits(x, "basilica"))
  return(x$fit$denovo_signatures)
}

