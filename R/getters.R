
get_exposure <- function(x) {
  stopifnot(inherits(x, "basilica"))
  return(x$exposure)
}

get_catalog_signatures <- function(x) {
  stopifnot(inherits(x, "basilica"))
  return(x$catalog_signatures)
}

get_denovo_signatures <- function(x) {
  stopifnot(inherits(x, "basilica"))
  return(x$denovo_signatures)
}
