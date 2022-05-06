
# constructor
init_object <- function(x) {
  #stopifnot("data.frame" %in% class(x))

  obj <- list(
    exposure=x[[1]],
    catalog_signatures=x[[2]],
    denovo_signatures=x[[3]]
    )

  #structure(res, class="basilica")
  class(obj) <- "basilica"

  return(obj)
}


