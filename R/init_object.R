
# constructor
init_object <- function(
    fit,
    x,
    groups,
    input_catalogue=NULL,
    reference_catalogue=basilica::COSMIC,
    k,
    lr,
    steps,
    phi,
    delta
    ) {

  obj = list()

  obj$input <- list(x=x, groups=groups, input_catalogue=input_catalogue)
  obj$reference_catalogue <- reference_catalogue
  obj$params <- list(k=k, lr=lr, steps=steps, phi=phi, delta=delta)
  obj$fit <- list(exposure=fit[[1]], catalogue_signatures=fit[[2]], denovo_signatures=fit[[3]])

  #structure(res, class="basilica")
  class(obj) <- "basilica"

  return(obj)
}


