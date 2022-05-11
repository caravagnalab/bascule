
# constructor
init_object <- function(
    fit,
    x,
    groups,
    input_catalog=NULL,
    reference_catalog=basilica::COSMIC,
    k,
    lr,
    steps,
    phi,
    delta
    ) {

  obj = list()

  obj$input <- list(x=x, groups=groups, input_catalog=input_catalog)
  obj$reference_catalog <- reference_catalog
  obj$params <- list(k=k, lr=lr, steps=steps, phi=phi, delta=delta)
  obj$fit <- list(exposure=fit[[1]], catalog_signatures=fit[[2]], denovo_signatures=fit[[3]])

  #structure(res, class="basilica")
  class(obj) <- "basilica"

  return(obj)
}


