#' creates an object of class \code(basilica}.
#'
#' @param x input mutational counts data (data.frame; rows as samples and columns as 96 mutational categories)
#' @param groups vector of discrete labels with one entry per sample, it defines the groups that will be considered by basilica
#' @param input_catalog input signature profiles, NULL by default
#' @param reference_catalog a catalog of reference signatures that basilica will use to compare input and de novo signatures (COSMIC catalogue by default)
#' @param k vector of possible number of de novo signatures to infer
#' @param lr stochastic variational inference learning rate
#' @param steps number of gradient steps
#' @param phi threshold to discard the signature based on its value in exposure matrix
#' @param delta threshold to consider inferred signature as COSMIC signature
#'
#' @return inferred exposure matrix, inferred signatures from reference catalogue and inferred de novo (not from reference catalogue) signatures
#'
#' @export
#'
#' @examples
fit <- function(
    x,
    groups=NULL,
    input_catalogue=NULL,
    reference_catalogue=basilica::COSMIC_catalogue,
    k=0:5,
    lr=0.05,
    steps=500,
    phi=0.05,
    delta=0.9
    ) {

  pybasilica <- reticulate::import("pybasilica")
  f <- pybasilica$pyfit(x, groups, input_catalogue, reference_catalogue, k, lr, steps, phi, delta)

  obj <- init_object(
    f,
    x,
    groups,
    input_catalogue,
    reference_catalogue,
    k,
    lr,
    steps,
    phi,
    delta
  )

  return(obj)
}




