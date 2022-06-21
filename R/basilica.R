
#' fit basilica model.
#'
#' @description fit the model and infer the underlying signatures and their contributions in the mutational catalogue counts.
#'
#' @param x input mutational counts data (data.frame; rows as samples and columns as 96 mutational categories)
#' @param reference_catalogue a catalog of reference signatures that basilica will use to compare input and de novo signatures (COSMIC catalogue by default)
#' @param k vector of possible number of de novo signatures to infer
#' @param lr stochastic variational inference learning rate
#' @param steps number of gradient steps
#' @param phi threshold to discard the signature based on its value in exposure matrix
#' @param delta threshold to consider inferred signature as COSMIC signature
#' @param groups vector of discrete labels with one entry per sample, it defines the groups that will be considered by basilica
#' @param input_catalogue input signature profiles, NULL by default
#'
#' @return inferred exposure matrix, inferred signatures from reference catalogue and inferred de novo (not from reference catalogue) signatures
#' @export fit
#'
#' @examples
fit <- function(
    x,
    reference_catalogue,
    k,
    lr,
    steps,
    phi,
    delta,
    groups=NULL,
    input_catalogue=NULL
    ) {

  counter <- 1
  while (TRUE) {

    obj <- basilica:::pyfit(
      x=x,
      k_list=k,
      lr=lr,
      n_steps=steps,
      groups=groups,
      input_catalogue=input_catalogue
      )

    if (is.null(input_catalogue)) {
      col_names <- colnames(x)
      input_catalogue = data.frame(matrix(nrow=0, ncol = length(col_names)))
      colnames(input_catalogue) = col_names
    }

    a <- basilica:::filter_fixed(x, obj$exposure, input_catalogue, phi)

    b <- basilica:::filter_denovo(obj$denovo_signatures, reference_catalogue, delta)

    if (nrow(a)==nrow(input_catalogue) & nrow(b)==0) {
      break
    }

    if (nrow(a)==0 & nrow(b)==0) {
      input_catalogue <- NULL
    } else {
      input_catalogue <- rbind(a, b)
    }

    counter <- counter + 1
  }

  if (nrow(input_catalogue)==0) {
    obj$catalogue_signatures <- NULL
  } else {
    obj$catalogue_signatures <- input_catalogue
  }
  #obj$reference_catalogue <- reference_catalogue
  #obj$phi <- phi
  #obj$delta <- delta

  # output
  #-------------------------------------:
  # exposure              --> data.frame
  # denovo_signatures     --> data.frame
  # bic                   --> numeric
  # losses                --> numeric
  # catalogue_signatures  --> data.frame

  return(obj)
}













