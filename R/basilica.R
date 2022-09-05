
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
    input_catalogue=NULL,
    lambda_rate=NULL,
    sigma=FALSE
    ) {

  # quality check --------------------------------------------------------------
  if (!identical(colnames(x), colnames(reference_catalogue))) {
    reference_catalogue = reference_catalogue[names(x)]
  }

  if (!is.null(input_catalogue)) {
    if (!identical(colnames(x), colnames(input_catalogue))) {
      input_catalogue = input_catalogue[names(x)]
    }
  }

  counter <- 1
  black_list <- c()
  while (TRUE) {

    cat("    iteration:", counter, '\n') # TEST

    obj <- basilica:::pyfit(
      x=x,
      k_list=k,
      lr=lr,
      n_steps=steps,
      groups=groups,
      input_catalogue=input_catalogue,
      lambda_rate=lambda_rate,
      sigma=sigma
      )

    # drop non-significant fixed signatures ------------------------------------
    a <- basilica:::filter.fixed(
      M = x,
      alpha = obj$exposure,
      beta_fixed = input_catalogue,
      phi = phi
    )
    remained_fixed <- a$remained_fixed                          # data.frame / NULL

    # TEST-----
    print('a$dropped_fixed')
    print(class(a$dropped_fixed))
    print(a$dropped_fixed)
    print('----------')
    print('black_list')
    print(class(black_list))
    print(black_list)
    # TEST-----

    if (!is.null(a$dropped_fixed)) {
      black_list <- union(black_list, rownames(a$dropped_fixed))  # character vector
    }

    # detect denovo signatures which are similar to reference signatures -------
    b <- basilica:::filter.denovo(
      reference_catalogue=reference_catalogue,
      beta_fixed=input_catalogue,
      beta_denovo=obj$denovo_signatures,
      black_list=black_list,
      delta=delta
    )
    new_fixed <- b$new_fixed  # data.frame / NULL (with labels from reference)
    reduced_denovo <- b$reduced_denovo  # data.frame / NULL (remaining denovo signatures)

    #TEST---------------------------------------------------------
    cat("        fixed          :", rownames(input_catalogue), '\n')
    cat("        remained fixed :", rownames(remained_fixed), '\n')
    cat("        black list     :", black_list, "\n\n")
    cat("        denovo         :", rownames(obj$denovo_signatures), '\n')
    cat("        new fixed      :", rownames(new_fixed), '\n')
    cat("        reduced denovo :", rownames(reduced_denovo), '\n')
    #TEST---------------------------------------------------------


    if (is.null(input_catalogue)) {
      col_names <- colnames(x)
      input_catalogue = data.frame(matrix(nrow=0, ncol = length(col_names)))
      colnames(input_catalogue) = col_names
    }

    if (is.null(new_fixed)) {
      col_names <- colnames(x)
      new_fixed = data.frame(matrix(nrow=0, ncol = length(col_names)))
      colnames(new_fixed) = col_names
    }

    if (is.null(remained_fixed)) {
      col_names <- colnames(x)
      remained_fixed = data.frame(matrix(nrow=0, ncol = length(col_names)))
      colnames(remained_fixed) = col_names
    }

    if (nrow(dplyr::setdiff(input_catalogue, remained_fixed))==0 & nrow(new_fixed)==0) {
      cat('        break loop\n')
      break
    }

    if (nrow(remained_fixed)==0 & nrow(new_fixed)==0) {
      input_catalogue <- NULL
    } else {
      input_catalogue <- rbind(remained_fixed, new_fixed)
    }

    counter <- counter + 1
    if (counter > 8) {
      cat('        limit reached! : 5\n')
      break
    }
    cat('    ------------------------------------\n')
  }

  if (nrow(input_catalogue)==0) {
    obj$catalogue_signatures <- NULL
  } else {
    obj$catalogue_signatures <- input_catalogue
  }

  # output ---> dtype: list
  #-------------------------------------:
  # exposure              --> data.frame
  # denovo_signatures     --> data.frame
  # bic                   --> numeric
  # losses                --> numeric
  # catalogue_signatures  --> data.frame

  return(obj)
}


