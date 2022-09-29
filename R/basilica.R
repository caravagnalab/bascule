

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
<<<<<<< HEAD
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
=======
fit <- function(x,
                reference_catalogue,
                k,
                lr = 0.05,
                steps = 500,
                phi,
                delta = 0.9,
                groups = NULL,
                input_catalogue = NULL,
                lambda_rate = NULL,
                sigma = FALSE)
{
  # cli::cli_h1("MUSICA - MUtational Signature Inference with a CAtalogue ")
  cli::cli_h1("Basilica - Bayesian signature learning with a catalogue")
  cat("\n")
>>>>>>> 3cb40b3527aabd136063f345b48acd13176212ed

  # First, sanitize inputs
  sanitized_inputs = sanitize_inputs(
    x = x,
    reference_catalogue = reference_catalogue,
    k = k,
    lr = lr,
    steps = steps,
    phi = phi,
    delta = delta,
    groups = NULL,
    input_catalogue = input_catalogue,
    lambda_rate = lambda_rate,
    sigma = sigma
  )

  x = sanitized_inputs$x
  reference_catalogue = sanitized_inputs$reference_catalogue
  input_catalogue = sanitized_inputs$input_catalogue

  # Report messages for the inputs
  cli::cli_alert_success(
    "Counts from {.field n = {nrow(x)}} samples."
  )

  # x %>% dplyr::as_tibble() %>% print()

  # Report messages for the reference
  cli::cli_alert_success(
    "Reference catalogue with {.field {nrow(reference_catalogue)}} signatures: {.field {head(rownames(reference_catalogue))}, ...}"
  )
  # reference_catalogue %>% dplyr::as_tibble() %>% print()

  # Report messages for the input
  if(!is.null(input_catalogue))
  {
    cli::cli_alert_success(
      "Detected input catalogue with {.field {nrow(reference_catalogue)}} signatures: {.field {head(rownames(reference_catalogue))}, ...}"
    )

    input_catalogue %>% dplyr::as_tibble() %>% print()
  }
  else
    cli::cli_alert_warning(
      "No input signatures given, detecting them from data"
    )

  cli::cli_h2("Iterative inference algorithm")

  counter <- 1
  black_list <- c()
  while (TRUE) {
    # cat("    iteration:", counter, '\n') # TEST

    cli::cli_h2("Bayesian NMF [step {.field {counter}}]")

    obj <- basilica:::pyfit(
      x = x,
      k_list = k,
      lr = lr,
      n_steps = steps,
      groups = groups,
      input_catalogue = input_catalogue,
      lambda_rate = lambda_rate,
      sigma = sigma
    )

    cli::cli_h2("Signatures matching")


    # drop non-significant fixed signatures ------------------------------------
    a <- basilica:::filter.fixed(
      M = x,
      alpha = obj$exposure,
      beta_fixed = input_catalogue,
      phi = phi
    )
    remained_fixed <-
      a$remained_fixed                          # data.frame / NULL

    # TEST-----
    cat('class(a$dropped_fixed):', class(a$dropped_fixed), '\n')
    cat('a$dropped_fixed:')
    print(a$dropped_fixed)
    cat('class(black_list):', class(black_list), '\n')
    cat('black_list:', black_list, '\n')
    # TEST-----

    if (!is.null(a$dropped_fixed)) {
      black_list <-
        union(black_list, rownames(a$dropped_fixed))  # character vector
    }

    # detect denovo signatures which are similar to reference signatures -------
    b <- basilica:::filter.denovo(
      reference_catalogue = reference_catalogue,
      beta_fixed = input_catalogue,
      beta_denovo = obj$denovo_signatures,
      black_list = black_list,
      delta = delta
    )
    new_fixed <-
      b$new_fixed  # data.frame / NULL (with labels from reference)
    reduced_denovo <-
      b$reduced_denovo  # data.frame / NULL (remaining denovo signatures)

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
      input_catalogue = data.frame(matrix(nrow = 0, ncol = length(col_names)))
      colnames(input_catalogue) = col_names
    }

    if (is.null(new_fixed)) {
      col_names <- colnames(x)
      new_fixed = data.frame(matrix(nrow = 0, ncol = length(col_names)))
      colnames(new_fixed) = col_names
    }

    if (is.null(remained_fixed)) {
      col_names <- colnames(x)
      remained_fixed = data.frame(matrix(nrow = 0, ncol = length(col_names)))
      colnames(remained_fixed) = col_names
    }

    if (nrow(dplyr::setdiff(input_catalogue, remained_fixed)) == 0 &
        nrow(new_fixed) == 0) {
      cat('        break loop\n')
      break
    }

    if (nrow(remained_fixed) == 0 & nrow(new_fixed) == 0) {
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

  if (nrow(input_catalogue) == 0) {
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


sanitize_inputs = function(x,
                           reference_catalogue,
                           k,
                           lr,
                           steps,
                           phi,
                           delta,
                           groups = NULL,
                           input_catalogue = NULL,
                           lambda_rate = NULL,
                           sigma = FALSE
                           )
{
  # Input counts
  if (!is.data.frame(x))
    cli::cli_abort("The count matrix should be a dataframe!")
  if (nrow(x) == 0)
    cli::cli_abort("The count matrix has no rows!")

  if (rownames(x) %>% is.null())
  {
    cm = paste0("Sample_", 1:nrow(x))
    rownames(x) = cm

    cli::cli_alert_warning("Sample names were missing, will use {.field {head(cm)}, ...}")
  }

  # Input reference_catalogue
  if (!is.data.frame(reference_catalogue))
    cli::cli_abort("The reference catalogue should be a dataframe!")
  if (nrow(reference_catalogue) == 0)
    cli::cli_abort("The reference catalogue has no rows!")

  if (rownames(reference_catalogue) %>% is.null)
    cli::cli_abort("The reference catalogue has no rownames (signature names)!")

  # What is there
  if (!all(colnames(x) %in% colnames(reference_catalogue)))
    cli::cli_abort("Some columns in the input counts miss from the reference catalogue")

  # ordering
  if (!all(colnames(x) == colnames(reference_catalogue)))
  {
    cli::cli_alert_warning("Reference signature and catalgous have different column orders, will re-order")

    reference_catalogue = reference_catalogue[names(x)]
  }

  # Input catalogue
  if (!is.null(input_catalogue))
  {
    if (!is.data.frame(input_catalogue))
      cli::cli_abort("The input catalogue should be a dataframe!")
    if (nrow(input_catalogue) == 0)
      cli::cli_abort("The input catalogue has no rows!")

    if (rownames(input_catalogue) %>% is.null)
      cli::cli_abort("The input catalogue has no rownames (signature names)!")

    # Ref inluded in input
    if (all(rownames(input_catalogue) %in% rownames(reference_catalogue)))
      cli::cli_abort("The input catalogue has signatures that are not in there reference!")

    # What is there
    if (!all(colnames(x) %in% colnames(input_catalogue)))
      cli::cli_abort("Some columns in the input counts miss from the input catalogue")

    # ordering
    if (!all(colnames(x) == colnames(input_catalogue)))
    {
      cli::cli_alert_warning("Reference catalgous and input counts have different column orders, will re-order")

      input_catalogue = input_catalogue[names(x)]
    }

    input_catalogue = input_catalogue%>% as.data.frame()
  }

  # Numerics
  if (!(is.numeric(k) |
        (!is.list(k) & all(sapply(k, is.numeric)))))
    cli::cli_abort("k must be a list of integers, or a single integer.")

  if(!(is.numeric(lr) & lr > 0)) cli::cli_abort("Invalid learning rate.")
  if(!(is.numeric(steps) & steps > 0)) cli::cli_abort("Invalid number of steps.")

  # TODO - complete
#   phi,
#   delta,
#   groups = NULL,
#   input_catalogue = NULL,
#   lambda_rate = NULL,
#   sigma = FALSE

  return(
    list(
      x = x %>% as.data.frame(),
      reference_catalogue = reference_catalogue %>% as.data.frame(),
      input_catalogue = input_catalogue
    )
  )
}
