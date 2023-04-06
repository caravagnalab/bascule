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
#' @param filt_pi threshold for COSMIC signature weight to be included as basis of the denovo signatures
#' @param groups vector of discrete labels with one entry per sample, it defines the groups that will be considered by basilica
#' @param input_catalogue input signature profiles, NULL by default
#' @param enforce_sparsity use Laplace prior over exposure weights (bayesian LASSO)
#' @param cohort
#' @param max_iterations
#' @param blacklist
#' @param lambda_rate
#' @param sigma
#'
#' @return inferred exposure matrix, inferred signatures from reference catalogue and inferred de novo (not from reference catalogue) signatures
#' @export fit
#'
#' @examples

fit <- function(x,
                k,
                py = NULL,
                reference_catalogue = basilica::COSMIC_catalogue,
                input_catalogue = basilica::COSMIC_catalogue["SBS1", ],
                filtered_cat = FALSE,
                cohort = "MyCohort",
                lr = 0.01,
                steps = 500,
                max_iterations = 20,
                blacklist = NULL,
                phi = 0.05,
                delta = 0.9,
                filt_pi =0.1,
                groups = NULL,
                lambda_rate = NULL,
                sigma = FALSE,
                CUDA = FALSE,
                compile = FALSE,
                enforce_sparsity = FALSE,
                store_parameters = FALSE,
                regularizer="cosine",
                reg_weight = 1,
                reg_bic = FALSE,
                cosine_by_subs = FALSE)

{

  sig_col = function(x) crayon::blue(x)

  fit = list()
  class(fit) = 'basilica_obj'

  # cli::cli_h1("MUSICA - MUtational Signature Inference with a CAtalogue ")
  cli::cli_h1("Basilica - Bayesian signature learning with a catalogue")
  cat("\n")

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
    sigma = sigma,
    blacklist = blacklist
  )

  x = sanitized_inputs$x
  reference_catalogue = sanitized_inputs$reference_catalogue
  input_catalogue = sanitized_inputs$input_catalogue

  # Report messages for the inputs

  cli::cli_alert("                       Input samples : {.field n = {nrow(x)}}")
  cli::cli_alert("                     Output clusters : {.field k = {k}}.")
  cli::cli_alert("                        Blacklist by : {.field {blacklist}}.")
  cli::cli_alert("                  Maximum iterations : {.field {max_iterations}}.")

  # Report messages for the reference
  cli::cli_alert(
    "Reference Catalogue Signatures (RCSs): {.field {rownames(reference_catalogue) %>% sig_col}} ({.field n = {nrow(reference_catalogue)}})"
  )

  # Report messages for the input
  if (!is.null(input_catalogue)) {
    cli::cli_alert(
      "    Input Catalogue Signatures (ICSs): {.field {head(rownames(input_catalogue)) %>% sig_col}, ...} ({.field n = {nrow(input_catalogue)}})"
    )
    input_catalogue %>% dplyr::as_tibble() %>% print()
    } else {
    cli::cli_alert_warning("    Input Catalogue Signatures (ICSs): none (no supervision).")
  }

  cli::cli_h2("Iterative Basilica inference algorithm")

  TIME = as.POSIXct(Sys.time(), format = "%H:%M:%S")

  counter <- 1
  black_list <- c()

  # Iterative objects
  BIC_trajectories = ICS_trajectories = DNS_trajectories =
    blacklist_trajectories = NULL

  # repeat untill convergence
  repeat {
    cli::cli_h2("Basilica step {.field {counter}}")

    ICS_trajectories = append(ICS_trajectories, list(input_catalogue))
    blacklist_trajectories = append(blacklist_trajectories, list(black_list))

    n_ICSs = ifelse(is.null(input_catalogue), 0, nrow(input_catalogue))
    cli::cli_h3("Bayesian NMF via SVI [{.field {steps}} steps, ICSs {.field k = {n_ICSs}}]")

    TIME_it = as.POSIXct(Sys.time(), format = "%H:%M:%S")

    ################### Bayesian NMF
    k_aux = k
    if(0 %in% k & n_ICSs == 0) {
      k = k[k != 0]
      cli::cli_alert_info("ICSs {.field k = {n_ICSs}}, removing {.field k = {0}} for this run, using {.field k = {k}} .")

      if(length(k) == 0) cli::cli_abort("Cannot proceed: ICSs {.field k = {n_ICSs}} and {.field k = {0}}.")
    }

    obj = pyfit(
      x = x,
      py = py,
      k_list = k,
      lr = lr,
      n_steps = steps,
      groups = groups,
      input_catalogue = input_catalogue,
      lambda_rate = lambda_rate,
      sigma = sigma,
      CUDA = CUDA,
      compile = compile,
      enforce_sparsity = enforce_sparsity,
      store_parameters = store_parameters,
      regularizer = regularizer,
      reg_weight = reg_weight,
      reg_bic = reg_bic
    )

    if (filtered_cat)
      obj$denovo_signatures = renormalize_denovo_thr(obj$denovo_signatures)

    k = k_aux

    TIME_it = difftime(as.POSIXct(Sys.time(), format = "%H:%M:%S"),
                       TIME_it, units = "mins") %>% round(2)

    loss = obj$losses[steps] %>% round(2)
    bic = obj$bic %>% round(2)
    BIC_trajectories = c(BIC_trajectories, bic)

    best_k = obj$exposure %>% ncol()
    fix_k = input_catalogue %>% nrow()
    if(fix_k %>% is.null) fix_k = 0

    best_k = best_k - fix_k

    cli::cli_alert_success(
      paste0(
        "NMF completed in {.field {TIME_it}} minutes [BIC {.field {bic}}, Loss {.field {loss}}, best DNSs {.field k = {best_k}} given ICSs {.field k = {n_ICSs}}]. "
      )
    )

    ################### Filter out small-exposure signatures
    cli::cli_h3("Checking ICSs exposure, \u03b1 > \u03A6 ({.field \u03A6 = {phi}}).")

    if (is.null(input_catalogue) || nrow(input_catalogue) == 0)
      cli::cli_alert_danger("No ICSs in this step were used.")


    if(!is.null(blacklist)) {

      # drop non-significant fixed signatures --------------------------------
      if(!is.null(blacklist) && blacklist == "TMB")

        if(!is.null(blacklist)) {
          if(blacklist == "TMB")
            a <- filter.fixed(
              M = x,
              alpha = obj$exposure,
              beta_fixed = input_catalogue,
              phi = phi)

          if(!is.null(blacklist) && blacklist == "freq")
            a = filter.fixed_minfreq(
              alpha = obj$exposure,
              beta_fixed = input_catalogue,
              phi = phi)

        } else {
          a = filter.fixed_nofilter( # fake function
            alpha = obj$exposure,
            beta_fixed = input_catalogue)
        }

    } else {
      a = filter.fixed_nofilter( # fake function
        alpha = obj$exposure,
        beta_fixed = input_catalogue)
    }


    remained_fixed <- a$remained_fixed # data.frame / NULL

    # Here we might have dropped something
    if (!is.null(a$dropped_fixed)) {
      n_dropped = a$dropped_fixed %>% nrow()
      w_dropped = a$dropped_fixed %>% rownames

      cli::cli_alert_warning(
        paste0(
          "{.field {n_dropped}} ICSs will be dropped ({.field {w_dropped %>% sig_col}}), and blacklisted."
        )
      )

      black_list <- union(black_list, w_dropped)  # character vector

      cli::cli_alert_warning(paste0(
        "The updated blacklist is: {.field {black_list %>% crayon::red()}}."
      ))
    } else {

      if (!is.null(input_catalogue)) {
        n_cat = input_catalogue %>% nrow()

        # Here we did NOT drop, and had some input signatures - we report it
        if (n_cat > 0)
          cli::cli_alert_success("No ICSs will be dropped.")
      }
    }

    # TEST
    # cat('class(a$dropped_fixed):', class(a$dropped_fixed), '\n')
    # cat('a$dropped_fixed:')
    # print(a$dropped_fixed)
    # cat('class(black_list):', class(black_list), '\n')
    # cat('black_list:', black_list, '\n')
    # TEST

    n_denovo = obj$denovo_signatures %>% nrow
    n_catalogue = reference_catalogue %>% nrow

    cli::cli_h3(
      "Comparing {.field {n_denovo}} DNSs to {.field {n_catalogue}} RCSs, imposing cosine similarlity above {.field \u0394 = {delta}}."
    )

    # detect denovo signatures which are similar to reference signatures -------
    if (cosine_by_subs)
      substitutions = get_contexts(obj) %>% dplyr::pull(subs) %>% unique() else
      substitutions = NULL

    obj$exposure = filter.denovo.phi(exposures = obj$exposure,
                                      denovo = obj$denovo_signatures %>% rownames(),
                                      phi = phi) # filter denovo if low exposure in all patients

    # check if reference are linear comb of denovo sigs co-occurring (exp>thr)
    b_denovo <- filter.denovo.QP(
      reference = rbind(input_catalogue, reference_catalogue),
      # beta_fixed = rbind(input_catalogue, reference_catalogue),
      beta_denovo = obj$denovo_signatures,
      thr_exposure = phi,
      exposures = obj$exposure,
      black_list = black_list,
      delta = delta,
      filt_pi = filt_pi,
      substitutions = substitutions)

    if (!is.null(b_denovo$new_fixed))
      print(b_denovo)

    denovo_filt = b_denovo$reduced_denovo

    ref = setdiff(reference_catalogue, rbind(input_catalogue, b_denovo$new_fixed))
    # check if denovo are linear comb of reference sigs
    b_reference <- filter.denovo.QP(
      reference = ref,
      # beta_fixed = rbind(input_catalogue, b_denovo$new_fixed),
      beta_denovo = denovo_filt,
      black_list = black_list,
      delta = delta,
      filt_pi = filt_pi,
      substitutions = substitutions)

    new_fixed <- rbind(b_denovo$new_fixed, b_reference$new_fixed) %>% unique()

    # if (!is.null(b_denovo$new_fixed))

    if (!is.null(new_fixed) && nrow(new_fixed) > 0) {
      n_denovo_denovo = b_reference$reduced_denovo %>% nrow()
      # cli::cli_alert_warning(
      #   "{.field {nrow(new_fixed)}}/{.field {n_denovo}} DNSs were found in the reference catalogue, and will become part of the ICSs!"
      # )
    } else {
      # cli::cli_alert_success("No DNSs were found in the reference catalogue!")
    }

    reduced_denovo <- b_reference$reduced_denovo  # data.frame / NULL (remaining denovo signatures)

    DNS_trajectories = append(DNS_trajectories, list(reduced_denovo))

    #TEST---------------------------------------------------------
    # cat("        fixed          :", rownames(input_catalogue), '\n')
    # cat("        remained fixed :", rownames(remained_fixed), '\n')
    # cat("        black list     :", black_list, "\n\n")
    # cat("        denovo         :", rownames(obj$denovo_signatures), '\n')
    # cat("        new fixed      :", rownames(new_fixed), '\n')
    # cat("        reduced denovo :", rownames(reduced_denovo), '\n')
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

    # Convergency test:
    # - we keep finding exactly the ICSs
    # - there are no DNSs that seem to actually be part of RCSs
    if (nrow(dplyr::setdiff(input_catalogue, remained_fixed)) == 0 &
        nrow(new_fixed) == 0) {
      cat('\n')
      cli::cli_alert_success("Converged - ICSs is stable and all DNSs are genuine.")

      break
    }

    if (nrow(remained_fixed) == 0 & nrow(new_fixed) == 0) {
      input_catalogue <- NULL
    } else {
      input_catalogue <- rbind(remained_fixed, new_fixed)

      cli::cli_alert_warning(
        paste0(
          "The updated ICSs contains now: {.field {input_catalogue %>% rownames %>% sig_col}}."
        )
      )
    }

    counter <- counter + 1
    if (counter > max_iterations) {

      cat('\n')
      cli::cli_alert_danger("Converged forced after {.value {crayon::red(counter)}} iterations.")

      break
    }
    # cat('    ------------------------------------\n')
  }
  # /end repeat untill convergence


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

  TIME = difftime(as.POSIXct(Sys.time(), format = "%H:%M:%S"),
                  TIME, units = "mins")
  TIME = TIME %>% round(2)

  cli::cli_h3("Basilica completed in {.field {TIME}} minutes and {.field {counter}} iterations.")

  if (length(BIC_trajectories) > 1) {
    cli::cli_h3("Fit statistics")

    # Print the BIC trajectory in +/- %
    BIC_trajectories_delta = sapply(1:(length(BIC_trajectories) - 1),
                                    function(i) {
                                      d = BIC_trajectories[i + 1] - BIC_trajectories[i]
                                      d_abs = abs(d)
                                      d_pro = (d_abs / BIC_trajectories[i])
                                      if (d > 0)
                                        d_pro * 100
                                      else-1 * d_pro  * 100
                                    })
    BIC_trajectories_delta = BIC_trajectories_delta %>% round(2)

    BIC_trajectories_delta = sapply(BIC_trajectories_delta,
                                    function(x) {
                                      if (x > 0)
                                        paste0(x, '%') %>% crayon::green()
                                      else
                                        paste0(x, '%') %>% crayon::red()
                                    })

    paste0(
      "> BIC ",
      BIC_trajectories[1],
      " : ",
      paste(BIC_trajectories_delta, collapse = " \u2192 ")
    ) %>% cat()
  }


  # Store in the output obj all relevant information
  fit$cohort = cohort

  fit$n_samples = x %>% nrow()

  #### CHECK ########
  # cat_expo = obj$denovo_signatures %>% rownames()
  # cat_expo = setdiff(obj$exposure %>% colnames, cat_expo)
  # fit$n_catalogue = obj$exposure[, cat_expo] %>% ncol()

  fit$n_denovo = ifelse(
    obj$denovo_signatures %>% is.null,
    0,
    obj$denovo_signatures %>% nrow
  )

  if (is.matrix(obj$alpha)) {
    obj$alpha = obj$alpha %>% as.data.frame()
    colnames(obj$alpha) = c(obj$catalogue_signatures %>% rownames(),
                            obj$denovo_signatures %>% rownames())
  }


  fit$input = list(
    counts = x,
    reference_catalogue = reference_catalogue,
    input_catalogue = ICS_trajectories[[1]]
    )

  fit$params = list(
    k = k,
    lr = lr,
    steps = steps,
    phi = phi,
    delta = delta,
    groups = groups,
    lambda_rate = lambda_rate,
    sigma = sigma
  )

  fit$iterations = list(
    BIC = BIC_trajectories,
    ICS = ICS_trajectories,
    DNS = DNS_trajectories,
    blacklist = blacklist_trajectories
  )

  fit$fit = obj


  return(fit)
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
                           sigma = FALSE,
                           blacklist)

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
    if (!all(rownames(input_catalogue) %in% rownames(reference_catalogue)))
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

    input_catalogue = input_catalogue %>% as.data.frame()
  }

  # Numerics
  if (!(is.numeric(k) |
        (!is.list(k) & all(sapply(k, is.numeric)))))
    cli::cli_abort("k must be a list of integers, or a single integer.")

  if (!(is.numeric(lr) &
        lr > 0))
    cli::cli_abort("Invalid learning rate.")
  if (!(is.numeric(steps) &
        steps > 0))
    cli::cli_abort("Invalid number of steps.")

  # TODO - complete
  #   phi,
  #   delta,
  #   groups = NULL,
  #   input_catalogue = NULL,
  #   lambda_rate = NULL,
  #   sigma = FALSE
  # blacklist
  # max_iterations (to be added yet...)


  return(
    list(
      x = x %>% as.data.frame(),
      reference_catalogue = reference_catalogue %>% as.data.frame(),
      input_catalogue = input_catalogue
    )
  )
}
