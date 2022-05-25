#' @import tibble
#' @import dplyr



#-------------------------------------------------------------------------------

generate_data <- function(ref_catalogue,
                          target_complexity,
                          input_complexity,
                          num_samples) {

  pybasilica <- reticulate::import("pybasilica")

  x <- pybasilica$input_generator(
    full_cosmic_df = ref_catalogue,
    target_complexity = target_complexity,
    input_complexity = input_complexity,
    num_samples = num_samples
  )

  obj = list()

  # counts (data.frame)
  obj$x <- reticulate::py_to_r(x$M)
  # expected exposure (data.frame)
  obj$exp_exposure <- reticulate::py_to_r(x$alpha)

  # expected catalogue signatures (data.frame / NULL)
  if (is.null(x$beta_fixed)) {
    obj$exp_fixed <- x$beta_fixed
  } else {
    obj$exp_fixed <- reticulate::py_to_r(x$beta_fixed)
  }
  # expected denovo signatures (data.frame / NULL)
  if (is.null(x$beta_denovo)) {
    obj$exp_denovo <- x$beta_denovo
  } else {
    obj$exp_denovo <- reticulate::py_to_r(x$beta_denovo)
  }
  # input catalogue (data.frame / NULL)
  if (is.null(x$beta_input)) {
    obj$input_catalogue <- x$beta_input
  } else {
    obj$input_catalogue <- reticulate::py_to_r(x$beta_input)
  }

  # reference catalogue (data.frame)
  obj$ref_catalogue <- reticulate::py_to_r(x$cosmic_df)

  obj$targetX <- target_complexity  # character
  obj$inputX <- input_complexity    # character
  obj$num_samples <- num_samples    # numeric

  return(obj)
}

#-------------------------------------------------------------------------------
sim_fit <- function(
    ref_catalogue,
    target_complexity,
    input_complexity,
    num_samples,
    k,
    lr,
    steps,
    phi,
    delta,
    iter
) {

  x <- generate_data(
    ref_catalogue,
    target_complexity,
    input_complexity,
    num_samples
  )

  pybasilica <- reticulate::import("pybasilica")
  fit <- pybasilica$pyfit(
    M=x$x,
    groups=NULL,
    input_catalogue=x$input_catalogue,
    reference_catalogue=x$ref_catalogue,
    k=k,
    lr=lr,
    steps=steps,
    phi=phi,
    delta=delta
  )

  x$inf_exposure = fit[[1]] # data.frame
  x$inf_fixed = fit[[2]]    # data.frame
  x$inf_denovo = fit[[3]]   # data.frame

  x$k = k         # integer
  x$lr = lr       # numeric
  x$steps = steps # numeric
  x$phi = phi     # numeric
  x$delta = delta # numeric
  x$iter = iter

  return(x)
}

#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------

#' Title
#'
#' @param targetX_list
#' @param inputX_list
#' @param num_samples_list
#' @param num_iters
#'
#' @import doParallel
#' @return
#' @export
#'
#' @examples
parFit <- function(
    targetX_list,
    inputX_list,
    num_samples_list,
    num_iters
    ) {

  my_params <- expand.grid(
    num_iters = num_iters,
    num_samples_list = num_samples_list,
    targetX_list = targetX_list,
    inputX_list = inputX_list
  )

  n.cores <- parallel::detectCores()
  my.cluster <- parallel::makeCluster(n.cores)  # create the cluster
  #print(my.cluster) # check cluster definition (optional)
  doParallel::registerDoParallel(cl = my.cluster) # register it to be used by %dopar%
  #foreach::getDoParRegistered() # check if it is registered (optional)
  #foreach::getDoParWorkers()  # how many workers are available? (optional)

  results <- foreach::foreach(i = 1:nrow(my_params)) %dopar% {
    library(reticulate)
    use_condaenv("pybasilica")

    basilica:::sim_fit(
      ref_catalogue = basilica::COSMIC_catalogue,
      target_complexity = as.character(my_params[i, ]$targetX_list),
      input_complexity = as.character(my_params[i, ]$inputX_list),
      num_samples = my_params[i, ]$num_samples_list,
      k=0:5,
      lr=0.05,
      steps=500,
      phi=0.05,
      delta=0.9,
      iter = my_params[i, ]$num_iters
    )
  }

  parallel::stopCluster(cl = my.cluster)

  return(results)

}


#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------






