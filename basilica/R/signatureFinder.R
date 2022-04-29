


#' Title
#'
#' @param x input mutational counts data (data.frame; rows as samples and columns as 96 mutational categories)
#' @param input_catalog input signature profiles, NULL by default
#' @param k vector of possible number of de novo signatures to infer
#' @param reference_catalog a catalog of reference signatures that basilica will use to compare input and de novo signatures
#' @param lr
#' @param steps_per_iter
#' @param fixedLimit threshold to discard the signature based on its value in exposure matrix
#' @param denovoLimit threshold to consider inferred signature as COSMIC signature
#'
#' @return
#' @export
#'
#' @examples
fit <- function(x, input_catalog=NULL, k=0:5, reference_catalog=basilica::COSMIC, lr=0.05, steps_per_iter=500, fixedLimit=0.05, denovoLimit=0.9) {

  # fit <- function(catalog, beta_input, k_list=0:5, beta_cosmic, fixedLimit=0.05, denovoLimit=0.9)
  x <- r_to_py(catalog)
  input_catalog <- r_to_py(input_catalog)
  #----------------------------- MUST BE CHANGED -------------------------------
  setwd("/home/azad/Documents/thesis/SigPhylo/PyBaSiLiCa")
  source_python("basilica.py")
  py_run_string("k = list(map(int, [0, 1, 2, 3, 4, 5]))")
  k <- py$k
  k <- r_to_py(k)
  #-----------------------------------------------------------------------------
  reference_catalog <- r_to_py(reference_catalog)
  lr <- r_to_py(lr)
  steps_per_iter <- r_to_py(steps_per_iter)
  fixedLimit <- r_to_py(fixedLimit)
  denovoLimit <- r_to_py(denovoLimit)

  # add groups as one of the arguments later
  # groups: vector of discrete labels with one entry per sample, it defines the groups that will be considered by basilica
  #output <- BaSiLiCa(x, groups, input_catalog, k, reference_catalog, fixedLimit, denovoLimit)
  output <- BaSiLiCa(x, input_catalog=NULL, k, reference_catalog, lr, steps_per_iter, fixedLimit, denovoLimit):

  alpha <- output[[1]]
  beta_fixed <- output[[2]]
  beta_denovo <- output[[3]]

  return(list(Alpha = alpha, Beta_Fixed = beta_fixed, Beta_Denovo = beta_denovo))
}






