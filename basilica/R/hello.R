# Hello, world!
#
# This is an example function named 'hello'
# which prints 'Hello, world!'.
#
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Install Package:           'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'

hello <- function() {
  print("Hello, world!")
}



#fit <- function(catalog, beta_input, k_list=0:5, beta_cosmic, fixedLimit=0.05, denovoLimit=0.9) {

fit <- function(x, groups, input_catalog=NULL, k=0:5, reference_catalog=basilica::COSMIC, fixedLimit=0.05, denovoLimit=0.9) {
  x <- r_to_py(catalog)
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
  fixedLimit <- r_to_py(fixedLimit)
  denovoLimit <- r_to_py(denovoLimit)

  output <- BaSiLiCa(x, groups, input_catalog, k, reference_catalog, fixedLimit, denovoLimit)

  alpha <- output[[1]]
  beta_fixed <- output[[2]]
  beta_denovo <- output[[3]]

  return(list(Alpha = alpha, Beta_Fixed = beta_fixed, Beta_Denovo = beta_denovo))

}






