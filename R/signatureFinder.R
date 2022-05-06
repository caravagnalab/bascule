#use_condaenv("pybasilica")


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
#' @return inferred exposure matrix, inferred COSMIC signatures and inferred de novo (not from referencecatalog) signatures
#'
#' @importFrom reticulate r_to_py
#' @import ggplot2
#' @import tidyr
#' @import data.table
#' @import gridExtra
#' @export
#'
#' @examples
fit <- function(x,
                input_catalog=NULL,
                k=0:5,
                reference_catalog=basilica::COSMIC,
                lr=0.05,
                steps_per_iter=500,
                fixedLimit=0.05,
                denovoLimit=0.9
                ) {

  pybasilica <- reticulate::import("pybasilica")

  #x <- reticulate::r_to_py(x)
  #input_catalog <- reticulate::r_to_py(input_catalog)
  #----------------------------- MUST BE CHANGED -------------------------------
  reticulate::py_run_string("k = list(map(int, [0, 1, 2, 3, 4, 5]))")
  k <- py$k
  k <- reticulate::r_to_py(k)
  reticulate::py_run_string("steps_per_iter = 500")
  steps_per_iter <- reticulate::r_to_py(py$steps_per_iter)
  #-----------------------------------------------------------------------------
  #reference_catalog <- reticulate::r_to_py(reference_catalog)
  #lr <- reticulate::r_to_py(lr)
  #fixedLimit <- reticulate::r_to_py(fixedLimit)
  #denovoLimit <- reticulate::r_to_py(denovoLimit)

  output <- pybasilica$pyfit(x, input_catalog, k, reference_catalog, lr, steps_per_iter, fixedLimit, denovoLimit)

  return(init_object(output))
}




