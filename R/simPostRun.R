

#' Title
#'
#' @return
#' @export
#'
#' @examples
init_tibble <- function() {

  obj <- tibble::tibble(
    x = list(),
    Input_Catalogue = list(),
    Ref_Catalogue = list(),

    Exp_Exposure = list(),
    Exp_Fixed = list(),
    Exp_Denovo = list(),

    Inf_Exposure = list(),
    Inf_Fixed = list(),
    Inf_Denovo = list(),

    TargetX = character(),
    InputX = character(),
    Num_Samples = numeric(),
    IterNum = numeric(),

    K = list(),
    Lr = numeric(),
    Steps = numeric(),
    Phi = numeric(),
    Delta = numeric()
  )

  return(obj)
}

#-------------------------------------------------------------------------------


#' Title
#'
#' @param x
#'
#' @import tibble
#'
#' @return
#' @export
#'
#' @examples
fill_tibble <- function(x) {

  data <- init_tibble()
  for (row in x) {
    data <- data %>% tibble::add_row(
      x = list(row$x),
      Input_Catalogue = list(row$input_catalogue),
      Ref_Catalogue = list(row$ref_catalogue),

      Exp_Exposure = list(row$exp_exposure),
      Exp_Fixed = list(row$exp_fixed),
      Exp_Denovo = list(row$exp_denovo),

      Inf_Exposure = list(row$inf_exposure),
      Inf_Fixed = list(row$inf_fixed),
      Inf_Denovo = list(row$inf_denovo),

      TargetX = row$targetX,
      InputX = row$inputX,
      Num_Samples = row$num_samples,
      IterNum = row$iter,

      K = list(row$k),
      Lr = row$lr,
      Steps = row$steps,
      Phi = row$phi,
      Delta = row$delta
    )
  }
  return(data)
}
