

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



#-------------------------------------------------------------------------------

visData <- function(data) {

  ACC <- c()
  MAE <- c()
  MSE <- c()
  fitness.rate <- c()

  data$TargetX <- factor(data$TargetX, levels = c("low", "medium", "high"))
  data$InputX <- factor(data$InputX, levels = c("low", "medium", "high"))
  data$Num_Samples <- factor(data$Num_Samples)

  n.cores <- parallel::detectCores()
  my.cluster <- parallel::makeCluster(n.cores)  # create the cluster
  #print(my.cluster) # check cluster definition (optional)
  doParallel::registerDoParallel(cl = my.cluster) # register it to be used by %dopar%
  #foreach::getDoParRegistered() # check if it is registered (optional)
  #foreach::getDoParWorkers()  # how many workers are available? (optional)

  results <- foreach::foreach(i = 1:nrow(data)) %dopar% {

    m <- data[i, ]$x[[1]]

    input_catalogue <- data[i, ]$Input_Catalogue[[1]]

    #exp_exposure <- data[i, ]$Exp_Exposure[[1]]
    inf_exposure <- data[i, ]$Inf_Exposure[[1]]

    exp_fixed <- data[i, ]$Exp_Fixed[[1]]
    inf_fixed <- data[i, ]$Inf_Fixed[[1]]

    #exp_denovo <- data[i, ]$Exp_Denovo[[1]]
    inf_denovo <- data[i, ]$Inf_Denovo[[1]]

    mr <- reconstruction_matrix(m, inf_exposure, inf_fixed, inf_denovo)

    ACC[length(ACC)+1] <- fixed.quality(input_catalogue, exp_fixed, inf_fixed)
    MAE[length(MAE)+1] <- compute.mae(m, mr)
    MSE[length(MSE)+1] <- compute.mse(m, mr)
    fitness.rate[length(fitness.rate)+1] <- fitness.quality(m, mr, 0.9)
  }

  parallel::stopCluster(cl = my.cluster)

  df <- data.frame(ACC = ACC, MAE = MAE, MSE = MSE, fitness = fitness.rate, TargetX = data$TargetX, InputX = data$InputX, n.samples = data$Num_Samples)
  return(df)
}





