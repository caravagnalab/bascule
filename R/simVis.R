

visData <- function(data) {

  ACC <- c()
  MAE <- c()
  MSE <- c()
  fitness.rate <- c()

  data$TargetX <- factor(data$TargetX, levels = c("low", "medium", "high"))
  data$InputX <- factor(data$InputX, levels = c("low", "medium", "high"))
  data$Num_Samples <- factor(data$Num_Samples)

  for (i in 1:nrow(data)) {

    m <- data[i, ]$x[[1]]

    input_catalogue <- data[i, ]$Input_Catalogue[[1]]

    #exp_exposure <- data[i, ]$Exp_Exposure[[1]]
    inf_exposure <- data[i, ]$Inf_Exposure[[1]]

    exp_fixed <- data[i, ]$Exp_Fixed[[1]]
    inf_fixed <- data[i, ]$Inf_Fixed[[1]]

    #exp_denovo <- data[i, ]$Exp_Denovo[[1]]
    inf_denovo <- data[i, ]$Inf_Denovo[[1]]

    #data[i, ]$IterNum
    #data[i, ]$Num_Samples
    #data[i, ]$TargetX
    #data[i, ]$InputX

    mr <- reconstruction_matrix(m, inf_exposure, inf_fixed, inf_denovo)

    ACC[length(ACC)+1] <- fixed.quality(input_catalogue, exp_fixed, inf_fixed)
    MAE[length(MAE)+1] <- compute.mae(m, mr)
    MSE[length(MSE)+1] <- compute.mse(m, mr)
    fitness.rate[length(fitness.rate)+1] <- fitness.quality(m, mr, 0.9)
  }

  df <- data.frame(ACC = ACC, MAE = MAE, MSE = MSE, fitness = fitness.rate, TargetX = data$TargetX, InputX = data$InputX, n.samples = data$Num_Samples)
  return(df)
}


#===============================================================================
# VISUALIZATION
#===============================================================================

plot_mae <- function(x) {
  ggplot(data = x, aes(x = n.samples, y = MAE, fill=n.samples)) +
    geom_boxplot() +
    facet_grid(InputX ~ TargetX , labeller = label_both, scales = 'fixed') +
    theme_fivethirtyeight() +
    xlab("No. of Samples") +
    ylab("Mean Average Error") +
    ggtitle("Evaluation of Various Categories of Sysnthetic Data")
    #scale_fill_lancet()
    #geom_hline(yintercept=c(0.5, 0.8, 1), linetype='dashed', color='red')
}

plot_mse <- function(x) {
  ggplot(data = x, aes(x = n.samples, y = MSE, fill=n.samples)) +
    geom_boxplot() +
    facet_grid(InputX ~ TargetX , labeller = label_both, scales = 'fixed') +
    theme_fivethirtyeight() +
    xlab("No. of Samples") +
    ylab("Mean Squared Error") +
    ggtitle("Evaluation of Various Categories of Sysnthetic Data")
    #scale_fill_lancet()
    #geom_hline(yintercept=c(0.5, 0.8, 1), linetype='dashed', color='red')
}


plot_acc <- function(x) {
  ggplot(data = x, aes(x = n.samples, y = ACC, fill=n.samples)) +
    geom_boxplot() +
    facet_grid(InputX ~ TargetX , labeller = label_both, scales = 'fixed') +
    theme_fivethirtyeight() +
    xlab("No. of Samples") +
    ylab("Accuracy") +
    ggtitle("Evaluation of Various Categories of Sysnthetic Data (From Catalogue)")
    #scale_fill_lancet()
    #geom_hline(yintercept=c(0.5, 0.8, 1), linetype='dashed', color='red')
}


