
#' @import ggthemes

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




