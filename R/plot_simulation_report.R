plot_simulation_report <- function(cohort) {

  plotlist <- lapply(1:nrow(cohort),FUN =  function(i) plot_simulation_report_aux(cohort[i,]))

  return(plotlist)

}


plot_simulation_report_aux <- function(x) {

  library(patchwork)

  fit <- x$fit[[1]]

  plot_signatures_fit <- plot_signatures(fit)
  sbs_real_fixed <- x$exp_fixed[[1]] %>% as.matrix() %>% reshape2::melt()
  sbs_real_fixed$Type = "Catalogue"
  sbs_real_denovo <- x$exp_denovo[[1]] %>% as.matrix() %>% reshape2::melt()
  sbs_real_denovo$Type = "De novo"
  sbs_real <- rbind(sbs_real_fixed, sbs_real_denovo)
  colnames(sbs_real) <- c("Signature", "Feature", "Value", "Type")
  plot_signatures_simulation <- plot_signatures_sbs(fit, sbs_real) + scale_fill_discrete() +
    ggplot2::geom_bar(
      ggplot2::aes(x = crazy_map, y = Value, fill = Signature, color = Type),
      stat="identity",
      position="identity") +
  scale_color_manual("Type",values = c("grey10", "grey60", "white"), breaks = c("Catalogue", "De novo", "")) +
    ggtitle(glue::glue("Real signatures ( n = {fit$n_samples})"))

  plot_exposure_fit <- plot_exposure(fit)
  exposure_real  <- x$exp_exposure[[1]] %>% as.matrix() %>% reshape2::melt() %>%
    mutate(Type = if_else(Var2  %in% unique(sbs_real_denovo$Var2) ,"De novo", "Catalogue"))
  colnames(exposure_real) <- c("Sample", "Signature", "Exposure", "Type")

  plot_exposure_real <- ggplot2::ggplot(
    data = exposure_real,
    ggplot2::aes(x=Sample, y=Exposure, fill=Signature)
  ) +
    ggplot2::geom_bar(stat = "identity") +
    my_ggplot_theme() +
    ggplot2::scale_y_continuous(labels=scales::percent) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 90, hjust = 1)
    ) +
    ggplot2::labs(
      title = paste0(fit$cohort, ' (n = ', fit$n_samples, ')')    )

  #fit_similarity <- plot_similarity_reference(fit)

  stats <- evaluate.cohort(x)
  ps_MAE <- simple_barplot(data.frame(value = stats$mae), "MAE")
  ps_MSE <- simple_barplot(data.frame(value = stats$mse), "MSE")
  ps_FA <- simple_barplot(data.frame(value = stats$fixed_acc), "Catalogue accuracy")
  ps_DA <- simple_barplot(data.frame(value = stats$denovo_ratio), "De novo accuracy")
  if(!is.null(stats$denovo_match[[1]])) {
    ps_DS <- ggplot(stats$denovo_match[[1]] %>% as.data.frame(), aes(x = match, y = similarity)) +
      geom_col(color = "grey30") + my_ggplot_theme() +
      ylab("") + xlab("") + ggtitle("De novo similarity") + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  } else {
    ps_DS <- ggplot()
  }

  ps_plot <- ps_MAE | ps_MSE | ps_FA | ps_DA | ps_DS

  final_plot <- ggpubr::ggarrange(plot_signatures_fit ,
                                  plot_signatures_simulation,
                                  plot_exposure_fit,
                                  plot_exposure_real,
                                  ps_plot, ncol = 1)

  return(final_plot)
}


simple_barplot <- function(df,title) {
  ggplot(df, aes(x = "", y = value)) + geom_col(color = "grey30") + my_ggplot_theme() + ylab("") + xlab("") + ggtitle (title)
}
