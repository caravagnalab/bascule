library(ggplot2)
library(ggrepel)
library(ggpubr)
library(reshape2)


remove.sparsity <- function(x, type) {
  x[x < 0.01] <- 0
  df <- data.frame(type = c(), name = c(), mean = c(), freq = c())
  for (i in 1:ncol(x)) {
    ls <- list(
      type = type, 
      name = colnames(x)[i], 
      mean = sum(x[,i]) / sum(x[,i] != 0), 
      freq = sum(x[,i] != 0)
    )
    df <- rbind(df, ls)
    df$type <- factor(df$type)
  }
  return(df)
}


# input
#   basilica_exposure --> wide dataframe
#   serena_exposure ----> wide dataframe
#   basilica_singles ---> unmapped basilica signatures (character vector)
#   serena_singles -----> unmapped serena signatures (character vector)
# output
#   ggplot object
plot_unexplained <- function(basilica_exposure, serena_exposure, basilica_singles, serena_singles) {
  
  df <- rbind(
    remove.sparsity(basilica_exposure[, basilica_singles], "basilica"), # unpaired
    remove.sparsity(serena_exposure[, serena_singles], "serena") # unpaired
  )
  
  df <- df %>% replace(is.na(.), 0)
  
  p <- ggplot(
    data = df, 
    aes(x = freq, y = mean, label=name)) + 
    geom_point(aes(color = factor(type)), size=2) + 
    ggtitle("Unexplained Signatures (exposure > 1%)") + 
    xlab("number of samples") + 
    ylab("mean exposure") + 
    #labs(fill='Program Type') + 
    geom_label_repel(
      aes(label = name),
      box.padding = 0.35, 
      point.padding = 0.5,
      segment.color = 'grey50', 
      max.overlaps = 25
    )
  return(p)
}

#-------------------------------------------------------------------------------

# input
#   basilica_exposure --> wide dataframe
#   serena_exposure ----> wide dataframe
#   basilica_serena ----> mapped signatures (basilica vs. serena)
# output
#   ggplot object
plot_exposure_comp <- function(basilica_exposure, serena_exposure, basilica_serena) {
  p_list <- list()
  for (i in 1:nrow(basilica_serena)) {
    df <- data.frame(
      b_pair = basilica_exposure[, basilica_serena[i,]$Basilica], 
      s_pair = serena_exposure[, basilica_serena[i,]$Serena]
    )
    p_list[[length(p_list)+1]] <- ggplot(df, aes(x=b_pair, y=s_pair)) + 
      geom_point() + 
      ggtitle("Exposure (Basilica vs. Serena)") +
      xlab(paste0(basilica_serena[i,]$Basilica, " exposure (basilica)")) +
      ylab(paste0(basilica_serena[i,]$Serena, " exposure (serena)"))
  }
  
  return(ggpubr::ggarrange(plotlist = p_list, common.legend = TRUE))
}


# ==============================================================================
# ==============================================================================
# ==============================================================================


# input
#   basilica_exposure --> wide dataframe
#   serena_exposure ----> wide dataframe
#   map ----------------> map.data function object
# output
#   ggplot object
plot_exposure_dist <- function(basilica_exposure, serena_exposure, map) {
  
  df <- data.frame(type = c(), pair = c(), exposure = c())
  for (i in 1:nrow(map$basilica_serena)) {
    nsamples <- nrow(serena_exposure)
    type <- c(rep("Baslica", nsamples), rep("Serena", nsamples))
    ls <- list(
      type = type, 
      pair = rep(paste0(map$basilica_serena$Basilica[i], "-", map$basilica_serena$Serena[i]), 2*nsamples), 
      exposure = c(
        basilica_exposure[, map$basilica_serena$Basilica[i]], 
        serena_exposure[, map$basilica_serena$Serena[i]]
      )
    )
    df <- rbind(df, ls)
  }
  
  p <- ggplot(df, aes(x = pair, y = exposure, color = type)) + 
    geom_violin()
  
  return(p)
}

#-------------------------------------------------------------------------------

# input : qc.linearCombination function output object
# output: ggplot object
plot.coef.heatmap <- function(x) {
  # Create the heatmap plot
  p <- ggplot(data = subset(x, coef > 0), aes(x = signature, y = denovos, fill = round(coef, 3), label = round(coef, 3))) +
    geom_tile(color = "white") + 
    scale_fill_gradient(low = "white", high = "red") + 
    geom_text(color = "black", size = 3) +  # Add text annotations
    #scale_fill_gradient(low = "white", high = "steelblue") +  # Choose your desired color gradient
    labs(fill='Coefficient') + 
    theme_minimal() + 
    theme(
      # remove the vertical grid lines
      panel.grid.major.x = element_blank(), 
      # explicitly set the horizontal lines (or they will disappear too)
      panel.grid.major.y = element_line( linewidth=.1, color="black" ), 
      #legend.position="none"
    ) + 
    labs(title = "Heatmap of Coefficient Values",
         x = "All Signatures (fixed + denovo)",
         y = "Denovo Signatures")
  
  return(p)
}





