
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

# input : basilica object
# output: ggplot (samples frequency in clusters)
plot.cluster.freq <- function(x) {
  df <- basilica:::get_cluster_assignments(x) # dataframe
  p <- ggplot(df, aes(x = clusters)) + 
    geom_bar(stat = "count", fill = "purple") +
    labs(x = "Clusters", y = "Number of samples", title = "Distribution of samples in clusters") + 
    theme_minimal()
  return(p)
}


#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

# input:
#   x           : basilica object
#   clusterName : cluster name
#   threshold   : exposure cut-off value to merge the corresponding signatures in low-exp signature
# output:
#   ggplot (signatures exposures samples-wise in specified cluster)
plot.cluster.exposure <- function(x, type, clusterName, threshold) {
  
  exposure <- basilica:::get_exposure(
    x, 
    types = type, 
    samples = basilica:::get_cluster_assignments(x, clusters = clusterName)$samples, 
    clusters = clusterName, 
    add_groups = TRUE, 
    matrix = TRUE
  )[[1]]
  
  all_signames <- exposure %>%dplyr::select(-c("clusters")) %>% apply(2, function(x) all(x < threshold))
  low_signames <- all_signames[all_signames == TRUE] %>% names()
  high_signames <- all_signames[all_signames == FALSE] %>% names()
  
  high_sigs <- exposure %>% 
    dplyr::select(all_of(high_signames))
  
  low_sigs <- exposure %>% 
    dplyr::select(all_of(low_signames)) %>% 
    apply(1, sum) %>% data.frame()
  colnames(low_sigs) <- "low-exp"
  
  df <- merge(high_sigs, low_sigs, by=0, all=TRUE) %>% data.frame(row.names = 1)
  df <- df %>% dplyr::rename_at('low.exp', ~'low-exp')
  
  df_long <- basilica:::wide_to_long(df, what = "exposures")
  
  colors <- basilica:::gen_palette(x=x, types=c(type), n=NULL)
  colors["low-exp"] <- "#030303"
  
  p <- ggplot(data = df_long, aes(x = samples, y = value, fill = sigs)) +
    geom_bar(stat = "identity", position = "stack") +
    labs(title = paste("Mutational Signature Exposures", clusterName, sep = " "),
         x = "Sample",
         y = "Exposure") +
    scale_fill_viridis_d() +
    #theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
    theme(axis.text.x = element_blank()) + 
    ggplot2::scale_fill_manual(values=colors)
  
  return(p)
}


#-------------------------------------------------------------------------------
#--------------------------------- SCORES --------------------------------------
#-------------------------------------------------------------------------------

# plot the final score for a cluster
plot.scores.aux1 <- function(df, types, clusterName) {
  
  df <- subset(df, type == types)
  
  val <- df %>% subset(cluster == clusterName) %>% dplyr::pull(score) %>% quantile(probs = c(0.9))
  # val <- quantile(df$score, probs = c(0.9))
  
  p <- ggplot(subset(df[c(1, 2, 6)], cluster == clusterName), aes(x = signature, y = score, color = cluster, group = cluster)) + 
    geom_line() +
    geom_point() +
    labs(x = "Mutational Signature", y = "Value", title = "Values by Mutational Signatures and Clusters") + 
    geom_hline(yintercept = val, linetype = "dashed", color = "black") + 
    theme_minimal()
  return(p)
}

#-------------------------------------------------------------------------------

# plot the individual scores for cluster
plot.scores.aux2 <- function(df, types, clusterName) {
  
  df <- subset(df, type == types)
  
  df_long <- pivot_longer(df[c(1,2,3,4,5)], cols = c(varRatio, activeRatio, mutRatio), names_to = "name", values_to = "value")
  p <- ggplot(subset(df_long, cluster == clusterName), aes(x = signature, y = value, color = name, group = name)) + 
    geom_line() + 
    geom_point() + 
    labs(x = "Mutational Signature", y = "Value", title = "Cluster Analysis") + 
    theme_minimal()
  return(p)
}

#-------------------------------------------------------------------------------
library(patchwork)
# all plots together
plot.cluster.score <- function(x, types, clusterName, threshold) {
  
  df <- compute.scores(x, types, threshold)
  
  p1 <- plot.scores.aux1(df, types, clusterName)
  p2 <- plot.scores.aux2(df, types, clusterName)
  
  p3 <- plot.cluster.exposure(x, types, clusterName, threshold)
  
  return((p1 | p2) / p3)
}

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

significant_signatures <- function(x, types, threshold) {
  
  df <- compute.scores(x, types, 0.05)
  
  aa <- tapply(df$score, df$cluster, function(x) quantile(x, probs = c(0.9)))
  df2 <- data.frame(matrix(ncol = 3, nrow = 0))
  colnames(df2) <- c("signature", "cluster", "score")
  
  for (i in 1:length(aa)) {
    #print(paste0( names(aa[i]), " - ", aa[[i]]) )
    
    df2 <- rbind(
      df2, 
      df %>% subset( cluster == names(aa[i]) & score > aa[[i]], select = c("signature", "cluster", "score") )
    )
  }
  return(df2)
}

#-------------------------------------------------------------------------------

# plots heatmap of scores
plot.score.heatmap <- function(x, types, threshold) {
  
  
  df2 <- significant_signatures(x, types, threshold)
  
  #df <- compute.scores(x = x, threshold = threshold)
  #aa <- tapply(df$score, df$cluster, function(x) quantile(x, probs = c(0.9)))
  #df2 <- data.frame(matrix(ncol = 3, nrow = 0))
  #colnames(df2) <- c("signature", "cluster", "score")
  
  #for (i in 1:length(aa)) {
    #print(paste0( names(aa[i]), " - ", aa[[i]]) )
  #  df2 <- rbind(
  #    df2, 
  #    df %>% subset( cluster == names(aa[i]) & score > aa[[i]], select = c("signature", "cluster", "score") )
  #    )
  #}
  
  p <- ggplot(data = df2, aes(x = cluster, y = signature, fill = round(score, 3), label = round(score, 3))) +
    geom_tile(color = "white") + 
    scale_fill_gradient(low = "grey", high = "red") + 
    geom_text(color = "black", size = 3) +  # Add text annotations
    #scale_fill_gradient(low = "white", high = "steelblue") +  # Choose your desired color gradient
    labs(fill='home made score') + 
    theme_minimal() + 
    theme(
      # remove the vertical grid lines
      panel.grid.major.x = element_blank(), 
      # explicitly set the horizontal lines (or they will disappear too)
      panel.grid.major.y = element_line( linewidth=.1, color="black" ), 
      #legend.position="none"
    ) + 
    labs(title = "signatures with score > 0.1 in each cluster",
         x = "Clusters",
         y = "Signatures")
  
  return(p)
  
}

plot.score.heatmap(x, "DBS", 0.05)

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------













