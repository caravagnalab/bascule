
# plots heatmap of scores
plot_cls_score_heatmap <- function(x, type, exposure_thr=0.05) {

  df = get_clusters_score(x, type, exposure_thr) %>% subset(significance == T)

  # df2 = significant_signatures(x, types, threshold)

  # df = compute.scores(x=x, threshold=threshold)
  # aa = tapply(df$score, df$cluster, function(x) quantile(x, probs=c(0.9)))
  # df2 = data.frame(matrix(ncol=3, nrow=0))
  # colnames(df2) = c("signature", "cluster", "score")

  # for (i in 1:length(aa)) {
  # print(paste0( names(aa[i]), " - ", aa[[i]]) )
  #   df2 = rbind(
  #     df2,
  #     df %>% subset( cluster == names(aa[i]) & score > aa[[i]], select=c("signature", "cluster", "score") )
  #     )
  # }

  p = ggplot(data=df, aes(x=cluster, y=signature, fill=round(score, 3), label=round(score, 3))) +
    geom_tile(color="white") +
    scale_fill_gradient(low="grey", high="dodgerblue1") +
    geom_text(color="black", size=3) +  # Add text annotations
    # scale_fill_gradient(low="white", high="steelblue") +  # Choose your desired color gradient
    labs(
      title="Significant signatures in clusters", 
      subtitle = "(based on clustering scores)", 
      x="Clusters",
      y="Signatures", 
      fill = "Score") + 
    theme_minimal() +
    theme(
      # remove the vertical grid lines
      panel.grid.major.x=element_blank(),
      # explicitly set the horizontal lines (or they will disappear too)
      panel.grid.major.y=element_line( linewidth=.1, color="black" ),
      #legend.position="none"
    ) + 
    theme_bw()

  return(p)
}

#===============================================================================

# input : basilica object
# output: ggplot (samples frequency in clusters)
plot_cluster_freq = function(x) {
  df = get_cluster_assignments(x) # dataframe
  p=ggplot(df, aes(x=clusters)) +
    geom_bar(stat="count", fill="purple") +
    labs(title="Distribution of samples in clusters", x="Clusters", y="Number of Samples") +
    theme_bw()
  return(p)
}

#===============================================================================


plot_cluster_scores <- function(x, context_type, cluster_label, final_score=TRUE) {
  
  if ( !(cluster_label %in% get_cluster_labels(x)) ) {
    warning("invalid cluster label!")
  }
  
  scores <- get_clusters_score(x, types = context_type, exposure_thr = 0.05, quantile_thr = 0.9)
  df = scores %>% subset(type == context_type & cluster == cluster_label)
  
  if (final_score) {
    # FINAL SCORE
    p = ggplot(df, aes(x=signature, y=score, color=cluster, group=cluster)) +
      geom_line() +
      geom_point() +
      geom_hline(yintercept=df$score_quantile, linetype="dashed", color="black") +
      labs(title=paste0("Signatures Score (Final) in Cluster: ", b), x="Mutational Signature", y="Score") +
      theme(axis.text.x=element_text(angle=90)) +
      theme_bw() + guides(color = guide_legend(title = "Clustering Scores"))
  } else if (!final_score) {
    # PARTIAL SCORES
    df_long = tidyr::pivot_longer(df[c(1,2,3,4,5)],
                                  cols=c(varRatio, activeRatio, mutRatio),
                                  names_to="score_title", values_to="value")
    p <- ggplot(subset(df_long, cluster == cluster_label),
                aes(x=signature, y=value, color=score_title, group=score_title)) +
      geom_line() +
      geom_point() +
      labs(
        title=paste0("Signatures Score (Ingredients) in Cluster: ", cluster_label), 
        x="Mutational Signature", 
        y="Value"
      ) +
      theme_bw()
  } else {warning("invalid final_score!")}
  
  return(p)
}





