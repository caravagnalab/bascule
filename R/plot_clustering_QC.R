
# plots heatmap of scores
plot_cls_score_heatmap = function(x, type, exposure_thr=0.05) {

  scores = get_clusters_score(x, type, exposure_thr) %>% subset(significance==T)

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

  p = ggplot(data=scores,
             aes(x=cluster, y=signature, fill=round(score, 3), label=round(score, 3))) +
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


## SAME AS plot_mixture_weights(x, empirical=TRUE)
# # input : basilica object
# # output: ggplot (samples frequency in clusters)
# plot_cluster_freq = function(x) {
#   df = get_cluster_assignments(x) # dataframe
#   ggplot(df, aes(x=clusters)) +
#     geom_bar(stat="count", fill="purple") +
#     labs(title="Distribution of samples in clusters", x="Clusters", y="Number of Samples") +
#     theme_bw()
# }


plot_cluster_scores = function(x, type, cluster_label, final_score=TRUE, exposure_thr=0.05, quantile_thr=0.9) {

  if (!(cluster_label %in% get_cluster_labels(x))) warning("invalid cluster label!")

  scores = get_clusters_score(x, types=type, exposure_thr=exposure_thr, quantile_thr=quantile_thr) %>%
    subset(type==type & cluster==cluster_label)

  if (final_score==TRUE) {
    # FINAL SCORE
    p = ggplot(scores, aes(x=signature, y=score, color=cluster, group=cluster)) +
      geom_line() +
      geom_point() +
      geom_hline(aes(yintercept=score_quantile), linetype="dashed", color="black") +
      labs(title=paste0("Signatures Score in Cluster ", cluster_label),
           x="Signature", y="Score") +
      theme_bw() + 
      theme(axis.text.x=element_text(angle = 45, vjust = 1, hjust = 1), legend.position = "none") +
      guides(color=guide_legend(title="Cluster"))

    return(p)

  } else if (final_score==FALSE) {
    # PARTIAL SCORES
    scores_long = scores %>%
      dplyr::select(signature, cluster, varRatio, activeRatio, mutRatio) %>%
      tidyr::pivot_longer(cols=c(varRatio, activeRatio, mutRatio),
                          names_to="score_title", values_to="value")
    p = ggplot(subset(scores_long, cluster==cluster_label),
               aes(x=signature, y=value, color=score_title, group=score_title)) +
      geom_line() +
      geom_point() +
      labs(title=paste0("Signatures Score in Cluster ", cluster_label, " (Partials)"),
           x="Signature", y="Score") +
      theme_bw() + 
      theme(axis.text.x=element_text(angle=45, vjust = 1, hjust = 1)) + 
      scale_color_manual(labels = c("Active Samples", "Mutational Burden", "Variance Stability"), values = c("red", "green", "blue")) + 
      guides(color=guide_legend(title="Partial Scores (relative)"))

    return(p)

  } else {cli::cli_abort("Argument `final_score` can be either {TRUE} or {FALSE}")}

}




