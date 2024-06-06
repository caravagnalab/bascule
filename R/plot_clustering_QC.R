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
    scale_fill_gradient(low="grey", high="red") +
    geom_text(color="black", size=3) +  # Add text annotations
    # scale_fill_gradient(low="white", high="steelblue") +  # Choose your desired color gradient
    labs(fill='clustering score') +
    theme_minimal() +
    theme(
      # remove the vertical grid lines
      panel.grid.major.x=element_blank(),
      # explicitly set the horizontal lines (or they will disappear too)
      panel.grid.major.y=element_line( linewidth=.1, color="black" ),
      #legend.position="none"
    ) +
    labs(title="significant signatures in each cluster",
         x="Clusters",
         y="Signatures")

  return(p)
}


# input : basilica object
# output: ggplot (samples frequency in clusters)
plot_cluster_freq = function(x) {
  df = get_cluster_assignments(x) # dataframe
  p=ggplot(df, aes(x=clusters)) +
    geom_bar(stat="count", fill="purple") +
    labs(x="Clusters", y="Number of samples", title="Distribution of samples in clusters") +
    theme_minimal()
  return(p)
}


# plot the final score for a cluster and single type
# a : scores data-frame
# b : single context type
# c : cluster name
plot_scores_aux1_1 = function(a, b, c) {

  # filter on context type and cluster name
  df = a %>% subset(type == b & cluster == c)

  # find the quantile value
  # val = df %>% dplyr::pull(score) %>% quantile(probs = c(0.9))
  # df1 =df[c(1,2,6)]

  p = ggplot(df, aes(x=signature, y=score, color=cluster, group=cluster)) +
    geom_line() +
    geom_point() +
    geom_hline(yintercept=df$score_quantile, linetype="dashed", color="black") +
    labs(x="Mutational Signature", y="Value", title=paste0("Signatures Score (Final) in Cluster: ", b)) +
    theme(axis.text.x=element_text(angle=90)) +
    theme_bw()

  return(p)
}


# plot the final score for a cluster and all context types
# a : scores data-frame
# b : all context types
# c : cluster name
plot_scores_aux1 = function(a, b, c) {

  plot_list = lapply(b, function(tid) plot_scores_aux1_1(a, tid, c))

  p = patchwork::wrap_plots(
    plot_list,
    nrow=length(plot_list),
    ncol=1,
    guides="collect",
    tag_level="keep",
    design="AAAA\nBBBB"
  ) & theme(axis.text.x=element_text(angle=90))

  return(p)
}


# plot the scores ingredients for single cluster
# a : scores data-frame
# b : single context type
# c : cluster name
plot_scores_aux2_1 = function(a, b, c) {

  df = a %>% subset(type == b & cluster == c)

  df_long = tidyr::pivot_longer(df[c(1,2,3,4,5)],
                                cols=c(varRatio, activeRatio, mutRatio),
                                names_to="name", values_to="value")
  p = ggplot(subset(df_long, cluster == c),
             aes(x=signature, y=value, color=name, group=name)) +
    geom_line() +
    geom_point() +
    labs(x="Mutational Signature", y="Value", title=paste0("Signatures Score (Ingredients) in Cluster: ", b)) +
    theme_minimal()
  return(p)
}


plot_scores_aux2 = function(a, b, c) {

  plot_list = lapply(b, function(tid) plot_scores_aux2_1(a, tid, c))

  p = patchwork::wrap_plots(
    plot_list,
    nrow=length(plot_list),
    ncol=1,
    guides="collect",
    tag_level="keep",
    design="AAAA\nBBBB",
    axes="collect"
  )

  return(p)
}


plot_cluster_scores = function(x, clusters=get_cluster_labels(x),
                               types=get_types(x), exposure_thr=0.05,
                               quantile_thr=0.9) {

  a = get_clusters_score(x, types, exposure_thr, quantile_thr)

  lapply(clusters, function(cluster_name) {
    p1 = plot_scores_aux1(a=a, b=types, c=cluster_name)
    p2 = plot_scores_aux2(a=a, b=types, c=cluster_name)
    p3 = plot_exposures(x, types, clusters=cluster_name, elim_thr=exposure_thr)

    return( (p1 | p2) / p3 )
  }) %>% setNames(clusters)
}

