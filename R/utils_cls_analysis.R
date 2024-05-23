


get_clusters_score_aux = function(x, type, exposure_thr, quantile_thr) {
  exposures = get_exposure(x, types=type, matrix=FALSE, add_groups=TRUE)[[type]] %>% subset(value > exposure_thr)
  df = data.frame(signature=c(), cluster=c(), varRatio=c(), activeRatio=c(), mutRatio=c(), score=c())

  for (cls in get_cluster_labels(x)) {
    for ( signature in ( exposures %>% subset(clusters == cls) %>% dplyr::pull(sigs) %>% unique ) ) {
      # specified signature exposure variance in one cluster / specified signature exposure variance in all clusters
      one_cluster_var = var(
        exposures %>% subset(sigs == signature & clusters == cls) %>% dplyr::pull(value) # numeric
      )
      all_clusters_var = var(
        exposures %>% subset(sigs == signature) %>% dplyr::pull(value) # numeric
      )

      if ( is.na(one_cluster_var) | is.na(all_clusters_var) ) {
        ratio_var = 0
      } else {
        ratio_var = (1 / ( 1 + exp( -(log(all_clusters_var / one_cluster_var)) ) ) ) # sigmoid(-log(ratio))
      }

      # samples with active signature / all samples
      num_one = exposures %>% subset(clusters == cls & sigs == signature, select=c("samples")) %>% unique() %>% nrow()
      num_all = exposures %>% subset(clusters == cls, select=c("samples")) %>% unique() %>% nrow()
      ratio_active = num_one / num_all

      # signature related mutations / all mutations
      input = get_input(
        x,
        types=type,
        samples=exposures %>% subset(clusters == cls & sigs == signature) %>% dplyr::pull(samples) %>% unique(),
        clusters=cls,
        matrix=TRUE,
        reconstructed=FALSE,
        add_groups=TRUE
      )[[type]]

      if (is.null(input)) {
        ratio_mut = 0
      } else {
        mut_all = sum(unlist(input)) # floor(rowSums(input) %>% sum)
        mut_one = (exposures %>% subset(clusters == cls & sigs == signature) %>%
                      dplyr::pull(value)) * rowSums(input)
        names(mut_one) = NULL
        mut_one = floor(mut_one %>% sum)
        ratio_mut = mut_one / mut_all
      }

      df = rbind(
        df, list(
          signature=signature,
          cluster=cls,
          varRatio=ratio_var,
          activeRatio=ratio_active,
          mutRatio=ratio_mut,
          score=ratio_var * ratio_active * ratio_mut
        )
      )
    }
  }
  
  df1 <- lapply(
    get_cluster_labels(x), # list of cluster labels
    function(cid) {
      df %>% 
        subset(cluster == cid) %>% 
        dplyr::mutate(
          score_quantile = df %>% 
            subset(cluster == cid) %>% 
            dplyr::pull(score) %>% 
            quantile(probs = c(quantile_thr))
        )
    }
  ) %>% do.call(rbind, .)
  
  df1 <- df1 %>% dplyr::mutate(significance = score >= score_quantile)

  return(df1)
}


get_clusters_score = function(x, types=get_types(x), exposure_thr=0.05, quantile_thr=0.9) {
  return(
    lapply(
      types,
      function(tid) get_clusters_score_aux(x, tid, exposure_thr, quantile_thr) %>%
        dplyr::mutate(type=tid)) %>% do.call(rbind, .) %>% dplyr::filter(!is.na(type))
  )
  
  #df <- lapply(
  #  cluster_labels, # list of cluster labels
  #  function(cid) {
  #    lapply(
  #      types, # list of context types
  #      function(tid) {
  #        scores %>% 
  #          subset(cluster == cid & type == tid) %>% 
  #          mutate(
  #            quantileValue = scores %>% 
  #              subset(cluster == cid & type == tid) %>% 
  #              dplyr::pull(score) %>% 
  #              quantile(probs = c(0.9))
  #          )
  #      }
  #    ) %>% do.call(rbind, .)
  #  }
  #) %>% do.call(rbind, .)
  
}


#significant.signatures = function(x, types, threshold) {
#
#  df = get_clusters_score(x, types, threshold)
#
#  q = tapply(df$score, df$cluster, function(x) quantile(x, probs=c(0.9)))
#
#  xx = data.frame(matrix(ncol=3, nrow=0))
#  colnames(xx) = c("signature", "cluster", "score")
#
#  for (i in 1:length(q)) {
#    xx = rbind(
#      xx,
#      df %>% subset( cluster == names(q[i]) & score > q[[i]], select=c("signature", "cluster", "score") )
#    )
#  }
#  rownames(xx)=seq(length=nrow(xx))
#
#  return(xx)
#}


# ==============================================================================
# ==============================================================================
# VISUALIZAION
# ==============================================================================
# ==============================================================================

# input : basilica object
# output: ggplot (samples frequency in clusters)
plot_cluster_freq <- function(x) {
  df <- basilica:::get_cluster_assignments(x) # dataframe
  p <- ggplot(df, aes(x = clusters)) + 
    geom_bar(stat = "count", fill = "purple") +
    labs(x = "Clusters", y = "Number of samples", title = "Distribution of samples in clusters") + 
    theme_minimal()
  return(p)
}

#-------------------------------------------------------------------------------

# plot the final score for a cluster and single type
# a : scores data-frame
# b : single context type
# c : cluster name
plot_scores_aux1_1 <- function(a, b, c) {
  
  # filter on context type and cluster name
  df <- a %>% subset(type == b & cluster == c)
  
  # find the quantile value
  #val <- df %>% dplyr::pull(score) %>% quantile(probs = c(0.9))
  #df1 <-df[c(1,2,6)]
  
  p <- ggplot(df, aes(x = signature, y = score, color = cluster, group = cluster)) + 
    geom_line() +
    geom_point() +
    labs(x = "Mutational Signature", y = "Value", title = paste0("Signatures Score (Final) in Cluster: ", b)) + 
    geom_hline(yintercept = df$score_quantile, linetype = "dashed", color = "black") + 
    theme_minimal()
  
  return(p)
}

#-------------------------------------------------------------------------------

# plot the final score for a cluster and all context types
# a : scores data-frame
# b : all context types
# c : cluster name
plot_scores_aux1 <- function(a, b, c) {
  
  plot_list <- lapply(b, function(tid) plot_scores_aux1_1(a, tid, c))
  
  p <- patchwork::wrap_plots(
    plot_list, 
    nrow = length(plot_list), 
    ncol = 1, 
    guides = "collect", 
    tag_level = "keep", 
    design = "AAAA
              BBBB",
    axes = "collect"
  )
  
  return(p)
}

#-------------------------------------------------------------------------------

# plot the scores ingredients for single cluster
# a : scores data-frame
# b : single context type
# c : cluster name
plot_scores_aux2_1 <- function(a, b, c) {
  
  df <- a %>% subset(type == b & cluster == c)
  
  df_long <- tidyr::pivot_longer(df[c(1,2,3,4,5)], cols = c(varRatio, activeRatio, mutRatio), names_to = "name", values_to = "value")
  p <- ggplot(subset(df_long, cluster == c), aes(x = signature, y = value, color = name, group = name)) + 
    geom_line() + 
    geom_point() + 
    labs(x = "Mutational Signature", y = "Value", title = paste0("Signatures Score (Ingredients) in Cluster: ", b)) + 
    theme_minimal()
  return(p)
}

#-------------------------------------------------------------------------------

plot_scores_aux2 <- function(a, b, c) {
  
  plot_list <- lapply(b, function(tid) plot_scores_aux2_1(a, tid, c))
  
  p <- patchwork::wrap_plots(
    plot_list, 
    nrow = length(plot_list), 
    ncol = 1, 
    guides = "collect", 
    tag_level = "keep", 
    design = "AAAA
              BBBB",
    axes = "collect"
  )
  
  return(p)
}

#-------------------------------------------------------------------------------

plot_cluster_scores <- function(x, types=get_types(x), clusterName, exposure_thr = 0.05, quantile_thr = 0.9) {
  
  a <- basilica:::get_clusters_score(x, types, exposure_thr, quantile_thr)
  
  p1 <- plot_scores_aux1(a, types, clusterName)
  p2 <- plot_scores_aux2(a, types, clusterName)
  
  p3 <- basilica:::plot_exposures(x, types, clusters = clusterName, elim_thr = exposure_thr)
  
  #p3 <- basilica::plot_exposures(x = x, types = stypes, clusters = clusterName)
  
  #p3 <- plot.cluster.exposure(x, types, clusterName, threshold)
  
  return( (p1 | p2) / p3 )
}


