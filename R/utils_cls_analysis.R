

get_clusters_score_aux <- function(x, type, threshold) {
  
  exposures <- basilica:::get_exposure(x, types = type, matrix = FALSE, add_groups = TRUE)[[type]] %>% subset(value > threshold)
  
  df <- data.frame(signature = c(), cluster = c(), varRatio = c(), activeRatio = c(), mutRatio = c(), score = c())
  
  for (cls in basilica:::get_cluster_labels(x)) {
    
    for ( signature in ( exposures %>% subset(clusters == cls) %>% dplyr::pull(sigs) %>% unique ) ) {
      
      #-------------------------------------------------------------------------
      # specified signature exposure variance in one cluster / specified signature exposure variance in all clusters
      one_cluster_var <- var(
        exposures %>% subset(sigs == signature & clusters == cls) %>% dplyr::pull(value) # numeric
      )
      all_clusters_var <- var(
        exposures %>% subset(sigs == signature) %>% dplyr::pull(value) # numeric
      )
      # ratio_var <- all_clusters_var / one_cluster_var
      
      if ( is.na(one_cluster_var) | is.na(all_clusters_var) ) {
        ratio_var <- 0
      } else {
        ratio_var <- (1 / ( 1 + exp( -(log(all_clusters_var / one_cluster_var)) ) ) ) # sigmoid(-log(ratio))
      }
      
      #-------------------------------------------------------------------------
      # samples with active signature / all samples
      num_one <- exposures %>% subset(clusters == cls & sigs == signature, select=c("samples")) %>% unique() %>% nrow
      num_all <- exposures %>% subset(clusters == cls, select=c("samples")) %>% unique() %>% nrow
      ratio_active <- num_one / num_all
      
      #-------------------------------------------------------------------------
      # signature related mutations / all mutations
      input <- basilica:::get_input(
        x, 
        types = type, 
        samples = exposures %>% subset(clusters == cls & sigs == signature) %>% dplyr::pull(samples) %>% unique(), 
        clusters = cls, 
        matrix = TRUE, 
        reconstructed = FALSE, 
        add_groups = TRUE
      )[[type]]
      
      if (is.null(input)) {
        ratio_mut <- 0
      } else {
        mut_all <- sum(unlist(input)) # floor(rowSums(input) %>% sum)
        mut_one <- (exposures %>% subset(clusters == cls & sigs == signature) %>% 
                      dplyr::pull(value)) * rowSums(input)
        names(mut_one) <- NULL
        mut_one <- floor(mut_one %>% sum)
        ratio_mut <- mut_one / mut_all
      }
      
      #-------------------------------------------------------------------------
      df <- rbind(
        df, list(
          signature = signature, 
          cluster = cls, 
          varRatio = ratio_var, 
          activeRatio = ratio_active, 
          mutRatio = ratio_mut, 
          score = ratio_var * ratio_active * ratio_mut
        )
      )
      
    }
    
  }
  return(df)
}

#-------------------------------------------------------------------------------


get_clusters_score <- function(x, types=get_types(x), threshold=0.05) {
  #if ( "SBS" %in% types ) {
  #  a <- compute.scores.aux(x = x, type = c("SBS"), threshold = 0.05)
  #} else a <- NULL
  #if ( "DBS" %in% types ) {
  #  b <- compute.scores.aux(x = x, type = c("DBS"), threshold = 0.05)
  #} else b <- NULL
  
  xx = lapply(
    types, 
    function(tid) get_clusters_score_aux(x = x, type = tid, threshold = threshold) %>% 
      dplyr::mutate(type=tid)) %>% do.call(rbind, .) %>% dplyr::filter(!is.na(type))
  
  #xx = rbind(a, b) %>%
  #  dplyr::mutate(type=match_type(types, signature)) %>%
  #  dplyr::filter(!is.na(type))# %>%
  #dplyr::rename(samples=cluster)
  
  return( xx )
}

#-------------------------------------------------------------------------------


significant.signatures <- function(x, types=get_types(x), threshold=0.05) {
  
  # TO-DO : check it -> for some high value of threshold may output error
  
  scores <- basilica:::get_clusters_score(x=x, types=types, threshold=threshold)
  
  q0 <- lapply(
    types, 
    function(tid) tapply( subset(scores, type==tid)$score, subset(scores, type==tid)$cluster, function(x) quantile(x, probs = c(0.9)))
  )
  names(q0) <- types
  
  q1 <- data.frame((sapply(q0, c)))
  q1["cluster"] <- q1 %>% rownames
  q2 <- q1 %>% tidyr::pivot_longer(cols=types, names_to='type', values_to='quantile')
  
  df <- scores %>% dplyr::right_join(q2, by=c("cluster","type"))
  #merge(df, q, by.x=c('cluster', 'type'), by.y=c('cluster', 'type'))
  
  xx <- df %>% subset(score >= quantile, select = c("signature", "cluster", "score", "type"))
  rownames(xx) = seq(length = nrow(xx))
  
  return(xx)
}


#significant.signatures <- function(x, types, threshold) {
#  
#  df <- get_clusters_score(x, types, threshold)
#  
#  q <- tapply(df$score, df$cluster, function(x) quantile(x, probs = c(0.9)))
#  
#  xx <- data.frame(matrix(ncol = 3, nrow = 0))
#  colnames(xx) <- c("signature", "cluster", "score")
#  
#  for (i in 1:length(q)) {
#    xx <- rbind(
#      xx, 
#      df %>% subset( cluster == names(q[i]) & score > q[[i]], select = c("signature", "cluster", "score") )
#    )
#  }
#  rownames(xx) = seq(length = nrow(xx))
#  
#  return(xx)
#}



# ==============================================================================
# VISUALIZAION
# ==============================================================================

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


