

#-------------------------------------------------------------------------------

compute.scores.aux <- function(x, type, threshold) {
  
  exposures <- basilica:::get_exposure(x, types = type, matrix = FALSE, add_groups = TRUE)[[1]] %>% subset(value > threshold)
  
  df <- data.frame(signature = c(), cluster = c(), varRatio = c(), activeRatio = c(), mutRatio = c(), score = c())
  
  for (cls in basilica:::get_cluster_labels(x)) {
    
    for ( signature in ( exposures %>% subset(clusters == cls) %>% dplyr::select(sigs) %>% unique() %>% dplyr::pull(sigs) ) ) {
      
      #-------------------------------------------------------------------------
      # signature exposure variance in one cluster / signature exposure variance in all clusters
      onecluster_var <- var(
        exposures %>% subset(sigs == signature & clusters == cls) %>% dplyr::pull(value) # numeric
      )
      allclusters_var <- var(
        exposures %>% subset(sigs == signature) %>% dplyr::pull(value) # numeric
      )
      # ratio_var <- allclusters_var / onecluster_var
      
      if ( is.na(onecluster_var) | is.na(allclusters_var) ) {
        ratio_var <- 0
      } else {
        ratio_var <- (1 / ( 1 + exp( -(log(allclusters_var / onecluster_var)) ) ) ) # sigmoid(-log(ratio))
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
        types = "SBS", 
        samples = exposures %>% subset(clusters == cls & sigs == signature) %>% dplyr::pull(samples) %>% unique(), 
        clusters = cls, 
        matrix = TRUE, 
        reconstructed = FALSE, 
        add_groups = TRUE
      )[[1]]
      
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
      
      #scaled_df <- as.data.frame( lapply( df[c(3,4,5,6)], function(x) (x - min(x)) / (max(x) - min(x)) ) )
      #scaled_df <- as.data.frame( lapply( df[c(3,4,5,6)], function(x) (x - mean(x)) / sd(x) ) )
      #df2 <- cbind(df[c(1,2)], scaled_df) #%>% mutate(value = var_value * num_value * mut_value)
      
    }
    
  }
  return(df)
}

#-------------------------------------------------------------------------------


#' Title
#'
#' @param x 
#' @param types 
#' @param threshold 
#'
#' @return
#' @export
#'
#' @examples
compute.scores = function(x, types=get_types(x), threshold=0.05) {
  if ( "SBS" %in% types ) {
    a <- compute.scores.aux(x = x, type = c("SBS"), threshold = 0.05)
  } else a <- NULL
  if ( "DBS" %in% types ) {
    b <- compute.scores.aux(x = x, type = c("DBS"), threshold = 0.05)
  } else b <- NULL
  
  xx = rbind(a, b) %>%
    dplyr::mutate(type=match_type(types, signature)) %>%
    dplyr::filter(!is.na(type))# %>%
    #dplyr::rename(samples=cluster)
  
  
  return( xx )
}





