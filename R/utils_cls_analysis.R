significant_signatures = function(x, types=get_types(x), exposure_thr=0.05) {
  if (!have_groups(x)) {
    cli::cli_alert_warning("The input data has not been clustered. Returning the input object.")
    return(x)
  }

  # TODO : check it -> for some high value of threshold may output error

  scores = get_clusters_score(x=x, types=types, threshold=exposure_thr)

  q0 = lapply(
    types,
    function(tid) tapply( subset(scores, type==tid)$score, subset(scores, type==tid)$cluster, function(x) quantile(x, probs=c(0.9)))
  ) %>% setNames(types)

  q1 = data.frame((sapply(q0, c)))
  q1["cluster"] = q1 %>% rownames
  q2 = q1 %>% tidyr::pivot_longer(cols=dplyr::all_of(types), names_to="type", values_to="quantile")

  df = scores %>% dplyr::right_join(q2, by=c("cluster","type"))
  # merge(df, q, by.x=c('cluster', 'type'), by.y=c('cluster', 'type'))

  return(
    df %>% dplyr::filter(score >= quantile) %>%
      dplyr::select(signature, cluster, score, type)
  )
}


get_clusters_score_aux = function(x, type, threshold) {
  exposures = get_exposure(x, types=type, matrix=FALSE, add_groups=TRUE)[[type]] %>% subset(value > threshold)
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

  return(df)
}


get_clusters_score = function(x, types=get_types(x), threshold=0.05) {
  return(
    lapply(
      types,
      function(tid) get_clusters_score_aux(x=x, type=tid, threshold=threshold) %>%
        dplyr::mutate(type=tid)) %>% do.call(rbind, .) %>% dplyr::filter(!is.na(type))
  )
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


