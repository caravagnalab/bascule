# Merge similar clusters #####
#' Function to merge similar clusters
#'
#' @description
#' This function will iteratively merge clusters that result similar in the centroid.
#' The merging will stop as soon as the cosine similarity between all pairs of clusters
#' is below `cutoff`
#'
#' @param x bascule object.
#' @param cutoff Minimum value of similarity to merge two clusters.
#'
#' @return bascule object.
#' @export merge_clusters
merge_clusters = function(x, cutoff=0.8) {
  if (!have_groups(x)) return(x)

  repeat {
    alpha_prior = empirical_centroids(x)
    if (nrow(alpha_prior) == 1) return(x)
    cosine_simil = lsa::cosine(t(alpha_prior)) %>% as.data.frame()
    cosine_simil[lower.tri(cosine_simil, diag=T)] = 0

    rownames(cosine_simil) = colnames(cosine_simil) = rownames(alpha_prior)

    merging = cosine_simil %>% tibble::rownames_to_column(var="gid1") %>%
      reshape2::melt(id="gid1", variable.name="gid2", value.name="cosine") %>%
      dplyr::mutate(gid1=as.character(gid1), gid2=as.character(gid2)) %>%
      dplyr::filter(cosine > cutoff) %>% dplyr::select(-cosine) %>%
      dplyr::group_by(gid1) %>%
      dplyr::summarise(cl_old=list(c(gid1,gid2) %>% unique())) %>% dplyr::ungroup() %>%
      dplyr::rename(cl_name=gid1)

    if (nrow(merging) == 0) return(x)

    grps = get_cluster_assignments(x)
    for (i in 1:nrow(merging)) {
      old_cl = dplyr::pull(merging[i,], cl_old)[[1]]
      new_cl = dplyr::pull(merging[i,], cl_name)
      grps = grps %>% dplyr::mutate(clusters=ifelse(clusters %in% old_cl,
                                                    new_cl, clusters))
    }
    x$clustering$centroids = alpha_prior %>%
      tibble::rownames_to_column(var="clusters") %>%
      reshape2::melt(id="clusters", variable.name="sigs", value.name="value") %>%
      tibble::as_tibble()
    x$clustering$clusters = grps
  }

  return(x)
}


empirical_centroids = function(x) {
  lapply(get_types(x), function(tid) {
    expos = get_exposure(x, add_groups=TRUE)[[tid]]
    expos %>% dplyr::group_by(clusters, sigs) %>%
      dplyr::reframe(centroid=mean(value)) %>%
      dplyr::ungroup() %>%
      tidyr::pivot_wider(values_from="centroid", names_from="sigs") %>%
      tibble::column_to_rownames(var="clusters")
  }) %>% do.call(cbind, .)
}


# Post-processing functions #####
compute_clustering_scores = function(x, types=get_types(x), threshold=0.05) {
  return(
    lapply(types, function(tid) {
      if (!tid %in% get_types(x))
        cli::cli_alert_warning("Variant type {tid} not present in the object. Available types are {paste0(get_types(x), collapse=', ')}")
      else
        compute_clustering_scores_aux(x=x, type=type, threshold=threshold)
    }) %>% do.call(rbind, .)
  )
}


compute_clustering_scores_aux = function(x, type, threshold) {
  exposures = get_exposure(x, types=type, matrix=FALSE, add_groups=TRUE)[[1]] %>% subset(value > threshold)

  df = data.frame(signature=c(), cluster=c(), varRatio=c(), activeRatio=c(), mutRatio=c(), score=c())

  for (cls in get_cluster_labels(x)) {
    for (signature in (exposures %>% subset(clusters == cls) %>% dplyr::select(sigs) %>% unique() %>% dplyr::pull(sigs))) {
      # signature exposure variance in one cluster / signature exposure variance in all clusters
      onecluster_var = var(
        exposures %>% subset(sigs == signature & clusters == cls) %>% dplyr::pull(value) # numeric
      )
      allclusters_var = var(
        exposures %>% subset(sigs == signature) %>% dplyr::pull(value) # numeric
      )
      # ratio_var = allclusters_var / onecluster_var

      if (is.na(onecluster_var) | is.na(allclusters_var)) {
        ratio_var = 0
      } else {
        ratio_var = (1 / ( 1 + exp( -(log(allclusters_var / onecluster_var)) ) ) ) # sigmoid(-log(ratio))
      }

      # samples with active signature / all samples
      num_one = exposures %>% subset(clusters == cls & sigs == signature, select=c("samples")) %>% unique() %>% nrow
      num_all = exposures %>% subset(clusters == cls, select=c("samples")) %>% unique() %>% nrow
      ratio_active = num_one / num_all

      # signature related mutations / all mutations
      input = get_input(
        x,
        types="SBS",
        samples=exposures %>% subset(clusters == cls & sigs == signature) %>% dplyr::pull(samples) %>% unique(),
        clusters=cls,
        matrix=TRUE,
        reconstructed=FALSE,
        add_groups=TRUE
      )[[1]]

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

      #scaled_df = as.data.frame( lapply( df[c(3,4,5,6)], function(x) (x - min(x)) / (max(x) - min(x)) ) )
      #scaled_df = as.data.frame( lapply( df[c(3,4,5,6)], function(x) (x - mean(x)) / sd(x) ) )
      #df2 = cbind(df[c(1,2)], scaled_df) #%>% mutate(value=var_value * num_value * mut_value)
    }
  }
  return(df)
}
