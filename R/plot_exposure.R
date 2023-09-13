#' Plot exposure matrix
#'
#' @description Creates bar plot of relative exposure matrix, where x-axis are
#'              samples and y-axis are their relative contribution.
#' @param x Basilica object
#' @param sample_name Logical. If set to \code{TRUE}, the sample IDs will be
#'                    reported on the x axis.
#' @param levels add
#' @param cls add
#' @param flip_coord add
#'
#' @return plot
#' @export plot_exposures

plot_exposures = function(x, sample_name=F, sigs_levels=NULL, cls=NULL,
                          plot_noise=FALSE, add_centroid=FALSE, flip_coord=FALSE,
                          muts=FALSE, sampleIDs=NULL, sort_by=NULL, centroids=F) {

  if (!have_groups(x)) add_centroid = centroids = FALSE
  if (is.null(sampleIDs)) sampleIDs = rownames(get_exposure(x))

  titlee = "Exposure"
  b = get_exposure(x, add_groups=TRUE)
  if (plot_noise) {
    titlee = "Exposure noise"
    b = x$fit$params$alpha_noise
  }

  if (muts) b = b * rowSums(x$fit$x)

  b = b[sampleIDs,]

  if (is.null(cls) && !have_color_palette(x)) {
    cls = gen_palette(ncol(b)) %>% setNames(colnames(b))
  } else if (is.null(cls) && have_color_palette(x))
    cls = get_color_palette(x)

  if (is.null(sigs_levels))
    sigs_levels = sort(colnames(b %>% dplyr::select(-dplyr::contains("groups"))))

  idcols = c("sample")
  if (have_groups(x)) idcols = c("sample","groups")

  if (!is.null(sort_by))
    sample_levels = b %>% as.data.frame() %>%
      tibble::rownames_to_column(var="sample") %>%
      reshape2::melt(cols=idcols, variable.name="Signature", value.name="alpha") %>%
      dplyr::filter(Signature==sort_by) %>%
      dplyr::arrange(desc(alpha)) %>% dplyr::pull(sample) else sample_levels = rownames(b)

  b = b %>% as.data.frame() %>%
    tibble::rownames_to_column(var="sample")

  if (add_centroid || centroids) {
    a_pr = get_centroids(x, normalize=T)
    # if ( all(grepl("G", rownames(a_pr))) )
    #   grps = paste0("G", unique(x$groups) %>% stringr::str_replace_all("G","") %>% sort()) else
    #     grps = unique(x$groups) %>% stringr::str_replace_all("G","") %>% sort()
    grps = unique(x$groups)

    a_pr = a_pr[as.character(grps),]
    rownames(a_pr) = paste0("G", grps %>% stringr::str_replace_all("G",""))

    a_pr$groups = "Exposure priors"

    if (add_centroid) {
      b = rbind(b, a_pr %>% as.data.frame() %>% tibble::rownames_to_column(var="sample"))
      sample_levels = c(sample_levels, rownames(a_pr))

    } else if (centroids) {
      a_pr = get_centroids(x, normalize=T)
      a_pr$groups = "Centroids"
      b = a_pr %>% as.data.frame() %>% tibble::rownames_to_column(var="sample")
      sample_levels = rownames(a_pr)
      sample_name = T
    }
  }

  b = b %>%
    reshape2::melt(id=idcols, variable.name="Signature", value.name="alpha") %>%
    dplyr::mutate(sample=factor(sample, levels=sample_levels))

  if (add_centroid) {
    centrs = dplyr::select(a_pr, -groups)[,rev(sigs_levels)]

    b = b %>% dplyr::mutate(groups=stringr::str_replace_all(groups,"G","")) %>%
      dplyr::full_join(centrs %>%
                         tibble::rownames_to_column(var="groups") %>%
                         dplyr::mutate(groups=stringr::str_replace_all(groups,"G","")) %>%
                         reshape2::melt(variable.name="Signature", value.name="alpha_centr") %>%
                         dplyr::group_by(groups) %>%
                         dplyr::mutate(alpha_centr=cumsum(alpha_centr)),
                         by=c("groups","Signature"))
  }


  p = plot_exposures_aux(expos=b, cls=cls, titlee=titlee, sigs_level=sigs_levels,
                         sample_name=sample_name, flip_coor=flip_coord)

  return(p)
}


plot_exposures_aux = function(expos, cls=NULL, titlee="", sigs_levels=NULL,
                              sample_name=FALSE, flip_coord=FALSE) {

  p = expos %>%
    ggplot(aes(x=factor(sample), y=alpha, fill=Signature)) +
    geom_bar(stat="identity", position="stack") +
    labs(title=titlee) +
    theme_bw() + theme(axis.text.x=element_text(angle=90)) +
    guides(fill=guide_legend(title="Signatures")) + ylab("") + xlab("")

  if (!is.null(cls))
    p = p + ggplot2::scale_fill_manual(values=cls) + ggplot2::scale_color_manual(values=cls)

  if ("alpha_centr" %in% colnames(expos))
    p = p + geom_hline(aes(yintercept=alpha_centr, color=Signature),
                       linetype="dashed") + guides(color="none")

  if (!sample_name)
    p = p + theme(axis.ticks.x=element_blank(), axis.text.x=element_blank()) + labs(x="")

  if (flip_coord)
    p = p + coord_flip()

  if ("groups" %in% colnames(expos))
    p = p + facet_grid(~ groups, scales="free_x", space="free_x")

  return(p)
}



plot_sigs_prevalence = function(x) {
  expos = get_exposure(x, long = T)

  expos %>%
    dplyr::group_by(Signature) %>%
    dplyr::arrange(Signature, Exposure) %>%
    dplyr::mutate(nnn=1, cs=cumsum(nnn)/x$n_samples) %>%
    ggplot(aes(x=Exposure, y=cs, color=Signature)) +
    geom_point(size=.5) +
    geom_line() +
    scale_color_manual(values=x$color_palette) +
    theme_bw() + ylab("Fraction of samples") +
    xlim(0,1) + ylim(0,1)

}



