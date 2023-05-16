#' Function to visualize the samples mutation profiles.
#'
#' @description creates bar plot of data profile for each sample,
#' where x-axis are 96 substitution bases and y-axis are their relative contribution.
#'
#' @param x Basilica object.
#' @param context Logical. If set to \code{TRUE}, the contexts on the x axis will be shown.
#' @param what Character among \code{["SBS","DBS","ID","CNV"]} reporting the type of signature.
#' @param sample_ids Array with the samples IDs to show.
#'
#' @return plot
#' @export plot_data

plot_data = function(x, sample_ids = NULL, what = "SBS", context = T) {

  a = x$input$counts

  a = a %>% dplyr::mutate(sbs=rownames(a)) %>% as_tibble() %>%
    reshape2::melt() %>%
    reformat_contexts(what=what)

  # library(CNAqc)

  if(!is.null(sample_ids)) a = a %>% filter(Var1 %in% samples_ids)

  all_plot = lapply(a$Var1 %>% unique(), function(var) {
    p = a %>% filter(Var1 == var) %>%
      ggplot() +
      geom_bar(aes(value, x = context, fill = "Indianred"), stat = 'identity') +
      facet_grid(~factor(substitution,levels = a$substitution %>% unique()), scales = 'free') +
      # CNAqc:::my_ggplot_theme()  +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 90)) +
      guides(fill = 'none')  +
      labs(y = "", title = var)

    if(!context)
      p = p + theme(axis.ticks.x = element_blank(), axis.text.x = element_blank()) + labs(x = "")

    p

  })

  all_plot
}


get_data = function(x, reconstructed=FALSE) {
  if(reconstructed)
    return((as.matrix(rowSums(x$fit$x) * get_exposure(x)) %*% as.matrix(get_signatures(x))) %>% as.data.frame())
  return(x$fit$x)
}



#' Function to plot the overall mutation profile
#'
#' @param x Basilica object.
#' @param sampleIDs List containing the sample IDs to consider in the counts.
#'                  If set to \code{NULL}, all samples will be considered.
#' @param by_sig Logical, if set to \code{TRUE} each context will report
#'               its signatures contribution.
#' @param reconstructed whether to plot the reconstructed or input matrix
#'
#' @return ggplot object
#' @export plot_mutations

plot_mutations = function(x, sampleIDs = NULL, by_sig = FALSE, reconstructed = FALSE) {
  if (is.null(sampleIDs)) sampleIDs = rownames(x$fit$x)

  if (by_sig) {
    cli::cli_alert_warning("Not implemented with mutations count by signature.
                           {.code by_sig} set to {.code FALSE}.")
    by_sig = FALSE
  }

  xx = get_data(x, reconstructed=reconstructed)
  xx_s = xx %>% tibble::rownames_to_column(var="sampleID") %>% dplyr::mutate(sig="")

  if (by_sig)
    xx_s = lapply(rownames(get_signatures(x)),
                  function(sname) {
                   ((as.matrix(x$fit$exposure[,sname], ncol=1) %*%
                       as.matrix(get_signatures(x)[sname,], nrow=1)) *
                      xx) %>%
                   dplyr::mutate(sig=sname) %>%
                   tibble::rownames_to_column(var="sampleID")
                  }
    ) %>% do.call(what=rbind, args=.)

  xx_s = xx_s %>%
    dplyr::filter(sampleID %in% sampleIDs) %>%
    reshape2::melt(variable.name="subs", value.name="mut_count") %>%
    dplyr::filter(mut_count>0) %>%
    tidyr::separate("subs", into=c("n1","subs.n2"), sep="[[]") %>%
    tidyr::separate("subs.n2", into=c("subs","n2"), sep="]") %>%
    dplyr::mutate(context=paste0(n1,"_",n2)) %>%

    dplyr::group_by(sig, context, subs) %>%
    dplyr::summarise(tot_muts=sum(mut_count)) %>%
    dplyr::ungroup()

  if (by_sig)
    return(
      xx_s %>%
        ggplot() +
        geom_histogram(aes(x=context, y=tot_muts, fill=sig), stat="identity", position="stack") +
        facet_grid(~subs) +
        ylab("Number of mutations") + xlab("") +
        theme_bw() + theme(axis.text.x=element_text(angle=90)) +
        scale_fill_manual(values=get_signature_colors(x))
    )

  return(
    xx_s %>%
      ggplot() +
      geom_histogram(aes(x=context, y=tot_muts), stat="identity") +
      facet_grid(~subs) +
      ylab("Number of mutations") + xlab("") +
      theme_bw() + theme(axis.text.x=element_text(angle=90))
  )

}




reformat_contexts = function(a, what) {
  if(what == 'SBS') {
    a = a %>% dplyr::rename(Var1 = sbs, Var2 = variable) %>%
      mutate(substitution=paste0(substr(start=3, stop=3, Var2), ">", substr(start=5, stop=5, Var2)),
             context=paste0(substr(start=1, stop=1, Var2), '_', substr(start=7, stop=7, Var2)))
  }

  if(what == "DBS") {
    a = a %>% dplyr::rename(Var1=sbs, Var2=variable) %>%
      mutate(substitution=paste0(substr(start=1, stop=2, Var2),">NN"),
             context=substr(start=4, stop=5, Var2))
  }

  if(what == "ID") {
    a = a %>% dplyr::rename(Var1=sbs, Var2=variable) %>%
      mutate(Var2=as.character(Var2),
             substitution=substr(start=1, stop=nchar(Var2) - 2, Var2),
             context=substr(start=nchar(Var2), stop=nchar(Var2), Var2))
  }

  if(what == "CNV") {
    a = a %>% dplyr::rename(Var1=sbs, Var2=variable) %>% rowwise() %>%
      mutate(Var2=as.character(Var2),
             substitution=paste0(str_split(Var2, pattern=":")[[1]][1], ":",
                                 str_split(Var2, pattern=":")[[1]][2]),
             context= paste0(str_split(Var2,pattern=":")[[1]][3]))
  }

  return(a)
}

