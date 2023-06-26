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
#' @import ggplot2
#'
#' @return plot
#' @export plot_data

plot_data = function(x, sample_ids = NULL, what = "SBS", context = T) {

  a = x$input$counts

  a = a %>% dplyr::mutate(sbs=rownames(a)) %>% as_tibble() %>%
    reshape2::melt() %>%
    reformat_contexts(what=what)

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
  if(reconstructed) {
    sigs = get_signatures(x)
    expos = get_exposure(x)
    sorted_names = get_signames(x) %>% sort
    return((as.matrix(rowSums(x$fit$x) * expos[,sorted_names]) %*% as.matrix(sigs[sorted_names,])) %>% as.data.frame())
  }
  return(x$fit$x)
}


get_epsilon = function(x) {
  if (have_epsilon(x))
    return(
      sapply(colnames(x$fit$eps_var), function(cc)
        lapply(rownames(x$fit$eps_var), function(nn)
          rhalfnorm(1, theta=sd2theta(x$fit$eps_var[nn, cc])) %>% round()) %>%
          setNames(rownames(x$fit$eps_var)) %>% unlist()
        ) %>% as_tibble()
      )

  tmp = get_data(x)
  tmp[tmp > 0] = 0
  return(tmp)
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

plot_mutations = function(x, sampleIDs=NULL, by_sig=FALSE,
                          reconstructed=FALSE, epsilon=FALSE,
                          what="SBS", cls=NULL) {
  if (is.null(sampleIDs)) sampleIDs = rownames(x$fit$x)
  if (is.null(cls) && have_color_palette(x)) cls = get_color_palette(x)
  if (is.null(cls) && !have_color_palette(x)) cls = get_signature_colors(x)

  xx = get_data(x, reconstructed=reconstructed)
  if (have_groups(x)) gid = x$groups else gid = 1

  x_eps = get_data(x); x_eps[x_eps > 0] = 0
  if (epsilon) x_eps = get_epsilon(x)
  x_eps = x_eps %>% tibble::rownames_to_column(var="sampleID") %>% dplyr::mutate(sbs="epsilon", groups=gid)

  xx_s = xx %>% tibble::rownames_to_column(var="sampleID") %>% dplyr::mutate(sbs="s1", groups=gid)

  if (by_sig) {
    rownm = rownames(get_exposure(x))
    xx_s = lapply(rownames(get_signatures(x)),
                  function(sname) {
                    ((as.matrix(x$fit$exposure[, sname], ncol=1) %*%
                        as.matrix(get_signatures(x)[sname,], nrow=1)) * rowSums(x$fit$x)) %>%
                      as.data.frame() %>%
                      dplyr::mutate(sbs=sname, groups=gid, sampleID=rownm)
                  }
    ) %>% do.call(what=rbind, args=.)
  }

  xx_s = xx_s %>%
    dplyr::filter(sampleID %in% sampleIDs) %>%
    reshape2::melt(id=c("sampleID","groups","sbs"), variable.name="variable", value.name="mut_count") %>%
    reformat_contexts(what=what) %>%
    dplyr::rename(sig=Var1) %>%
    # dplyr::filter(mut_count>0) %>%
    # tidyr::separate("subs", into=c("n1","subs.n2"), sep="[[]") %>%
    # tidyr::separate("subs.n2", into=c("subs","n2"), sep="]") %>%
    # dplyr::mutate(context=paste0(n1,"_",n2)) %>%

    dplyr::group_by(substitution, context, sig, groups) %>%
    dplyr::summarise(tot_muts=sum(mut_count)) %>%
    dplyr::ungroup() %>%

    dplyr::add_row(
      x_eps %>%
        dplyr::filter(sampleID %in% sampleIDs) %>%
        reshape2::melt(id=c("sampleID","groups","sbs"), variable.name="variable", value.name="tot_muts") %>%
        reformat_contexts(what=what) %>%
        dplyr::rename(sig=Var1) %>%
        dplyr::select(substitution, context, sig, groups, tot_muts)

        # tidyr::separate("subs", into=c("n1","subs.n2"), sep="[[]") %>%
        # tidyr::separate("subs.n2", into=c("subs","n2"), sep="]") %>%
        # dplyr::mutate(context=paste0(n1,"_",n2)) %>% dplyr::select(sig, context, subs, groups, tot_muts)
    )

  if (by_sig)
    p = xx_s %>%
      ggplot() +
      geom_histogram(aes(x=context, y=tot_muts, fill=sig), stat="identity", position="stack") +
      facet_grid(~substitution) +
      ylab("Number of mutations") + xlab("") +
      theme_bw() + theme(axis.text.x=element_text(angle=90)) +
      scale_fill_manual(values=c(cls, "epsilon"="gainsboro"))
  else
    p = xx_s %>%
      ggplot() +
      geom_histogram(aes(x=context, y=tot_muts, fill=sig), stat="identity") +
      facet_grid(~substitution) +
      ylab("Number of mutations") + xlab("") +
      theme_bw() + theme(axis.text.x=element_text(angle=90)) +
      scale_fill_manual(values=c("epsilon"="gainsboro", "s1"="grey5")) +
      guides(fill="none")

  if (have_groups(x)) return(p + facet_grid(groups~substitution))

  return(p)
}




reformat_contexts = function(a, what) {
  if(what == "SBS") {
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

