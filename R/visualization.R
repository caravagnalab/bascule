#' plot exposure matrix
#'
#' @description creates bar plot of relative exposure matrix, where x-axis are samples and y-axis are their relative contribution.
#' @param x basilica object
#'
#' @return
#' @export plot_exposure
#'
#' @examples
plot_exposure <- function(x) {
  alpha <- get_exposure(x, long = FALSE)
  alpha$Sample <- rownames(alpha)
  alpha_long <- tidyr::gather(alpha,
                       key="Signature",
                       value="Exposure",
                       c(-Sample)
  )

  ggplot(data = alpha_long, aes(x=Sample, y=Exposure, fill=Signature)) +
    geom_bar(stat = "identity") +
    theme_minimal() +
    ggtitle("Signatures Exposure in Samples") +
    scale_y_continuous(labels=scales::percent)
}


#' plot signatures
#'
#' @description creates bar plot of inferred signature profiles, where x-axis are 96 substitution bases and y-axis are their relative contribution.
#'
#' @param x basilica object
#' @param useRowNames using signature names from data.frame
#' @param xlabels axis label
#' @param denovoSignature if TRUE, plots inferred de-novo signatures, otherwise plots inferred catalogue signatures
#'
#' @return
#' @export plot_signatures
#'
#' @import ggplot2
#' @import data.table
#'
#' @examples
plot_signatures <- function( x, useRowNames = TRUE, xlabels = FALSE, denovoSignature = TRUE ) {

  if (denovoSignature==TRUE) {
    beta <- get_denovo_signatures(x)
  } else {
    beta <- get_catalogue_signatures(x)
  }

  # set names of the signatures
  if(!useRowNames) {
    rownames(beta) <- paste0("Signature ", 1:nrow(beta))
  }

  # separate context and alteration
  x <- data.table::as.data.table(reshape2::melt(as.matrix(beta),varnames=c("signature","cat")))
  x[, Context := paste0(substr(cat,1,1), ".", substr(cat, 7, 7)) ]
  x[, alt := paste0(substr(cat,3,3),">",substr(cat,5,5)) ]

  # make the ggplot2 object
  glist <- list()
  for(i in 1:nrow(beta)) {

    plt <- ggplot(x[signature==rownames(beta)[i]]) +
      geom_bar(aes(x=Context,y=value,fill=alt),stat="identity",position="identity") +
      facet_wrap(~alt,nrow=1,scales="free_x") +
      theme(axis.text.x=element_text(angle=90,hjust=1),panel.background=element_blank(),axis.line=element_line(colour="black")) +
      ggtitle(rownames(beta)[i]) + theme(legend.position="none") + ylab("Frequency of mutations")

    if(!xlabels) {
      plt <- plt + theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())
    }

    glist[[i]] <- plt

  }

  # make the final plot
  gridExtra::grid.arrange(grobs=glist,ncol=ceiling(nrow(beta)/3))

}


plot_cosine <- function(obj, limit) {
  reference <- obj$reference_catalogue
  denovo <- obj$denovo_signatures

  cosine_matrix <- data.frame( matrix( nrow = nrow(reference), ncol = nrow(denovo) ) )
  colnames(cosine_matrix) <- rownames(denovo)
  rownames(cosine_matrix) <- rownames(reference)

  for (i in rownames(reference)) {
    for (j in rownames(denovo)) {
      c <- cosine_sim(as.numeric(reference[i, ]), as.numeric(denovo[j, ]))
      cosine_matrix[i, j] <- c
    }
  }
  cosine_matrix <- cosine_matrix %>% dplyr::filter_all(dplyr::any_vars(. > limit))
  cosine_matrix <- tibble::rownames_to_column(cosine_matrix, var="ref")

  cos_long <- cosine_matrix %>% tidyr::pivot_longer(cols = -c(ref), names_to = 'denovo', values_to = 'cosine')

  p <- ggplot(cos_long, aes(x=ref, y=cosine)) +
    geom_point() +
    #facet_wrap(~denovo, ncol = 1) +
    facet_grid(denovo ~ .) +
    theme(axis.text.x=element_text(angle=-90, vjust=0.5, hjust=0)) +
    ggtitle("Similarity between Denovo and Catalogue Signatures") +
    xlab("Catalogue Signatures") +
    ylab("Cosine Similarity") +
    geom_hline(yintercept=c(limit), linetype='dashed', color='red')

  return(p)
}



