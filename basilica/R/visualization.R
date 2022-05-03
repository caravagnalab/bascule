library(tidyr)      # tibble() & gather()
library(ggplot2)    # ggplot()
library(data.table) # data.table()
library(gridExtra)  # grid.arrange


#' Title
#'
#' @param x basilica object
#'
#' @return plot exposure matrix
#' @export
#'
#' @examples
plot_exposure <- function(x) {
  alpha <- x$exposure
  alpha$Branch <- c(1:nrow(alpha))
  alpha_long <- gather(alpha,
                       key="Signature",
                       value="Exposure",
                       c(-Branch)
  )
  
  ggplot(data = alpha_long, aes(x=Branch, y=Exposure, fill=Signature)) +
    geom_bar(stat = "identity") +
    theme_minimal() +
    ggtitle("Signatures Exposure in Samples") +
    scale_y_continuous(labels=scales::percent)
}



#' Title
#'
#' @param x 
#' @param useRowNames 
#' @param xlabels 
#' @param denovoSignature 
#'
#' @return
#' @export
#'
#' @examples
plot_signatures <- function( x, useRowNames = FALSE, xlabels = FALSE, denovoSignature = TRUE ) {
  
  if (denovoSignature==TRUE) {
    beta <- x$denovo_signatures
  } else {
    beta <- x$catalog_signatures
  }
  
  # set names of the signatures
  if(!useRowNames) {
    rownames(beta) <- paste0("Signature ",1:nrow(beta))
  }
  
  # separate context and alteration
  x <- as.data.table(reshape2::melt(as.matrix(beta),varnames=c("signature","cat")))
  x[,Context:=paste0(substr(cat,1,1),".",substr(cat,7,7))]
  x[,alt:=paste0(substr(cat,3,3),">",substr(cat,5,5))]
  
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
  grid.arrange(grobs=glist,ncol=ceiling(nrow(beta)/3))
  
}



