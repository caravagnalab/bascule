


# plot exposure-----------------------------------------------------------------

plot.alpha <- function(exp_alpha, inf_alpha) {
  
  if (!identical(rownames(exp_alpha), rownames(inf_alpha))) {
    stop("expected and inferred exposure have no same sample names")
  }
  
  #rownames(exp_alpha) <- rownames(inf_alpha)  # just to be consistent, should be fixed later
  exp_alpha$sample <- rownames(exp_alpha)
  exp_alpha$sample <- factor(exp_alpha$sample, levels = as.character(1:length(exp_alpha$sample)))
  exp_alpha_long <- tidyr::gather(exp_alpha, key="signature", value="exposure", c(-sample))
  exp_alpha_long$type <- rep('expected', each=nrow(exp_alpha_long))
  
  #inf_alpha <- x$exposure[[1]]
  inf_alpha$sample <- rownames(inf_alpha)
  inf_alpha$sample <- factor(inf_alpha$sample, levels = as.character(1:length(inf_alpha$sample)))
  inf_alpha_long <- tidyr::gather(inf_alpha, key="signature", value="exposure", c(-sample))
  inf_alpha_long$type <- rep('inferred', each=nrow(inf_alpha_long))
  
  alpha <- rbind(exp_alpha_long, inf_alpha_long)
  
  # TEST -----------------------------------------------------------------------
  cols <- unique(alpha$signature)
  exp_denovo <- cols[stringr::str_detect(cols, '_D', negate = FALSE)]
  aa <- setdiff(cols, exp_denovo)
  fixed <- aa[stringr::str_detect(aa, 'SBS', negate = FALSE)]
  inf_denovo <- setdiff(aa, fixed)
  new_cols <- c(fixed, exp_denovo, inf_denovo)
  alpha$signature <- factor(alpha$signature, levels = new_cols)
  # TEST -----------------------------------------------------------------------
  
  plt <- ggplot(data = alpha, aes(x=sample, y=exposure, fill=signature)) +
    geom_bar(stat = "identity") +
    facet_grid(type ~ .) +
    theme_minimal() +
    ggtitle("Signatures exposure (Expected vs. Inferred)") +
    scale_fill_brewer(palette='Spectral')
    #ggsci::scale_fill_lancet()
    #scale_y_continuous(labels=scales::percent)
  
  return(plt)
}


# plot signatures --------------------------------------------------------------

plot.beta <- function(beta) {
  
  if (is.null(beta)) {
    p <- ggplot() +
      theme_void() +
      geom_text(aes(0,0,label='N/A')) +
      xlab(NULL) #optional, but safer in case another theme is applied later
  } else {
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
      
      #if(!xlabels) {
      plt <- plt + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank())
      #}
      
      glist[[i]] <- plt
    }
    
    # make the final plot
    #gridExtra::grid.arrange(grobs=glist,ncol=ceiling(nrow(beta)/3))
    
    p <- ggpubr::ggarrange(plotlist=glist, ncol = 1)
  }
  
  return(list(p))
}

# plot signatures cosine matrix ------------------------------------------------

plot.beta.cosine <- function(exp_denovo, inf_denovo) {
  
  if (is.null(exp_denovo) | is.null(inf_denovo)) {
    cplot <- ggplot() +
      theme_void() +
      geom_text(aes(0,0,label='N/A')) +
      xlab(NULL) #optional, but safer in case another theme is applied later
  } else {
    # expected vs inferred signatures cosine similarity matrix
    cos <- basilica:::cosine.matrix(inf_denovo, exp_denovo)
    cos1 <- tibble::rownames_to_column(cos, var = 'inferred_denovo')
    cos_long <- tidyr::gather(cos1, key="expected_denovo", value="cosine_similarity", c(-inferred_denovo))
    # plot data
    cplot <- ggplot(cos_long, aes(expected_denovo, inferred_denovo)) +
      geom_tile(aes(fill = cosine_similarity)) +
      geom_text(aes(label = round(cosine_similarity, 3))) +
      scale_fill_gradient(low = "white", high = "darkgreen") +
      ggtitle("Cosine similarity matrix (expected vs. inferred)") +
      xlab("Expected") +
      ylab("Inferred")
  }
  
  return(list(cplot))
}

#-------------------------------------------------------------------------------

multi.plot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

# ------------------------------------------------------------------------------

final.plot <- function(exp_alpha, inf_alpha, exp_denovo, inf_denovo) {
  
  alpha <- basilica:::plot.alpha(exp_alpha, inf_alpha)
  beta.cosine <- plot.beta.cosine(exp_denovo, inf_denovo)
  p1 <- list(ggpubr::ggarrange(plotlist=c(alpha, beta.cosine), ncol = 2))
  
  exp.beta <- plot.beta(exp_denovo)
  #ggpubr::annotate_figure(exp.beta, top = ggpubr::text_grob("Expected Signatures", color = "red", face = "bold", size = 14))
  inf.beta <- plot.beta(inf_denovo)
  #ggpubr::annotate_figure(inf.beta, top = ggpubr::text_grob("Inferred Signatures", color = "red", face = "bold", size = 14))
  p2 <- list(ggpubr::ggarrange(plotlist=c(exp.beta, inf.beta), ncol = 2))
  p <- c(p1, p2)
  
  return(basilica:::multi.plot(plotlist = p))
}




