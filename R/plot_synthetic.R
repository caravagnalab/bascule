

# plot exposure--------------------------------------------------------QC:almost

.plot.synthetic.alpha <- function(exp_alpha, inf_alpha) {

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


# plot signatures cosine matrix ------------------------------------------------

.plot.synthetic.beta.similarity <- function(exp_denovo, inf_denovo) {

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
    cplot <- ggplot(data = cos_long, aes(x=expected_denovo, y=inferred_denovo, fill=cosine_similarity)) +
      geom_tile() +
      geom_text(aes(label = round(cosine_similarity, 3))) +
      #scale_fill_gradient(low = "orange", high = "green") +
      scale_fill_gradient2(low = "orange",
                           mid = "white",
                           high = "purple") +
      ggtitle("Cosine similarity matrix (expected vs. inferred)") +
      xlab("Expected") +
      ylab("Inferred")
  }

  return(cplot)
}

#-------------------------------------------------------------------------------




