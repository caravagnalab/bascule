
get.one.sample <- function(
    synthetic,
    k,
    lr,
    steps,
    phi,
    delta
    ) {

  x <- synthetic$x[[1]]
  ref <- synthetic$ref_cat[[1]]
  input <- synthetic$input_cat[[1]]

  obj <- basilica:::fit(
    x=x,
    reference_catalogue=ref,
    k=k,
    lr=lr,
    steps=steps,
    phi=phi,
    delta=delta,
    groups=NULL,
    input_catalogue=input
  )

  obj$exposure <- list(obj$exposure)
  obj$denovo_signatures <- list(obj$denovo_signatures)
  obj$catalogue_signatures <- list(obj$catalogue_signatures)

  results <- c(synthetic, obj)

  return(results)
}


# plot exposure---------------------------------------------------------------

plot.alpha <- function(exp_alpha, inf_alpha) {
  #exp_alpha <- x$exp_exposure[[1]]
  rownames(exp_alpha) <- rownames(inf_alpha)  # just to be consistent, should be fixed later
  exp_alpha$sample <- rownames(exp_alpha)
  exp_alpha_long <- tidyr::gather(exp_alpha, key="signature", value="exposure", c(-sample))
  exp_alpha_long$type <- rep('expected', each=nrow(exp_alpha_long))

  #inf_alpha <- x$exposure[[1]]
  inf_alpha$sample <- rownames(inf_alpha)
  inf_alpha_long <- tidyr::gather(inf_alpha, key="signature", value="exposure", c(-sample))
  inf_alpha_long$type <- rep('inferred', each=nrow(inf_alpha_long))

  alpha <- rbind(exp_alpha_long, inf_alpha_long)

  plt <- ggplot(data = alpha, aes(x=sample, y=exposure, fill=signature)) +
    geom_bar(stat = "identity") +
    facet_grid(type ~ .) +
    theme_minimal() +
    ggtitle("Signatures exposure (Expected vs. Inferred)")
  #scale_y_continuous(labels=scales::percent)

  #glist <-
  #glist[[1]] <- plt

  return(list(plt))
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
      plt <- plt + theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())
      #}

      glist[[i]] <- plt
    }

    # make the final plot
    #gridExtra::grid.arrange(grobs=glist,ncol=ceiling(nrow(beta)/3))

    p <- ggpubr::ggarrange(plotlist=glist, ncol = 1)

  }

  return(list(p))
}

# plot signatures cosine matrix --------------------------------------------------------------


plot.beta.cosine <- function(exp_denovo, inf_denovo) {
  # expected vs inferred signatures cosine similarity matrix
  cos <- cosine_matrix(inf_denovo, exp_denovo)
  cos1 <- tibble::rownames_to_column(cos, var = 'inferred_denovo')
  cos_long <- tidyr::gather(cos1, key="expected_denovo", value="cosine_similarity", c(-inferred_denovo))
  # plot data
  cplot <- ggplot(cos_long, aes(expected_denovo, inferred_denovo)) +
    geom_tile(aes(fill = cosine_similarity)) +
    geom_text(aes(label = round(cosine_similarity, 3))) +
    scale_fill_gradient(low = "white", high = "darkgreen")

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
  inf.beta <- plot.beta(inf_denovo)
  p2 <- list(ggpubr::ggarrange(plotlist=c(exp.beta, inf.beta), ncol = 2))
  p <- c(p1, p2)

  return(basilica:::multi.plot(plotlist = p))
}


#-------------------------------------------------------------------------------

#' @import dplyr
fill_tibble <- function(x) {

  data <- tibble::tibble(
    x = list(),
    Input_Catalogue = list(),
    Ref_Catalogue = list(),

    Exp_Exposure = list(),
    Exp_Fixed = list(),
    Exp_Denovo = list(),

    Inf_Exposure = list(),
    Inf_Fixed = list(),
    Inf_Denovo = list(),

    TargetX = character(),
    InputX = character(),
    Num_Samples = numeric(),
    IterNum = numeric(),

    K = list(),
    Lr = numeric(),
    Steps = numeric(),
    Phi = numeric(),
    Delta = numeric()
  )

  for (row in x) {
    data <- data %>% tibble::add_row(
      x = list(row$x),
      Input_Catalogue = list(row$input_catalogue),
      Ref_Catalogue = list(row$ref_catalogue),

      Exp_Exposure = list(row$exp_exposure),
      Exp_Fixed = list(row$exp_fixed),
      Exp_Denovo = list(row$exp_denovo),

      Inf_Exposure = list(row$inf_exposure),
      Inf_Fixed = list(row$inf_fixed),
      Inf_Denovo = list(row$inf_denovo),

      TargetX = row$targetX,
      InputX = row$inputX,
      Num_Samples = row$num_samples,
      IterNum = row$iter,

      K = list(row$k),
      Lr = row$lr,
      Steps = row$steps,
      Phi = row$phi,
      Delta = row$delta
    )
  }
  return(data)
}

#-------------------------------------------------------------------------------

visData <- function(x) {

  df <- tibble::tibble(
    TargetX = character(),
    InputX = character(),
    N_Samples = numeric(),
    Fixed_TP = numeric(),
    Fixed_FP = numeric(),
    Fixed_TN = numeric(),
    Fixed_FN = numeric(),
    Fixed_TPR = numeric(),
    Fixed_Prec = numeric(),
    Fixed_Rec = numeric(),
    Fixed_Acc = numeric(),
    Denovo_Ratio = numeric(),
    Fitness = numeric(),
    MAE = numeric(),
    MSE = numeric()
  )

  for (i in 1:nrow(x)) {
    input_catalogue <- x[i, ]$Input_Catalogue[[1]]
    inferred_alpha <- x[i, ]$Inf_Exposure[[1]]
    inferred_fixed <- x[i, ]$Inf_Fixed[[1]]
    inferred_denovo <- x[i, ]$Inf_Denovo[[1]]
    expected_fixed <- x[i, ]$Exp_Fixed[[1]]
    expected_denovo <- x[i, ]$Exp_Denovo[[1]]
    m <- x[i, ]$x[[1]]
    mr <- reconstruction_matrix(m, inferred_alpha, inferred_fixed, inferred_denovo)

    # fixed signatures
    if (is.null(input_catalogue)) {inp <- c()} else {inp <- rownames(input_catalogue)}
    if (is.null(expected_fixed)) {exp <- c()} else {exp <- rownames(expected_fixed)}
    if (is.null(inferred_fixed)) {inf <- c()} else {inf <- rownames(inferred_fixed)}
    fixed_TP <- length(intersect(inf, exp))
    fixed_FP <- length(setdiff(inf, exp))
    fixed_TN <- length(setdiff(setdiff(inp, exp), inf))
    fixed_FN <- length(setdiff(exp, inf))
    if (length(exp)==0) {fixed_TPR <- fixed_TP+1} else {fixed_TPR <- fixed_TP / length(exp)}
    fixed_Prec <- fixed_TP / (fixed_TP + fixed_FP)
    fixed_Rec <- fixed_TP / (fixed_TP + fixed_FN)
    fixed_Acc <- (fixed_TP + fixed_TN) / (fixed_TP + fixed_TN + fixed_FP + fixed_FN)

    # denovo signatures
    if (is.null(expected_denovo)) {n_exp_denovo <- 0} else {n_exp_denovo <- nrow(expected_denovo)}
    if (is.null(inferred_denovo)) {n_inf_denovo <- 0} else {n_inf_denovo <- nrow(inferred_denovo)}
    denovo_Ratio <- (n_inf_denovo + 1) / (n_exp_denovo + 1)

    # goodness of fitness
    fitness <- fitness.quality(m, mr)
    mae <- compute.mae(m, mr)
    mse <- compute.mse(m, mr)


    df <- df %>% tibble::add_row(
      TargetX = x[i, ]$TargetX,
      InputX = x[i, ]$InputX,
      N_Samples = x[i, ]$Num_Samples,
      Fixed_TP = fixed_TP,
      Fixed_FP = fixed_FP,
      Fixed_TN = fixed_TN,
      Fixed_FN = fixed_FN,
      Fixed_TPR = fixed_TPR,
      Fixed_Prec = fixed_Prec,
      Fixed_Rec = fixed_Rec,
      Fixed_Acc = fixed_Acc,
      Denovo_Ratio = denovo_Ratio,
      Fitness = fitness,
      MAE = mae,
      MSE = mse
    )
  }

  df$TargetX <- factor(df$TargetX, levels = c("low", "medium", "high"))
  df$InputX <- factor(df$InputX, levels = c("low", "medium", "high"))

  return(df)
}


#-------------------------------------------------------------------------------

reconstruction_matrix <- function(m, alpha, fixed, denovo) {
  # all args are data.frame
  theta <- diag(rowSums(m)) # matrix
  alpha <- theta %*% as.matrix(alpha) # matrix
  beta <- as.matrix(rbind(fixed, denovo)) # matrix

  mr_matrix <- alpha %*% as.matrix(beta)
  mr <- round(as.data.frame(mr_matrix))
  return(mr)
}

#-------------------------------------------------------------------------------

fitness.quality <- function(m, mr) {
  total <- 0
  for (i in 1:nrow(m)) {
    cos <- cosine_sim(m[i,], mr[i, ])
    total <- total + cos
  }
  return(total / nrow(m))
}

#-------------------------------------------------------------------------------

compute.mae <- function(m , mr) {
  mae <- sum(abs(m - mr)) / (dim(m)[1] * dim(m)[2])
  return(mae)
}

#-------------------------------------------------------------------------------

compute.mse <- function(m , mr) {
  mse <- sum((m - mr)^2) / (dim(m)[1] * dim(m)[2])
  return(mse)
}

#-------------------------------------------------------------------------------

denovo.quality <- function(exp, inf) {

  if (length(exp)==0 | length(inf)==0) {
    return(list(matrix=0, ratio=0, cosine=0))
  } else {
    df <- data.frame(matrix(nrow = nrow(inf), ncol = nrow(exp)))
    colnames(df) <- rownames(exp)
    rownames(df) <- rownames(inf)

    for (i in 1:nrow(inf)) {
      inferred <- inf[i,]
      inferred_name <- rownames(inferred)
      maxScore <- 0
      bestMatch <- NULL
      for (j in 1:nrow(exp)) {
        target <- exp[j, ]
        target_name <- rownames(target)
        score <- cosine_sim(inferred, target)
        df[inferred_name, target_name] <- score
        '
        if (score > maxScore) {
          maxScore <- score
          bestMatch <- target_name
        }
        '
      }
    }
    '
    match_list <- colnames(df)[max.col(df, ties.method = "first")]
    s <- 0
    for (i in 1:length(match_list)) {
      s <- s + df[rownames(df)[i], match_list[i]]
    }
    cosine_score <- (s / length(match_list))
    '
  }

  a <- nrow(inf) / nrow(exp)
  b <- mean(apply(df, 1, max))
  return(list(matrix=df, ratio=a, cosine=b))
}

#-------------------------------------------------------------------------------






